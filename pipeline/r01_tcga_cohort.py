#!/usr/bin/env python3
"""
r01_tcga_cohort.py
Build the TCGA-GBM sample manifest and cohort flow (PREREGISTRATION.md, Section 3).

Uses file metadata, clinical records and mutation calls only; no expression
values are read.

Inputs
  data/downloads/gdc_meta/expression_files_meta.csv   (r00)
  data/downloads/gdc_meta/tcga_gbm_cases.json         (r00)
  data/downloads/gdc_maf/*.maf.gz, maf_files_meta.json (r00)
  data/downloads/cbioportal/*.json                     (r00)
  data/downloads/gdc_methylation/files_meta.json       (r00; file metadata only)
Outputs
  revision/manifests/tcga_samples.csv
  revision/results/cohort/tcga_cohort_flow.csv
  revision/results/cohort/idh_source_concordance.csv
  revision/source_data/FigS1_cohort_flow.csv
  revision/figures/FigS1_cohort_flow.png

Run from the repository root:
    python pipeline/r01_tcga_cohort.py
"""

import gzip
import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DL = Path("data/downloads")
MAN = Path("revision/manifests")
RES = Path("revision/results/cohort")
SRC = Path("revision/source_data")
FIG = Path("revision/figures")

BATCH = {"01": "BI", "41": "GDC2024"}


# ── expression file metadata ─────────────────────────────────────────────────

def load_files():
    f = pd.read_csv(DL / "gdc_meta/expression_files_meta.csv")
    f["center_code"] = f["aliquot"].str[-2:]
    f["batch"] = f["center_code"].map(BATCH).fillna("other_" + f["center_code"])
    f["sample_code"] = f["aliquot"].str[13:15]
    f["vial"] = f["aliquot"].str[15]
    f["portion"] = f["aliquot"].str[17:19].astype(int)
    f["plate"] = f["aliquot"].str[21:25]
    return f


def select_one_per_patient(prim):
    """Fixed rule: BI batch first, then lowest vial, portion, plate, full barcode."""
    order = prim.assign(batch_rank=(prim["batch"] != "BI").astype(int)).sort_values(
        ["patient", "batch_rank", "vial", "portion", "plate", "aliquot"])
    return set(order.groupby("patient").head(1)["aliquot"])


# ── clinical ─────────────────────────────────────────────────────────────────

def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def load_gdc_clinical():
    cases = json.load(open(DL / "gdc_meta/tcga_gbm_cases.json"))["hits"]
    rows = []
    for c in cases:
        demo = c.get("demographic") or {}
        diags = c.get("diagnoses") or []
        fups = c.get("follow_ups") or []
        prim = [d for d in diags if d.get("diagnosis_is_primary_disease") is True] or \
               [d for d in diags if _num(d.get("days_to_diagnosis")) == 0] or diags[:1]
        pd0 = prim[0] if prim else {}
        lfu = [_num(d.get("days_to_last_follow_up")) for d in diags] + \
              [_num(f.get("days_to_follow_up")) for f in fups]
        lfu = [x for x in lfu if not np.isnan(x)]
        age_days = _num(pd0.get("age_at_diagnosis"))
        rows.append({
            "patient": c["submitter_id"],
            "vital_status": demo.get("vital_status"),
            "days_to_death": _num(demo.get("days_to_death")),
            "days_to_last_follow_up": max(lfu) if lfu else np.nan,
            "age_at_diagnosis_years": age_days / 365.25 if not np.isnan(age_days) else _num(demo.get("age_at_index")),
            "age_source": "diagnosis.age_at_diagnosis" if not np.isnan(age_days) else "demographic.age_at_index",
            "sex": demo.get("sex_at_birth") or demo.get("gender"),
            "gdc_primary_diagnosis": pd0.get("primary_diagnosis"),
            "gdc_morphology": pd0.get("morphology"),
            "gdc_residual_disease": pd0.get("residual_disease"),
        })
    clin = pd.DataFrame(rows).drop_duplicates("patient")
    dead = clin["vital_status"] == "Dead"
    alive = clin["vital_status"] == "Alive"
    clin["os_days"] = np.where(dead, clin["days_to_death"], np.where(alive, clin["days_to_last_follow_up"], np.nan))
    clin["os_event"] = np.where(dead, 1, np.where(alive, 0, np.nan))
    clin["os_source"] = np.where(dead, "days_to_death", np.where(alive, "max last follow-up", "vital status missing"))
    return clin


def load_cbio(study, level):
    d = json.load(open(DL / f"cbioportal/{study}_{level}_clinical.json"))
    key = "patientId" if level == "patient" else "sampleId"
    out = defaultdict(dict)
    for r in d:
        out[r[key]][r["clinicalAttributeId"]] = r["value"]
    return pd.DataFrame.from_dict(out, orient="index")


def cbio_patient_level(df_sample):
    """Collapse sample-level annotations to patient, using primary tumor (-01) samples."""
    s = df_sample[df_sample.index.str[13:15] == "01"].copy()
    s["patient"] = s.index.str[:12]
    return s.groupby("patient").first()


def egfr_gistic():
    """EGFR GISTIC call per patient (primary tumor): PanCanAtlas, then Brennan 2013."""
    out = {}
    for study in ["gbm_tcga_pub2013", "gbm_tcga_pan_can_atlas_2018"]:  # later overrides earlier
        for r in json.load(open(DL / f"cbioportal/{study}_EGFR_gistic.json")):
            if r["sampleId"][13:15] == "01":
                out[r["sampleId"][:12]] = (int(r["value"]), study)
    return pd.DataFrame([(k, *v) for k, v in out.items()], columns=["patient", "egfr_gistic", "egfr_gistic_source"])


def methylation_platforms():
    """Methylation array platform(s) with an open beta-value file, primary tumor, per patient."""
    hits = json.load(open(DL / "gdc_methylation/files_meta.json"))["hits"]
    plat = defaultdict(set)
    for h in hits:
        if h["cases"][0]["samples"][0]["sample_type"] == "Primary Tumor":
            plat[h["cases"][0]["submitter_id"]].add(h["platform"].replace("Illumina Human Methylation ", "HM"))
    return pd.DataFrame({"patient": list(plat), "methylation_platforms": ["+".join(sorted(v)) for v in plat.values()]})


# ── IDH ──────────────────────────────────────────────────────────────────────

def idh_from_maf():
    meta = json.load(open(DL / "gdc_maf/maf_files_meta.json"))["hits"]
    status = {}
    for h in meta:
        tumor = [s for s in h["cases"][0]["samples"] if s["sample_type"] == "Primary Tumor"]
        if not tumor:
            continue
        pt = h["cases"][0]["submitter_id"]
        mut = []
        with gzip.open(DL / "gdc_maf" / h["file_name"], "rt") as f:
            header = None
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if header is None:
                    header = {k: i for i, k in enumerate(fields)}
                    continue
                gene, prot = fields[header["Hugo_Symbol"]], fields[header["HGVSp_Short"]]
                if (gene == "IDH1" and prot.startswith("p.R132")) or (gene == "IDH2" and prot.startswith("p.R172")):
                    mut.append(f"{gene} {prot}")
        prev = status.get(pt, [])
        status[pt] = prev + mut
    return pd.DataFrame({"patient": list(status),
                         "idh_maf": ["mutant" if v else "wildtype" for v in status.values()],
                         "idh_maf_variants": [";".join(sorted(set(v))) for v in status.values()]})


def main():
    for d in (MAN, RES, SRC, FIG):
        d.mkdir(parents=True, exist_ok=True)

    files = load_files()
    clin = load_gdc_clinical()
    idh_maf = idh_from_maf()

    cec = cbio_patient_level(load_cbio("lgggbm_tcga_pub", "sample")).join(
        load_cbio("lgggbm_tcga_pub", "patient"), rsuffix="_pt")
    bre = cbio_patient_level(load_cbio("gbm_tcga_pub2013", "sample")).join(
        load_cbio("gbm_tcga_pub2013", "patient"), rsuffix="_pt")
    pan = load_cbio("gbm_tcga_pan_can_atlas_2018", "patient")

    ext = pd.DataFrame(index=sorted(set(files["patient"])))
    ext.index.name = "patient"
    ext = ext.join(idh_maf.set_index("patient"))
    ext["idh_ceccarelli2016"] = cec["IDH_STATUS"].map({"WT": "wildtype", "Mutant": "mutant"}).reindex(ext.index)
    ext["idh_brennan2013"] = bre["IDH1_MUTATION"].map(lambda v: "wildtype" if v == "WT" else ("mutant" if isinstance(v, str) else np.nan)).reindex(ext.index)

    def resolve(r):
        if isinstance(r["idh_maf"], str):
            return r["idh_maf"], "GDC open MAF"
        if isinstance(r["idh_ceccarelli2016"], str):
            return r["idh_ceccarelli2016"], "Ceccarelli 2016"
        if isinstance(r["idh_brennan2013"], str):
            return r["idh_brennan2013"], "Brennan 2013"
        return "unknown", "none"
    ext[["idh_status", "idh_source"]] = ext.apply(lambda r: pd.Series(resolve(r)), axis=1)
    calls = ext[["idh_maf", "idh_ceccarelli2016", "idh_brennan2013"]]
    ext["idh_source_disagreement"] = calls.nunique(axis=1, dropna=True) > 1

    ext["mgmt_promoter"] = cec["MGMT_PROMOTER_STATUS"].reindex(ext.index).str.lower()
    ext["mgmt_promoter_source"] = np.where(ext["mgmt_promoter"].notna(), "Ceccarelli 2016", "")
    bre_m = bre["MGMT_STATUS"].reindex(ext.index).str.lower()
    fill = ext["mgmt_promoter"].isna() & bre_m.notna()
    ext.loc[fill, "mgmt_promoter"] = bre_m[fill]
    ext.loc[fill, "mgmt_promoter_source"] = "Brennan 2013"
    ext["kps"] = pd.to_numeric(cec.get("KARNOFSKY_PERFORMANCE_SCORE"), errors="coerce").reindex(ext.index)
    ext["gcimp"] = bre["G_CIMP_METHYLATION"].reindex(ext.index)
    ext["expression_subtype_brennan2013"] = bre["EXPRESSION_SUBTYPE"].reindex(ext.index)
    ext["histology_ceccarelli2016"] = cec["HISTOLOGICAL_DIAGNOSIS"].reindex(ext.index)
    ext["therapy_brennan2013"] = bre["THERAPY"].reindex(ext.index)
    ext["tmz"] = ext["therapy_brennan2013"].map(lambda v: ("yes" if "TMZ" in v else "no") if isinstance(v, str) else np.nan)
    ext["radiation"] = pan["RADIATION_THERAPY"].str.lower().reindex(ext.index)
    ext["tert_promoter_ceccarelli2016"] = cec["TERT_PROMOTER_STATUS"].reindex(ext.index)
    ext["chr7gain_chr10loss_ceccarelli2016"] = cec["CHR_7_GAIN_CHR_10_LOSS"].reindex(ext.index)
    ext = ext.join(egfr_gistic().set_index("patient"))
    ext["egfr_amplified"] = ext["egfr_gistic"].map(lambda v: np.nan if pd.isna(v) else ("yes" if v == 2 else "no"))
    ext = ext.join(methylation_platforms().set_index("patient"))

    # IDH source concordance (where more than one source exists)
    conc = []
    for a, b in [("idh_maf", "idh_ceccarelli2016"), ("idh_maf", "idh_brennan2013"), ("idh_ceccarelli2016", "idh_brennan2013")]:
        both = ext[[a, b]].dropna()
        conc.append({"source_a": a, "source_b": b, "n_both": len(both),
                     "agree": int((both[a] == both[b]).sum()), "disagree": int((both[a] != both[b]).sum()),
                     "disagree_patients": ";".join(both.index[both[a] != both[b]])})
    pd.DataFrame(conc).to_csv(RES / "idh_source_concordance.csv", index=False)

    # ── manifest ──
    m = files.merge(clin, on="patient", how="left").merge(ext.reset_index(), on="patient", how="left")
    selected = select_one_per_patient(m[m["sample_code"] == "01"])
    m["is_primary"] = m["sample_code"] == "01"
    m["selected_aliquot"] = m["aliquot"].isin(selected)
    m["os_evaluable"] = m["os_days"].notna() & (m["os_days"] > 0) & m["os_event"].notna()

    def reason(r):
        if not r["is_primary"]:
            return f"excluded: sample type {r['sample_type']}"
        if not r["selected_aliquot"]:
            return "excluded: additional aliquot of same patient"
        return "included"
    m["screen_status"] = m.apply(reason, axis=1)
    m["in_screen_sensitivity"] = m["screen_status"] == "included"
    m["in_screen_primary"] = m["in_screen_sensitivity"] & (m["idh_status"] == "wildtype")
    m["in_survival_sensitivity"] = m["in_screen_sensitivity"] & m["os_evaluable"]
    m["in_survival_primary"] = m["in_screen_primary"] & m["os_evaluable"]
    m.loc[m["in_screen_sensitivity"] & (m["idh_status"] != "wildtype"), "screen_status"] = \
        "sensitivity cohort only: IDH " + m["idh_status"]
    m["in_screen_idh_not_mutant"] = m["in_screen_sensitivity"] & (m["idh_status"] != "mutant")
    bb = m[m["is_primary"]].groupby("patient")["batch"].agg(lambda s: {"BI", "GDC2024"} <= set(s))
    m["patient_in_both_batches"] = m["patient"].map(bb).eq(True)
    for f in ["in_screen_primary", "in_screen_sensitivity", "in_screen_idh_not_mutant",
              "in_survival_primary", "in_survival_sensitivity"]:
        assert not m.loc[m[f], "patient"].duplicated().any(), f"duplicate patient in {f}"
    cols = ["aliquot", "patient", "sample_type", "batch", "center_code", "file_id", "file_name", "md5_gdc",
            "gdc_created", "workflow", "is_primary", "selected_aliquot", "screen_status",
            "in_screen_primary", "in_screen_sensitivity", "in_screen_idh_not_mutant",
            "in_survival_primary", "in_survival_sensitivity", "patient_in_both_batches",
            "vital_status", "os_days", "os_event", "os_source", "days_to_death", "days_to_last_follow_up",
            "age_at_diagnosis_years", "age_source", "sex", "kps", "idh_status", "idh_source", "idh_maf",
            "idh_maf_variants", "idh_ceccarelli2016", "idh_brennan2013", "idh_source_disagreement", "mgmt_promoter", "mgmt_promoter_source",
            "gcimp", "expression_subtype_brennan2013", "histology_ceccarelli2016", "gdc_primary_diagnosis",
            "gdc_morphology", "gdc_residual_disease", "therapy_brennan2013", "tmz", "radiation",
            "tert_promoter_ceccarelli2016", "chr7gain_chr10loss_ceccarelli2016", "egfr_gistic", "egfr_gistic_source",
            "egfr_amplified", "methylation_platforms"]
    m[cols].sort_values("aliquot").to_csv(MAN / "tcga_samples.csv", index=False)

    # ── flow ──
    sel = m[m["selected_aliquot"]]
    flow = [
        ("GDC STAR - Counts files (TCGA-GBM)", len(m), m["patient"].nunique()),
        ("Primary tumor (sample type 01)", int(m["is_primary"].sum()), m.loc[m["is_primary"], "patient"].nunique()),
        ("One aliquot per patient", len(sel), sel["patient"].nunique()),
        ("  of which BI batch", int((sel["batch"] == "BI").sum()), None),
        ("  of which GDC2024 batch", int((sel["batch"] == "GDC2024").sum()), None),
        ("  of which other batch", int((~sel["batch"].isin(["BI", "GDC2024"])).sum()), None),
        ("IDH-wildtype (primary cohort, screen)", int(sel["in_screen_primary"].sum()), None),
        ("IDH-mutant (sensitivity only)", int((sel["idh_status"] == "mutant").sum()), None),
        ("IDH unknown (sensitivity only)", int((sel["idh_status"] == "unknown").sum()), None),
        ("Survival evaluable, all primary (sensitivity)", int(sel["in_survival_sensitivity"].sum()), None),
        ("Survival evaluable, IDH-wildtype (primary)", int(sel["in_survival_primary"].sum()), None),
    ]
    fl = pd.DataFrame(flow, columns=["step", "n_samples", "n_patients"])
    fl["n_patients"] = fl["n_patients"].astype("Int64")
    fl.to_csv(RES / "tcga_cohort_flow.csv", index=False)
    fl.to_csv(SRC / "FigS1_cohort_flow.csv", index=False)
    excl = m.loc[~m["is_primary"], "sample_type"].value_counts()
    print(fl.to_string(index=False))
    print("Excluded non-primary:", excl.to_dict())
    print("IDH source among selected:", sel["idh_source"].value_counts().to_dict())
    print(pd.DataFrame(conc)[["source_a", "source_b", "n_both", "agree", "disagree"]].to_string(index=False))
    draw_flow(fl, excl, sel)


def draw_flow(fl, excl, sel):
    n = dict(zip(fl["step"], fl["n_samples"]))
    p = dict(zip(fl["step"], fl["n_patients"]))
    boxes = [
        (0.5, 0.92, f"TCGA-GBM STAR - Counts files (GDC)\n{n['GDC STAR - Counts files (TCGA-GBM)']} aliquots, "
                    f"{p['GDC STAR - Counts files (TCGA-GBM)']} patients"),
        (0.5, 0.74, f"Primary tumor (sample type 01)\n{n['Primary tumor (sample type 01)']} aliquots, "
                    f"{p['Primary tumor (sample type 01)']} patients"),
        (0.5, 0.56, f"One aliquot per patient\nn = {n['One aliquot per patient']} "
                    f"(BI {n['  of which BI batch']}, GDC2024 {n['  of which GDC2024 batch']}, "
                    f"other {n['  of which other batch']})"),
        (0.28, 0.36, f"Primary cohort: IDH-wildtype\nn = {n['IDH-wildtype (primary cohort, screen)']}"),
        (0.75, 0.36, f"Sensitivity cohort: all primary tumors\nn = {n['One aliquot per patient']}"),
        (0.28, 0.14, f"Survival evaluable\nn = {n['Survival evaluable, IDH-wildtype (primary)']}"),
        (0.75, 0.14, f"Survival evaluable\nn = {n['Survival evaluable, all primary (sensitivity)']}"),
    ]
    side = [
        (0.86, 0.83, "Excluded: " + ", ".join(f"{k} {v}" for k, v in excl.items())),
        (0.86, 0.65, f"Excluded: additional aliquots\nn = {n['Primary tumor (sample type 01)'] - n['One aliquot per patient']}"),
        (0.07, 0.46, f"IDH-mutant {n['IDH-mutant (sensitivity only)']}\nIDH unknown {n['IDH unknown (sensitivity only)']}\n"
                     "(sensitivity cohort only)"),
        (0.07, 0.25, f"No evaluable OS\nn = {n['IDH-wildtype (primary cohort, screen)'] - n['Survival evaluable, IDH-wildtype (primary)']}"),
        (0.97, 0.25, f"No evaluable OS\nn = {n['One aliquot per patient'] - n['Survival evaluable, all primary (sensitivity)']}"),
    ]
    fig, ax = plt.subplots(figsize=(8.5, 7.5))
    ax.axis("off")
    for x, y, t in boxes:
        ax.text(x, y, t, ha="center", va="center", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black", lw=0.8))
    for x, y, t in side:
        ax.text(x, y, t, ha="center", va="center", fontsize=8, color="0.25",
                bbox=dict(boxstyle="round,pad=0.4", fc="0.95", ec="0.6", lw=0.6))
    arrows = [((0.5, 0.88), (0.5, 0.78)), ((0.5, 0.70), (0.5, 0.60)),
              ((0.5, 0.52), (0.28, 0.40)), ((0.5, 0.52), (0.75, 0.40)),
              ((0.28, 0.32), (0.28, 0.18)), ((0.75, 0.32), (0.75, 0.18))]
    for a, b in arrows:
        ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.set_xlim(-0.05, 1.08)
    ax.set_ylim(0.05, 1.0)
    fig.savefig(FIG / "FigS1_cohort_flow.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
