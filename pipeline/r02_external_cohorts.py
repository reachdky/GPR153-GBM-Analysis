#!/usr/bin/env python3
"""
r02_external_cohorts.py
Build the external-cohort sample manifest (PREREGISTRATION.md, Section 7).

Selection uses sample annotations only. A duplicate-sample check correlates
whole expression profiles between samples with GPR153 and MGMT removed; no
GPR153 or MGMT value is read.

Inputs
  data/raw/geo/GSE13041/clinical.csv, data/raw/geo/GSE4290/clinical.csv
  data/raw/geo/GSE13041/expression.csv, data/raw/geo/GSE4290/expression.csv
  data/raw/geo/GSE13041_family.soft.gz (GPL96 annotation)
  data/downloads/cgga/CGGA.mRNAseq_{693,325}*.zip                         (r00)
Outputs
  revision/manifests/external_cohort_samples.csv
  revision/results/cohort/external_cohort_flow.csv
  revision/results/cohort/external_duplicate_check.csv

Run from the repository root:
    python pipeline/r02_external_cohorts.py
"""

import gzip
import io
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

RAW = Path("data/raw")
DL = Path("data/downloads")
MAN = Path("revision/manifests")
RES = Path("revision/results/cohort")
EXCLUDE_PROBES = {"221902_at", "64942_at", "204880_at", "2051_at", "2052_g_at"}  # GPR153, MGMT
EXCLUDE_GENES = {"GPR153", "MGMT"}
DUP_R = 0.98


def read_zip_table(path):
    z = zipfile.ZipFile(path)
    name = [x for x in z.namelist() if not x.startswith("__")][0]
    txt = z.read(name).decode().replace("\r\n", "\n").replace("\r", "\n")
    return pd.read_csv(io.StringIO(txt), sep="\t")


def duplicate_pairs(expr, cohort):
    """Pairs of samples with Spearman r >= DUP_R across all retained features."""
    x = expr.dropna().rank(axis=0)
    c = np.corrcoef(x.values.T)
    iu = np.triu_indices_from(c, 1)
    hit = c[iu] >= DUP_R
    off = c[iu]
    return pd.DataFrame({"cohort": cohort, "sample_a": x.columns[iu[0][hit]], "sample_b": x.columns[iu[1][hit]],
                         "spearman_r": c[iu][hit]}), float(np.median(off)), float(off.max())


def gpl96_check():
    rows = {}
    with gzip.open(RAW / "geo/GSE13041_family.soft.gz", "rt", errors="replace") as f:
        on = False
        for line in f:
            if line.startswith("^PLATFORM"):
                on = line.strip().endswith("GPL96")
            elif on and line.split("\t")[0] in EXCLUDE_PROBES:
                p = line.rstrip("\n").split("\t")
                rows[p[0]] = p[10] if len(p) > 10 else ""
    return rows


def geo_cohorts():
    out, dups, stats = [], [], []

    # GSE13041, GPL96 only
    k = pd.read_csv(RAW / "geo/GSE13041/clinical.csv", low_memory=False)
    k.columns = [c.replace("characteristics_ch1.", "") for c in k.columns]
    k["tumor_id"] = k["title"].str.replace(r"^Tumor_fresh-frozen_", "", regex=True)
    prior_cols = [c for c in k.columns if "prior to tumor resection" in c]
    k["prior_treatment"] = k[prior_cols].notna().any(axis=1)
    g = k[(k["platform_id"] == "GPL96") & (k["description"] == "Primary GBM")].sort_values("geo_accession")
    first = set(g.groupby("tumor_id").head(1)["geo_accession"])
    for _, r in k.iterrows():
        if r["platform_id"] != "GPL96":
            st, why = False, f"excluded: platform {r['platform_id']}"
        elif r["description"] != "Primary GBM":
            st, why = False, f"excluded: {r['description']}"
        elif r["geo_accession"] not in first:
            st, why = False, "excluded: repeat array of same tumor ID (later GSM)"
        else:
            st, why = True, "included"
        out.append({"cohort": "GSE13041_GPL96", "sample_id": r["geo_accession"], "title": r["title"],
                    "platform": r["platform_id"], "histology": r["description"],
                    "prs_type": r["description"], "grade": "", "idh_status": "not available (no IDH data deposited)",
                    "prior_treatment": bool(r["prior_treatment"]),
                    "os_days": pd.to_numeric(r.get("1.tts(days)"), errors="coerce"),
                    "os_event": {"DEAD": 1, "DECEASED": 1, "ALIVE": 0}.get(str(r.get("2.vital status")).upper(), np.nan),
                    "age": pd.to_numeric(r.get("3.age(years)"), errors="coerce"),
                    "subtype_deposited": r.get("4.hc"), "included": st, "reason": why,
                    "sex": {"M": "Male", "F": "Female"}.get(str(r.get("6.gender")).strip(), np.nan)})
    e = pd.read_csv(RAW / "geo/GSE13041/expression.csv", index_col=0)
    e = e.loc[~e.index.isin(EXCLUDE_PROBES), sorted(first)]
    d, med, mx = duplicate_pairs(e, "GSE13041_GPL96")
    dups.append(d)
    stats.append(("GSE13041_GPL96", e.shape[1], med, mx))

    # GSE4290, glioblastoma only
    k = pd.read_csv(RAW / "geo/GSE4290/clinical.csv")
    h = k["characteristics_ch1.0.Histopathological diagnostic"]
    for (_, r), hist in zip(k.iterrows(), h):
        st = hist == "glioblastoma, grade 4"
        out.append({"cohort": "GSE4290_GBM", "sample_id": r["geo_accession"], "title": r["title"],
                    "platform": r["platform_id"], "histology": hist, "prs_type": "", "grade": "",
                    "idh_status": "not available (no IDH data deposited)", "prior_treatment": np.nan, "os_days": np.nan,
                    "os_event": np.nan, "age": np.nan, "subtype_deposited": "",
                    "included": st, "reason": "included" if st else f"excluded: histology {hist}"})
    inc = [s for s, hist in zip(k["geo_accession"], h) if hist == "glioblastoma, grade 4"]
    e = pd.read_csv(RAW / "geo/GSE4290/expression.csv", index_col=0)
    e = e.loc[~e.index.isin(EXCLUDE_PROBES), inc]
    d, med, mx = duplicate_pairs(np.log2(e.clip(lower=1)), "GSE4290_GBM")
    dups.append(d)
    stats.append(("GSE4290_GBM", e.shape[1], med, mx))
    return out, dups, stats


def cgga_cohorts():
    out, dups, stats, exprs = [], [], [], {}
    for n in ["693", "325"]:
        c = read_zip_table(DL / f"cgga/CGGA.mRNAseq_{n}_clinical.20200506.txt.zip")
        for _, r in c.iterrows():
            if r["PRS_type"] != "Primary":
                st, why = False, f"excluded: {r['PRS_type']} tumor"
            elif r["Grade"] != "WHO IV":
                st, why = False, f"excluded: grade {r['Grade']}"
            elif r["IDH_mutation_status"] != "Wildtype":
                st, why = False, f"excluded: IDH {r['IDH_mutation_status'] if isinstance(r['IDH_mutation_status'], str) else 'unknown'}"
            else:
                st, why = True, "included"
            out.append({"cohort": f"CGGA_{n}", "sample_id": r["CGGA_ID"], "title": "", "platform": "RNA-seq (RSEM FPKM)",
                        "histology": r["Histology"], "prs_type": r["PRS_type"], "grade": r["Grade"],
                        "idh_status": r["IDH_mutation_status"], "prior_treatment": np.nan,
                        "os_days": r["OS"], "os_event": r["Censor (alive=0; dead=1)"], "age": r["Age"],
                        "sex": r["Gender"], "mgmt_promoter": r["MGMTp_methylation_status"],
                        "radiotherapy": r["Radio_status (treated=1;un-treated=0)"],
                        "tmz": r["Chemo_status (TMZ treated=1;un-treated=0)"],
                        "subtype_deposited": "", "included": st, "reason": why})
        e = read_zip_table(DL / f"cgga/CGGA.mRNAseq_{n}.RSEM-genes.20200506.txt.zip").set_index("Gene_Name")
        inc = [s for s, ok in zip(c["CGGA_ID"], [o["included"] for o in out[-len(c):]]) if ok]
        exprs[n] = np.log2(e.loc[~e.index.isin(EXCLUDE_GENES), inc] + 1)
        d, med, mx = duplicate_pairs(exprs[n], f"CGGA_{n}")
        dups.append(d)
        stats.append((f"CGGA_{n}", len(inc), med, mx))
    # cross-release check: same patient sequenced in both releases
    common = exprs["693"].index.intersection(exprs["325"].index)
    a, b = exprs["693"].loc[common].rank(), exprs["325"].loc[common].rank()
    c = np.corrcoef(np.hstack([a.values, b.values]).T)[: a.shape[1], a.shape[1]:]
    i, j = np.where(c >= DUP_R)
    dups.append(pd.DataFrame({"cohort": "CGGA_693_vs_325", "sample_a": a.columns[i], "sample_b": b.columns[j],
                              "spearman_r": c[i, j]}))
    stats.append(("CGGA_693_vs_325", c.size, float(np.median(c)), float(c.max())))
    return out, dups, stats


def main():
    MAN.mkdir(parents=True, exist_ok=True)
    RES.mkdir(parents=True, exist_ok=True)
    probes = gpl96_check()
    print("GPL96 annotation:", probes)

    g_out, g_d, g_s = geo_cohorts()
    c_out, c_d, c_s = cgga_cohorts()
    m = pd.DataFrame(g_out + c_out)
    dup = pd.concat(g_d + c_d, ignore_index=True)
    dup.to_csv(RES / "external_duplicate_check.csv", index=False)
    s = pd.DataFrame(g_s + c_s, columns=["cohort", "n_samples_or_pairs", "median_r", "max_r"])

    # duplicates within a cohort: keep the first sample ID, flag the other
    for _, r in dup[dup["cohort"] != "CGGA_693_vs_325"].iterrows():
        drop = max(r["sample_a"], r["sample_b"])
        sel = (m["cohort"] == r["cohort"]) & (m["sample_id"] == drop)
        m.loc[sel, ["included", "reason"]] = [False, f"excluded: duplicate profile of {min(r['sample_a'], r['sample_b'])}"]
    for _, r in dup[dup["cohort"] == "CGGA_693_vs_325"].iterrows():
        sel = (m["cohort"] == "CGGA_325") & (m["sample_id"] == r["sample_b"])
        m.loc[sel, ["included", "reason"]] = [False, f"excluded: duplicate profile of CGGA_693 {r['sample_a']}"]

    m.to_csv(MAN / "external_cohort_samples.csv", index=False)
    flow = (m.groupby(["cohort", "reason"]).size().rename("n").reset_index()
            .sort_values(["cohort", "reason"]))
    flow.to_csv(RES / "external_cohort_flow.csv", index=False)
    print(flow.to_string(index=False))
    print(s.to_string(index=False))
    print("duplicate pairs:", len(dup))
    if len(dup):
        print(dup.to_string(index=False))


if __name__ == "__main__":
    main()
