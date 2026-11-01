#!/usr/bin/env python3
"""
r07_survival_audit.py
Survival-data audit for every survival cohort (PREREGISTRATION.md, Section 6.6).
For each cohort: raw vital-status strings and their counts, the mapping to the
event indicator, the time column and unit, conversion to days, the handling of
missing time and time <= 0, and assertions that n and events recomputed here from
the raw source files match the pipeline outputs. Only clinical tables are read;
no expression value is read.

Inputs
  revision/manifests/tcga_samples.csv, external_cohort_samples.csv     (r01, r02)
  revision/results/cohort/survival_replication_candidates.csv          (r05)
  raw clinical files under data/raw and data/downloads
Outputs
  revision/results/cohort/survival_audit_strings.csv
  revision/results/cohort/survival_audit.csv

Run from the repository root:
    python pipeline/r07_survival_audit.py
"""

import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

DL = Path("data/downloads")
RAW = Path("data/raw")
MAN = Path("revision/manifests")
RES = Path("revision/results/cohort")
DAYS_PER_MONTH = 365.25 / 12


def zip_table(path):
    z = zipfile.ZipFile(path)
    name = [x for x in z.namelist() if not x.startswith("__")][0]
    return pd.read_csv(io.StringIO(z.read(name).decode().replace("\r\n", "\n").replace("\r", "\n")), sep="\t")


def audit(cohort, raw_status, mapping, time_days, time_col, unit, selected):
    """raw_status, time_days: Series over the selected cohort (after the cohort's selection rules)."""
    status = raw_status.astype(object).where(raw_status.notna(), "<missing>").astype(str)
    event = status.map(mapping)
    ok_time = time_days.notna() & (time_days > 0)
    usable = event.notna() & ok_time
    strings = [{"cohort": cohort, "raw_value": v, "n": int(n), "event": mapping.get(v, "unmapped (excluded)")}
               for v, n in status.value_counts().items()]
    summary = {"cohort": cohort, "selection": selected, "time_column": time_col, "unit_in_source": unit,
               "conversion_to_days": "x 30.4375" if unit == "months" else "none",
               "n_selected": len(status), "n_status_unmapped": int(event.isna().sum()),
               "n_time_missing": int(time_days.isna().sum()), "n_time_le_0": int((time_days <= 0).sum()),
               "time_le_0_rule": "excluded from survival analyses", "n_usable": int(usable.sum()),
               "events": int(event[usable].sum()), "censored": int((event[usable] == 0).sum()),
               "max_time_days": float(time_days[usable].max())}
    return strings, summary, usable, event


def tcga():
    m = pd.read_csv(MAN / "tcga_samples.csv")
    p = m[m["in_screen_primary"]]
    cases = {c["submitter_id"]: c for c in json.load(open(DL / "gdc_meta/tcga_gbm_cases.json"))["hits"]}
    raw = p["patient"].map(lambda x: (cases[x].get("demographic") or {}).get("vital_status"))
    assert (raw.fillna("NA").values == p["vital_status"].fillna("NA").values).all(), "TCGA vital status differs from GDC"
    s, summ, usable, ev = audit("TCGA-GBM", raw, {"Dead": 1, "Alive": 0}, p["os_days"],
                                "days_to_death (Dead) / max days_to_last_follow_up, days_to_follow_up (Alive)",
                                "days", "primary cohort (IDH-wildtype, sample type 01, one aliquot per patient)")
    assert summ["n_usable"] == int(p["in_survival_primary"].sum())
    assert summ["events"] == int(p.loc[p["in_survival_primary"], "os_event"].sum())
    return s, summ


def cgga(n, r05):
    c = zip_table(DL / f"cgga/CGGA.mRNAseq_{n}_clinical.20200506.txt.zip")
    c = c[(c.PRS_type == "Primary") & (c.Grade == "WHO IV") & (c.IDH_mutation_status == "Wildtype")]
    raw = c["Censor (alive=0; dead=1)"]
    s, summ, usable, ev = audit(f"CGGA mRNAseq_{n}", raw.map(lambda v: v if pd.isna(v) else str(int(v))),
                                {"1": 1, "0": 0}, c["OS"].astype(float), "OS",
                                "days (inferred; the source gives no unit)", "primary, WHO IV, IDH-wildtype")
    man = pd.read_csv(MAN / "external_cohort_samples.csv")
    g = man[(man.cohort == f"CGGA_{n}") & man.included]
    assert set(g.sample_id) == set(c.CGGA_ID)
    assert summ["n_usable"] == int((g.os_days.notna() & (g.os_days > 0) & g.os_event.notna()).sum())
    assert summ["events"] == int(g.loc[g.os_days > 0, "os_event"].sum())
    assert summ["n_usable"] == int(r05.loc[f"CGGA mRNAseq_{n}", "n_usable"])
    return s, summ


def gse13041(r05):
    k = pd.read_csv(RAW / "geo/GSE13041/clinical.csv", low_memory=False)
    man = pd.read_csv(MAN / "external_cohort_samples.csv")
    inc = man[(man.cohort == "GSE13041_GPL96") & man.included]
    k = k[k.geo_accession.isin(inc.sample_id)]
    raw = k["characteristics_ch1.2.vital status"]
    s, summ, usable, ev = audit("GSE13041 GPL96", raw, {"DEAD": 1, "DECEASED": 1, "ALIVE": 0, "dead": 1, "alive": 0},
                                pd.to_numeric(k["characteristics_ch1.1.tts(days)"], errors="coerce"), "1.tts(days)",
                                "days", "GPL96, description = Primary GBM, repeat and duplicate arrays excluded")
    assert summ["n_usable"] == int(r05.loc["GSE13041 GPL96 primary", "n_usable"])
    assert summ["events"] == int(inc.loc[inc.os_days > 0, "os_event"].sum())
    return s, summ


def rembrandt(r05):
    c = pd.read_csv(DL / "geo_replication/GSE108474_REMBRANDT_clinical.data.txt.gz", sep="\t", dtype=str)
    mp = pd.read_csv(DL / "geo_replication/GSE108474_REMBRANDT_biospecimen_mapping_GEO.txt.gz", sep="\t", dtype=str)
    import gzip
    with gzip.open(DL / "geo_replication/GSE108474_REMBRANDT_GeneExpression.txt.gz", "rt") as f:
        cols = set(f.readline().rstrip("\n").split("\t"))
    mapped = set(mp.loc[mp["BIOSPECIMEN_ID"].str.lstrip("X").isin(cols), "SUBJECT_ID"])
    g = c[(c.DISEASE_TYPE == "GBM") & c.SUBJECT_ID.isin(mapped)]
    t = pd.to_numeric(g.OVERALL_SURVIVAL_MONTHS, errors="coerce") * DAYS_PER_MONTH
    s, summ, usable, ev = audit("Rembrandt (GSE108474)", g.EVENT_OS, {"1": 1, "0": 0}, t, "OVERALL_SURVIVAL_MONTHS",
                                "months", "DISEASE_TYPE = GBM, with a mapped expression column")
    assert summ["n_usable"] == int(r05.loc["Rembrandt (GSE108474)", "n_usable"])
    return s, summ


def cptac(r05):
    meta = json.load(open(DL / "cptac/expression/files_meta.json"))["hits"]
    pts = {h["cases"][0]["submitter_id"] for h in meta if h["cases"][0]["samples"][0]["sample_type"] == "Primary Tumor"}
    rows = []
    for c in json.load(open(DL / "cptac/cptac_brain_cases.json"))["hits"]:
        diag = [x.get("primary_diagnosis") for x in c.get("diagnoses") or []]
        if c["submitter_id"] not in pts or "Glioblastoma" not in diag:
            continue
        d = c.get("demographic") or {}
        fu = [x.get("days_to_last_follow_up") for x in c.get("diagnoses") or []] + \
             [x.get("days_to_follow_up") for x in c.get("follow_ups") or []]
        fu = [v for v in fu if v is not None]
        vs = d.get("vital_status")
        t = d.get("days_to_death") if vs == "Dead" else (max(fu) if vs == "Alive" and fu else None)
        rows.append({"vital_status": vs, "time": t})
    g = pd.DataFrame(rows)
    s, summ, usable, ev = audit("CPTAC-3 GBM", g.vital_status, {"Dead": 1, "Alive": 0},
                                pd.to_numeric(g.time, errors="coerce"),
                                "days_to_death (Dead) / max days_to_last_follow_up, days_to_follow_up (Alive)", "days",
                                "primary_diagnosis Glioblastoma, primary tumor with open STAR counts")
    assert summ["n_usable"] == int(r05.loc["CPTAC-3 GBM", "n_usable"]), (summ["n_usable"], r05.loc["CPTAC-3 GBM", "n_usable"])
    return s, summ


def main():
    r05 = pd.read_csv(RES / "survival_replication_candidates.csv").set_index("cohort")
    parts = [tcga(), cgga("693", r05), cgga("325", r05), gse13041(r05), rembrandt(r05), cptac(r05)]
    strings = pd.DataFrame([x for s, _ in parts for x in s])
    summary = pd.DataFrame([summ for _, summ in parts])
    strings.to_csv(RES / "survival_audit_strings.csv", index=False)
    summary.to_csv(RES / "survival_audit.csv", index=False)
    pd.set_option("display.width", 250)
    print(strings.to_string(index=False))
    print(summary.drop(columns=["selection", "time_column"]).to_string(index=False))
    print("all survival assertions passed")


if __name__ == "__main__":
    main()
