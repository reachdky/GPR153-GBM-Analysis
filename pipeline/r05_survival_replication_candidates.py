#!/usr/bin/env python3
"""
r05_survival_replication_candidates.py
Check each candidate survival-replication cohort against the usability criteria
(PREREGISTRATION.md, Section 8) and log the reason for inclusion or exclusion.

Criteria (all required): glioblastoma; survival time (> 0) and vital status present;
a GPR153 probe or gene on the platform (from the annotation only); open access.
Only annotations, clinical tables, row identifiers and file metadata are read;
no expression value is read.

Outputs
  revision/results/cohort/survival_replication_candidates.csv

Run from the repository root:
    python pipeline/r05_survival_replication_candidates.py
"""

import gzip
import json
import zipfile
from pathlib import Path

import pandas as pd

DL = Path("data/downloads")
RAW = Path("data/raw")
RES = Path("revision/results/cohort")
GPR153_ENTREZ = "387509"
GPR153_PROBES = {"221902_at", "64942_at"}


def row_ids(path, sep="\t"):
    """First field of every row (gene or probe identifier)."""
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        next(f)
        return {line.split(sep, 1)[0].strip('"') for line in f}


def platform_ids(path, col):
    """GEO platform table (view=data); the table starts at the line beginning with 'ID<tab>'."""
    with open(path) as f:
        skip = next(i for i, line in enumerate(f) if line.startswith("ID\t"))
    return pd.read_csv(path, sep="\t", skiprows=skip, dtype=str, usecols=lambda c: c in ("ID", col))


def tcga():
    m = pd.read_csv("revision/manifests/tcga_samples.csv")
    n = int(m["in_survival_primary"].sum())
    return {"cohort": "TCGA-GBM (discovery)", "source": "GDC STAR - Counts, GENCODE v36", "glioblastoma": "yes",
            "survival_fields": "yes", "gpr153_on_platform": "yes (ENSG00000158292)", "open_access": "yes",
            "n_usable": n, "usable": True, "reason": "discovery cohort; included in meta-analysis with and without it"}


def cgga():
    out = []
    for n in ["693", "325"]:
        z = zipfile.ZipFile(DL / f"cgga/CGGA.mRNAseq_{n}_clinical.20200506.txt.zip")
        txt = z.read([x for x in z.namelist() if not x.startswith("__")][0]).decode().replace("\r\n", "\n").replace("\r", "\n")
        c = pd.read_csv(pd.io.common.StringIO(txt), sep="\t")
        ok = c[(c.PRS_type == "Primary") & (c.Grade == "WHO IV") & (c.IDH_mutation_status == "Wildtype")]
        ok = ok[(ok["OS"] > 0) & ok["Censor (alive=0; dead=1)"].notna()]
        z = zipfile.ZipFile(DL / f"cgga/CGGA.mRNAseq_{n}.RSEM-genes.20200506.txt.zip")
        name = [x for x in z.namelist() if not x.startswith("__")][0]
        genes = {l.split(b"\t", 1)[0].decode() for l in z.read(name).replace(b"\r\n", b"\n").replace(b"\r", b"\n").split(b"\n")[1:] if l}
        out.append({"cohort": f"CGGA mRNAseq_{n}", "source": "CGGA release 20200506, RSEM FPKM",
                    "glioblastoma": "yes (primary, WHO IV, IDH-wildtype)", "survival_fields": "yes (OS, Censor)",
                    "gpr153_on_platform": "yes" if "GPR153" in genes else "no", "open_access": "yes",
                    "n_usable": len(ok), "usable": "GPR153" in genes, "reason": "meets all criteria"})
    return out


def gse13041():
    m = pd.read_csv("revision/manifests/external_cohort_samples.csv")
    g = m[(m.cohort == "GSE13041_GPL96") & m.included]
    n = int(((g.os_days > 0) & g.os_event.notna()).sum())
    return {"cohort": "GSE13041 GPL96 primary", "source": "GEO GSE13041, GPL96 MAS5",
            "glioblastoma": "yes (Primary GBM)", "survival_fields": "yes (tts, vital status)",
            "gpr153_on_platform": "yes (221902_at, 64942_at)", "open_access": "yes", "n_usable": n, "usable": True,
            "reason": "meets all criteria; no IDH data"}


def gse4290():
    return {"cohort": "GSE4290 GBM", "source": "GEO GSE4290, GPL570", "glioblastoma": "yes",
            "survival_fields": "no", "gpr153_on_platform": "yes (221902_at, 64942_at)", "open_access": "yes",
            "n_usable": 0, "usable": False, "reason": "excluded: no survival data deposited"}


def characteristic_keys(path):
    """Sample characteristic field names from a series matrix header (stops before the data table)."""
    keys = set()
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("!Sample_characteristics_ch1"):
                keys |= {v.strip('"').split(":")[0] for v in line.rstrip("\n").split("\t")[1:] if v.strip('"')}
            if line.startswith("!series_matrix_table_begin"):
                break
    surv = any(k for k in keys if any(w in k.lower() for w in ("surviv", "vital", "death", "follow")))
    return keys, surv


def gse16011():
    keys, surv = characteristic_keys(DL / "geo_replication/GSE16011_series_matrix.txt.gz")
    ann = platform_ids(DL / "platforms/GPL8542_table.txt", "ORF")
    on = GPR153_ENTREZ in set(ann["ORF"])
    return {"cohort": "GSE16011 (Gravendeel 2009)", "source": "GEO GSE16011, GPL8542 (Brainarray ENTREZG v10)",
            "glioblastoma": "yes (subset)", "survival_fields": "yes" if surv else f"no (characteristics: {', '.join(sorted(keys))})",
            "gpr153_on_platform": "yes" if on else "no (Entrez 387509 not in GPL8542)", "open_access": "yes",
            "n_usable": 0, "usable": surv and on,
            "reason": "excluded: " + "; ".join(x for x, bad in [("GPR153 not on platform", not on),
                                                               ("no survival data in GEO record", not surv)] if bad)}


def gse68848():
    keys, surv = characteristic_keys(DL / "geo_replication/GSE68848_series_matrix.txt.gz")
    return {"cohort": "GSE68848 (Rembrandt subset)", "source": "GEO GSE68848, GPL570",
            "glioblastoma": "yes (subset)", "survival_fields": "yes" if surv else f"no (characteristics: {', '.join(sorted(keys))})",
            "gpr153_on_platform": "yes (221902_at, 64942_at)", "open_access": "yes", "n_usable": 0, "usable": surv,
            "reason": "meets all criteria" if surv else "excluded: no survival data in GEO record"}


def rembrandt():
    c = pd.read_csv(DL / "geo_replication/GSE108474_REMBRANDT_clinical.data.txt.gz", sep="\t", dtype=str)
    mp = pd.read_csv(DL / "geo_replication/GSE108474_REMBRANDT_biospecimen_mapping_GEO.txt.gz", sep="\t", dtype=str)
    mp["col"] = mp["BIOSPECIMEN_ID"].str.lstrip("X")
    with gzip.open(DL / "geo_replication/GSE108474_REMBRANDT_GeneExpression.txt.gz", "rt") as f:
        cols = set(f.readline().rstrip("\n").split("\t"))
    rows = row_ids(DL / "geo_replication/GSE108474_REMBRANDT_GeneExpression.txt.gz")
    g = c[c.DISEASE_TYPE == "GBM"]
    g = g[(pd.to_numeric(g.OVERALL_SURVIVAL_MONTHS, errors="coerce") > 0) & g.EVENT_OS.notna()]
    g = g[g.SUBJECT_ID.isin(mp.loc[mp["col"].isin(cols), "SUBJECT_ID"])]
    probes = sorted(GPR153_PROBES & rows)
    return {"cohort": "Rembrandt (GSE108474)", "source": "GEO GSE108474, GPL570, dChip-normalized",
            "glioblastoma": "yes (DISEASE_TYPE = GBM)", "survival_fields": "yes (OVERALL_SURVIVAL_MONTHS, EVENT_OS)",
            "gpr153_on_platform": f"yes ({', '.join(probes)})" if probes else "no", "open_access": "yes",
            "n_usable": len(g), "usable": bool(probes) and len(g) > 0,
            "reason": "meets all criteria; no IDH data; primary vs recurrent not recorded"}


def cptac():
    meta = json.load(open(DL / "cptac/expression/files_meta.json"))["hits"]
    pts = {h["cases"][0]["submitter_id"] for h in meta if h["cases"][0]["samples"][0]["sample_type"] == "Primary Tumor"}
    cases = json.load(open(DL / "cptac/cptac_brain_cases.json"))["hits"]
    n = 0
    for c in cases:
        if c["submitter_id"] not in pts:
            continue
        d = c.get("demographic") or {}
        diag = [x.get("primary_diagnosis") for x in c.get("diagnoses") or []]
        fu = [x.get("days_to_last_follow_up") for x in c.get("diagnoses") or []] + \
             [x.get("days_to_follow_up") for x in c.get("follow_ups") or []]
        fu = [v for v in fu if v is not None]
        vs = d.get("vital_status")
        t = d.get("days_to_death") if vs == "Dead" else (max(fu) if vs == "Alive" and fu else None)
        n += t is not None and t > 0 and "Glioblastoma" in diag
    return {"cohort": "CPTAC-3 GBM", "source": "GDC CPTAC-3, STAR - Counts, GENCODE v36",
            "glioblastoma": "yes (GDC primary_diagnosis = Glioblastoma)", "survival_fields": "yes (GDC vital status, days)",
            "gpr153_on_platform": "yes (same GDC pipeline as TCGA)", "open_access": "yes",
            "n_usable": n, "usable": n > 0, "reason": "meets all criteria; IDH from open masked MAFs"}


def main():
    rows = [tcga(), *cgga(), gse13041(), gse4290(), gse16011(), gse68848(), rembrandt(), cptac()]
    out = pd.DataFrame(rows)
    out.to_csv(RES / "survival_replication_candidates.csv", index=False)
    pd.set_option("display.width", 250)
    print(out[["cohort", "survival_fields", "gpr153_on_platform", "n_usable", "usable", "reason"]].to_string(index=False))


if __name__ == "__main__":
    main()
