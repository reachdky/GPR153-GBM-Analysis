#!/usr/bin/env python3
"""
a29_external_inputs.py
External expression matrices and clinical tables (PREREGISTRATION.md, Sections 7 and 8). Locked until the
analysis plan is registered. Whole matrices are kept (not only GPR153 and MGMT) for the array-distribution
diagnostics in Section 7.3. Values are stored as deposited; transformations are applied in a30 and a31.

  GSE13041   all 191 GPL96 arrays and all 27 GPL570 arrays (data/raw/geo/GSE13041/expression.csv), each with its
             own platform's probes (GPL8300 has no GPR153 probe set and is not loaded)
  GSE4290    all 180 arrays, with the manifest's inclusion flag (77 glioblastomas included; the others enter only
             the all-histology sensitivity analysis)
  CGGA       mRNAseq_693 and mRNAseq_325 RSEM FPKM, included samples (primary, WHO IV, IDH-wildtype)
  Rembrandt  GSE108474, DISEASE_TYPE = GBM subjects with a mapped expression column; clinical: OS months -> days
             (x 30.4375), EVENT_OS, AGE_RANGE, GENDER
  CPTAC-3    primary-tumor STAR - Counts (tpm_unstranded; one file per patient), md5-checked; primary_diagnosis
             Glioblastoma; survival as in r07; IDH from open masked MAFs (IDH1 R132 or IDH2 R172 in any tumor
             sample = mutant; MAF without either = wildtype; no MAF = unknown)

Inputs
  revision/manifests/external_cohort_samples.csv (r02); data/raw/geo/*; data/downloads/{cgga,geo_replication,cptac}
Outputs (revision/work/)
  ext_<cohort>.npz (X features x samples, feature, sample) and ext_<cohort>_clinical.csv for
  GSE13041 (GPL96), GSE13041_GPL570, GSE4290, CGGA_693, CGGA_325, Rembrandt, CPTAC

Run from the repository root:
    python pipeline/a29_external_inputs.py
"""

import gzip
import hashlib
import io
import json
import zipfile

import numpy as np
import pandas as pd

import common


def zip_table(path):
    z = zipfile.ZipFile(path)
    name = [x for x in z.namelist() if not x.startswith("__")][0]
    return pd.read_csv(io.StringIO(z.read(name).decode().replace("\r\n", "\n").replace("\r", "\n")), sep="\t")


def save(name, X, feature, sample, clinical):
    common.save_matrix(f"ext_{name}.npz", X=np.asarray(X, float), feature=np.asarray(feature, str),
                       sample=np.asarray(sample, str))
    clinical.to_csv(common.work_path(f"ext_{name}_clinical.csv"), index=False)
    print(f"ext_{name}: {np.shape(X)}")


def geo(man):
    for plat, name in (("GPL96", "GSE13041"), ("GPL570", "GSE13041_GPL570")):
        g = man[(man.cohort == "GSE13041_GPL96") & (man.platform == plat)]
        e = pd.read_csv(common.RAW / "geo/GSE13041/expression.csv", index_col=0, usecols=["ID_REF", *g.sample_id])
        e = e.dropna(how="all")
        save(name, e.values, e.index, e.columns, g)
    g = man[man.cohort == "GSE4290_GBM"]
    e = pd.read_csv(common.RAW / "geo/GSE4290/expression.csv", index_col=0, usecols=["ID_REF", *g.sample_id])
    save("GSE4290", e.values, e.index, e.columns, g)


def cgga(man):
    for n in ("693", "325"):
        g = man[(man.cohort == f"CGGA_{n}") & man.included].copy()
        g["sex"] = g["sex"].str.strip()
        e = zip_table(common.DL / f"cgga/CGGA.mRNAseq_{n}.RSEM-genes.20200506.txt.zip").set_index("Gene_Name")
        e = e[g.sample_id]
        save(f"CGGA_{n}", e.values, e.index, e.columns, g)


def rembrandt():
    c = pd.read_csv(common.DL / "geo_replication/GSE108474_REMBRANDT_clinical.data.txt.gz", sep="\t", dtype=str)
    mp = pd.read_csv(common.DL / "geo_replication/GSE108474_REMBRANDT_biospecimen_mapping_GEO.txt.gz", sep="\t", dtype=str)
    mp["col"] = mp["BIOSPECIMEN_ID"].str.lstrip("X")
    path = common.DL / "geo_replication/GSE108474_REMBRANDT_GeneExpression.txt.gz"
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    g = c[c.DISEASE_TYPE == "GBM"].merge(mp[mp["col"].isin(header)], on="SUBJECT_ID")
    g = g.drop_duplicates("SUBJECT_ID")
    e = pd.read_csv(path, sep="\t", header=None, skiprows=1, index_col=0)
    e.columns = header if len(header) == e.shape[1] else header[1:]
    e = e[g["col"]]
    ages = g["AGE_RANGE"].str.extract(r"(\d+)\s*-\s*(\d+)").astype(float)
    clin = pd.DataFrame({"sample_id": g["col"].values, "subject": g["SUBJECT_ID"].values,
                         "os_days": pd.to_numeric(g["OVERALL_SURVIVAL_MONTHS"], errors="coerce").values * common.DAYS_PER_MONTH,
                         "os_event": pd.to_numeric(g["EVENT_OS"], errors="coerce").values,
                         "age_range": g["AGE_RANGE"].values, "age": ages.mean(1).values,
                         "sex": g["GENDER"].str.strip().str.capitalize().replace({"": np.nan}).values})
    save("Rembrandt", e.values, e.index, e.columns, clin)


def cptac_idh():
    meta = json.load(open(common.DL / "cptac/maf/files_meta.json"))["hits"]
    status = {}
    for h in meta:
        if not any(s["sample_type"] == "Primary Tumor" for s in h["cases"][0]["samples"]):
            continue
        pt = h["cases"][0]["submitter_id"]
        mut = status.setdefault(pt, [])
        with gzip.open(common.DL / "cptac/maf" / h["file_id"] / h["file_name"], "rt") as f:
            header = None
            for line in f:
                if line.startswith("#"):
                    continue
                x = line.rstrip("\n").split("\t")
                if header is None:
                    header = {k: i for i, k in enumerate(x)}
                    continue
                gene, prot = x[header["Hugo_Symbol"]], x[header["HGVSp_Short"]]
                if (gene == "IDH1" and prot.startswith("p.R132")) or (gene == "IDH2" and prot.startswith("p.R172")):
                    mut.append(f"{gene} {prot}")
    return {k: ("mutant" if v else "wildtype") for k, v in status.items()}


def cptac():
    meta = json.load(open(common.DL / "cptac/expression/files_meta.json"))["hits"]
    files = {}
    for h in sorted(meta, key=lambda h: (h["cases"][0]["samples"][0]["submitter_id"], h["file_id"])):
        if h["cases"][0]["samples"][0]["sample_type"] == "Primary Tumor":
            files.setdefault(h["cases"][0]["submitter_id"], h)
    idh = cptac_idh()
    rows, cols, genes = [], [], None
    for c in json.load(open(common.DL / "cptac/cptac_brain_cases.json"))["hits"]:
        pt = c["submitter_id"]
        diag = c.get("diagnoses") or []
        if pt not in files or "Glioblastoma" not in [x.get("primary_diagnosis") for x in diag]:
            continue
        h = files[pt]
        path = common.DL / "cptac/expression" / h["file_id"] / h["file_name"]
        if hashlib.md5(path.read_bytes()).hexdigest() != h["md5sum"]:
            raise SystemExit(f"md5 mismatch: {path}")
        d = pd.read_csv(path, sep="\t", comment="#", usecols=["gene_id", "tpm_unstranded"])
        d = d[d.gene_id.str.startswith("ENSG") & ~d.gene_id.str.endswith("_PAR_Y")]
        if genes is None:
            genes = d.gene_id.str.replace(r"\.\d+$", "", regex=True).values
        cols.append(d.tpm_unstranded.values)
        dem = c.get("demographic") or {}
        fu = [x.get("days_to_last_follow_up") for x in diag] + [x.get("days_to_follow_up") for x in c.get("follow_ups") or []]
        fu = [v for v in fu if v is not None]
        vs = dem.get("vital_status")
        tt = dem.get("days_to_death") if vs == "Dead" else (max(fu) if vs == "Alive" and fu else None)
        age_d = next((x.get("age_at_diagnosis") for x in diag if x.get("age_at_diagnosis") is not None), None)
        rows.append({"sample_id": pt, "file_id": h["file_id"], "os_days": tt,
                     "os_event": {"Dead": 1, "Alive": 0}.get(vs), "vital_status": vs,
                     "age": age_d / 365.25 if age_d is not None else dem.get("age_at_index"),
                     "sex": (dem.get("sex_at_birth") or dem.get("gender") or "").capitalize() or None, "idh_status": idh.get(pt, "unknown")})
    clin = pd.DataFrame(rows)
    save("CPTAC", np.column_stack(cols), genes, clin.sample_id, clin)


def main():
    common.require_registration()
    man = common.external_manifest()
    geo(man)
    cgga(man)
    rembrandt()
    cptac()


if __name__ == "__main__":
    main()
