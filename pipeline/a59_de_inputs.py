#!/usr/bin/env python3
"""
a59_de_inputs.py
Count matrix for differential expression (PREREGISTRATION.md, Section 12). Locked until the analysis plan is
registered. UCSC Toil RSEM expected counts (GENCODE v23), stored as log2(x + 1), back-transformed as 2^x - 1 and
rounded to integers. Tumors: primary-cohort TCGA samples present in Toil, matched on the 15-character sample
barcode. Normals: GTEx Brain - Cortex and Brain - Frontal Cortex (BA9) present in the count matrix; primary set one
sample per donor (Brain - Cortex preferred, then lowest sample ID); sensitivity set all. Gene IDs version-stripped;
PAR_Y rows dropped.

Inputs
  data/downloads/toil/TcgaTargetGtex_gene_expected_count.gz, TcgaTargetGTEX_phenotype.txt.gz (r00)
  revision/manifests/tcga_samples.csv (r01)
Outputs (revision/work/)
  de_toil.npz (counts genes x samples, gene_id, sample), de_toil_samples.csv

Run from the repository root:
    python pipeline/a59_de_inputs.py
"""

import gzip

import numpy as np
import pandas as pd

import common

TISSUES = {"Brain - Cortex": 0, "Brain - Frontal Cortex (Ba9)": 1}


def main():
    common.require_registration()
    path = common.DL / "toil/TcgaTargetGtex_gene_expected_count.gz"
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    ph = pd.read_csv(common.DL / "toil/TcgaTargetGTEX_phenotype.txt.gz", sep="\t", encoding="latin1")
    prim = common.cohort_aliquots(common.tcga_manifest(), "in_screen_primary")
    tum = [s for s in prim["aliquot"].str[:15] if s in set(header)]
    nor = ph[ph["detailed_category"].isin(TISSUES) & ph["sample"].isin(header)].copy()
    nor["donor"] = nor["sample"].str.split("-").str[:2].str.join("-")
    nor["pref"] = nor["detailed_category"].map(TISSUES)
    first = nor.sort_values(["donor", "pref", "sample"]).groupby("donor").head(1)["sample"]
    samples = pd.concat([
        pd.DataFrame({"sample": tum, "group": "tumor", "tissue": "TCGA-GBM primary", "donor": [s[:12] for s in tum],
                      "normal_one_per_donor": False}),
        pd.DataFrame({"sample": nor["sample"], "group": "normal", "tissue": nor["detailed_category"],
                      "donor": nor["donor"], "normal_one_per_donor": nor["sample"].isin(first)})])
    e = pd.read_csv(path, sep="\t", usecols=["sample", *samples["sample"]], index_col=0)
    e = e[~e.index.str.endswith("_PAR_Y")]
    gid = e.index.str.replace(r"\.\d+$", "", regex=True)
    counts = np.rint(2 ** e[samples["sample"]].values - 1).astype(np.int64)
    common.save_matrix("de_toil.npz", counts=counts, gene_id=np.asarray(gid, str), sample=samples["sample"].to_numpy(str))
    samples.to_csv(common.work_path("de_toil_samples.csv"), index=False)
    print(f"de_toil.npz: {counts.shape}; tumors {len(tum)}, normals {len(nor)} ({len(first)} donors)")


if __name__ == "__main__":
    main()
