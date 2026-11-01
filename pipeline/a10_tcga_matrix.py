#!/usr/bin/env python3
"""
a10_tcga_matrix.py
Build the TCGA-GBM expression matrices (PREREGISTRATION.md, Section 4) from the 391 GDC STAR - Counts files:
tpm_unstranded, fpkm_unstranded and unstranded counts, genes keyed by version-stripped Ensembl ID (PAR_Y
duplicates dropped). Every file is md5-checked against the GDC md5 in the manifest. Locked until the analysis
plan is registered (common.require_registration).

Inputs
  revision/manifests/tcga_samples.csv                                              (r01)
  data/raw/tcga/expression/TCGA-GBM/Transcriptome_Profiling/Gene_Expression_Quantification/<file_id>/<file>
Outputs
  revision/work/tcga_expr.npz   tpm, fpkm, counts (genes x aliquots), gene_id, symbol, gene_type, aliquot
  revision/results/expression/tcga_gene_id_handling.csv

Run from the repository root:
    python pipeline/a10_tcga_matrix.py
"""

import hashlib

import numpy as np
import pandas as pd

import common

STAR_DIR = common.RAW / "tcga/expression/TCGA-GBM/Transcriptome_Profiling/Gene_Expression_Quantification"


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_star(path):
    d = pd.read_csv(path, sep="\t", comment="#", usecols=["gene_id", "gene_name", "gene_type", "unstranded",
                                                          "tpm_unstranded", "fpkm_unstranded"])
    return d[d["gene_id"].str.startswith("ENSG")]


def main():
    common.require_registration()
    m = common.tcga_manifest()
    genes, tpm, fpkm, cnt = None, [], [], []
    for r in m.itertuples():
        path = STAR_DIR / r.file_id / r.file_name
        if md5(path) != r.md5_gdc:
            raise SystemExit(f"md5 mismatch: {path}")
        d = read_star(path)
        if genes is None:
            genes = d[["gene_id", "gene_name", "gene_type"]].reset_index(drop=True)
        elif not (d["gene_id"].values == genes["gene_id"].values).all():
            raise SystemExit(f"gene order differs: {path}")
        tpm.append(d["tpm_unstranded"].to_numpy(float))
        fpkm.append(d["fpkm_unstranded"].to_numpy(float))
        cnt.append(d["unstranded"].to_numpy(np.int64))
    par = genes["gene_id"].str.endswith("_PAR_Y").to_numpy()
    stripped = genes["gene_id"].str.replace(r"\.\d+(_PAR_Y)?$", "", regex=True)
    keep = ~par
    assert stripped[keep].is_unique
    tpm, fpkm, cnt = (np.column_stack(a)[keep] for a in (tpm, fpkm, cnt))
    missing = int(np.isnan(tpm).any(1).sum())
    sym = genes.loc[keep, "gene_name"].fillna("")
    common.save_matrix("tcga_expr.npz", tpm=tpm, fpkm=fpkm, counts=cnt,
                       gene_id=stripped[keep].to_numpy(str), symbol=sym.to_numpy(str),
                       gene_type=genes.loc[keep, "gene_type"].fillna("").to_numpy(str),
                       aliquot=m["aliquot"].to_numpy(str))
    collisions = sym[sym.duplicated(keep=False) & (sym != "")]
    pd.DataFrame([
        {"step": "GENCODE v36 gene rows in STAR files", "n": len(genes)},
        {"step": "PAR_Y rows dropped", "n": int(par.sum())},
        {"step": "genes after version stripping", "n": int(keep.sum())},
        {"step": "genes with a missing value in any aliquot", "n": missing},
        {"step": "symbols shared by more than one Ensembl ID", "n": int(collisions.nunique())},
        {"step": "Ensembl IDs with a shared symbol", "n": len(collisions)},
        {"step": "aliquots", "n": tpm.shape[1]},
    ]).to_csv(common.outdir("expression") / "tcga_gene_id_handling.csv", index=False)
    print(f"tcga_expr.npz: {tpm.shape[0]} genes x {tpm.shape[1]} aliquots")


if __name__ == "__main__":
    main()
