#!/usr/bin/env python3
"""
a18_pca.py
PCA of the primary cohort (PREREGISTRATION.md, Section 13, S2 Fig): the 2,000 eligible genes with the largest
variance of log2(TPM + 1), each standardized; PCA (scikit-learn); k-means with k = 3, n_init = 50,
random_state = 42, on the same standardized matrix. Clusters are labeled Cluster 1-3 in order of size.

Inputs
  revision/work/tcga_expr.npz (a10), revision/manifests/tcga_samples.csv (r01),
  revision/results/screen/screen_gene_level.csv (a11)
Outputs (revision/results/pca/)
  pca_scores.csv, pca_variance.csv, pca_genes.csv

Run from the repository root:
    python pipeline/a18_pca.py
"""

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

import common


def main():
    common.set_seeds()
    common.require_registration()
    out = common.outdir("pca")
    P = common.PARAMS["pca"]
    t = common.TCGA()
    prim = common.cohort_aliquots(t.m, "in_screen_primary")
    gl = pd.read_csv(common.RES / "screen/screen_gene_level.csv")
    rows = np.array([t.row[g] for g in gl.loc[gl["eligible_20"], "gene_id"]])
    L = t.log2tpm(prim["aliquot"].values)[rows]
    top = rows[np.argsort(-L.var(1, ddof=1), kind="stable")[:P["n_genes"]]]
    X = t.log2tpm(prim["aliquot"].values)[top].T
    X = (X - X.mean(0)) / X.std(0, ddof=1)
    pca = PCA(n_components=10, random_state=common.SEED).fit(X)
    sc = pca.transform(X)
    km = KMeans(n_clusters=P["k"], n_init=P["n_init"], random_state=common.SEED).fit(X)
    sizes = pd.Series(km.labels_).value_counts()
    relabel = {lab: f"Cluster {i + 1}" for i, lab in enumerate(sizes.index)}
    d = pd.DataFrame(sc, columns=[f"PC{i + 1}" for i in range(10)])
    d.insert(0, "aliquot", prim["aliquot"].values)
    d.insert(1, "batch", prim["batch"].values)
    d["cluster"] = [relabel[c] for c in km.labels_]
    d.to_csv(out / "pca_scores.csv", index=False)
    pd.DataFrame({"component": [f"PC{i + 1}" for i in range(10)],
                  "variance_explained": pca.explained_variance_ratio_, "n": len(prim),
                  "genes": len(top)}).to_csv(out / "pca_variance.csv", index=False)
    pd.DataFrame({"gene_id": t.gene_id[top], "symbol": t.symbol[top]}).to_csv(out / "pca_genes.csv", index=False)
    print(pd.read_csv(out / "pca_variance.csv").head(3).to_string(index=False))


if __name__ == "__main__":
    main()
