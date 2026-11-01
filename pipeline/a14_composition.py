#!/usr/bin/env python3
"""
a14_composition.py
Composition adjustment (PREREGISTRATION.md, Section 5.2; decision condition 2).

ESTIMATE score: the tidyestimate 1.1.1 algorithm (estimate_score, is_affymetrix = FALSE) on log2(TPM + 1) of the
genes whose symbol is in tidyestimate's common_genes (filter_common_genes, find_alias = FALSE). Per sample:
average ranks, scaled by 10000 / number of genes, ordered decreasing (ties in input order), ssGSEA with weight
0.25 for the stromal and the immune set; ESTIMATE = stromal + immune. Gene sets are read from the package data
(data/downloads/estimate/extracted/, exported from gene_sets.rda and common_genes.rda).

Neuronal signature: mean z-score (across primary-cohort samples) of log2(TPM + 1) of RBFOX3, SNAP25, SYT1, STMN2,
GABRA1, SYN1, NEFL, GAD1, SLC17A7, CAMK2A. A symbol with more than one Ensembl ID uses the ID with the highest
mean expression.

For each score: GPR153-MGMT partial Spearman adjusted for batch and the score (k = 2; the score enters as ranks),
Bonett-Wright 95% CI; and the genome-wide screen repeated with the same adjustment (eligible genes, MGMT
excluded, BH), giving GPR153's q and rank.

Inputs
  revision/work/tcga_expr.npz (a10), revision/manifests/tcga_samples.csv (r01),
  revision/results/screen/screen_gene_level.csv (a11), data/downloads/estimate/extracted/*.csv
Outputs (revision/results/composition/)
  composition_scores.csv, composition_adjusted.csv, composition_screen_<score>.csv, neuronal_genes_used.csv

Run from the repository root:
    python pipeline/a14_composition.py
"""

import numpy as np
import pandas as pd

import common
from common import GPR153, MGMT

EST = common.DL / "estimate/extracted"


def estimate_scores(V, symbols, gene_sets):
    """tidyestimate::estimate_score on matrix V (genes x samples) with row names `symbols`."""
    R = common.rank_rows(V.T).T * 10000 / V.shape[0]
    out = {}
    for name, genes in gene_sets.items():
        common_g = set(genes) & set(symbols)
        hit_rows = np.isin(symbols, list(common_g))
        order = np.argsort(-R, axis=0, kind="stable")
        H = hit_rows[order].astype(float)
        val = np.take_along_axis(R, order, 0) ** 0.25
        reward = np.cumsum(H * val, 0) / (H * val).sum(0)
        penalty = np.cumsum(1 - H, 0) / (1 - H).sum(0)
        out[name] = (reward - penalty).sum(0)
        out[f"n_{name}_genes"] = len(common_g)
    return out


def load_estimate_sets():
    gs = pd.read_csv(EST / "gene_sets.csv")
    cg = pd.read_csv(EST / "common_genes.csv")
    return {c: gs[c].dropna().tolist() for c in gs.columns}, set(cg["hgnc_symbol"].dropna())


def neuronal_score(t, L):
    used, z = [], []
    for s in common.NEURONAL:
        rows = np.where(t.symbol == s)[0]
        if not len(rows):
            used.append({"symbol": s, "gene_id": "", "used": False})
            continue
        r = rows[np.argmax(L[rows].mean(1))]
        used.append({"symbol": s, "gene_id": t.gene_id[r], "used": True})
        z.append(common.standardize(L[r]))
    return np.mean(z, axis=0), pd.DataFrame(used)


def main():
    common.set_seeds()
    common.require_registration()
    out = common.outdir("composition")
    t = common.TCGA()
    prim = common.cohort_aliquots(t.m, "in_screen_primary")
    L = t.log2tpm(prim["aliquot"].values)
    batch = (prim["batch"].values == "GDC2024").astype(float)
    sets, cg = load_estimate_sets()
    keep = np.isin(t.symbol, list(cg))
    est = estimate_scores(L[keep], t.symbol[keep], {"stromal": sets["stromal_signature"],
                                                     "immune": sets["immune_signature"]})
    neu, used = neuronal_score(t, L)
    used.to_csv(out / "neuronal_genes_used.csv", index=False)
    scores = pd.DataFrame({"aliquot": prim["aliquot"], "patient": prim["patient"], "batch": prim["batch"],
                           "stromal": est["stromal"], "immune": est["immune"],
                           "estimate": est["stromal"] + est["immune"], "neuronal": neu})
    scores.to_csv(out / "composition_scores.csv", index=False)

    gl = pd.read_csv(common.RES / "screen/screen_gene_level.csv")
    elig = np.isin(t.gene_id, gl.loc[gl["eligible_20"], "gene_id"]) & (t.gene_id != MGMT)
    gi, mi = t.row[GPR153], t.row[MGMT]
    rows = []
    for name in ("estimate", "neuronal"):
        cov = np.column_stack([batch, scores[name].values])
        tested = elig & (L.std(1) > 0)
        rho = np.full(len(L), np.nan)
        p = np.full(len(L), np.nan)
        rho[tested], p[tested], df = common.partial_spearman(L[tested], L[mi], cov)
        q = np.full(len(L), np.nan)
        q[tested] = common.bh(p[tested])
        rank = common.signed_rank(rho, tested)
        pd.DataFrame({"gene_id": t.gene_id[tested], "symbol": t.symbol[tested], "rho": rho[tested],
                      "p": p[tested], "q": q[tested], "rank": rank[tested]}).to_csv(
            out / f"composition_screen_{name}.csv", index=False)
        lo, hi, _, _ = common.bw_ci(rho[gi], len(prim), 2)
        rows.append({"adjustment": f"batch + {name} score", "n": len(prim), "k": 2, "df": df,
                     "gpr153_rho": rho[gi], "gpr153_p": p[gi], "ci_low": lo, "ci_high": hi,
                     "gpr153_q_adjusted_screen": q[gi], "gpr153_rank_adjusted_screen": rank[gi],
                     "genes_tested": int(tested.sum()),
                     "score_genes": est.get("n_stromal_genes", "") if name == "estimate" else int(used["used"].sum()),
                     "condition2_pass": bool(rho[gi] > 0 and q[gi] < 0.05)})
    pd.DataFrame(rows).to_csv(out / "composition_adjusted.csv", index=False)
    print(pd.DataFrame(rows)[["adjustment", "gpr153_rho", "gpr153_q_adjusted_screen", "condition2_pass"]])


if __name__ == "__main__":
    main()
