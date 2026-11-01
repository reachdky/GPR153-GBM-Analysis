#!/usr/bin/env python3
"""
a40_enrichment.py
Pathway enrichment (PREREGISTRATION.md, Section 10): over-representation (gseapy enrich, local gene-set files)
and GSEA preranked (gseapy prerank); writes the inputs for CAMERA (a41_camera.R).

Ranking statistic: batch-adjusted partial Spearman of every eligible gene with GPR153 in the primary cohort
(GPR153 and MGMT excluded), BH. Gene sets use symbols: where a symbol has more than one eligible Ensembl ID, the ID
with the highest mean log2(TPM + 1) is kept; genes without a symbol are dropped.
ORA input lists: up to 200 genes with the largest positive rho among q < 0.05, and up to 200 with the most negative
rho among q < 0.05. Background: the eligible genes, restricted to the genes in the library; gene sets of 10-500
genes after that restriction; BH within each library. Libraries: Reactome_2022, KEGG_2021_Human,
GO_Biological_Process_2023. The hypergeometric p is also computed directly and reported beside gseapy's.
GSEA preranked: all eligible genes ranked by rho; MSigDB 2026.2.Hs Hallmark and Reactome; sizes 15-500;
1,000 gene-label permutations, seed 42; FDR per library. Targeted sets reported regardless of outcome:
HALLMARK_WNT_BETA_CATENIN_SIGNALING, HALLMARK_DNA_REPAIR.
Optional (last Tier 2 item, first to drop): the same ORA on up to 200 genes with the largest positive
batch-adjusted rho with MGMT among q < 0.05 (primary screen, a11; MGMT and GPR153 excluded), and GSEA preranked on
all eligible genes ranked by that rho; same libraries, sizes, background and seed.

Inputs
  revision/work/tcga_expr.npz (a10), revision/manifests/tcga_samples.csv (r01),
  revision/results/screen/screen_gene_level.csv (a11), data/downloads/genesets/*.gmt
Outputs (revision/results/enrichment/)
  gpr153_correlation.csv, ora_lists.csv, ora_<library>_<direction>.csv, gsea_<library>.csv, targeted_sets.csv,
  mgmt_ora_list.csv, mgmt_ora_<library>_positive.csv, mgmt_gsea_<library>.csv (optional analysis)
  revision/work/camera_expression.csv.gz, camera_covariates.csv (for a41)

Run from the repository root:
    python pipeline/a40_enrichment.py
"""

import gseapy
import numpy as np
import pandas as pd
from scipy import stats

import common
from common import GPR153, MGMT

GS = common.DL / "genesets"


def read_gmt(path):
    sets = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        sets[f[0]] = sorted({g.split(",")[0] for g in f[2:] if g})
    return sets


def ora(genes, universe, sets, P):
    """gseapy enrich on the restricted library; BH over every set tested (gseapy omits sets with no overlap,
    which enter with p = 1)."""
    u = set(universe)
    lib = {k: sorted(set(v) & u) for k, v in sets.items()}
    lib = {k: v for k, v in lib.items() if P["min_size"] <= len(v) <= P["max_size"]}
    bg = sorted(set().union(*lib.values()))
    q = sorted(set(genes) & set(bg))
    g = gseapy.enrich(gene_list=q, gene_sets=lib, background=bg, outdir=None, cutoff=1.0, no_plot=True).results
    N, n = len(bg), len(q)
    r = pd.DataFrame({"Term": list(lib), "set_size": [len(v) for v in lib.values()],
                      "overlap_n": [len(set(v) & set(q)) for v in lib.values()]})
    r["p_hypergeom_direct"] = [stats.hypergeom.sf(k - 1, N, m, n) for k, m in zip(r.overlap_n, r.set_size)]
    r = r.merge(g.drop(columns=["Gene_set"], errors="ignore"), on="Term", how="left")
    r["p"] = r["P-value"].fillna(1.0)
    r["q_bh"] = common.bh(r["p"].values)
    r["background_size"], r["list_size"], r["sets_tested"] = N, n, len(lib)
    return r.sort_values("p")


def mgmt_enrichment(gl, s, out):
    """Optional analysis: enrichment of the MGMT-ranked genes from the primary screen, same pipeline as above."""
    m = gl[gl["eligible_20"] & ~gl["gene_id"].isin([GPR153, MGMT])][["gene_id", "rho_adj", "q_adj"]]
    m = s[["gene_id", "symbol"]].merge(m, on="gene_id")
    P, G = common.PARAMS["ora"], common.PARAMS["gsea"]
    up = m[(m.q_adj < P["q"]) & (m.rho_adj > 0)].nlargest(common.PARAMS["mgmt_enrichment_max_genes"], "rho_adj")
    up.to_csv(out / "mgmt_ora_list.csv", index=False)
    for lib in P["libraries"]:
        if len(up):
            ora(up["symbol"], s["symbol"], read_gmt(GS / f"{lib}.gmt"), P).to_csv(out / f"mgmt_ora_{lib}_positive.csv", index=False)
    rnk = m.set_index("symbol")["rho_adj"].sort_values(ascending=False)
    for lib in G["libraries"]:
        gseapy.prerank(rnk=rnk, gene_sets=read_gmt(GS / f"{lib}.gmt"), min_size=G["min_size"], max_size=G["max_size"],
                       permutation_num=G["permutations"], seed=common.SEED, threads=4, outdir=None, no_plot=True,
                       verbose=False).res2d.to_csv(out / f"mgmt_gsea_{lib}.csv", index=False)
    print(f"MGMT-ranked ORA list: {len(up)} genes")


def main():
    common.set_seeds()
    common.require_registration()
    out = common.outdir("enrichment")
    t = common.TCGA()
    prim = common.cohort_aliquots(t.m, "in_screen_primary")
    gl = pd.read_csv(common.RES / "screen/screen_gene_level.csv")
    genes = gl.loc[gl["eligible_20"] & ~gl["gene_id"].isin([GPR153, MGMT]), "gene_id"].to_numpy()
    L = t.log2tpm(prim["aliquot"].values)
    rows = np.array([t.row[g] for g in genes])
    batch = (prim["batch"].values == "GDC2024").astype(float)[:, None]
    rho, p, _ = common.partial_spearman(L[rows], L[t.row[GPR153]], batch)
    c = pd.DataFrame({"gene_id": genes, "symbol": t.symbol[rows], "mean_log2tpm": L[rows].mean(1), "rho": rho, "p": p,
                      "q": common.bh(p)})
    c.to_csv(out / "gpr153_correlation.csv", index=False)
    s = c[c["symbol"] != ""].sort_values("mean_log2tpm", ascending=False).drop_duplicates("symbol")
    P = common.PARAMS["ora"]
    up = s[(s.q < P["q"]) & (s.rho > 0)].nlargest(P["max_genes"], "rho")
    down = s[(s.q < P["q"]) & (s.rho < 0)].nsmallest(P["max_genes"], "rho")
    pd.concat([up.assign(direction="positive"), down.assign(direction="negative")]).to_csv(out / "ora_lists.csv", index=False)
    for lib in P["libraries"]:
        sets = read_gmt(GS / f"{lib}.gmt")
        for lab, lst in (("positive", up), ("negative", down)):
            if len(lst):
                ora(lst["symbol"], s["symbol"], sets, P).to_csv(out / f"ora_{lib}_{lab}.csv", index=False)
    G = common.PARAMS["gsea"]
    rnk = s.set_index("symbol")["rho"].sort_values(ascending=False)
    targeted = []
    for lib in G["libraries"]:
        res = gseapy.prerank(rnk=rnk, gene_sets=read_gmt(GS / f"{lib}.gmt"), min_size=G["min_size"],
                             max_size=G["max_size"], permutation_num=G["permutations"], seed=common.SEED, threads=4,
                             outdir=None, no_plot=True, verbose=False).res2d
        res.to_csv(out / f"gsea_{lib}.csv", index=False)
        targeted.append(res[res["Term"].isin(common.PARAMS["targeted_sets"])].assign(library=lib, method="GSEA preranked"))
    pd.concat(targeted).to_csv(out / "targeted_sets.csv", index=False)
    mgmt_enrichment(gl, s, out)
    cam = gl.loc[gl["eligible_20"] & (gl["gene_id"] != GPR153), ["gene_id", "symbol"]]
    cam = cam.assign(mean=L[[t.row[g] for g in cam.gene_id]].mean(1))
    cam = cam[cam.symbol.notna() & (cam.symbol != "")].sort_values("mean", ascending=False).drop_duplicates("symbol")
    pd.DataFrame(L[[t.row[g] for g in cam.gene_id]], index=cam.symbol, columns=prim["aliquot"]).to_csv(
        common.work_path("camera_expression.csv.gz"))
    pd.DataFrame({"aliquot": prim["aliquot"], "gpr153": L[t.row[GPR153]], "batch_gdc": batch[:, 0]}).to_csv(
        common.work_path("camera_covariates.csv"), index=False)
    print(f"ORA lists: {len(up)} positive, {len(down)} negative; GSEA libraries: {len(G['libraries'])}")


if __name__ == "__main__":
    main()
