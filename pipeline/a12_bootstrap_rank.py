#!/usr/bin/env python3
"""
a12_bootstrap_rank.py
Rank stability of GPR153 in the primary screen (PREREGISTRATION.md, Section 5.1): 1,000 bootstrap resamples of
patients (seed 42). In each resample the eligible genes and MGMT are re-ranked, residualized on batch and
correlated (batch-adjusted partial Spearman, as in the primary run). The eligible gene set is fixed at the
primary analysis. Resamples run in parallel; each has its own seed from numpy SeedSequence(42), so results do
not depend on the number of workers.

Inputs
  revision/work/tcga_expr.npz (a10), revision/manifests/tcga_samples.csv (r01),
  revision/results/screen/screen_gene_level.csv (a11)
Outputs (revision/results/screen/)
  bootstrap_rank_draws.csv     per resample: GPR153 rho and rank, rho of ranks 1-10
  bootstrap_rank_summary.csv   proportion rank 1, proportion top 10, percentile 95% CI of rho, observed and
                               bootstrap median rho of ranks 1-10

Run from the repository root:
    python pipeline/a12_bootstrap_rank.py
"""

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import common
from common import GPR153, MGMT


def draw(seed, X, mg, batch, gpos):
    idx = np.random.default_rng(seed).choice(X.shape[1], X.shape[1], replace=True)
    rho, _, _ = common.partial_spearman(X[:, idx], mg[idx], batch[idx, None])
    rho = np.nan_to_num(rho, nan=-np.inf)  # a gene constant in a resample cannot rank
    top = np.sort(rho)[::-1][:10]
    rank = int((rho > rho[gpos]).sum() + 1)
    return [rho[gpos], rank, *top]


def main():
    common.set_seeds()
    common.require_registration()
    out = common.outdir("screen")
    t = common.TCGA()
    prim = common.cohort_aliquots(t.m, "in_screen_primary")
    gl = pd.read_csv(out / "screen_gene_level.csv")
    genes = gl.loc[gl["eligible_20"] & (gl["gene_id"] != MGMT), "gene_id"].to_numpy()
    rows = np.array([t.row[g] for g in genes])
    L = t.log2tpm(prim["aliquot"].values)
    X, mg = L[rows], L[t.row[MGMT]]
    batch = (prim["batch"].values == "GDC2024").astype(float)
    gpos = int(np.where(genes == GPR153)[0][0])
    B = common.PARAMS["n_bootstrap_rank"]
    seeds = np.random.SeedSequence(common.SEED).generate_state(B)
    res = Parallel(n_jobs=-1)(delayed(draw)(int(s), X, mg, batch, gpos) for s in seeds)
    d = pd.DataFrame(res, columns=["gpr153_rho", "gpr153_rank"] + [f"rho_rank{i}" for i in range(1, 11)])
    d.insert(0, "resample", np.arange(1, B + 1))
    d.to_csv(out / "bootstrap_rank_draws.csv", index=False)
    obs = gl.loc[gl["eligible_20"]].sort_values("rank_adj").head(10)
    lo, hi = np.percentile(d["gpr153_rho"], [2.5, 97.5])
    summ = {"resamples": B, "n": len(prim), "genes": len(genes),
            "prop_rank1": float((d["gpr153_rank"] == 1).mean()), "prop_top10": float((d["gpr153_rank"] <= 10).mean()),
            "rho_ci_low": lo, "rho_ci_high": hi, "rank_median": float(d["gpr153_rank"].median())}
    for i in range(1, 11):
        summ[f"observed_rho_rank{i}"] = float(obs["rho_adj"].iloc[i - 1])
        summ[f"observed_gene_rank{i}"] = obs["symbol"].iloc[i - 1]
        summ[f"bootstrap_median_rho_rank{i}"] = float(d[f"rho_rank{i}"].median())
    pd.DataFrame([summ]).to_csv(out / "bootstrap_rank_summary.csv", index=False)
    print({k: summ[k] for k in ("prop_rank1", "prop_top10", "rho_ci_low", "rho_ci_high")})


if __name__ == "__main__":
    main()
