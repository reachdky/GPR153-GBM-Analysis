#!/usr/bin/env python3
"""
a13_permutation_null.py
Family-wise permutation null for the primary screen (PREREGISTRATION.md, Section 5.4): 10,000 permutations
(seed 42) of MGMT within batch, with the full batch-adjusted screen over the eligible genes recomputed in each.

Because MGMT is permuted within batch, its ranks are permuted in the same way and its batch-residualized rank
vector is the permuted residual vector; the gene residuals do not change. Each permutation's screen is therefore
the product of the fixed, normalized gene-residual matrix with the permuted, normalized MGMT residuals, computed
in blocks of 500 permutations. This equals recomputing the screen in each permutation.

Reported (descriptive, not a pass condition): proportion of permutations with maximum rho >= observed maximum;
proportion in which the gene at GPR153's observed rank has rho >= GPR153's observed rho; the permutation
distribution of the maximum rho.

Inputs
  revision/work/tcga_expr.npz (a10), revision/manifests/tcga_samples.csv (r01),
  revision/results/screen/screen_gene_level.csv (a11)
Outputs (revision/results/screen/)
  permutation_null_draws.csv    per permutation: maximum rho, rho at GPR153's observed rank
  permutation_null_summary.csv

Run from the repository root:
    python pipeline/a13_permutation_null.py
"""

import numpy as np
import pandas as pd

import common
from common import GPR153, MGMT


def main():
    common.set_seeds()
    common.require_registration()
    out = common.outdir("screen")
    t = common.TCGA()
    prim = common.cohort_aliquots(t.m, "in_screen_primary")
    gl = pd.read_csv(out / "screen_gene_level.csv")
    genes = gl.loc[gl["eligible_20"] & (gl["gene_id"] != MGMT), "gene_id"].to_numpy()
    L = t.log2tpm(prim["aliquot"].values)
    batch = (prim["batch"].values == "GDC2024").astype(float)
    z = common.design(batch[:, None], len(prim))
    Y = common.residualize(common.rank_rows(L[[t.row[g] for g in genes]]), z)
    Y /= np.linalg.norm(Y, axis=1, keepdims=True)
    m = common.residualize(common.rank_rows(L[t.row[MGMT]]), z)
    m /= np.linalg.norm(m)
    obs = Y @ m
    gpos = int(np.where(genes == GPR153)[0][0])
    obs_rank = int((obs > obs[gpos]).sum() + 1)
    assert np.isclose(obs[gpos], gl.set_index("gene_id").loc[GPR153, "rho_adj"])

    rng = np.random.default_rng(common.SEED)
    groups = [np.where(batch == b)[0] for b in (0.0, 1.0)]
    P, blk = common.PARAMS["n_permutations"], common.PARAMS["permutation_block"]
    mx, at_rank = np.empty(P), np.empty(P)
    for s in range(0, P, blk):
        nb = min(blk, P - s)
        M = np.empty((len(m), nb))
        for j in range(nb):
            perm = np.arange(len(m))
            for g in groups:
                perm[g] = rng.permutation(g)
            M[:, j] = m[perm]
        R = Y @ M
        mx[s:s + nb] = R.max(0)
        at_rank[s:s + nb] = -np.partition(-R, obs_rank - 1, axis=0)[obs_rank - 1]
    pd.DataFrame({"permutation": np.arange(1, P + 1), "max_rho": mx, "rho_at_gpr153_rank": at_rank}).to_csv(
        out / "permutation_null_draws.csv", index=False)
    summ = {"permutations": P, "genes": len(genes), "n": len(prim), "observed_max_rho": float(obs.max()),
            "observed_gpr153_rho": float(obs[gpos]), "observed_gpr153_rank": obs_rank,
            "prop_max_ge_observed_max": float((mx >= obs.max()).mean()),
            "prop_rank_value_ge_gpr153": float((at_rank >= obs[gpos]).mean()),
            "max_rho_p50": float(np.median(mx)), "max_rho_p95": float(np.percentile(mx, 95)),
            "max_rho_p99": float(np.percentile(mx, 99))}
    pd.DataFrame([summ]).to_csv(out / "permutation_null_summary.csv", index=False)
    print(summ)


if __name__ == "__main__":
    main()
