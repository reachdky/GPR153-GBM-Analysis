#!/usr/bin/env python3
"""
T1-T4 (PREREGISTRATION.md, Section 13.1), on simulated data (seed 42).
  T1  partial Spearman with no covariate (k = 0) vs scipy.stats.spearmanr
  T2  partial Spearman with covariates (k = 1 and k = 3) vs pingouin.partial_corr(method="spearman")
  T2b partial Spearman with a continuous, skewed covariate (k = 2: batch and the covariate) vs pingouin; the
      covariate enters as ranks (pingouin ranks every column); the difference from entering it raw is reported
  T3  BH adjustment vs statsmodels multipletests(method="fdr_bh")
  T4  Fisher-z / Bonett-Wright CI vs hand calculation
Run from the repository root:  python tests/test_stats.py
"""

import sys

import numpy as np
import pandas as pd
import pingouin
from scipy import stats
from statsmodels.stats.multitest import multipletests

from _report import record, common

rng = np.random.default_rng(common.SEED)
ok = True

# T1
n, G = 120, 200
x = rng.normal(size=n)
y = 0.3 * x + rng.normal(size=(G, n))
y[:20] = np.round(y[:20], 1)  # ties
rho, p, df = common.partial_spearman(y, x)
ref = np.array([stats.spearmanr(x, yy) for yy in y])
d1 = max(np.max(np.abs(rho - ref[:, 0])), np.max(np.abs(p - ref[:, 1])))
ok &= record("T1", "partial Spearman, k = 0, vs scipy.stats.spearmanr", d1 <= 1e-10, "simulated",
             f"max abs difference {d1:.2e}; tolerance 1e-10; {G} genes, n = {n}, df = {df}")

# T2
batch = (rng.random(n) < 0.45).astype(float)
plat = (rng.random(n) < 0.3).astype(float)
cont = rng.normal(size=n) + batch
x2 = x + 0.5 * batch
y2 = y + batch + 0.3 * cont
worst = 0.0
for cov_names, cov in [(["batch"], batch[:, None]), (["batch", "plat", "cont"], np.column_stack([batch, plat, cont]))]:
    rho2, p2, df2 = common.partial_spearman(y2[:50], x2, cov)
    for g in range(50):
        d = pd.DataFrame({"x": x2, "y": y2[g], **{c: cov[:, i] for i, c in enumerate(cov_names)}})
        r = pingouin.partial_corr(d, x="x", y="y", covar=cov_names, method="spearman")
        worst = max(worst, abs(r["r"].iloc[0] - rho2[g]), abs(r["p-val"].iloc[0] - p2[g]))
        assert int(r["n"].iloc[0]) == n
ok &= record("T2", "partial Spearman, k = 1 and k = 3, vs pingouin.partial_corr", worst <= 1e-8, "simulated",
             f"max abs difference in rho and p {worst:.2e}; tolerance 1e-8; pingouin {pingouin.__version__}")

# T2b
skew = np.exp(rng.normal(size=n) + 0.8 * batch)  # continuous, right-skewed, batch-related (like a composition score)
x3 = x + 0.4 * np.log(skew)
y3 = y[:50] + 0.6 * np.log(skew) + batch
cov3 = np.column_stack([batch, skew])
rho3, p3, df3 = common.partial_spearman(y3, x3, cov3)
worst3 = 0.0
for g in range(50):
    d = pd.DataFrame({"x": x3, "y": y3[g], "batch": batch, "score": skew})
    r = pingouin.partial_corr(d, x="x", y="y", covar=["batch", "score"], method="spearman")
    worst3 = max(worst3, abs(r["r"].iloc[0] - rho3[g]), abs(r["p-val"].iloc[0] - p3[g]))
raw = common.residualize(common.rank_rows(y3), common.design(cov3, n))
rx = common.residualize(stats.rankdata(x3), common.design(cov3, n))
rho_raw = (raw @ rx) / (np.linalg.norm(raw, axis=1) * np.linalg.norm(rx))
ok &= record("T2b", "partial Spearman, continuous covariate (ranked), k = 2, vs pingouin.partial_corr", worst3 <= 1e-8,
             "simulated", f"max abs difference in rho and p {worst3:.2e}; tolerance 1e-8; df = {df3}; "
                          f"entering the covariate raw instead would change rho by up to {np.max(np.abs(rho_raw - rho3)):.3f}")

# T3
pv = np.concatenate([rng.random(500), rng.random(200) * 1e-3, np.repeat(0.01, 20), [np.nan] * 5])
q = common.bh(pv)
okp = ~np.isnan(pv)
qref = multipletests(pv[okp], method="fdr_bh")[1]
d3 = np.max(np.abs(q[okp] - qref))
ok &= record("T3", "BH vs statsmodels fdr_bh", d3 <= 1e-12 and np.isnan(q[~okp]).all(), "simulated",
             f"max abs difference {d3:.2e}; tolerance 1e-12; NaN inputs returned as NaN")

# T4
hand = {1: (0.309520, 0.104333, 0.1046, 0.4731), 3: (0.309520, 0.105437, 0.1025, 0.4747)}
fails = []
for k, (z, se, lo, hi) in hand.items():
    clo, chi, cz, cse = common.bw_ci(0.30, 100, k)
    got = (round(float(cz), 6), round(float(cse), 6), round(float(clo), 4), round(float(chi), 4))
    if got != (z, se, lo, hi):
        fails.append(f"k={k}: got {got}, expected {(z, se, lo, hi)}")
ok &= record("T4", "Fisher-z and Bonett-Wright CI vs hand calculation", not fails, "hand calculation",
             "; ".join(fails) or "rho 0.30, n 100: k = 1 and k = 3 match to the stated decimals")

sys.exit(0 if ok else 1)
