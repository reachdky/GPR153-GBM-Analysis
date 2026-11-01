#!/usr/bin/env python3
"""
T13 (PREREGISTRATION.md, Section 13.1): Cox recovery on simulated data (seed 42). 200 data sets of n = 252 with
age (normal, mean 58, SD 12; log HR 0.03 per year), sex (60% male; log HR 0.10), GPR153 on the log2 scale (normal,
SD 1.2; true HR 1.25 per SD), Weibull baseline hazard (shape 1.3), uniform censoring (about 17% censored, as in the
primary cohort) and times rounded to whole days (ties). GPR153 is standardized within each data set as in a20.
Model 3 (GPR153 + age + sex) is fitted with a20_cox.fit (lifelines, penalizer 0, Efron ties). Pass if the 95% CI
for GPR153 covers log(1.25) in 0.90-0.99 of the data sets.
Run from the repository root:  python tests/test_cox_recovery.py
"""

import sys

import numpy as np
import pandas as pd

from _report import record, common

sys.path.insert(0, "pipeline")
from a20_cox import MODELS, fit  # noqa: E402

rng = np.random.default_rng(common.SEED)
RUNS, n, HR = 200, 252, 1.25
cover, est, cens = [], [], []
for _ in range(RUNS):
    age = rng.normal(58, 12, n)
    male = (rng.random(n) < 0.6).astype(float)
    g = rng.normal(5, 1.2, n)
    lp = np.log(HR) * (g - 5) / 1.2 + 0.03 * (age - 58) + 0.10 * male
    tt = 400 * (-np.log(rng.random(n)) / np.exp(lp)) ** (1 / 1.3)
    cc = rng.uniform(0, 2600, n)
    d = pd.DataFrame({"time": np.maximum(np.round(np.minimum(tt, cc)), 1), "event": (tt <= cc).astype(int),
                      "gpr153": common.standardize(g), "age": age, "male": male})
    c, x = fit(d, MODELS["Model 3"])
    s = c.summary.loc["gpr153"]
    cover.append(s["coef lower 95%"] <= np.log(HR) <= s["coef upper 95%"])
    est.append(s["coef"])
    cens.append(1 - d["event"].mean())
cov = float(np.mean(cover))
ok = record("T13", "Cox Model 3 recovery: 95% CI coverage of the true HR per SD (200 data sets)", 0.90 <= cov <= 0.99,
            "simulated", f"coverage {cov:.3f} (required 0.90-0.99); true HR {HR}, mean estimated HR "
                         f"{np.exp(np.mean(est)):.3f}; mean censored fraction {np.mean(cens):.2f}; n = {n}")
sys.exit(0 if ok else 1)
