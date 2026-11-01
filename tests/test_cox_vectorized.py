#!/usr/bin/env python3
"""
T15 (PREREGISTRATION.md, Section 13.1): vectorized univariate Cox (common.cox_univariate, used by a22) vs R survival::coxph
(Efron ties) on simulated data with tied event times (seed 42); tolerance 1e-6 on coefficient and SE. The
difference between lifelines CoxPHFitter and R on the same data is reported for information (lifelines stops at a
looser convergence tolerance).
Run from the repository root:  python tests/test_cox_vectorized.py
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

from _report import record, common

rng = np.random.default_rng(common.SEED)
n, G = 200, 25
X = rng.normal(size=(G, n))
time = np.round(rng.exponential(300, n) * np.exp(-0.3 * X[0]), -1) + 10  # rounding creates ties
event = (rng.random(n) < 0.8).astype(int)
b, se, _, _ = common.cox_univariate(X, time, event)
with tempfile.TemporaryDirectory() as td:
    f = Path(td) / "d.csv"
    pd.DataFrame({"time": time, "event": event, **{f"x{g}": X[g] for g in range(G)}}).to_csv(f, index=False)
    code = (f"library(survival); d <- read.csv('{f}'); for (g in 0:{G - 1}) {{ m <- coxph(as.formula(paste0("
            f"'Surv(time, event) ~ x', g)), data = d, ties = 'efron'); cat(coef(m), sqrt(vcov(m)), '\\n') }}")
    out = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=True).stdout
r = np.array([list(map(float, line.split())) for line in out.strip().splitlines()])
d_r = max(np.abs(b - r[:, 0]).max(), np.abs(se - r[:, 1]).max())
ll = np.array([CoxPHFitter(penalizer=0.0).fit(pd.DataFrame({"t": time, "e": event, "x": X[g]}), "t", "e")
               .summary.loc["x", ["coef", "se(coef)"]].to_numpy(float) for g in range(G)])
d_l = max(np.abs(ll[:, 0] - r[:, 0]).max(), np.abs(ll[:, 1] - r[:, 1]).max())
ok = record("T15", "vectorized univariate Cox vs R coxph (Efron ties)", d_r <= 1e-6, "simulated",
            f"max abs difference in coef and SE vs R {d_r:.2e} (tolerance 1e-6); lifelines vs R {d_l:.2e}; "
            f"{G} genes, n = {n}, {n - len(np.unique(time))} tied times")
sys.exit(0 if ok else 1)
