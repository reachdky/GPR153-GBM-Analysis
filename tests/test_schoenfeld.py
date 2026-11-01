#!/usr/bin/env python3
"""
T5 (PREREGISTRATION.md, Section 13.1): Schoenfeld residual tests, lifelines proportional_hazard_test
(time_transform "km") vs R survival::cox.zph (transform "km") on R survival::lung (age, sex).
lifelines gives per-covariate tests only; R also gives the global test. If any p differs by more than 1e-3,
the pipeline reports R cox.zph as primary (a20_cox.py reads this result) and the difference is logged.
Run from the repository root:  python tests/test_schoenfeld.py
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

from _report import record, common

with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    subprocess.run(["Rscript", "-e", f"library(survival); d <- lung[, c('time','status','age','sex')]; "
                    f"d$event <- as.integer(d$status == 2); write.csv(d[, c('time','event','age','sex')], "
                    f"'{td}/lung.csv', row.names = FALSE)"], check=True)
    pd.DataFrame([{"input": td / "lung.csv", "prefix": td / "lung", "covariates": "age,sex", "strata": ""}]).to_csv(
        td / "jobs.csv", index=False)
    common.run_r("pipeline/a20_coxzph.R", td / "jobs.csv")
    r = pd.read_csv(td / "lung_zph.csv").set_index("term")
    d = pd.read_csv(td / "lung.csv")
cph = CoxPHFitter(penalizer=0.0).fit(d, "time", "event")
py = proportional_hazard_test(cph, d, time_transform="km").summary
diffs = {t: abs(float(py.loc[t, "p"]) - float(r.loc[t, "p"])) for t in ["age", "sex"]}
passed = max(diffs.values()) <= 1e-3
detail = "; ".join(f"{t}: lifelines p {py.loc[t, 'p']:.4f}, R p {r.loc[t, 'p']:.4f}" for t in diffs)
detail += f"; R global p {r.loc['GLOBAL', 'p']:.4f} (no lifelines global test)"
if not passed:
    detail += "; R cox.zph is used as primary for all PH tests"
record("T5", "Schoenfeld tests vs R cox.zph (lung)", passed, "R survival::lung", detail)
common.write_json(common.outdir("tests") / "t5_ph_source.json",
                  {"ph_primary": "lifelines" if passed else "R cox.zph", "differences": diffs})
sys.exit(0)
