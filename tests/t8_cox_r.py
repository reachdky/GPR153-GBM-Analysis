#!/usr/bin/env python3
"""
T8 (PREREGISTRATION.md, Section 13.1): primary Cox model (Model 3) re-run in R survival::coxph (Efron ties) on the
same data as a20; relative difference in HR, CI limits and p <= 1e-4. Reads expression values, so in real mode it
is locked until registration; in sim mode it runs on simulated data.
Run from the repository root:  python tests/t8_cox_r.py
"""

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from _report import record, common

common.require_registration()
d = pd.read_csv(common.work_path("survival_table_primary.csv"))[["time", "event", "gpr153", "age", "male"]].dropna()
py = pd.read_csv(common.RES / "survival/cox_models.csv")
py = py[(py.cohort == "primary") & (py.model == "Model 3") & (py.variant == "main") & (py.sample_set == "complete cases")].set_index("covariate")
with tempfile.TemporaryDirectory(dir=common.work_path("")) as td:
    td = Path(td)
    d.to_csv(td / "m3.csv", index=False)
    pd.DataFrame([{"input": td / "m3.csv", "prefix": td / "m3", "covariates": "gpr153,age,male", "strata": ""}]).to_csv(td / "jobs.csv", index=False)
    common.run_r("pipeline/a20_coxzph.R", td / "jobs.csv")
    r = pd.read_csv(td / "m3_coxph.csv").set_index("term")
worst = 0.0
for v in r.index:
    for a, b in (("hr", "hr"), ("hr_low", "hr_low"), ("hr_high", "hr_high"), ("p", "p")):
        worst = max(worst, abs(py.loc[v, a] - r.loc[v, b]) / abs(r.loc[v, b]))
ok = record("T8", "Cox Model 3 re-run in R survival::coxph", worst <= 1e-4 and len(d) == int(py.n.iloc[0]), common.MODE,
            f"n {len(d)}; max relative difference in HR, CI and p {worst:.2e}; tolerance 1e-4")
sys.exit(0 if ok else 1)
