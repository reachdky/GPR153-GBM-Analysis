#!/usr/bin/env python3
"""
T6 (PREREGISTRATION.md, Sections 6.5 and 13.1), on R data and clinical tables only.
  known answer   lifelines logrank_test on R survival::aml (Maintained vs Nonmaintained) vs survdiff
                 chi-square 3.3963887, p = 0.0653393
  event column   inverting the event indicator must change the result
  sample identity  common.survival_table (used by the Cox and KM scripts) equals the manifest survival columns
                 and flags for each TCGA analysis cohort, and the survival-data audit counts
Run from the repository root:  python tests/test_survival.py
"""

import io
import subprocess
import sys

import numpy as np
import pandas as pd
from lifelines.statistics import logrank_test

from _report import record, common

out = subprocess.run(["Rscript", "-e", "library(survival); write.csv(aml, stdout(), row.names = FALSE)"],
                     check=True, capture_output=True, text=True).stdout
aml = pd.read_csv(io.StringIO(out))
a, b = aml[aml.x == "Maintained"], aml[aml.x == "Nonmaintained"]
r = logrank_test(a.time, b.time, a.status, b.status)
known = abs(r.test_statistic - 3.3963887) < 1e-6 and abs(r.p_value - 0.0653393) < 1e-6
inv = logrank_test(a.time, b.time, 1 - a.status, 1 - b.status)
changed = abs(inv.test_statistic - r.test_statistic) > 1e-6

m = pd.read_csv("revision/manifests/tcga_samples.csv")
ident = []
for flag, sflag in [("in_screen_primary", "in_survival_primary"), ("in_screen_sensitivity", "in_survival_sensitivity")]:
    s = common.survival_table(m, flag)
    ref = m[m[sflag] & m["selected_aliquot"]].sort_values("patient")
    s = s.sort_values("patient")
    same = (list(s.patient) == list(ref.patient)
            and np.allclose(s.os_days.values, ref.os_days.values)
            and (s.os_event.values == ref.os_event.astype(int).values).all()
            and s.patient.is_unique)
    ident.append((flag, same, len(s), int(s.os_event.sum())))
aud = pd.read_csv("revision/results/cohort/survival_audit.csv").set_index("cohort").loc["TCGA-GBM"]
audit_ok = ident[0][2] == aud.n_usable and ident[0][3] == aud.events
passed = known and changed and all(x[1] for x in ident) and audit_ok
detail = (f"aml chi-square {r.test_statistic:.7f}, p {r.p_value:.7f}; inverted events chi-square "
          f"{inv.test_statistic:.3f}, p {inv.p_value:.3f}; "
          + "; ".join(f"{f}: identical={s}, n {n}, events {e}" for f, s, n, e in ident)
          + f"; audit n {aud.n_usable}, events {aud.events}")
record("T6", "log-rank known answer, event inversion, survival sample identity", passed,
       "R survival::aml; TCGA manifest (clinical columns only)", detail)
sys.exit(0 if passed else 1)
