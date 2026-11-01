#!/usr/bin/env python3
"""
T10 (PREREGISTRATION.md, Section 13.1): survival-data audit. Re-runs pipeline/r07_survival_audit.py (clinical
tables only), which asserts that n and events recomputed from the raw clinical sources equal the pipeline outputs,
and checks the usable n of every cohort against pipeline/r05's candidate table.
Run from the repository root:  python tests/test_survival_audit.py
"""

import subprocess
import sys

import pandas as pd

from _report import record

r = subprocess.run([sys.executable, "pipeline/r07_survival_audit.py"], capture_output=True, text=True,
                   env={**__import__("os").environ, "GBM_MODE": "real"})
a = pd.read_csv("revision/results/cohort/survival_audit.csv").set_index("cohort")
c = pd.read_csv("revision/results/cohort/survival_replication_candidates.csv").set_index("cohort")
pairs = {"TCGA-GBM": "TCGA-GBM (discovery)", "CGGA mRNAseq_693": "CGGA mRNAseq_693", "CGGA mRNAseq_325": "CGGA mRNAseq_325",
         "GSE13041 GPL96": "GSE13041 GPL96 primary", "Rembrandt (GSE108474)": "Rembrandt (GSE108474)", "CPTAC-3 GBM": "CPTAC-3 GBM"}
mism = [f"{k}: audit {a.loc[k, 'n_usable']} vs r05 {c.loc[v, 'n_usable']}" for k, v in pairs.items()
        if v in c.index and int(a.loc[k, "n_usable"]) != int(c.loc[v, "n_usable"])]
ok = r.returncode == 0 and "all survival assertions passed" in r.stdout and not mism
detail = "; ".join(f"{k} n {a.loc[k, 'n_usable']}, events {a.loc[k, 'events']}" for k in pairs)
record("T10", "survival-data audit", ok, "clinical tables", (detail + ("; " + "; ".join(mism) if mism else ""))
       if r.returncode == 0 else r.stderr[-300:])
sys.exit(0 if ok else 1)
