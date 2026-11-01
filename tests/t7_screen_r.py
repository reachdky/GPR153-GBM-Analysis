#!/usr/bin/env python3
"""
T7 (PREREGISTRATION.md, Section 13.1): primary screen re-run in R (tests/t7_screen.R: rank, lm residuals on batch,
cor) on the same eligible genes; max abs difference in rho <= 1e-8 and identical rank order. Reads expression
values, so in real mode it is locked until registration; in sim mode it runs on simulated data.
Run from the repository root:  python tests/t7_screen_r.py
"""

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from _report import record, common

common.require_registration()
t = common.TCGA()
prim = common.cohort_aliquots(t.m, "in_screen_primary")
gl = pd.read_csv(common.RES / "screen/screen_gene_level.csv")
el = gl[gl.eligible_20 & gl.rank_adj.notna()]
L = t.log2tpm(prim["aliquot"].values)
with tempfile.TemporaryDirectory(dir=common.work_path("")) as td:
    td = Path(td)
    rows = [t.row[g] for g in el.gene_id] + [t.row[common.MGMT]]
    pd.DataFrame(L[rows], index=list(el.gene_id) + ["MGMT_ROW"], columns=prim["aliquot"]).to_csv(td / "x.csv.gz")
    pd.DataFrame({"aliquot": prim["aliquot"], "batch_gdc": (prim["batch"] == "GDC2024").astype(int)}).to_csv(td / "b.csv", index=False)
    common.run_r("tests/t7_screen.R", td / "x.csv.gz", td / "b.csv", td / "r.csv")
    r = pd.read_csv(td / "r.csv").merge(el[["gene_id", "rho_adj", "rank_adj"]], on="gene_id")
d = float(np.max(np.abs(r.rho_r - r.rho_adj)))
same_order = (r.sort_values("rho_r", ascending=False, kind="stable").gene_id.values
              == r.sort_values("rank_adj").gene_id.values).all()
ok = record("T7", "primary screen re-run in R", d <= 1e-8 and same_order, common.MODE,
            f"{len(r)} genes; max abs difference in rho {d:.2e}; identical rank order {same_order}")
sys.exit(0 if ok else 1)
