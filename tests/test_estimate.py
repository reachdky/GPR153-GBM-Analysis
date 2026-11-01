#!/usr/bin/env python3
"""
T14 (PREREGISTRATION.md, Section 13.1): the Python ESTIMATE implementation in a14_composition.py vs tidyestimate 1.1.1's own
estimate_score() R code (sourced from the package tarball), on simulated data (seed 42).
Run from the repository root:  python tests/test_estimate.py
"""

import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from _report import record, common

sys.path.insert(0, "pipeline")
import a14_composition  # noqa: E402

rng = np.random.default_rng(common.SEED)
sets, cg = a14_composition.load_estimate_sets()
symbols = np.array(sorted(cg)[:3000] + sets["stromal_signature"][:100] + sets["immune_signature"][:100])
symbols = np.unique(symbols)
V = np.round(rng.gamma(2, 2, (len(symbols), 12)), 1)  # rounding creates ties
py = a14_composition.estimate_scores(V, symbols, {"stromal": sets["stromal_signature"],
                                                  "immune": sets["immune_signature"]})
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    with tarfile.open(common.DL / "estimate/tidyestimate_1.1.1.tar.gz") as tf:
        tf.extractall(td, filter="data")
    pd.DataFrame(V, columns=[f"s{i}" for i in range(12)]).assign(gene=symbols)[["gene"] + [f"s{i}" for i in range(12)]].to_csv(td / "v.csv", index=False)
    code = (f"suppressMessages({{library(dplyr)}}); load('{td}/tidyestimate/data/gene_sets.rda'); "
            f"src <- readLines('{td}/tidyestimate/R/estimate_score.R'); src <- gsub('tidyestimate::gene_sets', 'gene_sets', src, fixed = TRUE); "
            f"eval(parse(text = src)); d <- read.csv('{td}/v.csv'); s <- suppressMessages(estimate_score(d, is_affymetrix = FALSE)); "
            f"write.csv(s, '{td}/r.csv', row.names = FALSE)")
    subprocess.run(["Rscript", "-e", code], check=True, env=dict(__import__("os").environ, R_LIBS_USER=str(Path(".Rlib").resolve())))
    r = pd.read_csv(td / "r.csv")
d = max(np.max(np.abs(r["stromal"].values - py["stromal"])), np.max(np.abs(r["immune"].values - py["immune"])))
ok = record("T14", "ESTIMATE implementation vs tidyestimate estimate_score() R code", d <= 1e-8, "simulated",
            f"max abs difference {d:.2e}; tolerance 1e-8; {len(symbols)} genes, 12 samples, with ties")
sys.exit(0 if ok else 1)
