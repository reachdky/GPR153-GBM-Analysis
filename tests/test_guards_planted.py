#!/usr/bin/env python3
"""
T9b (PREREGISTRATION.md, Section 13.1): the T9 scan must fail on planted violations. The pipeline directory is
copied to a temporary directory, one forbidden pattern is planted in the copy, and T9's scan is run on it; this is
repeated for each pattern. The repository's own files are not changed. Passes only if the unmodified copy scans
clean and every planted pattern is reported in the planted file.
  1  `from scipy.stats import mannwhitneyu` appended to a survival script (a21_km.py)
  2  `from scipy.stats import pearsonr` appended to the screen script (a11_screen.py)
  3  a call with penalizer=0.1 appended to a Cox script (a20_cox.py)
  4  ax.set_title("p = 0.002") appended to a plotting script (f10_main_figures.py)
Run from the repository root:  python tests/test_guards_planted.py
"""

import shutil
import sys
import tempfile
from pathlib import Path

from _report import record
from test_guards import scan

PLANTS = [
    ("a21_km.py", "\nfrom scipy.stats import mannwhitneyu  # planted\n", "mannwhitneyu"),
    ("a11_screen.py", "\nfrom scipy.stats import pearsonr  # planted\n", "pearsonr"),
    ("a20_cox.py", "\n_planted = CoxPHFitter(penalizer=0.1)\n", "penalizer"),
    ("f10_main_figures.py", '\n\ndef _planted(ax):\n    ax.set_title("p = 0.002")\n', "numeric literal"),
]

results = []
with tempfile.TemporaryDirectory() as td:
    clean = Path(td) / "clean"
    shutil.copytree("pipeline", clean, ignore=shutil.ignore_patterns("__pycache__"))
    base = scan(clean)
    for f, code, key in PLANTS:
        d = Path(td) / f"plant_{f}"
        shutil.copytree(clean, d)
        (d / f).write_text((d / f).read_text() + code)
        found = [p for p in scan(d) if f in p and key in p]
        results.append((f, key, bool(found)))
ok = not base and all(r[2] for r in results)
record("T9b", "T9 scan fails on planted violations", ok, "source code (temporary copies)",
       f"unmodified copy: {len(base)} problems; " + "; ".join(f"{k} in {f}: {'detected' if x else 'NOT detected'}"
                                                             for f, k, x in results))
sys.exit(0 if ok else 1)
