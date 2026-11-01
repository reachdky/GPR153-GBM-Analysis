#!/usr/bin/env python3
"""
T9 (PREREGISTRATION.md, Section 13.1): source-code guard checks (static scan; reads no data).
  - no `mannwhitneyu` in survival scripts (pipeline/a2*.py, a31_*.py)
  - no `pearsonr` in the screen script (pipeline/a11_screen.py)
  - no `penalizer` argument other than 0 in any pipeline script
  - no typed-in statistic in a title, axis label, legend label or annotation string of a plotting script
    (pipeline/f1*.py, plotstyle.py): string constants passed to set_title, set_xlabel, set_ylabel, suptitle, text,
    annotate or a label= argument must not contain a decimal number or a comparison or equals sign followed by a
    number (e.g. "0.31", "p = 0.002", "n = 391", ">= 20"). Values formatted from data inside f-strings are allowed.
scan(root) is also used by T9b (tests/test_guards_planted.py) on planted copies.
Run from the repository root:  python tests/test_guards.py
"""

import ast
import re
import sys
from pathlib import Path

STAT = re.compile(r"\d+\.\d+|[=<>]\s*-?\d|\b[np]\s*[=<>]")
TEXT_FUNCS = {"set_title", "set_xlabel", "set_ylabel", "suptitle", "text", "annotate"}


def constants(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        yield node.value
    elif isinstance(node, ast.JoinedStr):
        for v in node.values:
            if isinstance(v, ast.Constant):
                yield v.value
    elif isinstance(node, ast.BinOp):
        yield from constants(node.left)
        yield from constants(node.right)


def scan(P):
    """Guard-check problems in the pipeline directory P (a list of strings; empty if none)."""
    problems = []
    for f in sorted(P.glob("a2*.py")) + sorted(P.glob("a31_*.py")):
        if "mannwhitneyu" in f.read_text():
            problems.append(f"{f}: mannwhitneyu in a survival script")
    if "pearsonr" in (P / "a11_screen.py").read_text():
        problems.append(f"{P / 'a11_screen.py'}: pearsonr in the screen script")
    for f in sorted(P.rglob("*.py")):
        for node in ast.walk(ast.parse(f.read_text())):
            if isinstance(node, ast.keyword) and node.arg == "penalizer":
                if not (isinstance(node.value, ast.Constant) and node.value.value == 0):
                    problems.append(f"{f}:{node.value.lineno}: penalizer other than 0")
    for f in sorted(P.glob("f1*.py")) + [P / "plotstyle.py"]:
        for node in ast.walk(ast.parse(f.read_text())):
            if not isinstance(node, ast.Call):
                continue
            name = getattr(node.func, "attr", getattr(node.func, "id", ""))
            args = list(node.args[:1]) if name in TEXT_FUNCS else []
            args += [k.value for k in node.keywords if k.arg == "label"]
            for a in args:
                for s in constants(a):
                    if STAT.search(s):
                        problems.append(f"{f}:{node.lineno}: numeric literal in text {s!r}")
    return problems


if __name__ == "__main__":
    from _report import record
    problems = scan(Path("pipeline"))
    ok = record("T9", "guard checks (source-code scan)", not problems, "source code",
                "; ".join(problems) or "no mannwhitneyu in survival scripts, no pearsonr in the screen, penalizer 0 only, "
                                        "no typed-in statistics in figure text")
    sys.exit(0 if ok else 1)
