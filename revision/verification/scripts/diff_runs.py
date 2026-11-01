"""Row-level diff of run 1 vs run 2 outputs. Usage: python diff_runs.py <run1 base> <run2 base> <out csv>
Each base holds results/, tables/, source_data/, figures/."""
import filecmp
import sys
from pathlib import Path

import pandas as pd

a, b, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
ID = ["cohort", "model", "variant", "sample_set", "covariate", "term", "test", "subset", "analysis", "quantity",
      "with_tcga", "gene", "split", "group", "part", "panel", "section", "item", "name", "run"]
rows = []
for d in ("results", "tables", "source_data", "figures"):
    fa = {p.relative_to(a / d) for p in (a / d).rglob("*") if p.is_file()}
    fb = {p.relative_to(b / d) for p in (b / d).rglob("*") if p.is_file()}
    for f in sorted(fa | fb):
        pa, pb = a / d / f, b / d / f
        if f not in fa or f not in fb:
            rows.append({"file": f"{d}/{f}", "change": "only in run 1" if f in fa else "only in run 2"})
            continue
        if filecmp.cmp(pa, pb, shallow=False):
            continue
        r = {"file": f"{d}/{f}", "change": "content differs"}
        if f.suffix == ".csv":
            x, y = pd.read_csv(pa, low_memory=False), pd.read_csv(pb, low_memory=False)
            r["rows_run1"], r["rows_run2"] = len(x), len(y)
            key = [c for c in ID if c in x.columns and c in y.columns]
            ka = x[key].astype(str).agg(" | ".join, axis=1) if key else pd.Series(dtype=str)
            kb = y[key].astype(str).agg(" | ".join, axis=1) if key else pd.Series(dtype=str)
            sa, sb = set(x.astype(str).agg("\x1f".join, axis=1)), set(y.astype(str).agg("\x1f".join, axis=1))
            ch_a = x[~x.astype(str).agg("\x1f".join, axis=1).isin(sb)]
            ch_b = y[~y.astype(str).agg("\x1f".join, axis=1).isin(sa)]
            r["rows_changed_or_removed"], r["rows_changed_or_added"] = len(ch_a), len(ch_b)
            if key:
                grp = [c for c in ("cohort", "analysis", "model", "test", "part", "panel", "quantity") if c in key][:2]
                g = (pd.concat([ch_a, ch_b])[grp].astype(str).agg(" / ".join, axis=1).value_counts() if grp else None)
                r["where"] = "; ".join(f"{k} ({v})" for k, v in g.items()) if g is not None else ""
            cols = [c for c in x.columns if c in y.columns and len(x) == len(y)
                    and not x[c].astype(str).reset_index(drop=True).equals(y[c].astype(str).reset_index(drop=True))]
            r["columns_differing"] = ", ".join(cols) if len(x) == len(y) else "row count differs"
        rows.append(r)
pd.DataFrame(rows).to_csv(out, index=False)
print(pd.DataFrame(rows)[["file", "change"] + [c for c in ("rows_run1", "rows_run2", "rows_changed_or_removed",
                                                             "rows_changed_or_added") if c in pd.DataFrame(rows)]].to_string())
