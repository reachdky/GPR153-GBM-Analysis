"""Rebuild tables, figures and source data from an archived results directory with the current figure/table code.
Usage (from the repo root): python rebuild_outputs.py <archived results dir> <output base dir>
Writes only under <output base dir>."""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, "pipeline")
import common  # noqa: E402

res, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
for d in ("results",):
    if (out / d).exists():
        shutil.rmtree(out / d)
shutil.copytree(res, out / "results")
shutil.copytree("revision/manifests", out / "manifests", dirs_exist_ok=True)
common.BASE, common.RES, common.MAN = out, out / "results", out / "manifests"
common.SRC, common.FIG, common.FIGMAP = out / "source_data", out / "figures", out / "figure_script_map.csv"

import f10_main_figures  # noqa: E402
import f11_supplement_figures  # noqa: E402
import f12_tables  # noqa: E402

f10_main_figures.R = f11_supplement_figures.R = common.RES
f12_tables.R, f12_tables.OUT = common.RES, out / "tables"
f10_main_figures.main()
f11_supplement_figures.main()
f12_tables.main()
