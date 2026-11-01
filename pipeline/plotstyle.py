"""
plotstyle.py
Shared figure style and output for the figure scripts: TIFF, 300 dpi, LZW compression, one source-data CSV per
figure, and a row in figure_script_map.csv. Figure text is filled from result files only.
"""

import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import common  # noqa: E402

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.5, "lines.linewidth": 1.5, "legend.frameon": False, "figure.dpi": 100,
})


def fmt_p(p):
    return f"{p:.2g}" if p >= 1e-3 else f"{p:.1e}"


def fmt_perm_q(q, permutations):
    """A permutation FDR of 0 is reported as below the resolution 1 / permutations."""
    return f"< {1 / permutations:g}" if q == 0 else fmt_p(q)


def short(term, n=40):
    """Shorten a gene-set name to n characters, keeping a trailing GO or Reactome identifier whole."""
    m = re.match(r"^(.*?)\s*(\(GO:\d+\)|R-HSA-\d+)$", term)
    name, tag = (m.group(1), m.group(2)) if m else (term, "")
    if len(name) > n:
        name = name[:n - 1].rstrip() + "…"
    return f"{name} {tag}".strip()


def legend_below(ax, ncol=2, y=-0.2):
    """Legend under the axes, clear of the data."""
    return ax.legend(loc="upper center", bbox_to_anchor=(0.5, y), ncol=ncol, fontsize=6, markerscale=1.5)


def save(fig, item, script, sources, frames=None):
    """Write <FIG>/<item>.tif and <SRC>/<item>.csv (concatenated frames, if given) and record the mapping."""
    common.FIG.mkdir(parents=True, exist_ok=True)
    common.SRC.mkdir(parents=True, exist_ok=True)
    path = common.FIG / f"{item.replace(' ', '_')}.tif"
    fig.savefig(path, dpi=300, format="tiff", pil_kwargs={"compression": "tiff_lzw"}, bbox_inches="tight")
    plt.close(fig)
    src = []
    if frames is not None:
        import pandas as pd
        sd = common.SRC / f"{item.replace(' ', '_')}.csv"
        pd.concat([f.assign(panel=k) for k, f in frames.items()], ignore_index=True).to_csv(sd, index=False)
        src.append(sd)
    common.record_figure(item, script, path, [*src, *sources])
    return path
