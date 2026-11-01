"""
plotstyle.py
Shared figure style and output for the figure scripts, to the PLOS ONE figure specifications
(journals.plos.org/plosone/s/figures, checked 2026-11-01): TIFF, LZW compression, 300 dpi, RGB without alpha
channel, flattened; width 2.63-7.5 in and height at most 8.75 in (789-2250 x at most 2625 px at 300 dpi); Arial,
Times or Symbol in 8-12 pt; main figures named Fig1.tif, Fig2.tif, ...; supporting figures S1_Fig.tif, ...
save() stops if a figure is outside these limits or any text is smaller than 8 pt. Arial is used when it is
installed; otherwise DejaVu Sans, and the font used is written into each TIFF (ImageDescription) and printed.
One source-data CSV per figure and a row in figure_script_map.csv. Figure text is filled from result files only.
"""

import io
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.text import Text  # noqa: E402
from PIL import Image  # noqa: E402

import common  # noqa: E402

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
DPI = 300
MIN_PT, MAX_PT = 8, 12
MIN_W_PX, MAX_W_PX, MAX_H_PX = 789, 2250, 2625


def _font():
    try:
        font_manager.findfont(font_manager.FontProperties(family="Arial"), fallback_to_default=False)
        return "Arial"
    except ValueError:
        return "DejaVu Sans"


FONT = _font()

plt.rcParams.update({
    "font.family": FONT, "font.size": MIN_PT, "axes.titlesize": 9, "axes.labelsize": MIN_PT,
    "xtick.labelsize": MIN_PT, "ytick.labelsize": MIN_PT, "legend.fontsize": MIN_PT, "figure.titlesize": 9,
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
    return ax.legend(loc="upper center", bbox_to_anchor=(0.5, y), ncol=ncol, markerscale=1.5)


def file_stem(item):
    """PLOS file names: "Fig 1" -> Fig1, "S1 Fig" -> S1_Fig."""
    return item.replace(" ", "") if item.startswith("Fig ") else item.replace(" ", "_")


def check_text(fig, item):
    small = sorted({f"{t.get_text()[:30]!r} {t.get_fontsize():g} pt" for t in fig.findobj(Text)
                    if t.get_visible() and t.get_text().strip() and not MIN_PT <= t.get_fontsize() <= MAX_PT})
    if small:
        raise SystemExit(f"{item}: text outside {MIN_PT}-{MAX_PT} pt: " + "; ".join(small))


def save(fig, item, script, sources, frames=None):
    """Write <FIG>/<stem>.tif and <SRC>/<stem>.csv (concatenated frames, if given) and record the mapping."""
    common.FIG.mkdir(parents=True, exist_ok=True)
    common.SRC.mkdir(parents=True, exist_ok=True)
    stem = file_stem(item)
    path = common.FIG / f"{stem}.tif"
    check_text(fig, item)
    buf = io.BytesIO()
    fig.savefig(buf, dpi=DPI, format="png", bbox_inches="tight", pad_inches=2 / 72, facecolor="white")
    plt.close(fig)
    im = Image.open(buf)
    rgb = Image.new("RGB", im.size, "white")
    rgb.paste(im, mask=im.getchannel("A") if im.mode == "RGBA" else None)
    w, h = rgb.size
    if not (MIN_W_PX <= w <= MAX_W_PX and h <= MAX_H_PX):
        raise SystemExit(f"{item}: {w} x {h} px at {DPI} dpi is outside {MIN_W_PX}-{MAX_W_PX} x at most {MAX_H_PX} px")
    rgb.save(path, format="TIFF", compression="tiff_lzw", dpi=(DPI, DPI),
             tiffinfo={270: f"{item}; font {FONT}; {script}"})
    src = []
    if frames is not None:
        import pandas as pd
        sd = common.SRC / f"{stem}.csv"
        pd.concat([f.assign(panel=k) for k, f in frames.items()], ignore_index=True).to_csv(sd, index=False)
        src.append(sd)
    common.record_figure(item, script, path, [*src, *sources])
    print(f"{item}: {path.name}, {w} x {h} px, font {FONT}")
    return path
