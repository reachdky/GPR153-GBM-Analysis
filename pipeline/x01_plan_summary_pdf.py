#!/usr/bin/env python3
"""
x01_plan_summary_pdf.py
Renders the summary (approval) page of revision/PREREGISTRATION.md to a one-page PDF (US Letter, 0.6-inch
margins, 9-point text) and checks that it fits: the script stops (exit 1) if the text runs past the bottom margin.
The page header gives the plan's title, version line and SHA-256. Reads only the plan file.

Outputs  revision/PREREGISTRATION_summary.pdf, revision/PREREGISTRATION_summary.png (preview)

Run from the repository root:
    python pipeline/x01_plan_summary_pdf.py
"""

import hashlib
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.font_manager import FontProperties  # noqa: E402

PLAN = Path("revision/PREREGISTRATION.md")
OUT = Path("revision/PREREGISTRATION_summary")
W, H, M = 8.5, 11.0, 0.6
SIZE, LEAD = 9.0, 1.32


def blocks(md):
    """Summary section -> list of (kind, indent, text): kind in heading, para, bullet, number."""
    i = md.index("## Summary")
    j = md.index("\n---", i)
    out = [("heading", 0, md[i:j].splitlines()[0].lstrip("# "))]
    for line in md[i:j].splitlines()[1:]:
        if not line.strip():
            out.append(("gap", 0, ""))
        elif re.match(r"^\s*- ", line):
            out.append(("bullet", (len(line) - len(line.lstrip())) // 2, line.strip()[2:]))
        elif re.match(r"^\s*\d+\. ", line):
            num, rest = line.strip().split(" ", 1)
            out.append(("number", 0, f"{num} {rest}"))
        elif line.startswith(" ") and out and out[-1][0] in ("bullet", "number", "para"):
            k, ind, t = out[-1]
            out[-1] = (k, ind, t + " " + line.strip())
        elif out and out[-1][0] == "para":
            out[-1] = ("para", 0, out[-1][2] + " " + line.strip())
        else:
            out.append(("para", 0, line.strip()))
    return out


def words(text):
    """Split into (word, bold) tokens; ** toggles bold; backticks are dropped."""
    toks, bold = [], False
    for part in re.split(r"(\*\*)", text.replace("`", "")):
        if part == "**":
            bold = not bold
            continue
        toks += [(w, bold) for w in part.split()]
    return toks


def main():
    md = PLAN.read_text()
    sha = hashlib.sha256(PLAN.read_bytes()).hexdigest()
    title = md.splitlines()[0].lstrip("# ")
    version = " ".join(md.split("\n\n")[1].splitlines()[:2])
    fig = plt.figure(figsize=(W, H))
    r = fig.canvas.get_renderer()
    fp = {False: FontProperties(family="DejaVu Sans", size=SIZE), True: FontProperties(family="DejaVu Sans", size=SIZE, weight="bold")}
    width = lambda w, prop: r.get_text_width_height_descent(w, prop, ismath=False)[0] / fig.dpi  # inches
    space = width(" ", fp[False])
    line_h = SIZE * LEAD / 72
    y = H - M

    def put(x, txt, prop):
        fig.text(x / W, y / H, txt, fontproperties=prop, va="top", ha="left")

    for txt, size, bold in ((title, 11, True), (version, 8, False), (f"SHA-256 of the plan file: {sha}", 7, False)):
        prop = FontProperties(family="DejaVu Sans", size=size, weight="bold" if bold else "normal")
        x = M
        for w in txt.split():
            if x + width(w, prop) > W - M and x > M:
                y -= size * LEAD / 72
                x = M
            put(x, w, prop)
            x += width(w + " ", prop)
        y -= size * LEAD / 72
    y -= line_h * 0.5
    right = W - M
    for kind, ind, text in blocks(md):
        if kind == "gap":
            y -= line_h * 0.45
            continue
        if kind == "heading":
            put(M, text, FontProperties(family="DejaVu Sans", size=10.5, weight="bold"))
            y -= 10.5 * LEAD / 72 + line_h * 0.2
            continue
        x0 = M + 0.18 * ind
        if kind == "bullet":
            put(x0 + 0.05, "•", fp[False])
            x0 += 0.2
        elif kind == "number":
            num, text = text.split(" ", 1)
            put(x0 + 0.02, num, fp[False])
            x0 += 0.22
        x = x0
        for w, b in words(text):
            wd = width(w, fp[b])
            if x + wd > right and x > x0:
                y -= line_h
                x = x0
            put(x, w, fp[b])
            x += wd + space
        y -= line_h
    fits = y >= M
    fig.savefig(OUT.with_suffix(".pdf"))
    fig.savefig(OUT.with_suffix(".png"), dpi=110)
    print(f"{OUT.with_suffix('.pdf')}: 1 page; text ends {y:.2f} in above the page bottom (margin {M} in); "
          f"{'fits' if fits else 'DOES NOT FIT'}")
    sys.exit(0 if fits else 1)


if __name__ == "__main__":
    main()
