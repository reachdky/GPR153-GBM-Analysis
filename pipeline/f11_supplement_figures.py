#!/usr/bin/env python3
"""
f11_supplement_figures.py
Supporting-information figures (PREREGISTRATION.md, Section 14). Reads result CSVs only. A figure whose results
are absent (an analysis dropped under Section 13.2) is skipped and listed.

  S1 Fig   mean log2(TPM + 1) of all non-constant genes, eligibility filter marked   results/screen/
  S2 Fig   PCA, k-means Cluster 1-3, batch shown                                       results/pca/
  S3 Fig   Kaplan-Meier, within-batch median split                                     results/survival/km_*
  S4 Fig   cross-batch reproducibility; separate post hoc panel for GPR153 and MGMT    results/screen/
  S5 Fig   scaled Schoenfeld residuals, Model 3                                        results/survival/
  S6 Fig   CAMERA, targeted and top sets                                               results/enrichment/camera_*
  S7 Fig   DESeq2 volcano (apeglm log2FC)                                              results/de/
  S8 Fig   survival forest: C by feature set; permutation importance                   results/forest/
Outputs: revision/figures/S<n>_Fig.tif, revision/source_data/S<n>_Fig.csv, revision/figure_script_map.csv

Run from the repository root:
    python pipeline/f11_supplement_figures.py
"""

import numpy as np
import pandas as pd

import common
from f10_main_figures import fig_km
from plotstyle import INK, INK2, SERIES, fmt_p, legend_below, plt, save, short

R = common.RES
ME = "pipeline/f11_supplement_figures.py"


def s1():
    gl = pd.read_csv(R / "screen/screen_gene_level.csv")
    u = pd.read_csv(R / "screen/screen_gene_universe.csv").set_index("step")["n"]
    fig, ax = plt.subplots(figsize=(4.5, 3))
    bins = np.linspace(0, gl.mean_log2tpm.max(), 60)
    ax.hist(gl.loc[~gl.eligible_20, "mean_log2tpm"], bins=bins, color=SERIES[1],
            label=f"not eligible ({int(u['non-constant in primary cohort']) - int(u['eligible, 20% filter (primary)'])})")
    ax.hist(gl.loc[gl.eligible_20, "mean_log2tpm"], bins=bins, color=SERIES[0],
            label=f"eligible, TPM >= {common.PARAMS['eligibility_tpm']:g} in >= "
                  f"{common.PARAMS['eligibility_fraction_primary']:.0%} ({u['eligible, 20% filter (primary)']})")
    ax.set_yscale("log")
    ax.set_xlabel(f"Mean log2(TPM + 1), primary cohort (n {int(gl.n.iloc[0])})")
    ax.set_ylabel("Genes")
    ax.legend(fontsize=6)
    save(fig, "S1 Fig", ME, [R / "screen/screen_gene_level.csv", R / "screen/screen_gene_universe.csv"],
         {"A": gl[["gene_id", "symbol", "mean_log2tpm", "frac_tpm_ge_1", "eligible_20"]]})


def s2():
    sc = pd.read_csv(R / "pca/pca_scores.csv")
    var = pd.read_csv(R / "pca/pca_variance.csv").set_index("component")
    fig, ax = plt.subplots(figsize=(4.5, 3.6))
    marks = dict(zip(sorted(sc.batch.unique()), ["o", "^"]))
    for k, (cl, g) in enumerate(sc.groupby("cluster")):
        for b, gb in g.groupby("batch"):
            ax.scatter(gb.PC1, gb.PC2, s=12, color=SERIES[k], marker=marks[b], edgecolor="white", lw=0.4,
                       label=f"{cl}, {b} (n {len(gb)})")
    ax.set_xlabel(f"PC1 ({var.loc['PC1', 'variance_explained']:.1%})")
    ax.set_ylabel(f"PC2 ({var.loc['PC2', 'variance_explained']:.1%})")
    ax.set_title(f"n {int(var.n.iloc[0])}, {int(var.genes.iloc[0])} most variable eligible genes", loc="left")
    legend_below(ax, ncol=2, y=-0.16)
    save(fig, "S2 Fig", ME, [R / "pca/pca_scores.csv", R / "pca/pca_variance.csv"], {"A": sc, "variance": var.reset_index()})


def s4():
    cb = pd.read_csv(R / "screen/crossbatch_reproducibility.csv")
    ph = pd.read_csv(R / "screen/posthoc_crossbatch_gpr153_mgmt.csv")
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3), gridspec_kw={"width_ratios": [1.4, 1]})
    ax[0].hist(cb.crossbatch_spearman, bins=60, color=SERIES[0])
    ax[0].set_xlabel("Cross-batch Spearman (BI vs GDC-2024 aliquot)")
    ax[0].set_ylabel("Genes")
    ax[0].set_title(f"A  {len(cb)} genes, {int(cb.n_patients.iloc[0])} patients (GPR153, MGMT excluded)", loc="left")
    ax[1].hist(cb.crossbatch_spearman, bins=60, color=GRID_C)
    for k, r in enumerate(ph.itertuples()):
        ax[1].axvline(r.crossbatch_spearman, color=SERIES[k + 1], lw=1.5,
                      label=f"{r.gene} {r.crossbatch_spearman:.2f} ({r.percentile_among_genes:.0f}th percentile)")
    ax[1].set_title("B  Post hoc: GPR153 and MGMT", loc="left")
    ax[1].set_xlabel("Cross-batch Spearman")
    legend_below(ax[1], ncol=1, y=-0.2)
    fig.tight_layout()
    save(fig, "S4 Fig", ME, [R / "screen/crossbatch_reproducibility.csv", R / "screen/posthoc_crossbatch_gpr153_mgmt.csv"],
         {"A": cb, "B_post_hoc": ph})


GRID_C = "#c3c2b7"


def s5():
    r = pd.read_csv(R / "survival/cox_model3_schoenfeld.csv")
    ph = pd.read_csv(R / "survival/cox_ph_tests.csv")
    ph = ph[(ph.cohort == "primary") & (ph.model == "Model 3") & (ph.variant == "main")].set_index("term")
    covs = [c for c in r.columns if c != "time"]
    fig, ax = plt.subplots(1, len(covs), figsize=(2.6 * len(covs), 2.6))
    for a, c in zip(np.atleast_1d(ax), covs):
        a.scatter(r.time / common.DAYS_PER_MONTH, r[c], s=6, color=SERIES[0], alpha=0.7)
        a.axhline(0, color=INK2, lw=0.8, ls="--")
        a.set_title(f"{c}: cox.zph p {fmt_p(ph.loc[c, 'p'])}", loc="left")
        a.set_xlabel("Months")
        a.set_ylabel("Scaled Schoenfeld residual")
    fig.suptitle(f"Model 3, global p {fmt_p(ph.loc['GLOBAL', 'p'])}", x=0.01, ha="left", fontsize=8)
    fig.tight_layout()
    save(fig, "S5 Fig", ME, [R / "survival/cox_model3_schoenfeld.csv", R / "survival/cox_ph_tests.csv"],
         {"A": r, "tests": ph.reset_index()})


def s6():
    files = sorted((R / "enrichment").glob("camera_*_igc_estimated.csv"))
    if not files:
        return "S6 Fig"
    c = pd.concat([pd.read_csv(f) for f in files])
    show = pd.concat([c.nsmallest(12, "PValue"), c[c.term.isin(common.PARAMS["targeted_sets"])]]).drop_duplicates("term")
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    y = np.arange(len(show))[::-1]
    cols = [SERIES[1] if t in common.PARAMS["targeted_sets"] else SERIES[0] for t in show.term]
    ax.barh(y, -np.log10(show.FDR.clip(lower=1e-300)), color=cols, height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{short(t, 45)} ({d})" for t, d in zip(show.term, show.Direction)], fontsize=5.5)
    ax.set_xlabel("-log10 FDR (CAMERA, inter-gene correlation estimated)")
    qmax = common.PARAMS["ora"]["q"]
    n_sig = int((c.FDR < qmax).sum())
    ax.set_title(f"{n_sig} of {len(c)} sets with FDR < {qmax:g}" + ("; smallest p and targeted sets shown" if not n_sig else ""),
                 loc="left")
    save(fig, "S6 Fig", ME, files, {"A": show})


def s7():
    f = R / "de/deseq2_primary.csv"
    if not f.exists() or not (R / "de/de_summary.csv").exists():
        return "S7 Fig"
    d = pd.read_csv(f).dropna(subset=["padj"])
    s = pd.read_csv(R / "de/de_summary.csv").set_index("comparison").loc["primary"]
    P = common.PARAMS["de"]
    sig = (d.padj < P["padj"]) & (d.log2FoldChange_apeglm.abs() >= P["abs_lfc"])
    fig, ax = plt.subplots(figsize=(4.5, 3.6))
    ax.scatter(d.loc[~sig, "log2FoldChange_apeglm"], -np.log10(d.loc[~sig, "padj"].clip(lower=1e-300)), s=2, color=GRID_C)
    ax.scatter(d.loc[sig, "log2FoldChange_apeglm"], -np.log10(d.loc[sig, "padj"].clip(lower=1e-300)), s=2, color=SERIES[0],
               label=f"significant ({int(sig.sum())})")
    g = d[d.gene_id == common.GPR153]
    if len(g):
        ax.scatter(g.log2FoldChange_apeglm, -np.log10(g.padj.clip(lower=1e-300)), s=24, color=SERIES[1], edgecolor="white",
                   lw=1, zorder=3, label=f"GPR153 (rank {s.gpr153_rank_padj:.0f} by padj)")
        gx, gy = g.log2FoldChange_apeglm.iloc[0], -np.log10(g.padj.clip(lower=1e-300).iloc[0])
        ax.annotate(f"GPR153\nlog2FC {gx:.2f}", (gx, gy), xytext=(-40, 40), textcoords="offset points", fontsize=6,
                    ha="right", arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6),
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))
    ax.set_xlabel("log2 fold change, apeglm (tumor vs normal)")
    ax.set_ylabel("-log10 adjusted p")
    ax.set_title(f"{int(s.tumors)} tumors vs {int(s.normals)} normals", loc="left")
    legend_below(ax, ncol=2, y=-0.16)
    save(fig, "S7 Fig", ME, [f, R / "de/de_summary.csv"], {"A": d})


def s8():
    f = R / "forest/forest_cv_folds.csv"
    if not f.exists():
        return "S8 Fig"
    cv = pd.read_csv(f)
    imp = pd.read_csv(R / "forest/forest_permutation_importance.csv")
    an = list(dict.fromkeys(cv.analysis))
    fig, ax = plt.subplots(1, 2, figsize=(7.5, 3.2))
    pos = 0
    ticks, labs = [], []
    for k, a in enumerate(an):
        for j, s in enumerate("ABC"):
            v = cv[(cv.analysis == a) & (cv.set == s)].c
            ax[0].boxplot(v, positions=[pos], widths=0.6, patch_artist=True,
                          boxprops=dict(facecolor=SERIES[j], color=INK2), medianprops=dict(color=INK), flierprops=dict(markersize=2))
            ticks.append(pos)
            labs.append(s)
            pos += 1
        ax[0].text(pos - 2, -0.13, a.replace("primary cohort, ", ""), fontsize=6, ha="center", va="top",
                   transform=ax[0].get_xaxis_transform())
        pos += 1
    ax[0].set_xticks(ticks)
    ax[0].set_xticklabels(labs)
    ax[0].set_ylabel("Harrell's C, held-out folds")
    ax[0].set_title("A  Feature sets A-C", loc="left")
    pi = imp[imp.analysis == an[0]].sort_values("mean")
    ax[1].barh(np.arange(len(pi)), pi["mean"], xerr=pi["std"], color=SERIES[0], height=0.6, error_kw=dict(lw=0.8, ecolor=INK2))
    ax[1].set_yticks(np.arange(len(pi)))
    ax[1].set_yticklabels(pi.feature, fontsize=6)
    ax[1].set_xlabel("Permutation importance (decrease in C)")
    ax[1].set_title(f"B  Set C, {an[0]}", loc="left")
    fig.tight_layout()
    save(fig, "S8 Fig", ME, [f, R / "forest/forest_permutation_importance.csv"], {"A": cv, "B": imp})


def main():
    skipped = [x for x in (s1(), s2(), fig_km("S3 Fig", "within-batch median", ME), s4(), s5(), s6(), s7(), s8())
               if isinstance(x, str)]
    print("supplementary figures written to", common.FIG, "; skipped (no results):", ", ".join(skipped) or "none")


if __name__ == "__main__":
    main()
