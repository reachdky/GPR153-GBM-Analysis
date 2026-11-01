#!/usr/bin/env python3
"""
f10_main_figures.py
Main-text figures (PREREGISTRATION.md, Section 14): Fig 1-6. Reads result CSVs only; every number shown in a
title, label, legend or annotation is taken from those files.

  Fig 1  cohort flow (TCGA and external cohorts)                 results/cohort/
  Fig 2  MGMT co-expression screen                               results/screen/
  Fig 3  Kaplan-Meier, pooled median split, with at-risk tables  results/survival/km_*
  Fig 4  GPR153 HR per SD across Cox Models 1-7                  results/survival/cox_models.csv
  Fig 5  external cohorts: co-expression and survival            results/external/
  Fig 6  enrichment: GSEA and ORA                                results/enrichment/
Outputs: revision/figures/Fig<n>.tif, revision/source_data/Fig<n>.csv, revision/figure_script_map.csv
Format and size limits: pipeline/plotstyle.py.

Run from the repository root:
    python pipeline/f10_main_figures.py
"""

import re

import numpy as np
import pandas as pd

import common
from plotstyle import INK, INK2, SERIES, fmt_p, fmt_perm_q, plt, save, short

R = common.RES
ME = "pipeline/f10_main_figures.py"


def fig1():
    tf = pd.read_csv(R / "cohort/tcga_cohort_flow.csv")
    ef = pd.read_csv(R / "cohort/external_cohort_flow.csv")
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 6.4), gridspec_kw={"width_ratios": [1, 1.1], "wspace": 0.05})
    for ax in (a, b):
        ax.axis("off")
    y = 1.0
    for r in tf.itertuples():
        txt = f"{r.step.strip()}\n{int(r.n_samples)} samples" + (f", {int(r.n_patients)} patients" if pd.notna(r.n_patients) else "")
        indent = 0.06 if r.step.startswith("  ") else 0.0
        a.text(indent, y, txt, color=INK, va="top",
               bbox=dict(boxstyle="round,pad=0.3", fc="#f0efec", ec=INK2, lw=0.5))
        y -= 1.0 / (len(tf) + 1)
    a.set_title("A  TCGA-GBM", loc="left")
    y = 1.0
    for cohort, g in ef.groupby("cohort", sort=False):
        inc = int(g.loc[g.reason == "included", "n"].sum())
        b.text(0, y, f"{cohort}: {int(g.n.sum())} samples, {inc} included", fontweight="bold", va="top")
        y -= 0.034
        for r in g[g.reason != "included"].itertuples():
            b.text(0.04, y, f"{re.sub(r'\bnan\b', 'missing', str(r.reason))}: {int(r.n)}", color=INK2, va="top")
            y -= 0.032
        y -= 0.015
    b.set_title("B  External co-expression cohorts", loc="left")
    save(fig, "Fig 1", ME, [R / "cohort/tcga_cohort_flow.csv", R / "cohort/external_cohort_flow.csv"],
         {"A": tf, "B": ef})


def fig2():
    gl = pd.read_csv(R / "screen/screen_gene_level.csv")
    s = pd.read_csv(R / "screen/screen_sensitivity.csv").set_index("run")
    v = pd.read_csv(R / "screen/screen_gpr153_mgmt_values.csv")
    bs = pd.read_csv(R / "screen/bootstrap_rank_draws.csv")
    bsum = pd.read_csv(R / "screen/bootstrap_rank_summary.csv").iloc[0]
    pn = pd.read_csv(R / "screen/permutation_null_draws.csv")
    psum = pd.read_csv(R / "screen/permutation_null_summary.csv").iloc[0]
    pr, po = s.loc["primary: batch-adjusted"], s.loc["pooled (unadjusted)"]
    fig, axd = plt.subplot_mosaic("AB;CD;EE", figsize=(7.0, 8.0))
    ax = [axd[k] for k in "ABCDE"]
    frames = {}
    for i, (tag, row, lab) in enumerate((("adj", pr, "Batch-adjusted"), ("pooled", po, "Pooled"))):
        e = gl[gl.eligible_20 & gl[f"rank_{tag}"].notna()].sort_values(f"rank_{tag}")
        ax[i].plot(e[f"rank_{tag}"], e[f"rho_{tag}"], color=SERIES[0], lw=1)
        g = gl[gl.symbol == "GPR153"].iloc[0]
        ax[i].scatter([g[f"rank_{tag}"]], [g[f"rho_{tag}"]], s=24, color=SERIES[1], zorder=3, edgecolor="white", lw=1,
                      label="GPR153")
        ax[i].set_xscale("log")
        ax[i].set_xlabel("Rank among eligible genes")
        ax[i].set_ylabel("Spearman rho with MGMT")
        ax[i].set_title(f"{'AB'[i]}  {lab}, n {int(row.n)}, {int(row.genes_tested)} genes\n"
                        f"GPR153 rank {g[f'rank_{tag}']:.0f}, rho {g[f'rho_{tag}']:.3f}, q {fmt_p(g[f'q_{tag}'])}",
                        loc="left")
        ax[i].legend(loc="lower left")
        frames["AB"[i]] = e[["gene_id", "symbol", f"rank_{tag}", f"rho_{tag}", f"q_{tag}"]]
    for k, (b, g) in enumerate(v.groupby("batch")):
        ax[2].scatter(g.gpr153_log2tpm, g.mgmt_log2tpm, s=9, color=SERIES[k], alpha=0.8, edgecolor="white", lw=0.4,
                      label=f"{b} (n {len(g)})")
    ax[2].set_xlabel("GPR153 log2(TPM + 1)")
    ax[2].set_ylabel("MGMT log2(TPM + 1)")
    ax[2].set_title(f"C  GPR153 vs MGMT, unadjusted values\nadjusted rho {pr.gpr153_rho:.3f}; "
                    f"pooled rho {po.gpr153_rho:.3f}", loc="left")
    lo_y, hi_y = v.mgmt_log2tpm.min(), v.mgmt_log2tpm.max()
    ax[2].set_ylim(lo_y - 0.05 * (hi_y - lo_y), hi_y + 0.3 * (hi_y - lo_y))
    ax[2].legend(loc="upper left", ncol=2, handletextpad=0.2, columnspacing=1.0)
    frames["C"] = v
    ax[3].hist(bs.gpr153_rank, bins=40, color=SERIES[0], edgecolor="white", lw=0.5)
    ax[3].set_xlabel("GPR153 rank in bootstrap resample")
    ax[3].set_ylabel("Resamples")
    ax[3].set_title(f"D  {int(bsum.resamples)} bootstrap resamples\nrank 1 in {bsum.prop_rank1:.1%}, "
                    f"top 10 in {bsum.prop_top10:.1%}", loc="left")
    frames["D"] = bs[["resample", "gpr153_rank", "gpr153_rho"]]
    ax[4].hist(pn.max_rho, bins=50, color=SERIES[2])
    ax[4].axvline(psum.observed_max_rho, color=INK, lw=1, ls="--", label=f"observed maximum {psum.observed_max_rho:.3f}")
    ax[4].axvline(psum.observed_gpr153_rho, color=SERIES[1], lw=1.5, label=f"GPR153 {psum.observed_gpr153_rho:.3f}")
    ax[4].set_xlabel("Maximum rho under permutation")
    ax[4].set_ylabel("Permutations")
    ax[4].set_title(f"E  {int(psum.permutations)} within-batch permutations", loc="left")
    ax[4].legend(loc="upper center")
    frames["E"] = pn
    fig.tight_layout()
    save(fig, "Fig 2", ME, [R / "screen"], frames)


def km_panel(ax, curves, groups, risk, tests, gene, split, title):
    for k, grp in enumerate(("high", "low")):
        c = curves[(curves.gene == gene) & (curves.split == split) & (curves.group == grp)]
        g = groups[(groups.gene == gene) & (groups.split == split) & (groups.group == grp)].iloc[0]
        t = c.time_days / common.DAYS_PER_MONTH
        ax.step(t, c.survival, where="post", color=SERIES[k], label=f"{grp} (n {g.n}, deaths {g.deaths})")
        ax.fill_between(t, c.ci_low, c.ci_high, step="post", color=SERIES[k], alpha=0.12, lw=0)
    p = tests[(tests.gene == gene) & (tests.split == split)].iloc[0].p
    ax.set_title(f"{title}  {gene}, log-rank p {fmt_p(p)}", loc="left")
    ax.set_xlabel("Months", labelpad=2)
    ax.set_ylabel("Overall survival")
    ax.legend(loc="upper right")
    ax.annotate("At risk", (0, 0), xytext=(-16, -40), textcoords="offset points", ha="right", color=INK,
                fontweight="bold", annotation_clip=False)
    rk = risk[(risk.gene == gene) & (risk.split == split)]
    months = sorted(rk.month.unique())
    ax.set_xlim(0, max(months))
    ax.set_xticks(months)
    for k, grp in enumerate(("high", "low")):
        r = rk[rk.group == grp].set_index("month").at_risk
        ax.annotate(grp, (0, 0), xytext=(-16, -52 - 11 * k), textcoords="offset points", ha="right", color=SERIES[k],
                    annotation_clip=False)
        for mo in months:
            ax.annotate(f"{r[mo]}", (mo, 0), xytext=(0, -52 - 11 * k), textcoords="offset points", ha="center",
                        color=SERIES[k], annotation_clip=False)


def fig_km(item, split, script=ME):
    curves = pd.read_csv(R / "survival/km_curves.csv")
    groups = pd.read_csv(R / "survival/km_groups.csv")
    risk = pd.read_csv(R / "survival/km_at_risk.csv")
    tests = pd.read_csv(R / "survival/km_logrank.csv")
    fig, ax = plt.subplots(2, 1, figsize=(6.0, 7.4))
    for a, gene, t in zip(ax, ("GPR153", "MGMT"), "AB"):
        km_panel(a, curves, groups, risk, tests, gene, split, t)
    fig.suptitle(f"High vs low: {split} split; shading 95% CI", x=0.01, ha="left")
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.75, bottom=0.14, left=0.14)
    save(fig, item, script, [R / "survival/km_curves.csv", R / "survival/km_groups.csv", R / "survival/km_at_risk.csv",
                             R / "survival/km_logrank.csv"],
         {"curves": curves[curves.split == split], "groups": groups[groups.split == split],
          "at_risk": risk[risk.split == split], "logrank": tests[tests.split == split]})


def forest(ax, labels, est, lo, hi, highlight=None, ref=1.0, xlabel="", log=True):
    y = np.arange(len(labels))[::-1]
    for i in range(len(labels)):
        c = SERIES[1] if highlight is not None and i == highlight else SERIES[0]
        ax.plot([lo[i], hi[i]], [y[i]] * 2, color=c, lw=1.5)
        ax.scatter([est[i]], [y[i]], color=c, s=24, zorder=3, edgecolor="white", lw=1)
    ax.axvline(ref, color=INK2, lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    if log:
        from matplotlib.ticker import FixedLocator, FormatStrFormatter, NullFormatter, NullLocator
        ax.set_xscale("log")
        lo_x, hi_x = min(np.min(lo), ref), max(np.max(hi), ref)
        for cand in ([0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0], [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]):
            cand = np.array(cand)
            ticks = cand[(cand >= lo_x * 0.98) & (cand <= hi_x * 1.02)]
            if len(ticks) >= 3:
                break
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(xlabel)
    ax.grid(axis="y", visible=False)


def fig4():
    c = pd.read_csv(R / "survival/cox_models.csv")
    c = c[(c.cohort == "primary") & (c.variant == "main") & (c.sample_set == "complete cases") & (c.covariate == "gpr153")]
    c = c.sort_values("model")
    labels = [f"{r.model} (n {r.n}, events {r.events})" for r in c.itertuples()]
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    forest(ax, labels, c.hr.values, c.hr_low.values, c.hr_high.values,
           highlight=list(c.model).index("Model 3"), xlabel="GPR153 hazard ratio per SD (95% CI)")
    ax.set_title("Cox Models 1-7, IDH-wildtype primary cohort\nModel 3 (primary) highlighted", loc="left")
    save(fig, "Fig 4", ME, [R / "survival/cox_models.csv"], {"A": c})


def fig5():
    ce = pd.read_csv(R / "external/external_coexpression.csv")
    em = pd.read_csv(R / "external/external_meta.csv")
    sr = pd.read_csv(R / "external/survival_replication_cox.csv")
    sm = pd.read_csv(R / "external/survival_replication_meta.csv")
    tc = pd.read_csv(R / "screen/screen_sensitivity.csv").set_index("run").loc["primary: batch-adjusted"]
    lo_t, hi_t, _, _ = common.bw_ci(tc.gpr153_rho, tc.n, 1)
    fig, ax = plt.subplots(2, 1, figsize=(6.4, 7.4))
    mm = em[em.analysis == "external cohorts + TCGA"].iloc[0]
    lab = [f"{r.cohort} (n {r.n})" for r in ce.itertuples()] + [f"TCGA-GBM, batch-adjusted (n {int(tc.n)})",
                                                                 f"Random effects, {int(mm.k)} cohorts incl. TCGA, "
                                                                 f"I2 {mm.I2:.0%}"]
    est = np.r_[ce.spearman_rho, tc.gpr153_rho, mm.mu_rho]
    lo = np.r_[ce.spearman_ci_low, lo_t, mm.ci_low_rho]
    hi = np.r_[ce.spearman_ci_high, hi_t, mm.ci_high_rho]
    forest(ax[0], lab, est, lo, hi, highlight=len(lab) - 1, ref=0, xlabel="Spearman rho, GPR153 vs MGMT (95% CI)", log=False)
    ax[0].set_title("A  Co-expression", loc="left")
    b = sr[(sr.model == "(b) GPR153 + age + sex") & ~((sr.cohort == "CPTAC-3 GBM") & (sr.subset == "all"))]
    mb = sm[(sm.model == "(b) GPR153 + age + sex") & sm.with_tcga].iloc[0]
    t3 = pd.read_csv(R / "survival/cox_models.csv")
    t3 = t3[(t3.cohort == "primary") & (t3.model == "Model 3") & (t3.variant == "main") & (t3.sample_set == "complete cases")
            & (t3.covariate == "gpr153")].iloc[0]
    lab = [f"{r.cohort}{', ' + r.subset if r.subset != 'all' else ''} (n {r.n})" for r in b.itertuples()] + [
        f"TCGA-GBM, Model 3 (n {t3.n})", f"Random effects, {int(mb.k)} cohorts incl. TCGA, I2 {mb.I2:.0%}"]
    forest(ax[1], lab, np.r_[b.hr, t3.hr, mb.mu_hr], np.r_[b.hr_low, t3.hr_low, mb.ci_low_hr],
           np.r_[b.hr_high, t3.hr_high, mb.ci_high_hr], highlight=len(lab) - 1, xlabel="GPR153 HR per SD (95% CI)")
    ax[1].set_title("B  Survival, GPR153 + age + sex (model b)", loc="left")
    fig.tight_layout()
    save(fig, "Fig 5", ME, [R / "external"], {"A": ce, "A_meta": em, "B": b, "B_meta": sm})


def fig6():
    G = common.PARAMS["gsea"]["libraries"]
    perm = common.PARAMS["gsea"]["permutations"]
    g = pd.concat([pd.read_csv(R / f"enrichment/gsea_{lib}.csv").assign(library=lib) for lib in G])
    g["FDR q-val"] = g["FDR q-val"].astype(float)
    top = g.sort_values("FDR q-val").groupby("library").head(8)
    tg = g[g.Term.isin(common.PARAMS["targeted_sets"])]
    show = pd.concat([top, tg]).drop_duplicates(["library", "Term"])
    o = pd.concat([pd.read_csv(R / f"enrichment/ora_{lib}_positive.csv").assign(library=lib)
                   for lib in common.PARAMS["ora"]["libraries"] if (R / f"enrichment/ora_{lib}_positive.csv").exists()])
    otop = o.sort_values("q_bh").head(12)
    fig, ax = plt.subplots(2, 1, figsize=(4.4, 8.0), gridspec_kw={"hspace": 0.3, "height_ratios": [len(show), len(otop)]})
    y = np.arange(len(show))[::-1]
    cols = [SERIES[1] if t in common.PARAMS["targeted_sets"] else SERIES[0] for t in show.Term]
    ax[0].barh(y, show.NES.astype(float), color=cols, height=0.7)
    ax[0].set_yticks(y)
    ax[0].set_yticklabels([f"{short(t, 36)} (FDR {fmt_perm_q(q, perm)})" for t, q in zip(show.Term, show["FDR q-val"])])
    ax[0].set_xlabel("Normalized enrichment score")
    ax[0].set_title("A  GSEA preranked\n(targeted sets in orange)", loc="left")
    y = np.arange(len(otop))[::-1]
    ax[1].barh(y, -np.log10(otop.q_bh.clip(lower=1e-300)), color=SERIES[2], height=0.7)
    ax[1].set_yticks(y)
    tag = {"Reactome_2022": "Reactome", "KEGG_2021_Human": "KEGG", "GO_Biological_Process_2023": "GO BP"}
    ax[1].set_yticklabels([f"{short(t, 36)} ({tag.get(lib, lib)})" for t, lib in zip(otop.Term, otop.library)])
    ax[1].set_xlabel("-log10 BH-adjusted p")
    qmax = common.PARAMS["ora"]["q"]
    n_sig = int((o.q_bh < qmax).sum())
    ax[1].set_title("B  Over-representation, positive list\n" + (f"{n_sig} sets with adjusted p < {qmax:g}" if n_sig
                    else f"no set with adjusted p < {qmax:g}; {len(otop)} smallest shown"), loc="left")
    save(fig, "Fig 6", ME, [R / "enrichment"], {"A": show, "B": otop})


def main():
    fig1()
    fig2()
    fig_km("Fig 3", "pooled median")
    fig4()
    fig5()
    fig6()
    print("Fig 1-6 written to", common.FIG)


if __name__ == "__main__":
    main()
