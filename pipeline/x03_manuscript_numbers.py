#!/usr/bin/env python3
"""
x03_manuscript_numbers.py
Write revision/numbers_for_manuscript.md: the values from the final result files that the revised manuscript
may cite, each with its source file, row key and column. Reads result files only; values are printed to four
significant digits (p-values and q-values to three), the full values are in the source files.

Inputs
  revision/results/**   (files named in each row)
Outputs
  revision/numbers_for_manuscript.md

Run from the repository root (after run_all.sh --analysis):
    python pipeline/x03_manuscript_numbers.py
"""

import pandas as pd

import common

R = common.RES
OUT = common.BASE / "numbers_for_manuscript.md"
ROWS = []


def fmt(x, col=""):
    if isinstance(x, str):
        return x
    if pd.isna(x):
        return "NA"
    if float(x).is_integer() and abs(x) < 1e7:
        return str(int(x))
    low = col.lower()
    digits = 3 if (low == "p" or low.startswith(("p_", "q", "p-")) or low.endswith(("_p", "_q", "p-val", "q-val"))
                   or low in ("pvalue", "padj", "fdr", "lr_p", "nb_p", "hk_p", "q_bh")) else 4
    return f"{x:.{digits}g}"


def row(path, key):
    d = pd.read_csv(R / path, low_memory=False)
    m = pd.Series(True, index=d.index)
    for k, v in key.items():
        m &= d[k].astype(str) == str(v)
    if m.sum() != 1:
        raise SystemExit(f"{path}: {key} matches {m.sum()} rows")
    return d[m].iloc[0]


def keystr(key):
    return "; ".join(f"{k} = {v}" for k, v in key.items()) if key else "(single row)"


def add(section, label, path, key, cols, value=None):
    """One output row. cols: one column, or (estimate, low, high) printed as 'est [low, high]'."""
    r = row(path, key)
    if value is None:
        if isinstance(cols, tuple):
            e, lo, hi = cols
            value = f"{fmt(r[e], e)} [{fmt(r[lo], lo)}, {fmt(r[hi], hi)}]"
        else:
            value = fmt(r[cols], cols)
    ROWS.append((section, label, value, f"results/{path}", keystr(key),
                 "; ".join(cols) if isinstance(cols, tuple) else cols))
    return r


def derived(section, label, value, path, rule, col):
    ROWS.append((section, label, value, f"results/{path}", rule, col))


def count(path, col, thr):
    d = pd.read_csv(R / path)
    return int((d[col] < thr).sum()), len(d)


def main():
    S = "1. Cohort flow"
    for step in ["GDC STAR - Counts files (TCGA-GBM)", "Primary tumor (sample type 01)", "One aliquot per patient",
                 "  of which BI batch", "  of which GDC2024 batch", "IDH-wildtype (primary cohort, screen)",
                 "IDH-mutant (sensitivity only)", "IDH unknown (sensitivity only)",
                 "Survival evaluable, all primary (sensitivity)", "Survival evaluable, IDH-wildtype (primary)"]:
        add(S, f"TCGA: {step.strip()} (aliquots)", "cohort/tcga_cohort_flow.csv", {"step": step}, "n_samples")
        if step == "GDC STAR - Counts files (TCGA-GBM)":
            add(S, f"TCGA: {step} (patients)", "cohort/tcga_cohort_flow.csv", {"step": step}, "n_patients")
    for c in ["TCGA-GBM", "CGGA mRNAseq_693", "CGGA mRNAseq_325", "GSE13041 GPL96", "Rembrandt (GSE108474)",
              "CPTAC-3 GBM"]:
        add(S, f"{c}: selected", "cohort/survival_audit.csv", {"cohort": c}, "n_selected")
        add(S, f"{c}: survival usable", "cohort/survival_audit.csv", {"cohort": c}, "n_usable")
        add(S, f"{c}: deaths", "cohort/survival_audit.csv", {"cohort": c}, "events")
    for c in ["GSE13041_GPL96", "GSE4290_GBM", "CGGA_693", "CGGA_325"]:
        add(S, f"{c}: included (co-expression)", "cohort/external_cohort_flow.csv", {"cohort": c, "reason": "included"},
            "n")

    S = "2. Screen (primary cohort, batch-adjusted partial Spearman, 20% filter)"
    g = {"gene_id": common.GPR153}
    add(S, "Genes tested", "screen/screen_sensitivity.csv", {"run": "primary: batch-adjusted"}, "genes_tested")
    add(S, "n (aliquots)", "screen/screen_gene_level.csv", g, "n")
    add(S, "GPR153 partial rho", "screen/screen_gene_level.csv", g, "rho_adj")
    add(S, "GPR153 p", "screen/screen_gene_level.csv", g, "p_adj")
    add(S, "GPR153 BH q", "screen/screen_gene_level.csv", g, "q_adj")
    add(S, "GPR153 rank (signed rho)", "screen/screen_gene_level.csv", g, "rank_adj")
    add(S, "GPR153 rank (|rho|)", "screen/screen_gene_level.csv", g, "abs_rank_adj")
    add(S, "Rank-1 gene", "decision/decision.csv", {}, "rank1_gene")
    add(S, "Rank-1 rho", "decision/decision.csv", {}, "rank1_rho")
    add(S, "GPR153 pooled (unadjusted) rho", "screen/screen_gene_level.csv", g, "rho_pooled")
    add(S, "GPR153 pooled q", "screen/screen_gene_level.csv", g, "q_pooled")
    add(S, "GPR153 pooled rank", "screen/screen_gene_level.csv", g, "rank_pooled")
    n05, n = count("screen/screen_gene_level.csv", "q_adj", 0.05)
    derived(S, "Genes with q_adj < 0.05 (of rows with a q_adj)", str(n05), "screen/screen_gene_level.csv",
            "count of rows with q_adj < 0.05", "q_adj")

    S = "3. Sensitivity analyses of the screen"
    for run in ["BI only (sign check)", "GDC2024 only (sign check)", "FPKM", "MGMT outliers removed (|z| > 3)",
                "GPR153 or MGMT outliers removed (single value)", "filter 10%", "filter 50%", "filter none",
                "IDH-not-mutant cohort", "sensitivity cohort", "excluding IDH source disagreement"]:
        k = {"run": run}
        add(S, f"{run}: n", "screen/screen_sensitivity.csv", k, "n")
        if run.endswith("(sign check)"):
            add(S, f"{run}: rho [95% CI]", "screen/screen_sensitivity.csv", k, ("gpr153_rho", "ci_low", "ci_high"))
        else:
            add(S, f"{run}: rho", "screen/screen_sensitivity.csv", k, "gpr153_rho")
        r = row("screen/screen_sensitivity.csv", k)
        if pd.notna(r["gpr153_q"]):
            add(S, f"{run}: q", "screen/screen_sensitivity.csv", k, "gpr153_q")
            add(S, f"{run}: rank", "screen/screen_sensitivity.csv", k, "gpr153_rank")
        else:
            add(S, f"{run}: p", "screen/screen_sensitivity.csv", k, "gpr153_p")
    for a in ["batch + estimate score", "batch + neuronal score"]:
        k = {"adjustment": a}
        add(S, f"Composition, {a}: rho [95% CI]", "composition/composition_adjusted.csv", k,
            ("gpr153_rho", "ci_low", "ci_high"))
        add(S, f"Composition, {a}: q (adjusted screen)", "composition/composition_adjusted.csv", k,
            "gpr153_q_adjusted_screen")
        add(S, f"Composition, {a}: rank (adjusted screen)", "composition/composition_adjusted.csv", k,
            "gpr153_rank_adjusted_screen")
    for gname in ["GPR153", "MGMT"]:
        k = {"gene": gname}
        add(S, f"Post hoc cross-batch reproducibility, {gname}: Spearman (n patients in both batches)",
            "screen/posthoc_crossbatch_gpr153_mgmt.csv", k, "crossbatch_spearman")
        add(S, f"Post hoc cross-batch reproducibility, {gname}: percentile among genes",
            "screen/posthoc_crossbatch_gpr153_mgmt.csv", k, "percentile_among_genes")
    add(S, "Post hoc cross-batch: patients in both batches", "screen/posthoc_crossbatch_gpr153_mgmt.csv",
        {"gene": "GPR153"}, "n_patients")

    S = "4. Pre-specified conditions and outcome branch"
    for c in range(1, 5):
        add(S, f"Condition {c}: result", "decision/decision_conditions.csv", {"condition": c}, "value")
        add(S, f"Condition {c}: met", "decision/decision_conditions.csv", {"condition": c}, "met")
    add(S, "Outcome branch", "decision/decision.csv", {}, "branch")
    add(S, "Failing conditions", "decision/decision.csv", {}, "failing_conditions")
    k = {}
    add(S, "Condition 4: n", "methylation/condition4.csv", k, "n")
    add(S, "Condition 4: methylation-adjusted partial rho [95% CI]", "methylation/condition4.csv", k,
        ("rho", "ci_low", "ci_high"))
    add(S, "Condition 4: p", "methylation/condition4.csv", k, "p")
    add(S, "Condition 4: covariates", "methylation/condition4.csv", k, "covariates")
    add(S, "STP27 vs Ceccarelli 2016 kappa", "methylation/methylation_concordance.csv", k, "kappa")
    add(S, "STP27 vs Ceccarelli 2016 overlap n", "methylation/methylation_concordance.csv", k,
        "overlap_with_ceccarelli")
    for a in ["bootstrap percentile 95% CI", "without platform (k = 2)", "binary STP27 status", "other measure"]:
        add(S, f"Condition 4 sensitivity, {a}: rho [95% CI]", "methylation/condition4_sensitivity.csv",
            {"analysis": a}, ("rho", "ci_low", "ci_high"))
        add(S, f"Condition 4 sensitivity, {a}: n", "methylation/condition4_sensitivity.csv", {"analysis": a}, "n")

    S = "5. Bootstrap rank stability (1000 resamples)"
    for c in ["resamples", "n", "genes", "prop_rank1", "prop_top10", "rank_median"]:
        add(S, c, "screen/bootstrap_rank_summary.csv", {}, c)
    r = row("screen/bootstrap_rank_summary.csv", {})
    add(S, "GPR153 rho, bootstrap 95% interval", "screen/bootstrap_rank_summary.csv", {}, "rho_ci_low",
        value=f"[{fmt(r.rho_ci_low)}, {fmt(r.rho_ci_high)}]")
    ROWS[-1] = ROWS[-1][:5] + ("rho_ci_low; rho_ci_high",)

    S = "6. Permutation null (10,000 permutations of GPR153 within batch)"
    for c in ["permutations", "genes", "n", "observed_max_rho", "observed_gpr153_rho", "observed_gpr153_rank",
              "prop_max_ge_observed_max", "prop_rank_value_ge_gpr153", "max_rho_p50", "max_rho_p95", "max_rho_p99"]:
        add(S, c, "screen/permutation_null_summary.csv", {}, c)

    S = "7. Cox Models 1-7 (primary cohort, complete cases, main specification)"
    cox = pd.read_csv(R / "survival/cox_models.csv")
    main_rows = cox[(cox.cohort == "primary") & (cox.variant == "main") & (cox.sample_set == "complete cases")]
    for _, r in main_rows.iterrows():
        k = {"cohort": "primary", "model": r.model, "variant": "main", "sample_set": "complete cases",
             "covariate": r.covariate}
        add(S, f"{r.model}, {r.covariate}: HR [95% CI]", "survival/cox_models.csv", k, ("hr", "hr_low", "hr_high"))
        add(S, f"{r.model}, {r.covariate}: p", "survival/cox_models.csv", k, "p")
        if r.covariate == "gpr153":
            add(S, f"{r.model}: n", "survival/cox_models.csv", k, "n")
            add(S, f"{r.model}: events", "survival/cox_models.csv", k, "events")
    lr = pd.read_csv(R / "survival/cox_lr_tests.csv")
    for t in lr[lr.cohort == "primary"].test:
        k = {"cohort": "primary", "test": t}
        add(S, f"Likelihood-ratio test, {t}: chi-square", "survival/cox_lr_tests.csv", k, "lr_stat")
        add(S, f"Likelihood-ratio test, {t}: p", "survival/cox_lr_tests.csv", k, "p")
        add(S, f"Likelihood-ratio test, {t}: n", "survival/cox_lr_tests.csv", k, "n")
    for m in ["Model 1", "Model 2", "Model 3"]:
        k = {"cohort": "primary", "model": m, "variant": "main", "sample_set": "Model 4 complete cases",
             "covariate": "gpr153"}
        add(S, f"{m} refitted on the Model 4 sample, gpr153: HR [95% CI]", "survival/cox_models.csv", k,
            ("hr", "hr_low", "hr_high"))
        add(S, f"{m} refitted on the Model 4 sample: n / events", "survival/cox_models.csv", k, "n",
            value=f"{fmt(row('survival/cox_models.csv', k)['n'])} / {fmt(row('survival/cox_models.csv', k)['events'])}")
        ROWS[-1] = ROWS[-1][:5] + ("n; events",)
    for m, v in [("Model 3", "stratified by batch"), ("Model 3", "batch covariate"),
                 ("Model 4", "+ temozolomide + radiotherapy (immortal-time bias)"),
                 ("Model 5", "+ temozolomide + radiotherapy (immortal-time bias)")]:
        k = {"cohort": "primary", "model": m, "variant": v, "sample_set": "complete cases", "covariate": "gpr153"}
        r = add(S, f"{m}, {v}, gpr153: HR [95% CI]", "survival/cox_models.csv", k, ("hr", "hr_low", "hr_high"))
        add(S, f"{m}, {v}, gpr153: p; n / events", "survival/cox_models.csv", k, "p",
            value=f"{fmt(r.p, 'p')}; {fmt(r.n)} / {fmt(r.events)}")
        ROWS[-1] = ROWS[-1][:5] + ("p; n; events",)
    for c in ["idh_not_mutant", "sensitivity"]:
        for m in ["Model 1", "Model 3", "Model 4", "Model 5", "Model 6", "Model 7"]:
            k = {"cohort": c, "model": m, "variant": "main", "sample_set": "complete cases", "covariate": "gpr153"}
            r = add(S, f"{c}, {m}, gpr153: HR [95% CI]", "survival/cox_models.csv", k, ("hr", "hr_low", "hr_high"))
            add(S, f"{c}, {m}, gpr153: p", "survival/cox_models.csv", k, "p")
            add(S, f"{c}, {m}: n / events", "survival/cox_models.csv", k, "n", value=f"{fmt(r.n)} / {fmt(r.events)}")
            ROWS[-1] = ROWS[-1][:5] + ("n; events",)
    add(S, "Model 3, restricted cubic spline for age vs linear: LR p", "survival/cox_age_spline.csv",
        {"model": "Model 3"}, "p")
    mi = pd.read_csv(R / "survival/cox_multiple_imputation.csv")
    for _, r in mi[mi.term == "gpr153"].iterrows():
        k = {"model": r.model, "term": "gpr153"}
        add(S, f"Multiple imputation (m = {r.m}), {r.model}, gpr153: HR [95% CI]",
            "survival/cox_multiple_imputation.csv", k, ("hr", "hr_low", "hr_high"))
        add(S, f"Multiple imputation, {r.model}, gpr153: p (n)", "survival/cox_multiple_imputation.csv", k, "p",
            value=f"{fmt(r.p, 'p')} (n = {fmt(r.n)})")
        ROWS[-1] = ROWS[-1][:5] + ("p; n",)
    k = {"analysis": "TCGA primary, Model 4 vs baseline"}
    for c in ["n", "events", "lr_p", "delta_c", "delta_c_optimism_corrected"]:
        add(S, f"Incremental value (Model 4 vs baseline): {c}", "survival/incremental_value.csv", k, c)
    for c in ["hr", "p", "genes", "percentile_by_p"]:
        add(S, f"Genome-wide univariable Cox, GPR153: {c}", "survival/cox_genomewide_gpr153.csv", {"gene": "GPR153"},
            c)

    S = "8. Kaplan-Meier and log-rank (primary cohort, survival evaluable)"
    for gname in ["GPR153", "MGMT"]:
        for sp in ["pooled median", "within-batch median"]:
            k = {"gene": gname, "split": sp}
            add(S, f"{gname}, {sp}: log-rank chi-square", "survival/km_logrank.csv", k, "chisq")
            add(S, f"{gname}, {sp}: log-rank p", "survival/km_logrank.csv", k, "p")
            for grp in ["high", "low"]:
                kg = {**k, "group": grp}
                r = row("survival/km_groups.csv", kg)
                add(S, f"{gname}, {sp}, {grp}: n / deaths", "survival/km_groups.csv", kg, "n",
                    value=f"{fmt(r.n)} / {fmt(r.deaths)}")
                ROWS[-1] = ROWS[-1][:5] + ("n; deaths",)
                add(S, f"{gname}, {sp}, {grp}: median OS, days [95% CI]", "survival/km_groups.csv", kg,
                    ("median_os_days", "median_os_ci_low", "median_os_ci_high"))
    add(S, "Median follow-up (reverse KM), days", "survival/km_logrank.csv",
        {"gene": "GPR153", "split": "pooled median"}, "median_followup_days_reverse_km")

    S = "9. Proportional-hazards diagnostics (R survival::cox.zph, transform km)"
    ph = pd.read_csv(R / "survival/cox_ph_tests.csv")
    for m in sorted(ph[(ph.cohort == "primary") & (ph.variant == "main")].model.unique()):
        for t in ["gpr153", "GLOBAL"]:
            k = {"cohort": "primary", "model": m, "variant": "main", "term": t}
            add(S, f"{m}, {t}: p", "survival/cox_ph_tests.csv", k, "p")

    S = "10. External cohorts: GPR153-MGMT co-expression"
    for c in ["GSE13041 GPL96 primary", "GSE4290 GBM", "CGGA mRNAseq_693", "CGGA mRNAseq_325"]:
        k = {"cohort": c}
        add(S, f"{c}: n", "external/external_coexpression.csv", k, "n")
        add(S, f"{c}: Spearman rho [95% CI]", "external/external_coexpression.csv", k,
            ("spearman_rho", "spearman_ci_low", "spearman_ci_high"))
        add(S, f"{c}: p", "external/external_coexpression.csv", k, "spearman_p")
        add(S, f"{c}: BH q", "external/external_coexpression.csv", k, "spearman_q_bh")
    for a in ["CGGA 693 + 325 pooled", "CGGA 693 + 325, adjusted for release (k = 1)"]:
        add(S, f"{a}: rho [95% CI] (n)", "external/cgga_pooled.csv", {"analysis": a},
            ("spearman_rho", "spearman_ci_low", "spearman_ci_high"))

    S = "11. External cohorts: survival replication"
    sr = pd.read_csv(R / "external/survival_replication_cox.csv")
    for _, r in sr.iterrows():
        k = {"cohort": r.cohort, "subset": r.subset, "model": r.model}
        add(S, f"{r.cohort} ({r.subset}), {r.model}: HR [95% CI]", "external/survival_replication_cox.csv", k,
            ("hr", "hr_low", "hr_high"))
        add(S, f"{r.cohort} ({r.subset}), {r.model}: p; n / events", "external/survival_replication_cox.csv", k,
            "p", value=f"{fmt(r.p, 'p')}; {fmt(r.n)} / {fmt(r.events)}")
        ROWS[-1] = ROWS[-1][:5] + ("p; n; events",)
    for a in ["CGGA mRNAseq_693: GPR153 + baseline vs baseline", "CGGA mRNAseq_325: GPR153 + baseline vs baseline"]:
        add(S, f"{a}: LR p", "external/survival_replication_incremental.csv", {"analysis": a}, "lr_p")
        add(S, f"{a}: optimism-corrected delta C", "external/survival_replication_incremental.csv",
            {"analysis": a}, "delta_c_optimism_corrected")

    S = "12. Meta-analyses (random effects, REML)"
    for a in ["external cohorts", "external cohorts + TCGA"]:
        k = {"analysis": a}
        add(S, f"Co-expression, {a}: pooled rho [95% CI]", "external/external_meta.csv", k,
            ("mu_rho", "ci_low_rho", "ci_high_rho"))
        add(S, f"Co-expression, {a}: p", "external/external_meta.csv", k, "p")
        add(S, f"Co-expression, {a}: Hartung-Knapp 95% CI (rho)", "external/external_meta.csv", k, "hk_p",
            value=f"[{fmt(row('external/external_meta.csv', k)['hk_ci_low_rho'])}, "
                  f"{fmt(row('external/external_meta.csv', k)['hk_ci_high_rho'])}]")
        ROWS[-1] = ROWS[-1][:5] + ("hk_ci_low_rho; hk_ci_high_rho",)
        add(S, f"Co-expression, {a}: I2", "external/external_meta.csv", k, "I2")
        add(S, f"Co-expression, {a}: k", "external/external_meta.csv", k, "k")
    for m in ["(a) GPR153", "(b) GPR153 + age + sex"]:
        for t in [False, True]:
            k = {"model": m, "with_tcga": t}
            lab = f"Survival {m}, {'with' if t else 'without'} TCGA"
            add(S, f"{lab}: pooled HR [95% CI]", "external/survival_replication_meta.csv", k,
                ("mu_hr", "ci_low_hr", "ci_high_hr"))
            add(S, f"{lab}: p", "external/survival_replication_meta.csv", k, "p")
            add(S, f"{lab}: Hartung-Knapp p", "external/survival_replication_meta.csv", k, "hk_p")
            add(S, f"{lab}: I2; k", "external/survival_replication_meta.csv", k, "I2",
                value=f"{fmt(row('external/survival_replication_meta.csv', k)['I2'])}; "
                      f"{fmt(row('external/survival_replication_meta.csv', k)['k'])}")
            ROWS[-1] = ROWS[-1][:5] + ("I2; k",)

    S = "13. Enrichment"
    q = common.PARAMS["ora"]["q"]
    for lib in common.PARAMS["ora"]["libraries"]:
        for d in ["positive", "negative"]:
            p = f"enrichment/ora_{lib}_{d}.csv"
            n05, n = count(p, "q_bh", q)
            derived(S, f"ORA {lib}, GPR153-{d} list: sets with q_bh < {q} (of {n} tested)", str(n05), p,
                    f"count of rows with q_bh < {q}", "q_bh")
        p = f"enrichment/mgmt_ora_{lib}_positive.csv"
        n05, n = count(p, "q_bh", q)
        derived(S, f"ORA {lib}, MGMT-positive list: sets with q_bh < {q} (of {n} tested)", str(n05), p,
                f"count of rows with q_bh < {q}", "q_bh")
    for lib in common.PARAMS["gsea"]["libraries"]:
        p = f"enrichment/gsea_{lib}.csv"
        n05, n = count(p, "FDR q-val", q)
        derived(S, f"GSEA preranked {lib}: sets with FDR < {q} (of {n})", str(n05), p,
                f"count of rows with FDR q-val < {q}", "FDR q-val")
        for igc in ["estimated", "0.01"]:
            p = f"enrichment/camera_{lib}_igc_{igc}.csv"
            n05, n = count(p, "FDR", q)
            derived(S, f"camera {lib}, inter-gene correlation {igc}: sets with FDR < {q} (of {n})", str(n05), p,
                    f"count of rows with FDR < {q}", "FDR")
    for t in common.PARAMS["targeted_sets"]:
        k = {"Term": t}
        add(S, f"Targeted GSEA, {t}: NES", "enrichment/targeted_sets.csv", k, "NES")
        add(S, f"Targeted GSEA, {t}: FDR", "enrichment/targeted_sets.csv", k, "FDR q-val")
        for igc in ["estimated", "0.01"]:
            add(S, f"Targeted camera ({igc}), {t}: FDR", "enrichment/camera_targeted_sets.csv",
                {"term": t, "inter_gene_cor": igc}, "FDR")

    S = "14. Survival forest (5-fold CV x 10 repeats; C-index)"
    for a in ["primary cohort", "primary cohort, complete cases", "sensitivity cohort"]:
        for s in ["A", "B", "C"]:
            k = {"analysis": a, "set": s}
            r = add(S, f"{a}, set {s}: mean C [2.5%, 97.5%]", "forest/forest_c_summary.csv", k,
                    ("c_mean", "c_p025", "c_p975"))
        add(S, f"{a}: n / events", "forest/forest_c_summary.csv", {"analysis": a, "set": "C"}, "n",
            value=f"{fmt(r.n)} / {fmt(r.events)}")
        ROWS[-1] = ROWS[-1][:5] + ("n; events",)
        k = {"analysis": a}
        add(S, f"{a}: C(set C) - C(set B) [95% CI]", "forest/forest_paired_difference.csv", k,
            ("mean_diff_C_minus_B", "nb_ci_low", "nb_ci_high"))
        add(S, f"{a}: Nadeau-Bengio p", "forest/forest_paired_difference.csv", k, "nb_p")

    S = "15. Differential expression (DESeq2, tumour vs normal)"
    for c in ["primary", "sensitivity", "tcga_normals"]:
        k = {"comparison": c}
        for col in ["tumors", "normals", "genes_tested", "significant", "gpr153_log2fc_apeglm", "gpr153_padj",
                    "gpr153_rank_padj"]:
            add(S, f"{c}: {col}", "de/de_summary.csv", k, col)
    for col in ["spearman_log2fc", "significance_agreement", "gpr153_log2fc_py", "gpr153_padj_py"]:
        add(S, f"PyDESeq2 vs DESeq2: {col}", "de/pydeseq2_vs_deseq2.csv", {}, col)

    S = "16. PCA (2000 most variable genes, primary cohort)"
    for pc in ["PC1", "PC2", "PC3"]:
        add(S, f"{pc} variance explained", "pca/pca_variance.csv", {"component": pc}, "variance_explained")
    add(S, "n", "pca/pca_variance.csv", {"component": "PC1"}, "n")
    sc = pd.read_csv(R / "pca/pca_scores.csv")
    tab = pd.crosstab(sc.cluster, sc.batch)
    for cl, r in tab.iterrows():
        derived(S, f"k-means {cl}: aliquots by batch",
                "; ".join(f"{b} {int(v)}" for b, v in r.items()), "pca/pca_scores.csv",
                f"count of rows with cluster = {cl}, by batch", "cluster; batch")

    md = ["# Numbers for the manuscript", "",
          "Generated by `pipeline/x03_manuscript_numbers.py` from the final result files (`revision/results/`;",
          "final run identified in `revision/verification/final_run.md`). Values are rounded to four significant",
          "digits (p and q to three); the full values are in the files. A row key lists the column values that",
          "select one row; \"(single row)\" means the file has one data row. \"count of rows\" entries are derived",
          "from the named file and column.", ""]
    for sec in dict.fromkeys(r[0] for r in ROWS):
        md += [f"## {sec}", "", "| Quantity | Value | File | Row key | Column |", "|---|---|---|---|---|"]
        for r in ROWS:
            if r[0] == sec:
                md.append("| " + " | ".join(str(x).replace("|", "\\|") for x in r[1:]) + " |")
        md.append("")
    OUT.write_text("\n".join(md))
    print(f"{OUT}: {len(ROWS)} values")


if __name__ == "__main__":
    main()
