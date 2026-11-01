#!/usr/bin/env python3
"""
a90_methods_record.py
Methods record (PREREGISTRATION.md, Section 13): methods_parameters.json (fixed parameters, seeds, package versions,
per-analysis n) and methods_concordance.csv (each Methods statement -> plan section -> script -> output file, with
whether the output exists). Under the scope rule (13.2) only analyses with a row here enter the manuscript; a row
whose output is missing marks an analysis that is not finished. The independent reproductions (a91) are appended
with both values: per quantity, except the per-gene screen checks (top 100, 500 random genes), which are summarized
by their largest difference (all values in results/independent/independent_check.csv).

Inputs   result files under revision/results/ (all analysis scripts)
Outputs  revision/methods_parameters.json, revision/methods_concordance.csv

Run from the repository root:
    python pipeline/a90_methods_record.py
"""

import importlib.metadata as md
import platform
import subprocess

import pandas as pd

import common
from f12_tables import per_analysis_n

CONCORDANCE = [
    ("TCGA cohort construction, one aliquot per patient, IDH status", "3, 3.1", "pipeline/r01_tcga_cohort.py", "manifests/tcga_samples.csv"),
    ("Survival-data audit", "6.7", "pipeline/r07_survival_audit.py", "results/cohort/survival_audit.csv"),
    ("Expression matrices from GDC STAR - Counts, md5-checked", "4", "pipeline/a10_tcga_matrix.py", "results/expression/tcga_gene_id_handling.csv"),
    ("Genome-wide batch-adjusted and pooled Spearman screen", "5", "pipeline/a11_screen.py", "results/screen/screen_gene_level.csv"),
    ("Sensitivity runs of the screen", "5.1", "pipeline/a11_screen.py", "results/screen/screen_sensitivity.csv"),
    ("Cross-batch reproducibility", "5.1", "pipeline/a11_screen.py", "results/screen/crossbatch_summary.csv"),
    ("Bootstrap rank stability", "5.1", "pipeline/a12_bootstrap_rank.py", "results/screen/bootstrap_rank_summary.csv"),
    ("Composition adjustment (ESTIMATE, neuronal score)", "5.2", "pipeline/a14_composition.py", "results/composition/composition_adjusted.csv"),
    ("MGMT-STP27, kappa, condition 4", "3.2, 5.3", "pipeline/a15_methylation.py", "results/methylation/condition4.csv"),
    ("Kappa population and 2x2 table (STP27 vs Ceccarelli)", "3.2", "pipeline/a15_methylation.py", "results/methylation/methylation_kappa_2x2.csv"),
    ("Permutation null", "5.4", "pipeline/a13_permutation_null.py", "results/screen/permutation_null_summary.csv"),
    ("Decision rule and outcome branch", "5.5, 5.6", "pipeline/a16_decision.py", "results/decision/decision.csv"),
    ("Post hoc cross-batch diagnostic", "5.7", "pipeline/a17_posthoc_crossbatch.py", "results/screen/posthoc_crossbatch_gpr153_mgmt.csv"),
    ("PCA and k-means", "13", "pipeline/a18_pca.py", "results/pca/pca_scores.csv"),
    ("Cox Models 1-7, LR tests, PH tests, age spline", "6.1, 6.2", "pipeline/a20_cox.py", "results/survival/cox_models.csv"),
    ("Kaplan-Meier, log-rank, at-risk, reverse KM follow-up", "6.4, 6.5", "pipeline/a21_km.py", "results/survival/km_logrank.csv"),
    ("Genome-wide univariate Cox and incremental value", "6.6", "pipeline/a22_cox_genomewide.py", "results/survival/incremental_value.csv"),
    ("Multiple imputation (mice)", "6.3", "pipeline/a23_mice.py", "results/survival/cox_multiple_imputation.csv"),
    ("External co-expression, GSE13041 runs, U1 criterion, meta-analysis, power", "7", "pipeline/a30_external_coexpression.py", "results/external/external_coexpression.csv"),
    ("GSE13041 by platform, cross-platform overlap check, random-effects pooling", "7.3", "pipeline/a30_external_coexpression.py", "results/external/gse13041_platforms.csv"),
    ("GSE4290 with all histologies", "7.3", "pipeline/a30_external_coexpression.py", "results/external/gse4290_all_histologies.csv"),
    ("Earlier same-session external values (not reproduced)", "7.3", "pipeline/a30_external_coexpression.py", "results/external/earlier_same_session_values.csv"),
    ("Survival replication and meta-analysis", "8", "pipeline/a31_survival_replication.py", "results/external/survival_replication_cox.csv"),
    ("Over-representation and GSEA preranked", "10", "pipeline/a40_enrichment.py", "results/enrichment/targeted_sets.csv"),
    ("Optional: ORA and GSEA of MGMT-ranked genes", "10", "pipeline/a40_enrichment.py", "results/enrichment/mgmt_ora_list.csv"),
    ("CAMERA", "10", "pipeline/a41_camera.py", "results/enrichment/camera_targeted_sets.csv"),
    ("Random survival forest", "11", "pipeline/a50_survival_forest.py", "results/forest/forest_paired_difference.csv"),
    ("Differential expression, DESeq2 with apeglm", "12", "pipeline/a60_deseq2.py", "results/de/de_summary.csv"),
    ("Differential expression, PyDESeq2 cross-check", "12", "pipeline/a61_pydeseq2.py", "results/de/pydeseq2_vs_deseq2.csv"),
    ("Verification tests T1-T15", "13.1", "tests/", "results/tests/test_results.csv"),
    ("Independent reproduction of screen, Model 3, KM / log-rank and meta-analyses", "13.2",
     "pipeline/a91_independent_check.py", "results/independent/independent_check_summary.csv"),
]
PER_GENE = ("screen, batch-adjusted: top 100", "screen, batch-adjusted: 500 random")


def r_versions():
    code = ('for (p in c("survival","DESeq2","apeglm","limma","mice")) cat(p, as.character(packageVersion(p)), "\\n"); '
            'cat("R", R.version$major, R.version$minor, "\\n"); cat("Bioconductor", as.character(BiocManager::version()), "\\n")')
    try:
        out = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=True,
                             env={**__import__("os").environ, "R_LIBS_USER": ".Rlib"}).stdout
        return dict(line.split(" ", 1) for line in out.strip().splitlines())
    except Exception as e:  # recorded, not fatal
        return {"error": str(e)}


def main():
    common.require_registration()
    py = {p: md.version(p) for p in ("numpy", "pandas", "scipy", "statsmodels", "scikit-learn", "lifelines", "gseapy",
                                     "pydeseq2", "pingouin", "matplotlib")}
    rec = {"mode": common.MODE, "parameters": common.PARAMS, "python": platform.python_version(),
           "os": platform.platform(), "python_packages": py,
           "forest_environment": ".venv_rsf (scikit-survival 0.28.0, numpy 2.x)", "r_packages": r_versions(),
           "data_release": "GDC Data Release 46.0 (2026-08-10)", "per_analysis_n": per_analysis_n().to_dict("records")}
    common.write_json(common.BASE / "methods_parameters.json", rec)
    rows = [{"statement": s, "plan_section": sec, "script": sc, "output": f"revision/{o}" if common.MODE == "real" else f"revision/sim/{o}",
             "output_exists": (common.BASE / o).exists() if not o.startswith("manifests") else (common.MAN / o.split("/", 1)[1]).exists()}
            for s, sec, sc, o in CONCORDANCE]
    ind = common.RES / "independent/independent_check.csv"
    if ind.exists():
        d = pd.read_csv(ind).rename(columns={"pass": "passed"})
        count = d.analysis.value_counts()
        per_gene = d[d.analysis.isin(PER_GENE)].sort_values("difference").groupby("analysis").tail(1)
        per_gene["quantity"] = [f"largest difference of {count[a]} values: {q}" for a, q in zip(per_gene.analysis, per_gene.quantity)]
        for r in pd.concat([d[~d.analysis.isin(PER_GENE)], per_gene]).itertuples():
            rows.append({"statement": f"Independent reproduction: {r.analysis}: {r.quantity}", "plan_section": "13.2",
                         "script": "pipeline/a91_independent_check.py", "output": f"{'revision' if common.MODE == 'real' else 'revision/sim'}/results/independent/independent_check.csv",
                         "output_exists": True, "pipeline_value": r.pipeline_value, "independent_value": r.independent_value,
                         "independent_implementation": r.independent_implementation, "criterion": r.criterion,
                         "tolerance": r.tolerance, "difference": r.difference, "pass": r.passed})
    pd.DataFrame(rows).to_csv(common.BASE / "methods_concordance.csv", index=False)
    d = pd.DataFrame(rows)
    main_rows = d[~d["statement"].str.startswith("Independent reproduction:")]
    print(main_rows[["plan_section", "script", "output_exists"]].to_string(index=False))
    print(f"independent reproduction rows: {len(d) - len(main_rows)}, failed: "
          f"{int((d.get('pass', pd.Series(dtype=bool)) == False).sum())}")


if __name__ == "__main__":
    main()
