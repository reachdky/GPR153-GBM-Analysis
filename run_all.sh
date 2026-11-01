#!/usr/bin/env bash
# Run the revision pipeline in order from the repository root (PREREGISTRATION.md, Section 13).
# Usage:
#   ./run_all.sh               cohort construction (clinical tables only) and the tests that need no expression data
#   ./run_all.sh --download    re-run pipeline/r00_download_external.py and pipeline/r10_download_raw_inputs.py
#                              first (md5-checked), then as above
#   ./run_all.sh --sim         as above, then the full analysis chain on simulated data (revision/sim/), incl. T7, T8
#                              and the independent reproduction (a91)
#   ./run_all.sh --analysis    as above, then the full analysis chain on real data. Every step that reads expression
#                              or methylation values stops unless the analysis plan is registered (pipeline/common.py)
# Stops on the first failure.
set -euo pipefail
cd "$(dirname "$0")"
source .venv/bin/activate
export PYTHONHASHSEED=42
export R_LIBS_USER="$PWD/.Rlib"
MODE="${1:-}"

if [[ "$MODE" == "--download" ]]; then
    python pipeline/r00_download_external.py
    python pipeline/r10_download_raw_inputs.py
fi

# Cohort construction and annotation (no GPR153 or MGMT expression values)
python pipeline/r01_tcga_cohort.py
python pipeline/r02_external_cohorts.py
python pipeline/r03_batch_comparison.py
python pipeline/r04_tcga_annotation_tables.py
python pipeline/r05_survival_replication_candidates.py
python pipeline/r06_download_manifest.py
python pipeline/r07_survival_audit.py
python pipeline/r08_cgga_clinical_match.py
python pipeline/r09_resection_field_check.py

# Verification tests on simulated data, R data sets, clinical tables or source code (Section 13.1)
GBM_MODE=real python tests/test_stats.py            # T1-T4, T2b
GBM_MODE=real python tests/test_schoenfeld.py       # T5 (records the PH-test source used by a20)
GBM_MODE=real python tests/test_survival.py         # T6
GBM_MODE=real python tests/test_guards.py           # T9
GBM_MODE=real python tests/test_guards_planted.py   # T9b
GBM_MODE=real python tests/test_survival_audit.py   # T10
GBM_MODE=real python tests/test_screen_simulation.py  # T11, T12
GBM_MODE=real python tests/test_cox_recovery.py     # T13
GBM_MODE=real python tests/test_estimate.py         # T14
GBM_MODE=real python tests/test_cox_vectorized.py   # T15

analysis_chain() {
    # Loaders (real mode: locked until registration; sim mode: replaced by pipeline/sim/s00_make_sim_data.py)
    if [[ "$GBM_MODE" == "real" ]]; then
        python pipeline/a10_tcga_matrix.py
        python pipeline/a19_methylation_inputs.py
        python pipeline/a29_external_inputs.py
        python pipeline/a59_de_inputs.py
    else
        python pipeline/sim/s00_make_sim_data.py
    fi
    # Section 5 and 3.2
    python pipeline/a11_screen.py
    python pipeline/a12_bootstrap_rank.py
    python pipeline/a13_permutation_null.py
    python pipeline/a14_composition.py
    python pipeline/a15_methylation.py
    python pipeline/a16_decision.py
    python pipeline/a17_posthoc_crossbatch.py
    python pipeline/a18_pca.py
    # Section 6
    python pipeline/a20_cox.py
    python pipeline/a21_km.py
    python pipeline/a22_cox_genomewide.py
    python pipeline/a23_mice.py
    # Sections 7 and 8
    python pipeline/a30_external_coexpression.py
    python pipeline/a31_survival_replication.py
    # Section 10
    python pipeline/a40_enrichment.py
    python pipeline/a41_camera.py
    # Section 11 (scikit-survival needs numpy 2; separate environment)
    .venv_rsf/bin/python pipeline/a50_survival_forest.py
    # Section 12
    python pipeline/a60_deseq2.py
    python pipeline/a61_pydeseq2.py
    # Tests that read expression values (Section 13.1)
    python tests/t7_screen_r.py                      # T7
    python tests/t8_cox_r.py                         # T8
    # Independent reproduction of the numbers that may enter the manuscript (Section 13.2)
    python pipeline/a91_independent_check.py
    # Figures, tables, methods record (Sections 13, 14)
    python pipeline/f10_main_figures.py
    python pipeline/f11_supplement_figures.py
    python pipeline/f12_tables.py
    python pipeline/a90_methods_record.py
    # Claims ledger: revised value for each submitted claim (real data only; reads results)
    if [[ "$GBM_MODE" == "real" ]]; then
        python pipeline/x02_claims_ledger.py
    fi
}

case "$MODE" in
    --sim)      export GBM_MODE=sim;  analysis_chain ;;
    --analysis) export GBM_MODE=real; analysis_chain ;;
esac

echo "run_all.sh: completed"
