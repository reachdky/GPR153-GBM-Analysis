# GPR153 expression, MGMT co-expression and survival in IDH-wildtype glioblastoma: revised analysis

Code and outputs for the revised analysis of PLOS ONE submission PONE-D-26-34494. The analysis follows a
time-stamped analysis plan, registered on OSF before any revised GPR153 or MGMT value was computed:

- Plan: `revision/PREREGISTRATION.md` (version 1.5), SHA-256 in `revision/PREREGISTRATION.sha256`
- Registration: https://osf.io/qpd2f/ (2026-10-27 23:44:45 UTC; record in `revision/PREREGISTRATION_REGISTRATION.md`)
- Departures from the plan, with reasons: `revision/DEVIATIONS.md`

Code versions:

- Original submission: tag `v1-submitted` (commit `e0d6bb5`, 2026-07-06). Its scripts are kept unchanged in
  `legacy/` and are not used by the revised analysis; see `legacy/README.md`. The README of that version is
  available at the tag.
- Revised analysis: the commit that adds `pipeline/` and this README; its full hash is recorded in the
  submission package (`revision/submission_package/REVISED_COMMIT.txt`).

## Layout

| Path | Contents |
|---|---|
| `pipeline/` | Revised analysis scripts (`r*` cohort construction, `a*` analysis, `f*` figures and tables, `x*` reporting), shared code in `common.py` and `plotstyle.py` |
| `pipeline/sim/` | Simulated inputs for testing the chain without real expression values |
| `tests/` | Verification tests T1-T15 (plan Section 13.1) |
| `run_all.sh` | Runs everything in order |
| `revision/manifests/` | Sample manifests (TCGA, external cohorts), download log with md5, GDC query records |
| `revision/results/` | All analysis outputs, untruncated |
| `revision/figures/`, `revision/tables/`, `revision/source_data/` | Figures (TIFF, 300 dpi), tables (CSV), one source-data CSV per figure |
| `revision/figure_script_map.csv` | Each figure and table -> script -> output file -> input files |
| `revision/methods_concordance.csv` | Each Methods statement -> plan section -> script -> output; independent re-implementation values |
| `revision/methods_parameters.json` | Fixed parameters, seeds, package versions, per-analysis n |
| `revision/environment/` | Python and R package versions, OS, install logs |
| `revision/verification/` | Final-run identification, run-to-run comparisons |
| `revision/numbers_for_manuscript.md`, `claims_not_reproduced.md`, `data_availability_draft.md`, `editor_checklist.md` | Manuscript-support documents |
| `legacy/` | Scripts as posted at initial submission (superseded) |

## Scripts

Paths are relative to the repository root. `work/` means `revision/work/` (intermediate files, not distributed);
`results/` means `revision/results/`. "Locked" scripts stop unless the analysis plan is registered (hash and OSF
URL checked by `common.require_registration`).

### Downloads and cohort construction (no GPR153 or MGMT expression values read)

| Script | Reads | Writes | Feeds |
|---|---|---|---|
| `r00_download_external.py` | Internet: GDC API, cBioPortal, CGGA, UCSC Toil, GEO, MSigDB, Enrichr, Europe PMC, DepMap, GDSC | `data/downloads/**`, `data/downloads/DOWNLOAD_LOG.csv` (URL, UTC time, bytes, md5) | all |
| `r01_tcga_cohort.py` | `data/downloads/gdc_meta/` (expression file metadata, cases), `data/downloads/gdc_maf/*.maf.gz`, `data/downloads/cbioportal/*.json`, `data/downloads/gdc_methylation/files_meta.json` | `revision/manifests/tcga_samples.csv`, `results/cohort/tcga_cohort_flow.csv`, `results/cohort/idh_source_concordance.csv`, `revision/figures/FigS1_cohort_flow.png`, `revision/source_data/FigS1_cohort_flow.csv` | Fig 1, S1 Table |
| `r02_external_cohorts.py` | `data/raw/geo/GSE13041/`, `data/raw/geo/GSE4290/` (`clinical.csv`, `expression.csv`), `data/raw/geo/GSE13041_family.soft.gz`, `data/downloads/cgga/CGGA.mRNAseq_{693,325}*.zip` | `revision/manifests/external_cohort_samples.csv`, `results/cohort/external_cohort_flow.csv`, `results/cohort/external_duplicate_check.csv` | Fig 1, S2 Table |
| `r03_batch_comparison.py` | `revision/manifests/tcga_samples.csv` | `results/cohort/tcga_batch_comparison.csv`, `tcga_batch_os.csv` | Table 1 |
| `r04_tcga_annotation_tables.py` | `revision/manifests/tcga_samples.csv` | `results/cohort/tcga_annotation_tables.csv`, `tcga_cohort_overlap.csv` | — |
| `r05_survival_replication_candidates.py` | `revision/manifests/tcga_samples.csv`; clinical and expression files of each candidate cohort under `data/downloads/` and `data/raw/` (e.g. `cgga/CGGA.mRNAseq_{693,325}*.zip`) | `results/cohort/survival_replication_candidates.csv` | Table 5 |
| `r06_download_manifest.py` | `data/downloads/DOWNLOAD_LOG.csv`, GDC query records | `revision/manifests/DOWNLOAD_LOG.csv`, `download_dates.md`, `gdc_queries/*.json` | Data Availability |
| `r07_survival_audit.py` | both manifests, `results/cohort/survival_replication_candidates.csv`, raw clinical files under `data/raw` and `data/downloads` | `results/cohort/survival_audit_strings.csv`, `survival_audit.csv` | Table 5, S5 Table |
| `r08_cgga_clinical_match.py` | CGGA clinical tables (`data/downloads/cgga/`) | `results/cohort/cgga_clinical_match.csv` | — |
| `r09_resection_field_check.py` | `revision/manifests/tcga_samples.csv`; `data/downloads/tcga_pub/TCGA-CDR-SupplementalTableS1.xlsx`, `clinical_PANCAN_patient_with_followup.tsv`; `data/downloads/gdc_bcr_xml/*/*clinical*.xml`; GDC diagnoses | `results/cohort/resection_field_check.csv` | — |
| `r10_download_raw_inputs.py` | Internet: GDC API (`/files`, `/status`, `/data`), GEO FTP; `revision/manifests/tcga_samples.csv` (file IDs, GDC md5), `revision/manifests/raw_inputs_expected_md5.csv` | `data/raw/tcga/expression/TCGA-GBM/Transcriptome_Profiling/Gene_Expression_Quantification/<file_id>/<file>`, `data/raw/geo/{GSE13041,GSE4290}_family.soft.gz`, `data/raw/geo/{GSE13041,GSE4290}/{expression,clinical}.csv` (only when the md5 matches); `revision/manifests/gdc_queries/tcga_star__files_meta.json`, `revision/manifests/raw_inputs_manifest.csv` | inputs of a10, r02, r07, a29 |

### Loaders (locked)

| Script | Reads | Writes |
|---|---|---|
| `a10_tcga_matrix.py` | `revision/manifests/tcga_samples.csv`; 391 GDC STAR - Counts files under `data/raw/tcga/expression/TCGA-GBM/Transcriptome_Profiling/Gene_Expression_Quantification/<file_id>/` (each md5-checked against GDC) | `work/tcga_expr.npz`, `results/expression/tcga_gene_id_handling.csv` |
| `a19_methylation_inputs.py` | `data/downloads/gdc_methylation/` (SeSAMe level-3 beta files) | `work/methylation_inputs.csv` |
| `a29_external_inputs.py` | `revision/manifests/external_cohort_samples.csv`, `data/raw/geo/*`, `data/downloads/{cgga,geo_replication,cptac}` | `work/ext_<cohort>.npz`, `work/ext_<cohort>_clinical.csv` (GSE13041 GPL96 and GPL570, GSE4290, CGGA 693 and 325, Rembrandt, CPTAC-3) |
| `a59_de_inputs.py` | `data/downloads/toil/TcgaTargetGtex_gene_expected_count.gz`, `TcgaTargetGTEX_phenotype.txt.gz`, `revision/manifests/tcga_samples.csv` | `work/de_toil.npz`, `work/de_toil_samples.csv` |

### Analysis

| Script | Reads | Writes (`results/`) | Feeds |
|---|---|---|---|
| `a11_screen.py` | `work/tcga_expr.npz`, `tcga_samples.csv` | `screen/screen_gene_level.csv`, `screen_sensitivity.csv`, `screen_gene_universe.csv`, `crossbatch_*.csv` | Fig 2, S1 Fig, S4 Fig, Table 2, S3 Table, S4 Table |
| `a12_bootstrap_rank.py` | `work/tcga_expr.npz`, `tcga_samples.csv`, `screen/screen_gene_level.csv` | `screen/bootstrap_rank_draws.csv`, `bootstrap_rank_summary.csv` | Fig 2 |
| `a13_permutation_null.py` | as a12 | `screen/permutation_null_draws.csv`, `permutation_null_summary.csv` | Fig 2 |
| `a14_composition.py` | `work/tcga_expr.npz`, `tcga_samples.csv`, `screen/screen_gene_level.csv`, `data/downloads/estimate/extracted/*.csv` | `composition/composition_scores.csv`, `composition_adjusted.csv`, `composition_screen_<score>.csv`, `neuronal_genes_used.csv` | Table 2 |
| `a15_methylation.py` | `work/tcga_expr.npz`, `work/methylation_inputs.csv`, `tcga_samples.csv` | `methylation/stp27_calls.csv`, `methylation_concordance.csv`, `methylation_kappa_2x2.csv`, `condition4.csv`, `condition4_sensitivity.csv`, `methylation_descriptive.csv`; `work/methylation_measure.csv` | Table 2, S4 Table |
| `a16_decision.py` | `screen/screen_sensitivity.csv`, `composition/composition_adjusted.csv`, `methylation/condition4.csv` | `decision/decision_conditions.csv`, `decision.csv` | Table 2 |
| `a17_posthoc_crossbatch.py` | `work/tcga_expr.npz`, `tcga_samples.csv`, `screen/crossbatch_reproducibility.csv` | `screen/posthoc_crossbatch_gpr153_mgmt.csv` | S4 Fig |
| `a18_pca.py` | `work/tcga_expr.npz`, `tcga_samples.csv`, `screen/screen_gene_level.csv` | `pca/pca_scores.csv`, `pca_variance.csv`, `pca_genes.csv` | S2 Fig |
| `a20_cox.py` (+ `a20_coxzph.R`) | `work/tcga_expr.npz`, `work/methylation_measure.csv`, `tcga_samples.csv` | `survival/cox_models.csv`, `cox_lr_tests.csv`, `cox_ph_tests.csv`, `cox_age_spline.csv`, `cox_model3_schoenfeld.csv`; `work/survival_table_<cohort>.csv`, `work/mi_input.csv` | Fig 4, Fig 5B, S5 Fig, Table 3, S4 Table |
| `a21_km.py` | `work/survival_table_primary.csv` | `survival/km_groups.csv`, `km_logrank.csv`, `km_at_risk.csv`, `km_curves.csv` | Fig 3, S3 Fig |
| `a22_cox_genomewide.py` | `work/tcga_expr.npz`, `work/survival_table_primary.csv`, `screen/screen_gene_level.csv` | `survival/cox_genomewide.csv`, `cox_genomewide_gpr153.csv`, `incremental_value.csv` | — |
| `a23_mice.py` (+ `a23_mice.R`) | `work/mi_input.csv` | `survival/cox_multiple_imputation.csv` | S7 Table |
| `a30_external_coexpression.py` | `work/ext_*.npz` and `_clinical.csv`, `screen/screen_sensitivity.csv`; original validation CSVs in `~/gpr153_project_orig/data/results/validation/` (read only, for the "earlier values" rows) | `external/external_coexpression.csv`, `external_probe_combinations.csv`, `gse13041_runs.csv`, `gse13041_u1_criterion.csv`, `cgga_pooled.csv`, `external_meta.csv`, `external_power.csv`, `external_discordance.csv`, `external_composition.csv`, `gse13041_platforms.csv`, `gse13041_platform_overlap.csv`, `gse4290_all_histologies.csv`, `earlier_same_session_values.csv` | Fig 5A, Table 4, S4 Table |
| `a31_survival_replication.py` | `work/ext_*.npz` and `_clinical.csv`, `survival/cox_models.csv`, `cohort/survival_audit.csv` | `external/survival_replication_cox.csv`, `survival_replication_meta.csv`, `survival_replication_incremental.csv` | Fig 5B, Table 5, S4 Table |
| `a40_enrichment.py` | `work/tcga_expr.npz`, `tcga_samples.csv`, `screen/screen_gene_level.csv`, `data/downloads/genesets/*.gmt` | `enrichment/gpr153_correlation.csv`, `ora_lists.csv`, `ora_<library>_<direction>.csv`, `gsea_<library>.csv`, `targeted_sets.csv`, `mgmt_ora_*.csv`, `mgmt_gsea_*.csv`; `work/camera_expression.csv.gz`, `camera_covariates.csv` | Fig 6, Table 6, S6 Table |
| `a41_camera.py` (+ `a41_camera.R`) | `work/camera_expression.csv.gz`, `camera_covariates.csv`, `data/downloads/genesets/` | `enrichment/camera_<library>_igc_<estimated or 0.01>.csv`, `camera_targeted_sets.csv` | S6 Fig, S6 Table |
| `a50_survival_forest.py` (runs in `.venv_rsf`) | `work/survival_table_{primary,sensitivity}.csv` | `forest/forest_cv_folds.csv`, `forest_c_summary.csv`, `forest_paired_difference.csv`, `forest_permutation_importance.csv` | S8 Fig |
| `a60_deseq2.py` (+ `a60_deseq2.R`) | `work/de_toil.npz`, `de_toil_samples.csv`, `work/tcga_expr.npz` | `de/deseq2_<comparison>.csv`, `de_summary.csv`, `de_gene_id_handling.csv`; `work/de_primary_counts.csv.gz`, `de_primary_samples.csv` | S7 Fig, S8 Table |
| `a61_pydeseq2.py` | `work/de_primary_counts.csv.gz`, `de_primary_samples.csv`, `de/deseq2_primary.csv` | `de/pydeseq2_primary.csv`, `pydeseq2_vs_deseq2.csv` | — |
| `a91_independent_check.py` (+ `a91_independent.R`) | `work/tcga_expr.npz`, `work/survival_table_primary.csv`, `results/screen/`, `survival/`, `external/` | `independent/independent_check.csv`, `independent_check_summary.csv` (stops on any failure) | `methods_concordance.csv` |
| `a90_methods_record.py` | all result files | `revision/methods_parameters.json`, `revision/methods_concordance.csv` | — |

### Figures, tables, reporting (read result files only)

| Script | Reads | Writes | Produces |
|---|---|---|---|
| `f10_main_figures.py` | `results/cohort/`, `screen/`, `survival/`, `external/`, `enrichment/` | `revision/figures/Fig_<n>.tif`, `revision/source_data/Fig_<n>.csv`, rows in `figure_script_map.csv` | Fig 1-6 |
| `f11_supplement_figures.py` | `results/screen/`, `pca/`, `survival/`, `enrichment/`, `de/`, `forest/` | `revision/figures/S<n>_Fig.tif`, `revision/source_data/S<n>_Fig.csv`, rows in `figure_script_map.csv` | S1-S8 Fig |
| `f12_tables.py` | result files listed per table in `figure_script_map.csv`; both manifests | `revision/tables/<item>.csv`, rows in `figure_script_map.csv` | Tables 1-6, S1-S8 Table |
| `x02_claims_ledger.py` | `results/claims_trace_v1.csv` and the result files it cites | `results/claims_trace_v2.csv`, `revision/tables/S9_Table.csv`, "Revised value" rows in `methods_concordance.csv` | S9 Table |
| `x03_manuscript_numbers.py` | the result files named in each output row | `revision/numbers_for_manuscript.md` (value, file, row key, column) | — (not run by `run_all.sh`) |
| `x04_claims_not_reproduced.py` | `revision/tables/S9_Table.csv` | `revision/claims_not_reproduced.md` | — (not run by `run_all.sh`) |
| `x01_plan_summary_pdf.py` | `revision/PREREGISTRATION.md` | `revision/PREREGISTRATION_summary.pdf`, `.png` | — (not run by `run_all.sh`) |
| `sim/s00_make_sim_data.py` | both manifests, `results/cohort/` | `revision/sim/**` (simulated stand-ins for every loader output) | — (`--sim` only) |

### Tests (`tests/`, plan Section 13.1)

| File | Test |
|---|---|
| `test_stats.py` | T1-T4, T2b: partial Spearman vs scipy and pingouin; BH vs statsmodels; Fisher z and Bonett-Wright CI |
| `test_schoenfeld.py` | T5: lifelines vs R `cox.zph` (fails; R `cox.zph` is used, `DEVIATIONS.md` 2026-10-26) |
| `test_survival.py` | T6: log-rank known answer, event inversion, sample identity |
| `t7_screen_r.py` (+ `t7_screen.R`), `t8_cox_r.py` | T7, T8: screen and Model 3 re-run in R (after the analysis chain) |
| `test_guards.py`, `test_guards_planted.py` | T9, T9b: source-code guards and planted-violation check |
| `test_survival_audit.py` | T10: survival-data audit |
| `test_screen_simulation.py` | T11, T12: null calibration and planted effect |
| `test_cox_recovery.py` | T13: Cox CI coverage |
| `test_estimate.py` | T14: ESTIMATE implementation vs tidyestimate 1.1.1 |
| `test_cox_vectorized.py` | T15: vectorized Cox vs R `coxph` |
| `_report.py` | Shared recorder: one row per test in `revision/results/tests/test_results.csv` |

## Reproducing the analysis

### Environment (exact versions in `revision/environment/`)

- Ubuntu 24.04.5 LTS, Python 3.12.3, R 4.6.1
- `.venv` (analysis): `revision/environment/pip_freeze_2026-10-28.txt`
- `.venv_rsf` (survival forest; scikit-survival 0.28.0 needs numpy >= 2): `revision/environment/pip_freeze_venv_rsf_2026-10-28.txt`
- R packages in `.Rlib` (Bioconductor 3.23): survival 3.8-12, metafor 5.2-1, limma 3.68.5, DESeq2 1.52.0,
  apeglm 1.34.0, mice 3.19.0, dplyr 1.2.1 (T14); full list `revision/environment/R_packages_2026-10-28.csv`,
  session `revision/environment/R_sessionInfo_2026-10-28.txt`

```bash
python3.12 -m venv .venv && .venv/bin/pip install -r revision/environment/pip_freeze_2026-10-28.txt
python3.12 -m venv .venv_rsf && .venv_rsf/bin/pip install -r revision/environment/pip_freeze_venv_rsf_2026-10-28.txt
mkdir -p .Rlib && R_LIBS_USER=$PWD/.Rlib Rscript -e 'install.packages(c("BiocManager","metafor","mice","dplyr"), repos="https://cloud.r-project.org"); BiocManager::install(version="3.23"); BiocManager::install(c("limma","DESeq2","apeglm"))'
```

### Data

All inputs are open access; `data/` is not distributed and is rebuilt by two scripts.

```bash
source .venv/bin/activate
python pipeline/r00_download_external.py     # GDC metadata, MAFs, methylation, CPTAC, cBioPortal, CGGA, UCSC Toil,
                                             # GEO, MSigDB, Enrichr, ESTIMATE, DepMap, GDSC, CELLxGENE; URL, time, md5
python pipeline/r10_download_raw_inputs.py   # 391 TCGA-GBM STAR - Counts files (GDC) and GSE13041, GSE4290 (GEO)
```

`r10` checks every file against its expected md5 (GDC md5 in `revision/manifests/tcga_samples.csv`; the md5 of
the GEO files read by the final run in `revision/manifests/raw_inputs_expected_md5.csv`) and keeps a file only if
it matches; the outcome per file is in `revision/manifests/raw_inputs_manifest.csv`. `a10` checks the STAR md5s
again before reading them. GDC data release, query parameters and retrieval times are in
`revision/manifests/gdc_queries/` and `revision/manifests/download_dates.md`.

No controlled-access data are used. TCGA mutation calls are the open-access masked MAFs.

### Run

```bash
./run_all.sh --download    # r00 and r10 first, then cohort construction and tests
./run_all.sh               # cohort construction and tests that need no expression data
./run_all.sh --sim         # full chain on simulated data (revision/sim/)
./run_all.sh --analysis    # full chain on real data (requires the registered plan)
```

`run_all.sh --analysis` runs, in order: r01-r09; tests T1-T6, T9, T9b, T10-T15; loaders a10, a19, a29, a59;
analysis a11-a18, a20-a23, a30-a31, a40-a41, a50, a60-a61; tests T7, T8; independent reproduction a91;
figures and tables f10-f12; methods record a90; claims ledger x02. It stops on the first failure. Afterwards,
`python pipeline/x03_manuscript_numbers.py` and `python pipeline/x04_claims_not_reproduced.py` write the two
summary documents. Run time on
22 cores: about 34 minutes. Seeds are fixed (42; 43 for the independent check), and repeated runs give
byte-identical files under `revision/results/` apart from run timestamps and timings in
`results/tests/test_results.csv`; in `methods_concordance.csv` the computed relative-difference column can
differ in the last printed digit.
