# Editor checklist: PONE-D-26-34494

Each requirement of the decision letter mapped to its supporting file or to an outstanding task. Paths are
relative to the repository root. "Outstanding" means not done; manuscript text, the literature review and every
submission-system action are outstanding unless stated otherwise. The analysis supports and constrains the
manuscript changes; it does not make them.

Status key: **Provided** (file exists and covers the item) · **Partly provided** (file exists; an author step
remains) · **Outstanding** (author action, not started here).

## A. Materials required with the revised submission

| # | Material | Status | Supporting file or outstanding task |
|---|---|---|---|
| 1 | Clean revised manuscript | Outstanding | Author: write after reviewing `revision/numbers_for_manuscript.md` and `revision/claims_not_reproduced.md` |
| 2 | Tracked-changes version | Outstanding | Author |
| 3 | Point-by-point response | Outstanding | Author. Sources: `revision/tables/S9_Table.csv`, `revision/claims_not_reproduced.md`, `revision/DEVIATIONS.md`, this checklist |
| 4 | Frozen code archive or permanent release | Partly provided | Original submission: tag `v1-submitted` (commit `e0d6bb5`). Revised analysis: commit in `revision/submission_package/REVISED_COMMIT.txt`. Outstanding: archival release with DOI (Zenodo or OSF) |
| 5 | Exact software environment | Provided | `revision/environment/` (`pip_freeze_2026-10-28.txt`, `pip_freeze_venv_rsf_2026-10-28.txt`, `R_sessionInfo_2026-10-28.txt`, `R_packages_2026-10-28.csv`, `os_info.txt`, `python_version.txt`) |
| 6 | Complete TCGA and external-cohort manifests | Provided | `revision/manifests/tcga_samples.csv` (S1 Table), `revision/manifests/external_cohort_samples.csv` (S2 Table) |
| 7 | Full eligible-gene correlation table | Provided | `revision/results/screen/screen_gene_level.csv` (all 57,660 non-constant genes, eligibility flags, 26,370 tested); S3 Table |
| 8 | Complete Cox, KM, DE, ML, enrichment, validation outputs | Provided | `revision/results/survival/`, `de/`, `forest/`, `enrichment/`, `external/` (untruncated) |
| 9 | Source data for every figure and table | Provided | `revision/source_data/` (one CSV per figure), `revision/tables/` (each table is its own source data) |
| 10 | Figure-to-script map | Provided | `revision/figure_script_map.csv` (Fig 1-6, S1-S8 Fig, Tables 1-6, S1-S9 Table) |
| 11 | Corrected Data Availability statement | Partly provided | `revision/data_availability_draft.md` (draft). Outstanding: Zenodo DOI; author edit; submission-form update |
| 12 | Updated PubMed search summary with final search date | Outstanding | `revision/literature_search_protocol.md` is a protocol only; no search has been run |
| 13 | AI-use declaration | Outstanding | Author: write the declaration (tool, task, verification) for every task listed in Section 14.9 |
| 14 | Statement of analyses or claims not reproduced and removed | Partly provided | `revision/claims_not_reproduced.md` (generated from S9 Table). Outstanding: author decision on which claims are removed and the statement text |

## B. Journal requirements

| # | Requirement | Status | Supporting file or outstanding task |
|---|---|---|---|
| 1 | PLOS ONE style and file naming | Partly provided | Figures and tables use PLOS names (`Fig_1.tif`, `S1_Fig.tif`, `S1_Table.csv`). Outstanding: manuscript |
| 2 | PLOS LaTeX template | Outstanding | Author |
| 3 | Data sharing plan in the submission form | Outstanding | Author; all inputs are open access (`revision/data_availability_draft.md`) |
| 4 | Identical title on form and manuscript | Outstanding | Author |
| 5 | Cite Figures 6, 7 and 8 in the text | Outstanding | Author. The revised figure set is Fig 1-6 (`revision/figure_script_map.csv`); figure numbering in the text must follow it |
| 6 | Captions for Figures 6, 7 and 8 | Outstanding | Author (captions for the revised Fig 1-6 and S1-S8 Fig) |
| 7 | Unreadable "Supplemental Fig 1.zip" | Outstanding | Author: replace with the individual TIFF files in `revision/figures/` |
| 8 | Evaluate reviewer-recommended citations | Outstanding | Author |

## C. Essential revisions

### 1. One complete and reproducible version

| Item | Status | Supporting file or outstanding task |
|---|---|---|
| 1.1 Exact commit or dated archival release | Partly provided | Revised commit hash in `revision/submission_package/REVISED_COMMIT.txt`; `v1-submitted` = `e0d6bb5`. Outstanding: Zenodo/OSF archive |
| 1.2 Environment file with exact versions | Provided | `revision/environment/` |
| 1.3 Master script from downloaded data | Provided | `run_all.sh` (`--download`, `--analysis`); commands in `README.md` |
| 1.4 TCGA manifest with IDs, specimen type, inclusion, survival, covariates | Provided | `revision/manifests/tcga_samples.csv` |
| 1.5 One primary-tumor sample per patient | Provided | `revision/results/cohort/tcga_cohort_flow.csv`; `selected_aliquot` and `is_primary` in `tcga_samples.csv` |
| 1.6 GDC query parameters, file and workflow types, annotation, units, access dates | Provided | `revision/manifests/gdc_queries/*.json` (incl. `tcga_star__files_meta.json`), `download_dates.md`, `DOWNLOAD_LOG.csv`, `raw_inputs_manifest.csv`; GENCODE v36 and units in `revision/methods_parameters.json` |
| 1.7 All intermediate outputs | Provided | `revision/results/` |
| 1.8 Figure-to-script map | Provided | `revision/figure_script_map.csv` |
| README does not match scripts | Provided | `README.md` (every current script: reads, writes, figure/table); `legacy/README.md` |
| Inconsistent access dates (Nov 12, 2026) | Partly provided | Actual dates in `revision/manifests/download_dates.md`, `raw_inputs_manifest.csv` (claim C25 in S9 Table). Outstanding: manuscript text |

### 2. Genome-wide MGMT correlation screen

| Item | Status | Supporting file or outstanding task |
|---|---|---|
| 2.1-2.2 Gene universe and filters defined in advance | Provided | `revision/PREREGISTRATION.md` Section 5 (registered https://osf.io/qpd2f/); `revision/results/screen/screen_gene_universe.csv`; `results/expression/tcga_gene_id_handling.csv` |
| 2.3-2.4 Spearman vs every eligible gene, BH | Provided | `revision/results/screen/screen_gene_level.csv` |
| 2.5 Complete ranked table with required columns | Provided | S3 Table (`revision/tables/S3_Table.csv`) |
| 2.6 Meaning of "rank first" | Partly provided | Signed and absolute ranks in `screen_gene_level.csv` (`rank_adj`, `abs_rank_adj`). Outstanding: manuscript text |
| 2.7 TPM/FPKM and outlier sensitivity values | Provided | `revision/results/screen/screen_sensitivity.csv` |
| 2.8 Scatterplot from ranking data | Provided | Fig 2 (`revision/figures/Fig_2.tif`; `revision/source_data/Fig_2.csv`) |
| 2.9 Revise or remove claims not reproduced | Partly provided | `revision/claims_not_reproduced.md` (C01, C02, C03). Outstanding: manuscript text |

### 3. Disease cohort (IDH-wildtype)

| Item | Status | Supporting file or outstanding task |
|---|---|---|
| 3.1-3.2 Diagnosis, IDH status and other diagnostic variables | Provided | `revision/manifests/tcga_samples.csv`; `revision/results/cohort/idh_source_concordance.csv`, `tcga_annotation_tables.csv`; S1 Table |
| 3.3 Primary analysis in IDH-wildtype | Provided | Primary cohort (252 aliquots); `results/cohort/tcga_cohort_flow.csv` |
| 3.4 Sensitivity analysis excluding IDH-mutant | Provided | IDH-not-mutant and sensitivity cohorts in `screen_sensitivity.csv`, `survival/cox_models.csv` |
| 3.5 No IDH1 expression as status | Provided | IDH from MAF and publication calls (`idh_source` in `tcga_samples.csv`) |
| 3.6 Exclusions at each stage | Provided | `results/cohort/tcga_cohort_flow.csv`, `external_cohort_flow.csv`; Fig 1 |
| 3.7 Current terminology throughout | Outstanding | Author (manuscript) |

### 4. Survival analysis

| Item | Status | Supporting file or outstanding task |
|---|---|---|
| KM 4.1-4.3 (time and event, log-rank, n, deaths, censored, follow-up, medians with CI) | Provided | `revision/results/survival/km_groups.csv`, `km_logrank.csv`, `km_curves.csv` |
| KM 4.4 At-risk table under every KM figure | Provided | Fig 3, S3 Fig (`km_at_risk.csv`) |
| KM 4.5-4.6 Combination of death/follow-up times; missing or implausible values | Provided | `revision/results/cohort/survival_audit.csv`, `survival_audit_strings.csv`. Outstanding: Methods text |
| KM 4.7-4.8 Median split exploratory; Cox primary | Partly provided | Plan Section 6 (Cox primary). Outstanding: manuscript text |
| Cox 1 Penalized or not | Provided | Penalizer 0 (`revision/methods_parameters.json`, T9 guard). Outstanding: Methods text |
| Cox 2 Complete output for every covariate | Provided | `revision/results/survival/cox_models.csv`; Table 3, S4 Table |
| Cox 3 PH diagnostics | Provided | `revision/results/survival/cox_ph_tests.csv`, `cox_model3_schoenfeld.csv`; S5 Fig |
| Cox 4 No "independent prognostic biomarker" from age-only model | Outstanding | Author (manuscript wording) |
| Cox 5 Clinical/molecular covariates | Provided | Models 4-7 (KPS, methylation, MGMT, treatment variant); extent-of-resection availability in `results/cohort/resection_field_check.csv` |
| Cox 6 State unavailable covariates | Outstanding | Author; sources `resection_field_check.csv`, `methods_parameters.json` (per-analysis n) |
| Cox 7 Nested models with and without MGMT and methylation | Provided | Models 4-7 and `cox_lr_tests.csv` |
| Cox 8 Remove "suppressor variable" | Outstanding | Author (manuscript) |

### 5. Differential expression

| Item | Status | Supporting file or outstanding task |
|---|---|---|
| 5.1 Count-scale inputs | Provided | DESeq2 on Toil expected counts (`pipeline/a59_de_inputs.py`, `a60_deseq2.py`) |
| 5.2 TCGA vs GTEx confounding | Partly provided | Comparison with TCGA normals (`de/deseq2_tcga_normals.csv`). Outstanding: limitation text |
| 5.3 Number and origin of normals | Provided | `revision/results/de/de_summary.csv` (`normal_tissue`) |
| 5.4 Identifier harmonization | Provided | `revision/results/de/de_gene_id_handling.csv` |
| 5.5-5.7 Adjusted p, effect size with uncertainty, complete table | Provided | `revision/results/de/deseq2_*.csv`; S8 Table |
| 5.8 Volcano plot from final output | Provided | S7 Fig |
| 5.9 Reconcile gene counts | Provided | `de_summary.csv` (`genes`, `genes_tested`) |
| 5.10 Remove "most significantly upregulated" | Outstanding | Author (claim C22) |

### 6. Random-forest analysis

| Item | Status | Supporting file or outstanding task |
|---|---|---|
| 6.1-6.5 Outcome and censoring, predictors, repeated CV, distributions, baselines | Provided | Random survival forest (`pipeline/a50_survival_forest.py`), `revision/results/forest/forest_cv_folds.csv`, `forest_c_summary.csv`, `forest_paired_difference.csv` |
| 6.6 Class balance and confusion matrix | Not applicable to the revised method | Survival forest on censored time (C-index); no classification. Outstanding: state this in Methods |
| 6.7 Permutation importance | Provided | `revision/results/forest/forest_permutation_importance.csv`; S8 Fig |
| 6.8-6.10 Descriptive importance; "among prespecified features"; AUC 0.634 | Outstanding | Author; claims C12, C13 |

### 7. Pathway enrichment

| Item | Status | Supporting file or outstanding task |
|---|---|---|
| 7.1-7.4 ORA and GSEA, input lists, background, ranked list and statistic | Provided | `revision/results/enrichment/ora_lists.csv`, `ora_*.csv` (background size per row), `gpr153_correlation.csv` (ranking), `gsea_*.csv`, `camera_*.csv` |
| 7.5 Database names, versions, access dates | Provided | `revision/methods_parameters.json`; `revision/manifests/download_dates.md` |
| 7.6 Complete output files | Provided | `revision/results/enrichment/` |
| 7.7-7.8 Adjusted p in text, tables, figures | Partly provided | Fig 6, S6 Fig, Table 6, S6 Table use adjusted values. Outstanding: text |
| 7.9 Derivation of list; direction | Provided | `ora_lists.csv` (`direction`) |
| 7.10-7.11 Instability, circularity; mechanistic wording | Outstanding | Author |

### 8. External validation

| Item | Status | Supporting file or outstanding task |
|---|---|---|
| 8.1 Accession, platform, selection, classification, n | Provided | `revision/manifests/external_cohort_samples.csv`, `results/cohort/external_cohort_flow.csv`; S2 Table |
| 8.2 Probe-to-gene mapping rule | Provided | `revision/results/external/external_probe_combinations.csv`; plan Section 7 |
| 8.3-8.4 Same statistic; CI, p, adjusted p | Provided | `revision/results/external/external_coexpression.csv`; Table 4 |
| 8.5 Independence of the two GSE13041 analyses | Provided | `gse13041_runs.csv`, `gse13041_platform_overlap.csv`, `results/cohort/external_duplicate_check.csv` |
| 8.6 IDH and classification | Provided | `external_cohort_flow.csv`; `external_composition.csv` |
| 8.7-8.8 Negative results; remove "consistently reproduced" | Outstanding | Author; values in `external_meta.csv` |

### 9. Figures

| Item | Status | Supporting file or outstanding task |
|---|---|---|
| 9.1-9.4 Figure content and labels | Provided | All figures regenerated by `pipeline/f10_main_figures.py`, `f11_supplement_figures.py` |
| 9.5-9.6 README descriptions of 16 and 17 | Provided | `legacy/README.md`; `README.md` |

### 10-13. Literature, Introduction, Discussion, title/Abstract/Conclusion

| Item | Status | Supporting file or outstanding task |
|---|---|---|
| 10 Updated PubMed search and citation audit | Outstanding | `revision/literature_search_protocol.md` (protocol only) |
| 11 Introduction | Outstanding | Author |
| 12 Discussion | Outstanding | Author; values in `revision/numbers_for_manuscript.md` |
| 13 Title, Abstract, Conclusion after analyses | Outstanding | Author; outcome branch d (`revision/results/decision/decision.csv`) |

### 14. Reporting and transparency

| Item | Status | Supporting file or outstanding task |
|---|---|---|
| 14.1 Dates, titles, figure numbers, gene lists, counts, versions | Partly provided | Values in `numbers_for_manuscript.md`, `methods_parameters.json`. Outstanding: manuscript |
| 14.2 Identical title everywhere | Outstanding | Author (the repository README title must also match the final title) |
| 14.3-14.4 Abbreviations, grammar | Outstanding | Author |
| 14.5-14.6 Exact p-values; CIs for primary estimates | Partly provided | All in result files and `numbers_for_manuscript.md`. Outstanding: manuscript |
| 14.7 Data Availability with archived code version and result files | Partly provided | `revision/data_availability_draft.md`. Outstanding: Zenodo DOI |
| 14.8 Acknowledged contributor | Outstanding | Author |
| 14.9 Generative-AI and automated-tool declaration | Outstanding | Author |
| Biostatistics consultation (closing paragraph) | Outstanding | Author |

## D. PLOS administrative issue

| Item | Status | Supporting file or outstanding task |
|---|---|---|
| AI use in the decision letter; contact the journal office | Outstanding | Addressed to the journal office; author to follow any instructions received |
