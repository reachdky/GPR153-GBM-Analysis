# Audit of the analysis against the decision letter (2026-11-01)

Each analytical requirement of the decision letter (PONE-D-26-34494) checked against the code and the final
outputs, not only the presence of a file. "Verified" means the code path and the output were inspected and agree
with the requirement; the check made is stated. Values refer to the final run (`final_run.md`). Manuscript-text
requirements are in `revision/editor_checklist.md`.

## Corrections made as a result of this audit

| Finding | Correction | Rerun |
|---|---|---|
| Proportional-hazards tests were missing for Models 1-3 refitted on the Model 4 complete-case sample (plan 6.2 requires a PH test for every reported fit; letter Section 4, Cox 3) | `pipeline/a20_cox.py` runs R `cox.zph` for those fits; `cox_ph_tests.csv` gains `sample_set` | Run 6 (DEVIATIONS 2026-11-01) |
| Figures did not meet the PLOS specifications (width, alpha channel, text below 8 pt) and had readability defects | `pipeline/plotstyle.py`, `f10`, `f11` (DEVIATIONS 2026-11-01) | Run 6 |
| New result from the added PH tests (not a correction) | On the Model 4 sample, Models 2-3 violate PH globally (p 0.0035-0.0079 across the three cohorts), driven by age (and sex in Model 3); GPR153 term p 0.16-0.28 (Model 2) and 0.34-0.54 (Model 3) (corrected 2026-11-02; the earlier text gave 0.34-0.54 for both). Reported in `final_run.md`; to be stated in the manuscript | — |
| `numbers_for_manuscript.md` labelled the meta-analyses "REML"; the code and the plan use DerSimonian-Laird | Header corrected in `pipeline/x03_manuscript_numbers.py` | Document regenerated |

## Section 1: reproducible version

| Item | Check | Result |
|---|---|---|
| 1.1 Commit / archive | Revised commit and tag; archive prepared (`revision/archive_release/`) | Commit verified locally; push and DOI pending account action |
| 1.2 Environment | Both pip freezes resolve with exact pins (`pip install --dry-run`, 60 and 29 packages); R versions recorded; README pins them | Verified |
| 1.3 Master script from downloaded data | `run_all.sh --download` (r00, r10) then `--analysis`; r10 fresh retrieval reproduced every input byte-for-byte | Verified |
| 1.4 TCGA manifest fields | `tcga_samples.csv` has patient, aliquot, sample type, inclusion flags, OS days, event, age, sex, KPS, IDH (sources), MGMT, treatment, EGFR | Verified |
| 1.5 One sample per patient | 284 selected aliquots, 0 duplicate patients (`selected_aliquot`); primary cohort 252 | Verified by count |
| 1.6 Query parameters, workflow, annotation, units, dates | `gdc_queries/*.json`; STAR - Counts, GENCODE v36, log2(TPM + 1); `provenance.md` | Verified |
| 1.7 Intermediate outputs | `revision/results/` (110 files) | Verified |
| 1.8 Figure-to-script map | Every figure and table regenerated from results by the mapped script and compared by SHA-256 | Verified |

## Section 2: genome-wide screen

| Item | Check | Result |
|---|---|---|
| 2.1-2.2 Gene universe defined in advance; filters | Plan Section 4-5 registered before the run; `screen_gene_universe.csv`: 60,616 genes, 57,660 non-constant, 26,371 eligible (TPM >= 1 in >= 20%), 26,370 tested (MGMT excluded); PAR_Y and version handling in `tcga_gene_id_handling.csv` | Verified |
| 2.3 Spearman with every eligible gene, same samples | `common.partial_spearman`, batch-adjusted, n 252 for every gene; T1-T2 vs scipy/pingouin; T7 vs R | Verified |
| 2.4 BH over the complete set | `a11_screen.py` applies `common.bh` over all eligible p-values | Verified |
| 2.5 Ranked table columns | S3 Table: gene_id, symbol, n, rho, p, q, signed rank, absolute rank (adjusted and pooled) | Verified |
| 2.6 Meaning of "rank first" | Signed rank (`rank_adj`) and absolute rank (`abs_rank_adj`) both reported | Output verified; wording is manuscript |
| 2.7 TPM/FPKM and outlier values | `screen_sensitivity.csv` (FPKM, MGMT |z| > 3 outliers: 0 removed; single-value outlier; filters) | Verified |
| 2.8 Scatterplot from ranking data | Fig 2C plots `screen_gpr153_mgmt_values.csv` (the 252 aliquots of the screen); title states the plotted values are unadjusted | Verified |

## Section 3: cohort

| Item | Check | Result |
|---|---|---|
| 3.1-3.2 Classification; IDH1/IDH2 | IDH from open MAFs (IDH1 R132 and IDH2 R172; 21 IDH1 variants, 0 IDH2), Ceccarelli 2016, Brennan 2013; concordance in `idh_source_concordance.csv` | Verified |
| 3.3 Primary analysis in IDH-wildtype | Primary cohort = IDH-wildtype (252; survival 249) | Verified |
| 3.4 Sensitivity excluding IDH-mutant | IDH-not-mutant and sensitivity cohorts in the screen and Cox models | Verified |
| 3.5 No IDH1 expression as status | IDH status from mutation calls only (`r01_tcga_cohort.py`) | Verified |
| 3.6 Exclusions by stage | `tcga_cohort_flow.csv`, `external_cohort_flow.csv`, Fig 1 | Verified |

## Section 4: survival

| Item | Check | Result |
|---|---|---|
| KM 1-2 Time and event; log-rank | `a21_km.py` uses lifelines `logrank_test` with durations and events; T6 known-answer test; a91 reproduces in R | Verified |
| KM 3 n, deaths, censored, median follow-up, medians with CI, at-risk | `km_groups.csv`, `km_logrank.csv`, `km_at_risk.csv`. Median follow-up (reverse KM) is not reached in the MGMT-high groups (written as `inf`) | Verified; `inf` = not reached, to be stated |
| KM 4 At-risk table under every KM figure | Fig 3 and S3 Fig, rows labelled by group, from `km_at_risk.csv` | Verified (visual) |
| KM 5-6 Death/follow-up combination; missing/implausible | `r01` (days_to_death if dead, else max follow-up); `survival_audit.csv` (time missing, time <= 0 excluded) | Verified |
| Cox 1 Penalization | `penalizer = 0` everywhere (T9 guard); `methods_parameters.json` | Verified |
| Cox 2 Full output per covariate | `cox_models.csv`: coef, SE, HR, CI, p, n, events for every covariate of every fit | Verified |
| Cox 3 PH diagnostics | R `cox.zph` (km) per term and global for every fit, incl. the Model 4-sample refits after the correction above; T5 failed comparison with lifelines, R used as prescribed | Verified after correction |
| Cox 5 Clinical/molecular covariates | Age, sex, KPS, MGMT methylation (STP27), MGMT expression, IDH (sensitivity cohort), batch, treatment (immortal-time sensitivity only) | Verified; extent of resection unavailable (`resection_field_check.csv`) |
| Cox 7 Nested with/without MGMT expression and methylation | Models 4-7, LR tests on identical samples (`cox_lr_tests.csv`) | Verified |

## Section 5: differential expression

| Item | Check | Result |
|---|---|---|
| 5.1 Count inputs | `a59_de_inputs.py` back-transforms Toil log2(x + 1) expected counts and rounds; DESeq2 `~ condition`, apeglm | Verified |
| 5.2 TCGA-GTEx confounding | Same-pipeline (Toil) comparison; separate TCGA-normals comparison (n 5) | Verified; limitation is manuscript |
| 5.3 Normals | `de_summary.csv`: 151 GTEx (Brain - Cortex 102, Frontal Cortex BA9 49) | Verified |
| 5.4 Identifiers | `de_gene_id_handling.csv` | Verified |
| 5.5-5.7 Adjusted p, effect size and SE, full table | `deseq2_*.csv` (padj, log2FC, lfcSE, apeglm); PyDESeq2 cross-check | Verified |
| 5.8 Volcano | S7 Fig from `deseq2_primary.csv`; -log10 cap stated | Verified |
| 5.9 Gene counts | `de_summary.csv`: 60,498 genes, 37,987 tested | Verified |

## Section 6: random forest

| Item | Check | Result |
|---|---|---|
| 6.1 Outcome, censoring | Random survival forest on (time, event); censoring retained | Verified |
| 6.2 Predictors | Feature sets fixed in plan Section 11 | Verified |
| 6.3-6.4 Repeated CV, distribution | 10 x 5-fold stratified by event; C per fold, mean and percentiles | Verified |
| 6.5 Baselines | Set A (age + sex; + KPS, methylation in complete cases); Nadeau-Bengio paired test C vs B | Verified |
| 6.6 Class balance, confusion matrix | Not applicable: no classification; n and events reported per analysis | To be stated in Methods |
| 6.7 Held-out permutation importance | `permutation_importance` on test folds, 10 repeats, SD over 50 folds | Verified |

## Section 7: enrichment

| Item | Check | Result |
|---|---|---|
| 7.1-7.4 ORA and GSEA; lists; background; ranking statistic | `ora_lists.csv`, `mgmt_ora_list.csv`; background = eligible genes in each library (`background_size`); GSEA ranked by batch-adjusted partial Spearman with GPR153 (`gpr153_correlation.csv`) | Verified |
| 7.5 Versions and dates | MSigDB 2026.2.Hs, Enrichr libraries, retrieval dates in `DOWNLOAD_LOG.csv` | Verified |
| 7.6 Complete output | `results/enrichment/` (all sets, untruncated) | Verified |
| 7.7-7.8 Adjusted values in figures and tables | Fig 6, S6 Fig, Table 6 use BH `q_bh`, GSEA FDR, CAMERA FDR | Verified |

## Section 8: external validation

| Item | Check | Result |
|---|---|---|
| 8.1 Accession, platform, selection, n | `external_cohort_samples.csv`, `external_cohort_flow.csv` | Verified |
| 8.2 Probe rule | Highest-mean probe among 221902_at and 64942_at; all combinations in `external_probe_combinations.csv` | Verified |
| 8.3-8.4 Same statistic; CI, p, adjusted p | Spearman, Bonett-Wright CI, BH across cohorts (`spearman_q_bh`) | Verified |
| 8.5 GSE13041 independence | GPL96 primary arrays only; duplicate profile removed; `gse13041_platform_overlap.csv` | Verified |
| 8.6 IDH | CGGA IDH-wildtype only; CPTAC IDH from MAFs; GSE13041, GSE4290, Rembrandt have no IDH data | Verified; limitation is manuscript |
| Meta-analysis method | DerSimonian-Laird with Hartung-Knapp sensitivity (`common.meta_random`), reproduced with R metafor (a91) | Verified |

## Section 14 (reporting values)

| Item | Check | Result |
|---|---|---|
| 14.5-14.6 Exact p and CIs | Every estimate in `numbers_for_manuscript.md` has its CI and source | Verified |
