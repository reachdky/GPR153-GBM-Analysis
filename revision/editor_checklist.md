# Editor checklist: PONE-D-26-34494

Every requirement of the decision letter, mapped to its evidence and status. Analytical evidence was checked
against code and outputs (`revision/verification/revision_letter_audit.md`), not only for file presence. Paths
are relative to the repository root.

**Status key**
- **Done**: completed and verified; evidence given.
- **Done, account action pending**: prepared and verified; a step needs the authors' account (Zenodo).
- **Drafted (author files)**: written in the revised manuscript or submission documents, which are kept with the author and are not in this repository until publication.
- **Not supplied**: required, and not yet provided by the author.
- **Next phase: submission**: an action in Editorial Manager or the submission form.
- **Not applicable**: the item does not apply to the revised method (reason given).
- **Not pursued**: an author decision; does not block completion.
- **Journal-office matter**: addressed to the journal, not to the authors; no author action without instruction.

## A. Materials required with the revised submission

| # | Material | Status | Evidence / next step |
|---|---|---|---|
| 1 | Clean revised manuscript | Drafted (author files) | Sources: `revision/numbers_for_manuscript.md`, `revision/claims_not_reproduced.md`, `revision/literature/literature_search_summary.md` |
| 2 | Tracked-changes version | Drafted (author files) | — |
| 3 | Point-by-point response | Drafted (author files) | Sources: this checklist, `revision/verification/revision_letter_audit.md`, `revision/tables/S9_Table.csv`, `revision/DEVIATIONS.md` |
| 4 | Frozen code archive / permanent release | Done, account action pending | Analysis version `55e41b3` (final run 7), public on GitHub; later commits licence, documentation and archive build script only; tag `v1-submitted` = `e0d6bb5`; archive files in `revision/archive_release/` (archive commit in `REVISED_COMMIT.txt`). Zenodo DOI 10.5281/zenodo.23872434 reserved; record not yet published (authors' account) |
| 5 | Exact software environment | Done | `revision/environment/` (pip freezes resolve exactly on PyPI, checked 2026-11-01; R versions); `README.md` install commands with exact versions |
| 6 | TCGA and external-cohort manifests | Done | `revision/manifests/tcga_samples.csv` (S1 Table), `external_cohort_samples.csv` (S2 Table) |
| 7 | Eligible-gene correlation table | Done | `revision/results/screen/screen_gene_level.csv`; S3 Table |
| 8 | Cox, KM, DE, ML, enrichment, validation outputs | Done | `revision/results/{survival,de,forest,enrichment,external}/` |
| 9 | Source data for every figure and table | Done | `revision/source_data/Fig1.csv`-`Fig6.csv`, `S1_Fig.csv`-`S8_Fig.csv`; tables in `revision/tables/` |
| 10 | Figure-to-script map | Done | `revision/figure_script_map.csv`; every entry regenerated from results and compared by SHA-256 (`revision/verification/final_outputs_regenerated_sha256.csv`) |
| 11 | Corrected Data Availability statement | Drafted (author files) | `revision/data_availability_draft.md`; reserved Zenodo DOI stated as not yet published |
| 12 | Updated PubMed search summary with final date | Done (evidence); statement is next phase | `revision/literature/literature_search_summary.md` (searches run 2026-11-01; queries, counts, PMIDs) |
| 13 | AI-use declaration | **Not supplied** | The author's declaration of tools, tasks and verification has not been provided |
| 14 | Statement of analyses or claims not reproduced and removed | Drafted (author files) | `revision/claims_not_reproduced.md` (evidence); statement lists the 10 removed and 15 changed claims |

## B. Journal requirements

| # | Requirement | Status | Evidence / next step |
|---|---|---|---|
| 1 | PLOS style and file naming | Done for figures and tables; manuscript next phase | `Fig1.tif`-`Fig6.tif`, `S1_Fig.tif`-`S8_Fig.tif`, `S1_Table.csv`-`S9_Table.csv` |
| 2 | PLOS LaTeX template | Drafted (author files) | — |
| 3 | Data sharing plan on the submission form | Next phase: submission | All data open access (`revision/manifests/provenance.md`) |
| 4 | Identical title on form and manuscript | Next phase: submission | — |
| 5 | Cite Figures 6, 7, 8 in the text | Drafted (author files) | The revised set is Fig 1-6 and S1-S8 Fig; numbering in the text follows it |
| 6 | Captions for Figures 6, 7, 8 | Drafted (author files) | Captions for the revised figures |
| 7 | Unreadable "Supplemental Fig 1.zip" | Next phase: submission | Upload individual files from `revision/figures/` |
| 8 | Evaluate reviewer-recommended citations | Done (verification); citing is next phase | All 8 verified (`revision/literature/citation_audit.csv`, R1-R8) |

## C. Essential revisions

### 1. One complete and reproducible version

| Item | Status | Evidence |
|---|---|---|
| 1.1 Exact commit or archival release (Zenodo/OSF preferred) | Done, account action pending | Analysis version `55e41b3`; `revision/archive_release/`; Zenodo DOI 10.5281/zenodo.23872434 reserved, not yet published |
| 1.2 Environment file with exact versions | Done | `revision/environment/`; dry-run resolution check in `README.md` |
| 1.3 Master script from downloaded data | Done | `run_all.sh --download` / `--analysis`; r10 fresh retrieval byte-identical (`revision/manifests/raw_inputs_manifest.csv`) |
| 1.4 TCGA manifest content | Done | `revision/manifests/tcga_samples.csv` |
| 1.5 One primary sample per patient | Done | `selected_aliquot` (284 patients, no duplicates); `results/cohort/tcga_cohort_flow.csv` |
| 1.6 GDC query parameters, file and workflow types, annotation, units, dates | Done | `revision/manifests/gdc_queries/`, `provenance.md`, `methods_parameters.json` |
| 1.7 Intermediate outputs | Done | `revision/results/` |
| 1.8 Figure-to-script map | Done | `revision/figure_script_map.csv` |
| README does not match scripts | Done | `README.md` (every script: reads, writes, figure/table); `legacy/README.md` |
| Access dates (Nov 12, 2026) | Done (evidence); text next phase | Documented vs inferred dates in `revision/manifests/provenance.md` |

### 2. Genome-wide screen

| Item | Status | Evidence |
|---|---|---|
| 2.1 Gene universe defined first | Done | Registered plan (https://osf.io/qpd2f/); `screen_gene_universe.csv` |
| 2.2 Filters for low expression, zeros, duplicates, IDs | Done | `screen_gene_universe.csv`, `results/expression/tcga_gene_id_handling.csv` |
| 2.3 Spearman, same samples, every gene | Done | `screen_gene_level.csv` (n 252 for every gene); T1, T2, T7 |
| 2.4 BH over all p-values | Done | `a11_screen.py`; `q_adj` |
| 2.5 Ranked table with required columns | Done | S3 Table |
| 2.6 Meaning of "rank first" | Drafted (author files) | Signed and absolute ranks both in S3 Table |
| 2.7 TPM/FPKM and outlier values | Done | `screen_sensitivity.csv` |
| 2.8 Scatterplot from ranking data | Done | Fig 2C |
| 2.9 Revise or remove unreproduced claims | Drafted (author files) | `claims_not_reproduced.md` (C01-C03) |

### 3. Disease cohort

| Item | Status | Evidence |
|---|---|---|
| 3.1 Diagnosis and classification in the source | Done (evidence); text next phase | `tcga_samples.csv` (`histology_ceccarelli2016`, `gdc_primary_diagnosis`) |
| 3.2 IDH1/IDH2 and diagnostic variables | Done | `tcga_samples.csv`; `idh_source_concordance.csv` |
| 3.3 Primary analysis in IDH-wildtype | Done | Primary cohort |
| 3.4 Sensitivity excluding IDH-mutant | Done | IDH-not-mutant and sensitivity cohorts |
| 3.5 No IDH1 expression as status | Done | Mutation calls only |
| 3.6 Exclusions per stage | Done | `tcga_cohort_flow.csv`, `external_cohort_flow.csv`, Fig 1 |
| 3.7 Current terminology | Drafted (author files) | — |

### 4. Survival

| Item | Status | Evidence |
|---|---|---|
| KM 1-2 Time and event; real log-rank | Done | `a21_km.py`; T6; a91 |
| KM 3 n, deaths, censored, follow-up, medians with CI, at risk | Done | `km_groups.csv`, `km_logrank.csv`, `km_at_risk.csv` (median follow-up not reached in MGMT-high groups, written `inf`) |
| KM 4 At-risk table under every KM figure | Done | Fig 3, S3 Fig |
| KM 5-6 Time combination; missing/implausible values | Done (evidence); Methods next phase | `survival_audit.csv` |
| KM 7-8 Median split exploratory; Cox primary | Drafted (author files) | Plan Section 6 |
| Cox 1 Penalization stated | Done (penalizer 0, T9); Methods next phase | `methods_parameters.json` |
| Cox 2 Complete output | Done | `cox_models.csv`; Table 3; S4 Table |
| Cox 3 PH diagnostics | Done | `cox_ph_tests.csv` (R `cox.zph` for every fit; T5 failed comparison with lifelines, R used as prescribed); S5 Fig |
| Cox 4 No "independent prognostic biomarker" | Drafted (author files) | — |
| Cox 5 Clinical and molecular covariates | Done | Models 4-7; extent of resection unavailable (`resection_field_check.csv`) |
| Cox 6 State unavailable covariates | Drafted (author files) | Evidence above |
| Cox 7 Nested models with/without MGMT expression and methylation | Done | Models 4-7, `cox_lr_tests.csv` |
| Cox 8 Remove "suppressor variable" | Drafted (author files) | — |

### 5. Differential expression

| Item | Status | Evidence |
|---|---|---|
| Approach | Done: reperformed with DESeq2 on counts; no ComBat-seq | `a59`-`a61`; the manuscript must drop the ComBat-seq claim (next phase) |
| 5.1 Count-scale input | Done | Toil expected counts, back-transformed and rounded |
| 5.2 TCGA-GTEx confounding | Done (same-pipeline data; TCGA-normals check); limitation next phase | `de_summary.csv` |
| 5.3 Normals: number and origin | Done | `de_summary.csv` |
| 5.4 Identifier harmonization | Done | `de_gene_id_handling.csv` |
| 5.5-5.7 Adjusted p, effect and uncertainty, full table | Done | `deseq2_*.csv`; S8 Table |
| 5.8 Volcano from final output | Done | S7 Fig |
| 5.9 Reconcile gene counts | Done | `de_summary.csv` |
| 5.10 Remove "most significantly upregulated" | Drafted (author files) | C22 |

### 6. Random forest

| Item | Status | Evidence |
|---|---|---|
| Censoring-aware model | Done | Random survival forest (`a50_survival_forest.py`) |
| 6.1-6.5 Outcome, predictors, repeated CV, distributions, baselines | Done | `forest_*.csv`; S8 Fig |
| 6.6 Class balance, confusion matrix | Not applicable | No classification; n and events reported |
| 6.7 Held-out permutation importance | Done | `forest_permutation_importance.csv` |
| 6.8-6.10 Descriptive importance; wording; AUC 0.634 | Drafted (author files) | C12, C13 |

### 7. Enrichment

| Item | Status | Evidence |
|---|---|---|
| 7.1-7.4 Methods, lists, background, ranked list | Done | `results/enrichment/` |
| 7.5 Database versions and dates | Done | `methods_parameters.json`, `DOWNLOAD_LOG.csv` |
| 7.6 Complete outputs | Done | `results/enrichment/` |
| 7.7-7.8 Adjusted p in tables and figures | Done; text next phase | Fig 6, S6 Fig, Table 6 |
| 7.9 List derivation and direction | Done | `ora_lists.csv` |
| 7.10-7.11 Circularity; cautious wording | Drafted (author files) | — |

### 8. External validation

| Item | Status | Evidence |
|---|---|---|
| 8.1-8.5 | Done | `external_cohort_samples.csv`, `external_probe_combinations.csv`, `external_coexpression.csv`, `gse13041_*.csv` |
| 8.6 IDH and classification | Done where data exist | CGGA, CPTAC; others lack IDH (limitation, next phase) |
| 8.7-8.8 Negative results; "mixed" wording | Drafted (author files) | `external_meta.csv` |

### 9. Figures

| Item | Status | Evidence |
|---|---|---|
| 9.1-9.6 Mismatched figures and README | Done | All figures regenerated by `f10`/`f11` from results; `README.md`, `legacy/README.md` |
| Figure requirements (legend match, sample sets, labels, resolution, CIs, at-risk, no manual curves) | Done | PLOS specifications enforced in `plotstyle.py` (300 dpi, RGB, size limits, 8-12 pt, Arial, checked per text element); source data per figure. Legends are next phase |
| PCA clusters neutral | Done | S2 Fig "Cluster 1-3" |
| NAAS check | Not performed | Main figures should be run through NAAS before upload (submission phase) |

### 10. Literature

| Item | Status | Evidence |
|---|---|---|
| PubMed search through revision date | Done | `revision/literature/` (2026-11-01) |
| Fève 2014 and other named references | Done (verified) | `citation_audit.csv` R1-R8 |
| Line-by-line citation audit | Done | All references checked against Crossref and PubMed; statements about GPR153 checked against abstract or full text; Lee 2010 could not be accessed, so the Gli1-target statement was removed (`revision/literature/literature_search_summary.md`, Section 6) |
| Avoid novelty claims; give search date | Drafted (author files) | Search found no GPR153 + glioblastoma/glioma/MGMT record in PubMed |
| Search for GPR153 modulators | Done | No agonist, antagonist or ligand found (searches 11-12) |

### 11-13. Introduction, Discussion, title/Abstract/Conclusion

| Item | Status | Evidence |
|---|---|---|
| 11.1-11.7 Introduction content | Drafted (author files) | Literature evidence file |
| 12.1-12.10 Discussion and limitations | Drafted (author files) | `numbers_for_manuscript.md` |
| 13 Title, Abstract, Conclusion | Drafted (author files) | Outcome branch d (`results/decision/decision.csv`) |

### 14. Reporting and transparency

| Item | Status | Evidence |
|---|---|---|
| 14.1 Dates, titles, numbers, versions | Drafted (author files) | `numbers_for_manuscript.md`, `provenance.md` |
| 14.2 Identical title everywhere | Drafted (author files); submission form pending | — |
| 14.3-14.4 Abbreviations, grammar | Drafted (author files) | — |
| 14.5-14.6 Exact p, CIs | Done (values); text next phase | `numbers_for_manuscript.md` |
| 14.7 Data Availability with code version | Drafted (author files) | `data_availability_draft.md` (analysis version, archive commit, reserved DOI, licences) |
| 14.8 Acknowledged contributor | Done | Acknowledgment removed at the author's instruction |
| 14.9 AI and automated-tool declaration | **Not supplied** | Not provided by the author |
| Biostatistical consultation | Not pursued | Author decision; stated in the response letter |

## D. PLOS administrative note

| Item | Status |
|---|---|
| AI use in preparing the decision letter; contact with the journal office | Journal-office matter; no author action unless instructed |
