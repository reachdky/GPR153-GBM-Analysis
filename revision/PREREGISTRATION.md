# Analysis plan: GPR153 expression, MGMT co-expression and survival in glioblastoma (revised analysis)

Version 1.5, 2026-10-27. Manuscript: PLOS ONE PONE-D-26-34494 (revision; revised
PLOS deadline 2026-11-01).
The registration record (OSF link, timestamp, DOI) is kept in
`revision/PREREGISTRATION_REGISTRATION.md`, outside this file. The SHA-256 of
this exact file, as uploaded, is recorded in `revision/PREREGISTRATION.sha256`.

---

## Summary (approval page)

**Status.** Exploratory. The authors have seen earlier results on overlapping
data. This plan fixes definitions and decision rules before revised values are
computed. The manuscript calls it a time-stamped analysis plan, not a
preregistration.

**Primary hypothesis.** In IDH-wildtype primary glioblastoma, GPR153 is the gene
with the strongest positive co-expression with MGMT genome-wide.

**Primary analysis.**
- TCGA-GBM primary tumors, one aliquot per patient, IDH-wildtype: n = 252
  (BI 146, GDC-2024 106).
- Batch-adjusted Spearman of MGMT with every eligible gene (TPM >= 1 in >= 20%),
  from log2(TPM + 1); BH correction. The pooled Spearman is reported with equal
  prominence.

**Pass rule** (all four required; branches in 5.6):
1. GPR153 is rank 1 with q < 0.05.
2. GPR153-MGMT is positive with q < 0.05 after adjusting for ESTIMATE score, and
   separately for a neuronal signature.
3. GPR153-MGMT is positive in the BI-only and the GDC-2024-only subsets (sign
   check only).
4. GPR153-MGMT adjusted for MGMT promoter methylation, platform and batch is
   positive with a 95% CI excluding 0.

**Expected outcome.** The original code's Pearson top-100 MGMT list did not
contain GPR153 (listed genes r 0.425-0.548; GPR153 r 0.3645). The authors
therefore expect condition 1 to fail and the outcome to be branch (b) or (d).

**Primary survival estimate.** GPR153 per SD in a Cox model with age and sex
(Model 3), unpenalized, in the IDH-wildtype primary cohort; described at most as
an age- and sex-adjusted association.

**Tiers** (Section 15).
- Tier 1: cohort, screen, survival, external co-expression, enrichment,
  reproducibility deliverables.
- Tier 2: decision-rule adjustments, permutation null, survival replication,
  multiple imputation, CAMERA, checklists; optional MGMT-ranked enrichment.
- Tier 3: differential expression, survival forest.

Drop order: MGMT-ranked enrichment, DE, forest, multiple imputation, CAMERA.
Each dropped analysis is logged in DEVIATIONS.md. Scope rule and independent
reproduction rule in 13.2.

**Seen before this plan** (Section 1):
- submitted manuscript values
- re-runs of the original scripts on the original 391 aliquots
- earlier external values (GSE13041, GSE4290, CGGA), including same-session
  all-platform and GPL570 GSE13041 values and an all-histology GSE4290 value
- earlier enrichment
- the original code's own Pearson top-100 MGMT list (GPR153 not in it)
- clinical and molecular annotations
- sample-pair duplicate checks
- three raw MGMT intensities in GSE16011 (excluded)

The original-code re-runs computed GPR153-MGMT statistics on all 391 aliquots,
which include the revised cohort's samples. No value has been computed by the
revised pipeline.

---

## 1. Status and disclosure

This is a time-stamped analysis plan for a re-analysis, not a confirmatory
study, and the manuscript refers to it that way (not as a preregistration).
The hypothesis arose from an earlier exploratory analysis whose code could not
be reproduced: the genome-wide MGMT screen reported in the submission (GPR153
rank 1, BH-adjusted) has no code or output. The only MGMT screen in the original
code is a Pearson top-100 list without multiple-testing correction, and GPR153
is not in that list (table below).

The authors have already seen results from the earlier analysis of overlapping
data (TCGA-GBM, GEO GSE13041 and GSE4290, CGGA). All analyses are therefore
exploratory. The plan fixes definitions, inclusion rules and decision rules
before the revised values are computed.

Earlier results: manuscript-reported values beside the values traced from the
original code and outputs. Full list: `revision/results/claims_trace_v1.csv`.
The claims ledger (Section 13) adds the revised value.

Paths in the last column are relative to the original project directory
(`~/gpr153_project_orig/`, read-only). The scripts under `scripts/analysis/` and
`scripts/cox_regression_v2.py` are byte-identical to the copies in `legacy/` of
this repository; `scripts/cox_regression.py` is not in `legacy/`. Script line
numbers refer to those files. In output files, line 1 is the header.

| Analysis | Manuscript-reported | Traced value | Source of the traced value (file: line) |
|---|---|---|---|
| Genome-wide screen | GPR153 rank 1 of 59,569 genes; BH q < 0.05 | No screen with BH adjustment exists. The original code's only MGMT screen is a Pearson top-100 list with no multiple-testing correction; GPR153 is not in it (listed genes r 0.425 to 0.548; GPR153 r 0.3645) | `scripts/analysis/09_correlation_analysis.py`: `find_correlated_genes()` lines 69-106 (Pearson r at line 88, keeps \|r\| >= 0.3 at line 90, sorts by \|r\| at line 103, keeps the first 100 at line 106). GPR153 list: call at line 148, saved at line 157 to `data/results/correlation/gpr153_correlated_genes.csv`. MGMT list: call at line 167, saved at line 176 to `data/results/correlation/mgmt_correlated_genes.csv` (lines 2-101; no GPR153 row) |
| GPR153-MGMT correlation, TCGA | "Spearman" r = 0.364, p = 9.96e-14, n = 391 | Pearson r = 0.3645, p = 9.96e-14; Spearman rho = 0.328, p = 2.84e-11 | `scripts/analysis/09_correlation_analysis.py`: call at line 127; Pearson at line 54, Spearman at line 57 (both computed by the script); saved at line 138 to `data/results/correlation/gpr153_mgmt_correlation.csv` line 2. Input: `data/processed/combined_gbm_normal_expression.csv` (391 aliquots = 372 primary, 14 recurrent, 5 solid-tissue normal; 293 patients) |
| GPR153 survival, median split | log-rank p = 0.002 | Log-rank p = 0.0795 (n = 287). The script's own test: Mann-Whitney U p = 0.0666 | Mann-Whitney U on survival days, ignoring censoring: `scripts/analysis/10_survival_analysis.py` line 121; saved at line 205 to `data/results/survival/gpr153_survival_stats.csv` line 2. Log-rank computed during tracing (lifelines `logrank_test`; no script in the original project) on `data/results/survival/gpr153_survival_data.csv` (written at line 210; lines 2-288). The reported 0.002 equals the MGMT log-rank value (0.00219). `scripts/cox_regression.py` line 225 contains "Univariate log-rank p=0.002 (pre-existing canonical result) is unchanged."; `scripts/cox_regression_v2.py` lines 205-206 print "Univariate log-rank p=0.002 remains the primary survival finding" |
| MGMT survival, median split | log-rank p = 0.014 | Mann-Whitney U p = 0.0144; log-rank p = 0.00219 (n = 287) | Mann-Whitney U: `scripts/analysis/10_survival_analysis.py` line 121; saved at line 236 to `data/results/survival/mgmt_survival_stats.csv` line 2. Log-rank computed during tracing on `data/results/survival/mgmt_survival_data.csv` (written at line 240) |
| Cox, GPR153 per SD | univariate HR ≈ 1.08, p = 0.169; age-adjusted HR 1.160 [1.018-1.322], p = 0.026; PH met (p = 0.77) | univariate HR 1.095, p = 0.169; age-adjusted HR 1.160, p = 0.026 (n = 286, penalizer 0.1); with MGMT and age HR 1.109, p = 0.137; Schoenfeld GPR153 p = 0.766, age p = 0.034 | `scripts/cox_regression_v2.py`: GPR153 = mean of all aliquots per patient, including recurrent and normal (line 37); z-scores (line 64); Models A and B (lines 72-73); `CoxPHFitter(penalizer=0.1)` (line 91); Schoenfeld test with `time_transform='rank'` (line 116; p-values printed to the console, no file); results saved at line 128 to `data/results/survival/cox_v2_results.csv` lines 2-4. Model with MGMT: `scripts/cox_regression.py` (penalizer 0.1 at line 105; saved at line 149 to `data/results/survival/cox_multivariate_results.csv`, GPR153 in model C at line 5) |
| Classifier | AUC = 0.634 (80/20 split); GPR153 first by Gini importance (0.31), above TERT | AUC 0.506 (70/30 split, 365-day label); Gini: TP53 0.233, GPR153 0.193; TERT not in model; 0.634 not found in any output | `scripts/analysis/14_predictive_modeling.py`: 365-day label (line 161), `test_size=0.3` (line 168), Gini importance (line 100); saved at line 195 to `data/results/ml/rf_feature_importance.csv` lines 2-3, and at line 233 to `data/results/ml/model_comparison.csv` line 2 |
| External, GSE13041 | +0.119 (n = 191); +0.116 (n = 163) | Pearson +0.119 (n = 191, p = 0.101); +0.116 (n = 163, p = 0.142). Same session, all platforms pooled: -0.285 (n = 218 = 191 GPL96 + 27 GPL570 arrays); GPL570 only: -0.360 (n = 27); "Primary": -0.421 (n = 11, the number of GPL570 "Primary GBM" arrays) | `data/results/validation/complete_validation_picture.csv` lines 3 (n = 191) and 4 (n = 163); `data/results/validation/cross_dataset_validation.csv` line 4 (n = 218); `data/results/validation/cross_dataset_validation_filtered.csv` lines 4 (n = 27) and 5 (n = 11). Files dated 2026-03-16; generating code not on the machine; they store N, r and p only, no sample list |
| External, GSE4290 | -0.134 (n = 77) | Pearson -0.134 (n = 77, p = 0.246). Same session, all histologies: -0.081 (n = 180, p = 0.280) | `data/results/validation/cross_dataset_validation_filtered.csv` line 3 (n = 77); `data/results/validation/cross_dataset_validation.csv` line 3 (n = 180) (2026-03-16; code not on the machine) |
| External, CGGA | -0.293 (n = 14, "Grade IV") | Pearson -0.293 (n = 14, p = 0.310) | `data/results/validation/complete_validation_picture.csv` line 2 ("CGGA (G4 only)"; code not on the machine). Recomputed during tracing (no script): CGGA_IDH_A cohort (IDH-mutant astrocytoma) with Grade_2021 = G4, untransformed RSEM |
| Combined external n | 445 | 191 + 163 + 77 + 14 | Sum of the four N values above. The 163 arrays are a subset of the 191 (Section 7.1), so 445 counts them twice |
| Enrichment | Reactome Signaling By WNT adj. p = 0.028 | adj. p = 0.0280, combined score 54.12 | `scripts/analysis/11_pathway_analysis.py`: input `data/results/correlation/gpr153_correlated_genes.csv` (line 23; the GPR153 Pearson list from `09_correlation_analysis.py`); Enrichr query, default background (line 197); saved at line 209 to `data/results/pathways/reactome_pathways.csv` line 2 |
| GSEA | WNT ranks 1st, DNA repair 4th | No GSEA output exists | none |
| Differential expression | DESeq2 + ComBat-seq; GPR153 "most significantly upregulated", log2FC ≈ 3.8 | Student's t-test (equal variances) on the expression values; log2 of the ratio of group means; GPR153 log2FC 1.13 vs 40 GTEx normals | `scripts/analysis/12_differential_expression.py`: input line 23, normals line 36, group means lines 63-64, ratio line 69, `stats.ttest_ind` with default equal variances line 73; saved at line 161 to `data/results/differential/differential_expression_all.csv` line 14482 (GPR153) |
| PCA | PC1 23.4%, PC2 14.7%, n = 115, with subtype labels | PC1 70.8%, PC2 5.3%, n = 391; labels assigned by a fixed mapping | `scripts/analysis/13_clustering.py`: PCA line 55; saved at line 139 to `data/results/ml/pca_variance.csv` lines 2-3. Subtype names: fixed mapping of cluster numbers in `scripts/export_dashboard_data.py` line 112 |

The following were all chosen after these results were seen:
- the six genes in Section 11
- the 20% expression filter in Section 4
- MGMT as the screening target
- GPL96 primary arrays only for GSE13041 (the same-session outputs also held
  the all-platform value, -0.285, n = 218)
- glioblastoma only for GSE4290 (the same-session outputs also held the
  all-histology value, -0.081, n = 180)

The two external choices are therefore paired with prespecified sensitivity
analyses (Section 7.3): GSE13041 per platform, pooled only after a
cross-platform overlap check, and GSE4290 with all histologies. All are
reported in Table 4.

### 1.1 Work completed before registration

The original-code re-runs computed GPR153-MGMT statistics on all 391 aliquots,
which include the revised cohort's samples. No value has been computed by the
revised pipeline. The one other exception is marked "Exception".

- Downloads of all sources in Section 2, with URLs, dates and md5 checksums.
  Installation of R 4.6.1 packages (DESeq2 1.52.0, limma 3.68.5, apeglm 1.34.0,
  mice 3.19.0; Bioconductor 3.23).
- GDC file metadata for the 391 TCGA expression files (sample type, batch,
  workflow, md5 checks against GDC). GDC metadata for the TCGA-GBM methylation
  files (platform per patient).
- Clinical data, IDH calls from MAFs, and cBioPortal annotations:
  - Ceccarelli 2016: histology, MGMT promoter status, TERT promoter status,
    chr7 gain / chr10 loss.
  - EGFR GISTIC calls.
  - TCGA cohort construction, cohort flow, batch comparison, and annotation
    tables (`revision/results/cohort/`).
- External-cohort selection from sample annotations. Survival-replication
  candidates checked from clinical tables, platform annotation tables and row
  identifiers only (`revision/results/cohort/survival_replication_candidates.csv`).
- Survival-data audit of every survival cohort from clinical tables only
  (Section 6.7).
- CGGA 325 vs 693 clinical-record match (Section 7.2).
- Extent-of-resection field check in GDC, TCGA-CDR, PanCanAtlas clinical data and
  the BCR clinical XML (field names and non-missing counts only; Section 3.1).
- Whole-profile duplicate checks in GSE13041 GPL96 (primary-only, and all 191
  arrays), GSE4290 and CGGA. These are sample-by-sample Spearman correlations
  across all probes or genes except the GPR153 and MGMT probes and genes. Only
  sample-pair correlations were inspected.
- Re-runs of the original scripts on the original data (391 aliquots), to trace
  the submitted values (Section 1 table).
- A known-answer check of the log-rank implementation on the R `survival::aml`
  data (Section 6.5).
- Verification of the MGMT-STP27 model from its publication and of the probe
  annotation (Section 3.2).
- Cell-line and drug metadata: DepMap Model.csv and OmicsProfiles.csv; the
  DepMap expression file's header and row identifiers; GDSC compound lists, and
  dose-response fields for temozolomide (concentration range, IC50 relative to
  it).
- Single-cell metadata: GBmap cell annotations, donor and stage fields, and the
  gene index. No expression values were read.
- Gene-set files: set names and counts only. ESTIMATE stromal and immune gene
  sets and the common-gene list exported from the tidyestimate 1.1.1 package data
  to `data/downloads/estimate/extracted/` (2026-10-26).
- Python packages installed 2026-10-26: pingouin 0.5.5, gseapy 1.3.1, pydeseq2
  0.5.2 (in `.venv`), and scikit-survival 0.28.0 in a separate `.venv_rsf`
  (it requires numpy 2; the main environment keeps numpy 1.26.4).
- The analysis scripts for Sections 3.2, 5-8 and 10-14 (`pipeline/a10`-`a91`,
  `f10`-`f12`) were written and run end to end on simulated expression and
  methylation values (`pipeline/sim/s00_make_sim_data.py`, seed 42; outputs in
  `revision/sim/`, not results). The simulation copies the real clinical tables
  (sample IDs, cohort flags, survival, covariates); it reads no expression or
  methylation value. Every script that reads real expression or methylation
  values stops until the registration record has an OSF URL and this file's
  SHA-256 matches `revision/PREREGISTRATION.sha256`.
- Verification tests T1-T6, T2b, T9, T9b and T10-T15 (Section 13.1) run on
  simulated data, R data sets, clinical tables or source code; T7, T8 and the
  independent reproduction script (13.2, `pipeline/a91_independent_check.py`)
  run on simulated data only. T5 found that lifelines and R `cox.zph` differ
  (lung data: age p 0.723 vs 0.647), so R `cox.zph` is the source of all
  proportional-hazards tests, as T5 prescribes.
- R package metafor 5.2.1 installed 2026-10-27 (independent meta-analysis
  check, 13.2). The lifelines convergence precision was set to 1e-12 (6.1)
  after a check on simulated data: at the default, Model 3 differed from R
  `coxph` by a relative 3.6e-5; at 1e-12, by 5e-12.
- GSE13041 platform tables in the GEO family SOFT file (annotation only):
  GPL8300 has no GPR153 probe set (MGMT: 2051_at, 2052_g_at); GPL96 and GPL570
  have 221902_at, 64942_at and 204880_at.
- GSE13041 clinical table: the characteristics field "gender" (M/F) was not
  read into the manifest before 2026-10-27; it is now (all 162 primary GPL96
  arrays have a value). Cross-platform overlap check from clinical fields
  (7.3): no tumor ID (sample title) and no clinical record (age, sex, survival
  time, vital status) is shared between the 191 GPL96 and the 27 GPL570
  arrays.
- On 2026-10-27, the original outputs
  `data/results/validation/cross_dataset_validation.csv` and
  `cross_dataset_validation_filtered.csv` were read again for the
  same-session values in the Section 1 table.
- On 2026-10-26, while adding line numbers to the Section 1 table, the original
  output `data/results/correlation/mgmt_correlated_genes.csv` (original 391
  aliquots) was read: GPR153 is not among its 100 genes.
- Exception: during a header check of the GSE16011 series matrix on 2026-10-25,
  the first three raw values of the MGMT probe set 4255_at were printed. No
  statistic was computed. GSE16011 is excluded from all analyses (Section 8).

## 2. Data sources (all open access)

Retrieval dates per file are taken from `data/downloads/DOWNLOAD_LOG.csv`
(UTC). The log is copied to `revision/manifests/DOWNLOAD_LOG.csv` and lists
URL, size and md5 for every file.

**Downloaded in December 2025** (by the original project; dates are file
modification times; not re-downloaded):

| Files | Date |
|---|---|
| TCGA-GBM STAR - Counts files, 391 aliquots (`augmented_star_gene_counts.tsv`, GENCODE v36) | 2025-12-24. All 391 md5s match GDC Data Release 46.0 metadata (checked 2026-10-24) |
| GEO GSE13041 and GSE4290 series data and family SOFT files (GSE13041: GPL96, GPL570, GPL8300; GSE4290: GPL570) | 2025-12-25 |

**Downloaded in October 2026** (`pipeline/r00_download_external.py`). All GDC
metadata, MAFs, clinical records, clinical XML and methylation files are from
GDC Data Release 46.0 (2026-08-10).

| File (under data/downloads/) | Content | Retrieved (UTC) | Note |
|---|---|---|---|
| `cbioportal/gbm_tcga_pan_can_atlas_2018_EGFR_gistic.json` | cBioPortal gbm_tcga_pan_can_atlas_2018 EGFR GISTIC calls | 2026-10-25 |  |
| `cbioportal/gbm_tcga_pan_can_atlas_2018_attributes.json` | cBioPortal gbm_tcga_pan_can_atlas_2018 clinical attribute definitions | 2026-10-24 |  |
| `cbioportal/gbm_tcga_pan_can_atlas_2018_patient_clinical.json` | cBioPortal gbm_tcga_pan_can_atlas_2018 patient clinical data | 2026-10-24 |  |
| `cbioportal/gbm_tcga_pan_can_atlas_2018_sample_clinical.json` | cBioPortal gbm_tcga_pan_can_atlas_2018 sample clinical data | 2026-10-24 |  |
| `cbioportal/gbm_tcga_pub2013_EGFR_gistic.json` | cBioPortal gbm_tcga_pub2013 EGFR GISTIC calls | 2026-10-25 |  |
| `cbioportal/gbm_tcga_pub2013_attributes.json` | cBioPortal gbm_tcga_pub2013 clinical attribute definitions | 2026-10-24 |  |
| `cbioportal/gbm_tcga_pub2013_patient_clinical.json` | cBioPortal gbm_tcga_pub2013 patient clinical data | 2026-10-24 |  |
| `cbioportal/gbm_tcga_pub2013_sample_clinical.json` | cBioPortal gbm_tcga_pub2013 sample clinical data | 2026-10-24 |  |
| `cbioportal/lgggbm_tcga_pub_attributes.json` | cBioPortal lgggbm_tcga_pub clinical attribute definitions | 2026-10-24 |  |
| `cbioportal/lgggbm_tcga_pub_patient_clinical.json` | cBioPortal lgggbm_tcga_pub patient clinical data | 2026-10-24 |  |
| `cbioportal/lgggbm_tcga_pub_sample_clinical.json` | cBioPortal lgggbm_tcga_pub sample clinical data | 2026-10-24 |  |
| `cellxgene/core_GBmap_861acfd8-25f0-418b-a445-aa96da232827.h5ad` | Core GBmap H5AD, dataset version 861acfd8-25f0-418b-a445-aa96da232827 | 2026-10-25 | not used (Section 9 removed in v1.5) |
| `cellxgene/gbmap_collection.json` | GBmap collection metadata | 2026-10-25 | not used (Section 9 removed in v1.5) |
| `cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt.zip` | CGGA mRNAseq_325 RSEM FPKM | 2026-10-24 |  |
| `cgga/CGGA.mRNAseq_325_clinical.20200506.txt.zip` | CGGA mRNAseq_325 clinical | 2026-10-24 |  |
| `cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt.zip` | CGGA mRNAseq_693 RSEM FPKM | 2026-10-24 |  |
| `cgga/CGGA.mRNAseq_693_clinical.20200506.txt.zip` | CGGA mRNAseq_693 clinical | 2026-10-24 |  |
| `cptac/cptac_brain_cases.json` | GDC CPTAC-3 brain case clinical records | 2026-10-25 | 220 records; Data Release 46.0 - August 10, 2026 |
| `cptac/expression/files_meta.json` | GDC query and metadata: CPTAC-3 brain STAR counts | 2026-10-25 | 238 records; Data Release 46.0 - August 10, 2026 |
| `cptac/maf/files_meta.json` | GDC query and metadata: CPTAC-3 brain open masked MAFs | 2026-10-25 | 269 records; Data Release 46.0 - August 10, 2026 |
| `depmap/Model.csv` | DepMap 24Q4 Public: Model.csv | 2026-10-25 | not used (Section 9 removed in v1.5) |
| `depmap/OmicsExpressionProteinCodingGenesTPMLogp1.csv` | DepMap 24Q4 Public: OmicsExpressionProteinCodingGenesTPMLogp1.csv | 2026-10-25 | not used (Section 9 removed in v1.5) |
| `depmap/OmicsProfiles.csv` | DepMap 24Q4 Public: OmicsProfiles.csv | 2026-10-25 | not used (Section 9 removed in v1.5) |
| `depmap/README.txt` | DepMap 24Q4 Public: README.txt | 2026-10-25 | not used (Section 9 removed in v1.5) |
| `depmap/figshare_article.json` | DepMap figshare record: DepMap 24Q4 Public, 10.25452/figshare.plus.27993248.v1 | 2026-10-25 | not used (Section 9 removed in v1.5) |
| `estimate/tidyestimate_1.1.1.tar.gz` | ESTIMATE stromal/immune gene sets (tidyestimate 1.1.1, CRAN) | 2026-10-24 |  |
| `gdc_bcr_xml/files_meta.json` | GDC query and metadata: TCGA-GBM clinical supplements (BCR XML) | 2026-10-26 | 599 records; Data Release 46.0 - August 10, 2026 |
| `gdc_maf/maf_files_meta.json` | GDC query and metadata for TCGA-GBM open masked MAFs | 2026-10-24 | 464 records; Data Release 46.0 - August 10, 2026 |
| `gdc_meta/expression_files_meta.json` | GDC metadata for the TCGA-GBM STAR count files on disk | 2026-10-24 | 391 records; Data Release 46.0 - August 10, 2026 |
| `gdc_meta/tcga_gbm_cases.json` | GDC TCGA-GBM case clinical records | 2026-10-24 | 617 records; Data Release 46.0 - August 10, 2026 |
| `gdc_methylation/files_meta.json` | GDC query and metadata: TCGA-GBM methylation beta values (27K, 450K) | 2026-10-25 | 450 records; Data Release 46.0 - August 10, 2026 |
| `gdsc/Cell_Lines_Details.xlsx` | GDSC cell line details, release 8.5 | 2026-10-25 | not used (Section 9 removed in v1.5) |
| `gdsc/GDSC1_fitted_dose_response_27Oct23.xlsx` | GDSC1 fitted dose response (LN_IC50, AUC), release 8.5 | 2026-10-25 | not used (Section 9 removed in v1.5) |
| `gdsc/GDSC2_fitted_dose_response_27Oct23.xlsx` | GDSC2 fitted dose response (LN_IC50, AUC), release 8.5 | 2026-10-25 | not used (Section 9 removed in v1.5) |
| `gdsc/screened_compounds_rel_8.5.csv` | GDSC release 8.5 screened compounds | 2026-10-25 | not used (Section 9 removed in v1.5) |
| `genesets/GO_Biological_Process_2023.gmt` | Enrichr library GO_Biological_Process_2023 | 2026-10-25 |  |
| `genesets/KEGG_2021_Human.gmt` | Enrichr library KEGG_2021_Human | 2026-10-25 |  |
| `genesets/Reactome_2022.gmt` | Enrichr library Reactome_2022 | 2026-10-25 |  |
| `genesets/c2.cp.reactome.v2026.2.Hs.symbols.gmt` | MSigDB 2026.2.Hs Reactome gene sets | 2026-10-25 |  |
| `genesets/h.all.v2026.2.Hs.symbols.gmt` | MSigDB 2026.2.Hs Hallmark gene sets | 2026-10-25 |  |
| `geo_replication/GSE108474_REMBRANDT_GeneExpression.txt.gz` | Rembrandt GPL570 expression | 2026-10-25 |  |
| `geo_replication/GSE108474_REMBRANDT_biospecimen_mapping_GEO.txt.gz` | Rembrandt subject to GSM mapping | 2026-10-25 |  |
| `geo_replication/GSE108474_REMBRANDT_clinical.data.txt.gz` | Rembrandt clinical data | 2026-10-25 |  |
| `geo_replication/GSE108474_REMBRANDT_data.dictionary.xlsx` | Rembrandt clinical data dictionary | 2026-10-25 |  |
| `geo_replication/GSE16011_series_matrix.txt.gz` | GSE16011 series matrix | 2026-10-25 |  |
| `geo_replication/GSE68848_series_matrix.txt.gz` | GSE68848 (Rembrandt subset) series matrix | 2026-10-25 |  |
| `geo_sc/GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx` | GSE131928 (Neftel 2019) cell-to-tumor annotation | 2026-10-25 | not used (Section 9 removed in v1.5) |
| `literature/PMC3444709.xml` | Bady 2012 Acta Neuropathol 124:547, full text | 2026-10-25 |  |
| `literature/PMC3444709_M1.gif` | Bady 2012 M-value equation | 2026-10-25 |  |
| `literature/PMC3444709_M2.gif` | Bady 2012 MGMT-STP27 model equation | 2026-10-25 |  |
| `literature/PMC4079675.xml` | Bady 2012 erratum (Acta Neuropathol 2013) | 2026-10-25 |  |
| `platforms/GPL13534_table.txt` | Illumina HumanMethylation450 annotation | 2026-10-25 |  |
| `platforms/GPL570_table.txt` | Affymetrix HG-U133 Plus 2.0 annotation | 2026-10-25 |  |
| `platforms/GPL8490_table.txt` | Illumina HumanMethylation27 annotation | 2026-10-25 |  |
| `platforms/GPL8542_table.txt` | GSE16011 platform (Brainarray HGU133Plus2_Hs_ENTREZG v10) annotation | 2026-10-25 |  |
| `tcga_pub/TCGA-CDR-SupplementalTableS1.xlsx` | TCGA Pan-Cancer Clinical Data Resource (Liu 2018, Cell 173:400), Table S1 | 2026-10-26 |  |
| `tcga_pub/TCGA_GBM_dnameth_call_categories.txt` | Brennan 2013 publication page: DNA methylation call categories | 2026-10-26 |  |
| `tcga_pub/clinical_PANCAN_patient_with_followup.tsv` | PanCanAtlas merged BCR clinical table, patient with follow-up | 2026-10-26 |  |
| `tcga_pub/lgggbm_2016_AWG_clinical_supplementary_info` | Ceccarelli 2016 publication page: AWG Clinical Supplementary Info | 2026-10-26 |  |
| `tcga_pub/lgggbm_2016_GBM_clinical_data` | Ceccarelli 2016 publication page: GBM Clinical data | 2026-10-26 |  |
| `toil/TcgaTargetGTEX_phenotype.txt.gz` | UCSC Toil sample phenotype (study, sample type, GTEx tissue) | 2026-10-24 |  |
| `toil/TcgaTargetGtex_gene_expected_count.gz` | UCSC Toil RSEM expected counts, log2(x+1), GENCODE v23 | 2026-10-24 |  |
| `toil/gencode.v23.annotation.gene.probemap` | UCSC Toil gene ID map, GENCODE v23 | 2026-10-24 |  |

Files fetched through a GDC query are listed by their query-metadata file. The
note gives the number of records and the GDC data release. Each downloaded file
is md5-checked against GDC at download:
- TCGA-GBM open MAFs: 464 files
- TCGA-GBM methylation beta files: 450
- TCGA-GBM clinical supplements (BCR XML): 599
- CPTAC-3 brain STAR counts: 238
- CPTAC-3 brain open MAFs: 269

Query parameters are saved as JSON in `revision/manifests/gdc_queries/`.

## 3. TCGA cohort

Step-wise inclusion, with counts reported at each step (Fig 1 and
`revision/manifests/tcga_samples.csv`):

1. All TCGA-GBM STAR - Counts files on disk (391 aliquots, 293 patients).
2. Sample type 01 (primary tumor); 372 aliquots, 284 patients. Recurrent (02,
   n = 14) and solid-tissue normal (11, n = 5) samples excluded.
3. One aliquot per patient, chosen by a fixed rule:
   (a) prefer the original Broad Institute batch (aliquot barcode center code 01)
   over the 2024 GDC addition (center code 41) when a patient has both;
   (b) then lowest vial letter, lowest portion number, lowest plate, and
   lexicographic order of the full aliquot barcode. Result: 284 patients
   (BI 155, GDC-2024 129).
4. Overall survival: as defined in the survival-data audit (Section 6.7).
   Patients with unmapped vital status, missing time, or time <= 0 are excluded
   from survival analyses but not from the screen.
5. IDH status (Section 3.1). Three cohorts:
   - **Primary cohort**: IDH-wildtype, n = 252 (survival-evaluable 249, 207
     deaths).
   - **IDH-not-mutant cohort**: wildtype plus unknown, n = 264, excluding only
     known IDH-mutant tumors (sensitivity).
   - **Sensitivity cohort**: all primary tumors, n = 284 (survival-evaluable 281).

An automated check stops the pipeline if any patient ID occurs more than once
in any analysis cohort (`pipeline/r01_tcga_cohort.py`).

Batch: each aliquot is labeled BI (center code 01) or GDC-2024 (center code
41). Both batches use GDC STAR - Counts, GENCODE v36. In the primary cohort,
the batches do not differ significantly on age, sex, KPS, MGMT promoter status,
treatment or OS (log-rank p = 0.24;
`revision/results/cohort/tcga_batch_comparison.csv`).

**TCGA manifest** (`revision/manifests/tcga_samples.csv`) columns:
- identifiers: aliquot, patient ID, sample ID, specimen type, batch
- file details: file ID and md5
- inclusion status and reason, with flags for every analysis cohort
- survival: time, event and source
- clinical: age, sex, KPS
- IDH (each source, and disagreement), MGMT promoter call and source
- molecular: methylation array platform, EGFR GISTIC, chr7/chr10, TERT promoter,
  G-CIMP, expression subtype
- GDC primary_diagnosis, morphology, residual_disease
- treatment

### 3.1 IDH status and diagnostic variables

- Source 1: GDC open-access masked MAFs. A patient is **IDH-mutant** if any
  tumor aliquot carries IDH1 R132 or IDH2 R172; **IDH-wildtype** if a MAF exists
  and contains neither; **unknown** if no MAF exists.
- Source 2, used only when Source 1 is unknown: Ceccarelli 2016 IDH_STATUS
  (cBioPortal lgggbm_tcga_pub), then Brennan 2013 IDH1_MUTATION
  (gbm_tcga_pub2013).
- Priority when sources disagree: MAF, then Ceccarelli 2016, then Brennan 2013.
  The disagreement is recorded in the manifest. Before registration, 2 patients
  had a disagreement:
  - TCGA-06-0178 (MAF wildtype vs Ceccarelli mutant; in the primary cohort)
  - TCGA-06-0156 (Ceccarelli mutant vs Brennan wildtype; IDH-mutant, sensitivity
    cohort only)

  A sensitivity analysis excludes patients with any disagreement.
- Unknown remains a separate tier. IDH1 expression is not used as a proxy for
  mutation status.
- Diagnosis, tabulated by batch for each cohort
  (`revision/results/cohort/tcga_annotation_tables.csv`). Primary cohort:
  - GDC primary_diagnosis "Glioblastoma", morphology 9440/3: 251; 1 missing
  - Ceccarelli histology "glioblastoma": 251; 1 missing
- Molecular markers reported where available. Primary cohort, n = 252:
  - EGFR amplification (GISTIC = 2; PanCanAtlas, then Brennan 2013): 118
    amplified, 128 not, 6 missing.
  - chr7 gain / chr10 loss (Ceccarelli 2016): 176 combined, 69 not, 7 missing.
  - TERT promoter (Ceccarelli 2016): 29 mutant, 223 missing. Ceccarelli
    reports only mutant calls for this group, so "missing" cannot be read as
    wildtype.
- **Extent of resection: not available in any TCGA-GBM source checked.**
  Checked before registration (`pipeline/r09_resection_field_check.py`,
  `revision/results/cohort/resection_field_check.csv`; field names matched on
  resect / residual / extent / surg / margin):
  - GDC diagnoses `residual_disease`: 0 of 252 primary-cohort patients have a
    value.
  - TCGA Pan-Cancer Clinical Data Resource (Liu 2018, Table S1; 251 of 252
    patients present): `residual_tumor` and `margin_status` empty for all 251.
  - PanCanAtlas merged BCR clinical table (251 present): 35 matching fields,
    including `residual_tumor`, `tumor_residual_disease` and
    `surgical_procedure_name`, all empty for all 251.
  - Original BCR clinical XML (GDC, 599 files; 251 patients present): the
    TCGA-GBM clinical schema has no extent-of-resection element. The only
    matching elements are `nte:new_tumor_event_additional_surgery_procedure`
    (54 patients) and its day count (56). These describe surgery for a new tumor
    event after initial treatment, not the initial resection, and are not used.

  It is stated as a limitation and cannot be adjusted for.
- Other covariates:
  - KPS (Ceccarelli 2016)
  - G-CIMP and expression subtype (Brennan 2013)
  - temozolomide (Brennan 2013 THERAPY)
  - radiotherapy (PanCanAtlas RADIATION_THERAPY)
  - sex and age (GDC). Age = age_at_diagnosis / 365.25, with age_at_index
    where missing. Last follow-up = the maximum of the diagnosis and follow-up
    records.

### 3.2 MGMT promoter methylation

**Measures, verified from source files:**

| Item | Finding | Source |
|---|---|---|
| STP27 probes | cg12434587 and cg12981137. The original abstract misspelled the first as "cg1243587"; the erratum corrects it | Bady 2012 full text; erratum |
| Coefficients and cutoff | logit(y) = 4.3215 + 0.5271 × M(cg12434587) + 0.9265 × M(cg12981137); methylated if y >= 0.358. Equation given as an image in the article | Bady 2012, equation image M2 |
| Original M-value definition | log2((max(signal B, 0) + 1) / (max(signal A, 0) + 1)), from GenomeStudio-normalized intensities | Bady 2012, equation image M1 |
| Probe annotation, 27K (GPL8490) | Both probes present; gene MGMT (GeneID 4255); chr10:131155199 and 131155565 (build 36) | GPL8490 table |
| Probe annotation, 450K (GPL13534) | Both probes present; MGMT; chr10:131265209 (TSS1500) and 131265575 (body); CpG island chr10:131264948-131265710 (build 37); Methyl27_Loci = TRUE | GPL13534 table |
| TCGA-GBM methylation files | Beta values (SeSAMe): 27K, 295 files, 287 patients; 450K, 155 files, 142 patients. Masked IDAT intensities also open: 27K, 590 files; 450K, 310 files | GDC Release 46.0 file metadata |
| Ceccarelli 2016 call | Table S1 (GDC publication page, "AWG Clinical Supplementary Info", sheet S1A) has a column "MGMT promoter status" with no method note | Ceccarelli 2016 Table S1 |

**Computation.**
- MGMT-STP27 (Bady et al. 2012, Acta Neuropathol 124:547, PMC3444709; erratum
  PMC4079675) is applied to M-values computed from GDC SeSAMe beta values as
  M = log2(beta / (1 - beta)), with beta clipped to [0.001, 0.999]. The output
  is the linear predictor (logit) and the binary status.
- If a patient has both arrays, the 450K file is used.
- **STP27 population**: primary-cohort patients whose beta file has both STP27
  probes non-missing. Not yet counted (Appendix A, U2). Upper bound 202:
  - 127 patients have 27K data only
  - 73 have 450K only
  - 2 have both
  - 50 have no methylation file
- Platform by batch in the primary cohort:
  - BI: 27K 71, 450K 47, both 1, none 27
  - GDC-2024: 27K 56, 450K 26, both 1, none 23

  Platform and batch are therefore not collinear.
- **Concordance**: Cohen's kappa between STP27 binary status and the Ceccarelli
  2016 call. **Kappa population**: primary-cohort patients with both an STP27
  call and a Ceccarelli call. Its n and the 2 × 2 table (STP27 methylated /
  unmethylated × Ceccarelli methylated / unmethylated) are reported
  (`revision/results/methylation/methylation_kappa_2x2.csv`). Computed after
  registration: the n requires reading the STP27 CpG values (Appendix A, U2).
  Upper bound 197, the primary-cohort patients with a Ceccarelli call, all of
  whom have a methylation file.
- **Binary status for descriptive analyses**, in fixed order: STP27 status;
  else the Ceccarelli 2016 call; else Brennan 2013 MGMT_STATUS. The number of
  patients under each source is reported. The mean beta of a single available
  STP27 CpG is reported as a continuous descriptive only (no published cutoff).

**Deviation from the original model.** GDC provides SeSAMe-processed beta
values, not GenomeStudio intensities. STP27 applied to beta-derived M-values
approximates the published input (Appendix A, U3). The kappa above and the
platform covariate address this.

**Condition 4 (decision rule).** Fixed before registration:
- **Statistic**: partial Spearman of GPR153 with MGMT expression, adjusted for
  batch, array platform (450K vs 27K) and the methylation measure; k = 3. The
  STP27 linear predictor, like every covariate, enters as ranks (Section 5,
  step 2).
- **Primary methylation measure**: the STP27 linear predictor (continuous), in
  the STP27 population, if both hold:
  - the STP27 population has n >= 100
  - kappa(STP27, Ceccarelli) >= 0.60
- **Fallback**, if either fails: the Ceccarelli 2016 binary call is the
  methylation measure, in primary-cohort patients with a Ceccarelli call (197
  before registration, all with a methylation file). k = 3 (batch, platform,
  call).
- If the fallback population also has n < 100, condition 4 is "not evaluable",
  which counts as not met (branch c, Section 5.6).
- The measure actually used and the reason are reported.
- **CI**: Fisher z with Bonett-Wright SE = sqrt((1 + rho²/2) / (n - 3 - k)), k = 3.
  Condition 4 is met if the partial rho is positive and the 95% CI excludes 0.
- **Sensitivity checks** (do not change the decision):
  - bootstrap percentile 95% CI (1,000 resamples of patients, seed 42,
    residualization repeated in each resample)
  - without the platform covariate (k = 2)
  - binary STP27 status instead of the linear predictor
  - the other measure (Ceccarelli if STP27 was primary, and the reverse)

**Other methylation analyses** (descriptive; primary cohort; batch-adjusted
partial Spearman throughout):
1. GPR153-MGMT expression correlation within methylated and within unmethylated
   patients (binary status).
2. GPR153 expression vs methylation:
   - Mann-Whitney U by binary status
   - batch-adjusted Spearman with the STP27 linear predictor

## 4. Expression data and gene universe

- Units: log2(TPM + 1) from tpm_unstranded, GENCODE v36. Genes are keyed by
  version-stripped Ensembl ID. Symbols are reported but not used as keys;
  symbol collisions are counted and listed. PAR_Y duplicate IDs are dropped.
- Gene universe: all genes with non-constant expression in the primary cohort.
- **Eligible genes**: TPM >= 1 in at least 20% of primary-cohort samples.
  Sensitivity filters: >= 10%, >= 50%, and none. Whether GPR153 and MGMT pass
  each filter is reported.
- Any gene with a missing value is excluded and counted.
- **Gene-level output** covers all non-constant genes. Each row carries an
  eligibility flag for each filter; it is not limited to eligible genes.
- A gene-universe reconciliation table reports counts at each step for each
  analysis:
  - GDC v36 total, non-constant, eligible
  - Toil v23 overlap, DE tested
  - GSEA and CAMERA input

## 5. Primary analysis: genome-wide MGMT co-expression screen

The editor asked for Spearman correlation. Two Spearman statistics are computed
for every gene and reported with **equal prominence** in the abstract, Results,
tables and Fig 2:
- batch-adjusted (partial) Spearman
- pooled (unadjusted) Spearman

Pass/fail and the outcome branch (5.5, 5.6) use the batch-adjusted result. If
GPR153's rank-1 status differs between the two, the abstract states both.

**Computation of the batch-adjusted statistic:**
1. For each gene and for MGMT, log2(TPM + 1) values are ranked across all
   primary-cohort samples (average ranks for ties).
2. Each rank vector is regressed by OLS on an intercept and a batch indicator
   (1 = GDC-2024, 0 = BI); k = 1. A third batch, if present, gets its own
   indicator and k increases. **Covariates enter as ranks**: every covariate
   column (here batch; in 5.2 the ESTIMATE score or the neuronal score; in 3.2
   platform, the STP27 linear predictor or a binary call) is ranked across the
   samples (average ranks) before the regression. This matches
   `pingouin.partial_corr(method="spearman")`, which ranks x, y and all
   covariates (tests T2 and T2b). Ranking a binary indicator only rescales it,
   so binary covariates give the same residuals ranked or not; continuous
   covariates do not.
3. The partial Spearman rho is the Pearson correlation of the MGMT residuals
   with the gene residuals.
4. t = rho × sqrt(df / (1 - rho²)), df = n - 2 - k; two-sided p from the t
   distribution.

**Primary run:**
- primary cohort, log2(TPM + 1)
- eligible genes (20% filter), MGMT excluded
- BH over eligible genes

Every other run is a sensitivity analysis.

**Output for every non-constant gene:**
- Ensembl ID, symbol, eligibility flags, n
- batch-adjusted rho, p, q, signed rank (rank 1 = largest positive rho) and
  absolute-rho rank
- the same columns for the pooled Spearman

Fig 2 is drawn from this table, the same one that gives the ranking.

**Definition of "rank 1":** the largest positive batch-adjusted rho among
eligible genes, MGMT excluded. Genes with exactly equal rho (to 1e-12) get
consecutive ranks in gene-table order; a tie involving GPR153 is reported.

### 5.1 Sensitivity analyses (one table: rho, p, q and rank for each)

- Pooled Spearman (reported with equal prominence; see above).
- **BI batch only; GDC-2024 batch only.** These are decision condition 3, a
  sign-only check. For each batch the GPR153-MGMT Spearman rho is reported with
  a Bonett-Wright 95% CI (k = 0) and GPR153's within-batch rank. The two
  batches come from the same cohort and pipeline. The manuscript does not
  describe their agreement as replication.
- FPKM (fpkm_unstranded) instead of TPM.
- **Outlier removal.** For the screen: samples with |z| > 3 for MGMT
  log2(TPM + 1) are removed (z within the primary cohort), and the full screen
  is rerun. GPR153 values do not affect which samples enter the screen. For
  the single GPR153-MGMT value: samples with |z| > 3 for either GPR153 or MGMT
  are removed.
- Expression filters at 10%, 50% and none.
- IDH-not-mutant cohort (n = 264); sensitivity cohort (n = 284); primary cohort
  excluding patients with IDH source disagreement.
- **Rank stability.** 1,000 bootstrap resamples of patients (seed 42). Report:
  - the proportion with GPR153 at rank 1 and in the top 10
  - the percentile 95% CI of its rho
  - the rho values of ranks 1-10, to show the gap
- **Cross-batch reproducibility.** 79 patients in the primary cohort have
  primary-tumor aliquots in both batches (85 across all primary tumors). For
  them, the per-gene correlation between batches of log2(TPM + 1) is computed.
  MGMT and GPR153 are excluded from this prespecified computation and its
  outputs. Their own cross-batch values are a separate post hoc diagnostic
  (5.7).

### 5.2 Composition adjustment

GPR153 is enriched in brain tissue and may track non-tumor neural content. Two
adjusted GPR153-MGMT correlations are computed, each a partial Spearman that
also adjusts for batch:
- **ESTIMATE score**: single-sample GSEA of the ESTIMATE stromal and immune
  gene sets (tidyestimate 1.1.1) on log2(TPM + 1).
- **Neuronal signature score**: mean z-score of RBFOX3, SNAP25, SYT1, STMN2,
  GABRA1, SYN1, NEFL, GAD1, SLC17A7 and CAMK2A.

Each score enters the partial Spearman as ranks, with the batch indicator
(Section 5, step 2; k = 2). The genome-wide screen is repeated with each
adjustment, and GPR153's rank is reported.

### 5.3 Methylation adjustment

Condition 4, as specified in Section 3.2.

### 5.4 Family-wise permutation null

10,000 permutations (seed 42). In each, MGMT values are shuffled within batch
and the full screen is recomputed (eligible genes, batch-adjusted). Reported:
- the proportion of permutations whose maximum rho is >= the observed maximum
- the proportion in which the gene at GPR153's observed rank has rho >= GPR153's
  observed rho
- the permutation distribution of the maximum rho

This is descriptive; it is not a pass condition.

### 5.5 Decision rule (fixed before the screen is run)

**Pass** if all four hold:
1. GPR153 is rank 1 with q < 0.05 in the primary analysis (batch-adjusted).
2. The GPR153-MGMT partial correlation is positive with q < 0.05 after each
   composition adjustment in 5.2. Each q comes from the corresponding adjusted
   genome-wide screen.
3. GPR153-MGMT rho is positive in both the BI-only and the GDC-2024-only
   analyses (sign only; 5.1).
4. The methylation-adjusted partial correlation is positive with a 95% CI
   excluding 0 (Section 3.2, condition 4).

### 5.6 Outcome branches and wording

Exactly one branch applies. It is assigned in this order:
1. If GPR153 is not positive with q < 0.05 in the primary analysis → **(d)**.
2. Otherwise, if condition 1 fails → **(b)**.
3. Otherwise, if any of conditions 2-4 fails → **(c)**.
4. Otherwise → **(a)**.

All four conditions are evaluated and reported in every branch.

| Branch | Definition | Reported in addition |
|---|---|---|
| (a) | All four conditions pass | rho, q, pooled rho and rank, condition 2-4 values |
| (b) | Condition 1 fails, but GPR153 is positive with q < 0.05 | GPR153's rank; the rho gap to rank 1 (rho of rank 1 minus rho of GPR153) and the name of the rank-1 gene; adjusted results (conditions 2-4) are reported but do not change the branch |
| (c) | Condition 1 passes, but one or more of conditions 2-4 fails | each failing condition, by name, with its value |
| (d) | GPR153 is not positive with q < 0.05 | GPR153's rho, q and rank |

**Wording.** Bracketed items are filled from the results. The title never
contains "biomarker", "driver", "predicts" or "independent".

| Branch | Title | Abstract (result sentence) | Conclusion |
|---|---|---|---|
| (a) | "GPR153 is the gene most strongly co-expressed with MGMT in IDH-wildtype glioblastoma: an exploratory re-analysis" | "In a time-stamped exploratory analysis of [n] IDH-wildtype primary tumors, GPR153 had the largest batch-adjusted Spearman correlation with MGMT among [N] genes (rho = [x], q = [y]; pooled rho = [x'], rank [r'])." The same paragraph states that the association remained after adjustment for tumor composition and MGMT promoter methylation, and had the same sign in both sequencing batches | "GPR153 co-expression with MGMT is a correlational finding in TCGA, with [external result] in external cohorts. It warrants experimental study; it does not establish a functional link" |
| (b) | "GPR153 and MGMT co-expression in IDH-wildtype glioblastoma: an exploratory re-analysis" | "GPR153 was positively correlated with MGMT (rho = [x], q = [y]) but ranked [r] of [N] genes (rank 1: [gene], rho = [x1])." | "GPR153 is one of many genes positively co-expressed with MGMT; the data do not single it out" |
| (c) | As (b) | "GPR153 had the largest batch-adjusted correlation with MGMT (rho = [x], q = [y]), but this did not hold after [failing condition(s)]." | "The ranking of GPR153 depends on [batch / tumor composition / MGMT promoter methylation] and is not robust" |
| (d) | "GPR153 expression, MGMT co-expression and survival in IDH-wildtype glioblastoma: an exploratory re-analysis" | "GPR153 was not positively correlated with MGMT at q < 0.05 (rho = [x], q = [y], rank [r] of [N])." | "The earlier report of GPR153-MGMT co-expression was not reproduced in this re-analysis" |

In every branch:
- The abstract and Methods state that the analysis followed a time-stamped
  analysis plan, and that the hypothesis arose from an earlier exploratory
  analysis whose code could not be reproduced.
- The abstract reports the pooled and batch-adjusted results side by side.
- Downstream analyses are described as exploratory characterization of GPR153,
  not as follow-up of a screen winner, unless the branch is (a).

### 5.7 Post hoc diagnostic (labeled post hoc)

After registration, GPR153's and MGMT's own cross-batch correlations are
reported in the 79 overlapping IDH-wildtype primary-cohort patients:
- Spearman correlation between the BI and the GDC-2024 aliquot, for each gene
- their percentile among all genes in 5.1

This is labeled "post hoc" in the text, tables and figure legends. It is not
part of the decision rule and does not change the branch.

## 6. Survival analyses (TCGA)

### 6.1 Primary survival estimate

- **Model 3**: GPR153 per SD of log2(TPM + 1), with age and sex.
  - unpenalized Cox model (lifelines CoxPHFitter, penalizer = 0, Efron ties,
    Newton-Raphson convergence precision 1e-12; the same setting is used for
    every lifelines Cox fit)
  - IDH-wildtype primary cohort: 249 patients, 207 deaths
- Reported: HR, 95% CI, p, n, events; proportional-hazards test per covariate
  and global.
- Described at most as an "age- and sex-adjusted association". It is not called
  an "independent prognostic biomarker".
- **Result wording** (bracketed items filled from the result; the direction is
  "higher" if HR > 1 and "lower" if HR < 1):
  - 95% CI excludes 1: "In an age- and sex-adjusted Cox model, higher GPR153
    expression was associated with a [higher/lower] hazard of death (HR [x] per
    SD, 95% CI [l]-[u], p = [p]; n = [n], [e] deaths)."
  - 95% CI includes 1: "In an age- and sex-adjusted Cox model, GPR153 expression
    was not clearly associated with the hazard of death (HR [x] per SD, 95% CI
    [l]-[u], p = [p]; n = [n], [e] deaths)."
  - In both cases the estimate is called an "age- and sex-adjusted
    association", and no causal or predictive term is used.
- All other models (6.2) and all other cohorts (sensitivity cohorts, external
  cohorts in Section 8) are secondary.

### 6.2 Cox models

All models are unpenalized (penalizer = 0), with GPR153 per SD of log2(TPM + 1).
The numbering below is used everywhere in the plan and the manuscript.

| Model | Covariates | Sample | Role |
|---|---|---|---|
| Model 1 | GPR153 | survival primary cohort (249) | secondary |
| Model 2 | GPR153 + age | survival primary cohort | secondary |
| **Model 3** | GPR153 + age + sex | survival primary cohort | **primary** |
| Model 4 | GPR153 + age + sex + KPS + MGMT methylation | complete cases | secondary |
| Model 5 | Model 4 + MGMT expression | Model 4 complete cases | secondary |
| Model 6 | GPR153 + age + sex + MGMT expression | survival primary cohort | secondary |
| Model 7 | GPR153 + age + sex + MGMT methylation | methylation subset | secondary |

- MGMT methylation in Models 4, 5 and 7 is the measure selected for condition 4
  (Section 3.2).
- Nested comparisons use likelihood-ratio tests on identical samples. Models
  1-3 are refitted on the Model 4 complete-case set, so that changes in the
  GPR153 estimate can be separated from changes in sample.
- Batch (sensitivity): Models 1-7 refitted (a) stratified by batch and (b) with
  batch as a covariate.
- Treatment (sensitivity only): Models 4 and 5 with temozolomide and
  radiotherapy added. Treatment is recorded after diagnosis and requires
  survival long enough to receive it. These models are therefore subject to
  immortal-time bias and are not used for inference about GPR153.
- Sensitivity cohorts: Models 1-7 refitted in the IDH-not-mutant and the
  sensitivity cohorts, with IDH status added as a covariate in the sensitivity
  cohort.
- Reported for every covariate: coefficient, HR, 95% CI, p, n, events.
  - Proportional hazards: Schoenfeld residual test per covariate and global.
  - Age linearity: restricted cubic spline (3 knots) vs linear, LR test.
- Interpretation: no "suppressor variable" interpretation of any change in
  estimates between models, anywhere in the manuscript or response.
- Extent of resection is not available (Section 3.1) and cannot be adjusted for.

### 6.3 Missing data

- Primary: complete-case analysis for Models 4, 5 and 7.
- Sensitivity: multiple imputation of KPS and MGMT methylation status (R mice
  3.19.0):
  - m = 20 imputations, chained equations, 20 iterations, seed 42
  - predictive mean matching for KPS; logistic regression for methylation status
  - imputation model: all model covariates, GPR153, MGMT expression, batch, the
    event indicator and the Nelson-Aalen cumulative hazard
  - estimates pooled by Rubin's rules

### 6.4 Kaplan-Meier (visualization)

- **Displayed figure (Fig 3)**: pooled median split (one cutoff for the primary
  cohort), for GPR153 (A) and MGMT (B). This is the only KM figure in the main
  text.
- **Supplement (S3 Fig)**: within-batch median split (each patient compared with
  the median of their own batch).
- Reported per group, for both versions:
  - n, deaths, censored
  - median OS with 95% CI
  - log-rank p
  - numbers at risk at 0, 6, 12, 18, 24 and 36 months
  - median follow-up by reverse Kaplan-Meier
- Each KM figure draws its at-risk table under the curves, from the fitted
  objects.
- The median split is not presented as a clinical cutoff.

### 6.5 Log-rank implementation and tests

The log-rank test uses lifelines `logrank_test` with durations and event
indicators. Automated tests (`tests/test_survival.py`) run in `run_all.sh`:
- **Known answer.** On R `survival::aml` (Maintained vs Nonmaintained),
  survdiff gives chi-square 3.3963887, p = 0.0653393. Before registration,
  lifelines reproduced it: 3.3963887, p = 0.0653393.
- **Event column matters.** Inverting the event indicator on the same data must
  change the result (before registration: chi-square 0.304, p = 0.581).
- **Sample identity.** The pipeline's survival table equals the manifest
  survival columns for the analysis cohort.

Further verification tests are listed in Section 13.1.

### 6.6 Selection context and incremental value

- Univariate Cox (per SD) for every eligible gene in the primary cohort; GPR153's
  percentile by p and by |log HR| is reported.
- Incremental value of GPR153 over the clinical baseline (Model 4 without
  GPR153: age, sex, KPS, MGMT methylation; complete cases):
  - likelihood-ratio test, Model 4 vs baseline
  - change in Harrell's C
  - bootstrap optimism-corrected C for both models (1,000 resamples, seed 42)

### 6.7 Survival-data audit

For every survival cohort, `pipeline/r07_survival_audit.py` (run in
`run_all.sh`) reads the raw clinical source and:
- tabulates the raw vital-status strings
- maps them to the event indicator
- converts time to days
- applies the time <= 0 rule
- asserts that n and events recomputed from the raw source equal the pipeline's
  outputs. The pipeline stops on any mismatch.

Outputs: `revision/results/cohort/survival_audit_strings.csv` and
`survival_audit.csv`. Values before registration:

**Raw strings and mapping** (within each cohort's selection):

| Cohort | Raw vital-status value (n) | Mapping to event |
|---|---|---|
| TCGA-GBM primary | "Dead" (208), "Alive" (42), "Not Reported" (1), missing (1) | Dead = 1, Alive = 0; others excluded |
| CGGA mRNAseq_693 | "1" (92), "0" (13), missing (4) | column "Censor (alive=0; dead=1)": 1 = death, 0 = censored |
| CGGA mRNAseq_325 | "1" (64), "0" (10) | as above |
| GSE13041 GPL96 primary | "DECEASED" (149), "ALIVE" (13) | DECEASED (or DEAD) = 1, ALIVE = 0 |
| Rembrandt (GSE108474) | "1" (187), "0" (23), missing (10) | EVENT_OS: 1 = death, 0 = censored (Appendix A, U12) |
| CPTAC-3 GBM | "Dead" (141), "Alive" (50) | Dead = 1, Alive = 0 |

**Time, units and counts:**

| Cohort | Time field | Unit in source | Conversion | Selected | Time missing | Time <= 0 | Usable n | Events | Censored |
|---|---|---|---|---|---|---|---|---|---|
| TCGA-GBM primary | days_to_death (Dead); max of days_to_last_follow_up, days_to_follow_up (Alive) | days | none | 252 | 3 | 0 | 249 | 207 | 42 |
| CGGA mRNAseq_693 | OS | days (inferred; Appendix A, U11) | none | 109 | 4 | 0 | 105 | 92 | 13 |
| CGGA mRNAseq_325 | OS | days (inferred; U11) | none | 74 | 0 | 0 | 74 | 64 | 10 |
| GSE13041 GPL96 primary | 1.tts(days) | days | none | 162 | 0 | 0 | 162 | 149 | 13 |
| Rembrandt (GSE108474) | OVERALL_SURVIVAL_MONTHS | months | × 30.4375 | 220 | 10 | 0 | 210 | 187 | 23 |
| CPTAC-3 GBM | as TCGA | days | none | 191 | 0 | 3 | 188 | 141 | 47 |

**Rules:**
- Time <= 0, missing time, or an unmapped status: excluded from survival
  analyses in every cohort.
- All analyses use days. KM axes are shown in months (days / 30.4375).
- The audit found one CPTAC case with negative follow-up (C3L-06625, Alive, -8
  days). It had been counted by the candidate check, which now applies the
  time > 0 rule (CPTAC usable n 189 → 188).

## 7. External co-expression cohorts

| Cohort | Selection | Expression |
|---|---|---|
| GSE13041 | GPL96 arrays with description "Primary GBM"; whole-profile duplicates excluded (7.1) | MAS5 as deposited, log2 |
| GSE4290 | Histology = glioblastoma, grade 4 (n = 77); all 180 arrays in a sensitivity analysis (7.3) | MAS5 as deposited, log2 |
| CGGA mRNAseq_693 | Primary, WHO grade IV, IDH-wildtype (n = 109) | log2(FPKM + 1) |
| CGGA mRNAseq_325 | Primary, WHO grade IV, IDH-wildtype (n = 74) | log2(FPKM + 1) |

The two CGGA releases are analyzed as separate cohorts (different sequencing
batches) and pooled as a sensitivity analysis (7.2).

### 7.1 GSE13041 reconciliation

GPL96 holds 191 arrays: 163 with description "Primary GBM" and 28 with
"Recurrent GBM". The duplicate rule (Spearman r >= 0.98 across all probes
except the GPR153 and MGMT probes; the later GSM accession is excluded) was
fixed before the check and is not changed. It finds two pairs among the 191:

| Pair | r | Group | Effect |
|---|---|---|---|
| GSM326825 / GSM326826 (GBM_00PA / GBM_00QP) | 0.983 | both Primary GBM | GSM326826 excluded: primary 163 → **162** |
| GSM326662 / GSM326670 (GBM_1342 / GBM_2158) | 0.981 | both Recurrent GBM | outside the primary set; GSM326670 excluded only in the all-191 run: 191 → 189 |

So the earlier report (2 pairs, 189 of 191) and this plan (1 pair, 163 → 162)
are the same check applied to different sets.

**The earlier n = 163.** The output `complete_validation_picture.csv` stores
only N, r and p ("GSE13041 (GPL96 pri)", n = 163), and the generating code is
not on the machine. The only metadata rule that gives 163 and matches the label
is description == "Primary GBM". It is a subset of the 191 and not independent
of them.

**Criterion for confirming the earlier values** (computed after
registration):
- Candidate combinations: each GPL96 probe set annotated to GPR153 (221902_at,
  64942_at; GPL96 annotation) × the MGMT probe set (204880_at), each on the
  deposited scale and on log2 (4 combinations).
- Confirmed only if one combination meets both:
  - Pearson r on the 163 "Primary GBM" arrays reproduces +0.116
  - on all 191 GPL96 arrays it reproduces +0.119

  "Reproduces" means it rounds to the reported value at the reported precision
  (3 decimals).
- The combination that does so is reported, and all four are tabulated.
- Otherwise the claims ledger records "not reproducible" for both values.

**Runs** (Spearman primary; Pearson reported for bridging):
- (a) primary: 162 arrays
- (b) duplicate retained: 163 arrays
- (c) all 191 GPL96 arrays, no duplicate exclusion: bridges to the earlier
  +0.119
- (d) all GPL96 arrays after duplicate exclusion: 189
- (e) with and without samples treated before resection

### 7.2 CGGA 325 vs 693 overlap

- **Expression profiles**: no whole-profile duplicates (Spearman r >= 0.98)
  within either release or between them. The releases share no sample IDs.
- **Clinical records**: a pair is a match if age, sex, grade, histology, OS time
  and censor status are all equal, with text fields stripped of surrounding
  whitespace (`pipeline/r08_cgga_clinical_match.py`,
  `revision/results/cohort/cgga_clinical_match.csv`).
  - Before registration: 0 matched pairs among all 325 × 693 records, and 0
    among the selected samples.
  - For information, 5 pairs agree on age, sex and OS only. All 5 differ in
    grade, histology or censor status, and none has both members selected.
- **Rule**: a matched pair (by profile or clinical record) is analyzed once. The
  mRNAseq_693 record is kept and the mRNAseq_325 record is excluded, in both the
  separate and the pooled analyses.

### 7.3 Statistics

- Probes, from the platform annotation tables:
  - GPR153: 221902_at and 64942_at (both on GPL96 and GPL570)
  - MGMT: 204880_at (GPL96, GPL570)
  - Primary rule: the GPR153 probe set with the highest mean expression in the
    cohort. Every probe combination is reported as a sensitivity analysis.
- Spearman rho (primary) and Pearson r (sensitivity). 95% CI via Fisher z with
  the Bonett-Wright SE, sqrt((1 + rho²/2) / (n - 3)). BH across cohorts.
- **Meta-analysis**, with and without TCGA:
  - random effects, DerSimonian-Laird on Fisher z, with tau² and I²
  - sensitivity: Hartung-Knapp(-Sidik-Jonkman) adjustment of the pooled CI
  - tau² uncertainty: with 4-5 cohorts, tau² and I² are imprecise. A
    Q-profile 95% CI for tau² and a 95% prediction interval are reported, and
    heterogeneity is described as poorly estimated.
- Power: for each cohort, power to detect rho = 0.12, 0.15, 0.20 and the TCGA
  estimate at alpha = 0.05, two-sided.
- **Platform and histology sensitivity analyses** (prespecified because the
  GSE13041 and GSE4290 selections were made after earlier values were seen;
  Section 1). They do not enter the meta-analysis above. All are reported in
  Table 4.
  - GSE13041 per platform, for the platforms with both a GPR153 and an MGMT
    probe set: GPL96 (all 191 arrays; the 162 primary arrays) and GPL570 (all
    27 arrays; the 11 "Primary GBM" arrays). GPL8300 (49 arrays) has no GPR153
    probe set. Each platform uses its own highest-mean GPR153 probe set and the
    log2 rule above.
  - Cross-platform overlap check, before pooling: (1) a shared tumor ID (sample
    title); (2) an identical clinical record (age, sex, survival time, vital
    status); (3) whole-profile Spearman r >= 0.98 between a GPL96 and a GPL570
    array on the probe sets common to both platforms, GPR153 and MGMT probe
    sets excluded. Checks (1) and (2) were run before registration from the
    clinical table: no shared tumor ID and no identical clinical record
    (Section 1.1). Check (3) is run after registration. The result of all three
    is reported.
  - The two platforms are pooled by DerSimonian-Laird random effects on Fisher
    z (as above), for all arrays and for primary arrays, only if all three
    checks find no overlap. Otherwise only the per-platform values are
    reported.
  - GSE4290 with all histologies: all 180 arrays (77 glioblastoma, 38
    oligodendroglioma grade 2, 23 non-tumor, 19 astrocytoma grade 3, 12
    oligodendroglioma grade 3, 7 astrocytoma grade 2, 4 missing histology).
  - The earlier same-session values (Section 1 table: GSE13041 all platforms
    -0.285, n = 218; GPL570 -0.360, n = 27; "Primary" -0.421, n = 11; GSE4290
    all histologies -0.081, n = 180) are listed beside these results, read
    from the original output files, and labeled "reported in same-session
    outputs; generating code unavailable; not reproduced".
- Discordance diagnostics:
  - GPR153 probe expression and detection relative to the array distribution
  - sample composition (histology, subtype, IDH where available)
- GSE13041 and GSE4290 deposit no IDH data, so IDH-mutant tumors cannot be
  excluded. This is recorded in the manifest and stated as a limitation.
- Manifest of every sample, with inclusion or exclusion reason:
  `revision/manifests/external_cohort_samples.csv`.

## 8. Survival replication

**Usable cohort**, all required:
- glioblastoma
- survival time > 0 and vital status present (audit, Section 6.7)
- a GPR153 probe or gene on the platform (checked from the annotation only)
- open access

Every cohort checked is logged with its reason
(`revision/results/cohort/survival_replication_candidates.csv`):

| Cohort | Survival fields | GPR153 on platform | Usable n | Decision |
|---|---|---|---|---|
| TCGA-GBM (discovery) | yes | yes (ENSG00000158292) | 249 | included: discovery cohort; included in meta-analysis with and without it |
| CGGA mRNAseq_693 | yes (OS, Censor) | yes | 105 | included: meets all criteria |
| CGGA mRNAseq_325 | yes (OS, Censor) | yes | 74 | included: meets all criteria |
| GSE13041 GPL96 primary | yes (tts, vital status) | yes (221902_at, 64942_at) | 162 | included: meets all criteria; no IDH data |
| GSE4290 GBM | no | yes (221902_at, 64942_at) | 0 | excluded: no survival data deposited |
| GSE16011 (Gravendeel 2009) | no (characteristics: age at diagnosis, gender, histology, tissue) | no (Entrez 387509 not in GPL8542) | 0 | excluded: GPR153 not on platform; no survival data in GEO record |
| GSE68848 (Rembrandt subset) | no (characteristics: disease state, organism part, sample, tumor grading) | yes (221902_at, 64942_at) | 0 | excluded: no survival data in GEO record |
| Rembrandt (GSE108474) | yes (OVERALL_SURVIVAL_MONTHS, EVENT_OS) | yes (221902_at, 64942_at) | 210 | included: meets all criteria; no IDH data; primary vs recurrent not recorded |
| CPTAC-3 GBM | yes (GDC vital status, days) | yes (same GDC pipeline as TCGA) | 188 | included: meets all criteria; IDH from open masked MAFs |

Usable n is the number with survival time > 0, mapped vital status and
expression after the cohort's selection rules:
- TCGA: primary cohort.
- CGGA: primary, WHO IV, IDH-wildtype.
- GSE13041: the 162 primary arrays.
- Rembrandt: DISEASE_TYPE = GBM, with a mapped expression column.
- CPTAC: primary_diagnosis Glioblastoma, primary tumor with open STAR counts.
  Of the 200 CPTAC primary-tumor RNA patients, 184 are IDH-wildtype and 16
  IDH-mutant by open MAFs.

**Analysis:**
- GPR153 per SD in a Cox model in each usable cohort: (a) alone; (b) with age
  and sex where available, matching Model 3.
- Pooled by random-effects meta-analysis (DerSimonian-Laird on log HR) with tau²
  and I², with and without TCGA. Hartung-Knapp sensitivity, Q-profile CI for
  tau² and prediction interval as in 7.3.
- Expression units: as in Section 7. Rembrandt (GSE108474) values are
  dChip-normalized as deposited (GEO data-processing field). They are
  log2-transformed if the matrix maximum exceeds 100; otherwise they are used as
  deposited. CPTAC uses log2(TPM + 1).
- Rembrandt probe for GPR153: as in Section 7.3 (highest-mean probe set among
  221902_at and 64942_at).
- CPTAC IDH status comes from open MAFs (same rule as TCGA). The primary
  analysis is IDH-wildtype; sensitivity uses all.
- Rembrandt has no IDH data and no primary/recurrent field; stated as a
  limitation.
- Sex in GSE13041 comes from the characteristics field "gender" (M/F; present
  for all 162 primary arrays).
- Incremental value (LR test, change in C, optimism-corrected C) is computed in
  TCGA (Section 6.6). It is also computed in CGGA, with baseline age, sex and
  MGMT methylation (CGGA has no KPS).

## 9. (Removed in version 1.5)

The cell-line and single-cell analyses were removed before registration. No
script for them was written and no result was computed. Section numbers are
unchanged. Composition confounding is therefore addressed only by the
adjustments in 5.2 (Appendix A, U13).

## 10. Pathway enrichment

- **Input lists**, from the primary-cohort batch-adjusted Spearman correlation
  with GPR153 (MGMT and GPR153 excluded):
  - up to 200 genes with the largest positive rho among those with q < 0.05
  - separately, up to 200 with the most negative rho and q < 0.05

  If fewer than 200 pass, all passing genes are used and the number reported.
- **Over-representation analysis** (gseapy enrich, run locally):
  - background: the eligible genes, restricted to the intersection of
    background and gene sets
  - libraries: Reactome_2022, KEGG_2021_Human, GO_Biological_Process_2023
    (Enrichr library files, downloaded and logged)
  - gene-set size 10-500; BH within each library
- **GSEA preranked** (gseapy prerank):
  - all eligible genes, ranked by batch-adjusted rho with GPR153
  - MSigDB 2026.2.Hs Hallmark and Reactome (c2.cp.reactome); size 15-500
  - 1,000 permutations, seed 42; FDR per library
  - gseapy prerank permutes gene labels. It therefore ignores inter-gene
    correlation, and its p-values are anti-conservative for correlated gene
    sets.
- **CAMERA** (limma 3.68.5; supplement only, 13.2), correlation-aware:
  - design ~ GPR153 + batch on log2(TPM + 1) of eligible genes; test of the
    GPR153 coefficient
  - inter-gene correlation estimated per set (inter.gene.cor = NA); sensitivity
    with the limma 3.68.5 default preset (inter.gene.cor = 0.01)
  - GPR153 is removed from the expression matrix before testing
  - same gene sets and size limits; FDR per library
- **Targeted tests**, reported under GSEA (and CAMERA if kept) regardless of
  outcome: HALLMARK_WNT_BETA_CATENIN_SIGNALING and HALLMARK_DNA_REPAIR.
- **Optional, MGMT-ranked enrichment** (Tier 2, last item; the first analysis
  dropped under 13.2; supplement only):
  - ORA, as above, on up to 200 genes with the largest positive batch-adjusted
    rho with MGMT among q < 0.05 in the primary screen (Section 5; MGMT and
    GPR153 excluded)
  - GSEA preranked, as above, on all eligible genes ranked by batch-adjusted
    rho with MGMT (GSEA uses the full ranking, not a top-200 list)
  - same libraries, size limits, background, permutations and seed
- Enrichment is described as descriptive and potentially circular, not as
  confirmation.

## 11. Random survival forest (Tier 3; supplement only)

- scikit-survival RandomSurvivalForest. Hyperparameters fixed in advance:
  n_estimators = 500, min_samples_leaf = 15, max_features = "sqrt",
  random_state = 42.
- Feature sets (labeled A-C, distinct from Cox Models 1-7):
  - (A) age + sex
  - (B) A + MGMT, EGFR, PTEN, IDH1, TERT expression
  - (C) B + GPR153

  In the complete-case subset, KPS and MGMT methylation are added to all sets.
  In the sensitivity cohort, IDH status is added.
- Validation: 10 repeats of 5-fold CV stratified by event, on identical folds
  for every set. Harrell's C on held-out folds.
- Reported:
  - the distribution of C across repeats
  - the paired difference C(C) - C(B), with a Nadeau-Bengio corrected interval
  - permutation importance on held-out folds
- No train/test split (70/30 or 80/20) is used. The Methods state that the
  earlier classifier (binary 1-year label) was replaced by a survival forest.
  Classification metrics (confusion matrix, class balance, AUC) do not apply
  to a survival forest, and the Methods say so.
- The six genes were selected after the earlier results were seen. They are
  never called "prespecified".
- TRIPOD checklist completed if this analysis is kept.

## 12. Differential expression (Tier 3; supplement only)

- **Tools:**
  - Primary: R DESeq2 1.52.0 on integer counts, design ~ condition, log2FC
    shrinkage by apeglm 1.34.0.
  - Cross-check: PyDESeq2. All three R packages (plus limma) were installed and
    load before registration.
- No ComBat-seq: tissue source and batch are fully confounded, so batch
  correction is not identifiable.
- **Data:** Toil RSEM expected counts (GENCODE v23), stored as log2(x + 1);
  back-transformed as 2^x - 1 and rounded to integers.
- **Tumors:** primary-cohort TCGA tumors present in Toil, matched on the
  15-character sample barcode. 144 of 252, all BI batch; the GDC-2024 aliquots
  are not in Toil.
- **Normals:**
  - GTEx Brain - Cortex (105 samples) and Brain - Frontal Cortex (BA9) (101
    samples), from 151 donors
  - 51 donors contribute both tissues. Primary uses one sample per donor
    (Brain - Cortex preferred), giving 151 normals. Sensitivity uses all 206.
- **Thresholds:** BH-adjusted p < 0.05 and |apeglm-shrunken log2FC| >= 1.
- **Gene IDs:**
  - Ensembl IDs are version-stripped.
  - Genes duplicated after stripping (PAR_Y) are dropped.
  - Genes matched between v23 and v36 by stripped ID; the overlap is reported.
  - Symbols are not used for matching.
- **Sensitivity:** TCGA primary tumors vs the 5 TCGA solid-tissue normals (GDC
  STAR unstranded counts).
- **Reported:**
  - the complete table: baseMean, log2FC, lfcSE, shrunken log2FC, p, padj
  - GPR153's rank by padj and by shrunken log2FC
  - the number and tissue site of normals
  - the number of tumors
  - the gene-ID handling counts
- "Most significantly upregulated" is claimed only if exactly supported.

## 13. Reproducibility and reporting

- **Environment:**
  - `revision/environment/`: pip lockfile with exact versions, R
    `sessionInfo()` and package versions, Python and OS versions
  - `run_all.sh` runs every pipeline step in order from the repository root,
    with all tests in 13.1 and the survival audit (6.7). It stops on the first
    failure.
- **Provenance:**
  - GDC queries (filters, fields, release) are saved as JSON next to each
    download and copied to `revision/manifests/gdc_queries/`.
  - Automated checks: unique patient IDs per analysis cohort; md5 of every input
    against the download log.
- **Methods concordance:**
  - The pipeline writes `revision/methods_parameters.json`: filters, seeds, CV
    design, library and package versions, and per-analysis n.
  - It also writes `revision/methods_concordance.csv`: each Methods statement →
    script → output file.
  - The revised Methods are drafted only from these two files.
- **Complete outputs**, untruncated: KM, Cox, DE, forest, enrichment (ORA,
  GSEA, CAMERA), validation, meta-analyses, permutation null, bootstrap, and the
  independent reproduction (13.2).
- **Figures:**
  - Plotting scripts read only result CSVs; no typed-in numbers. n, rho and p
    in titles and legends are filled from data.
  - KM curves are drawn from fitted objects. Enrichment plots show adjusted p.
    PCA uses the full primary cohort with k-means clusters labeled Cluster 1-3
    (no subtype names).
  - Format: TIFF or EPS at 300-600 dpi, numbered consecutively, each cited in
    the text and captioned, no zip files.
  - Each figure has a source-data CSV (`revision/source_data/`) and an entry in
    `revision/figure_script_map.csv`.
- **Tables:**
  - a per-analysis n table, with final n and events for every analysis
  - a replacement table for the response letter: for every analysis the editor
    flagged, what the original code did, what the manuscript said, and what the
    new pipeline does
- **Claims ledger** (`revision/results/claims_ledger.csv`), neutral wording.
  Columns: claim, location, reported value, original-code value, revised value.
- **Checklists:** REMARK and STROBE (and TRIPOD if the forest is kept), saved
  under `revision/checklists/`.
- **Literature:** searches and the citation audit follow
  `revision/literature_search_protocol.md`. AI tool use is logged in
  `revision/ai_tool_log.md`.
- **Data release:** a Zenodo release of derived, non-identifiable tables. These
  are the manifests, all result CSVs and source data, without raw expression
  files or controlled-access data. It is created only after the authors
  approve.
- **Analysis-plan integrity:** the SHA-256 of the uploaded copy of this file is
  recorded in `revision/PREREGISTRATION.sha256`. The OSF link, timestamp and DOI
  are recorded in `revision/PREREGISTRATION_REGISTRATION.md`, so this file does
  not change after upload.
- **Other:**
  - PCA of the primary cohort (top 2,000 most variable eligible genes,
    standardized); k-means k = 3, n_init = 50, random_state = 42.
  - Seeds: 42 for numpy, random, scikit-learn, scikit-survival, gseapy and R.
  - All p-values two-sided. Secondary analyses are not corrected across
    analyses and are labeled exploratory.

### 13.1 Verification tests (run in `run_all.sh`)

Each test fails the run if its tolerance is not met. Tests on simulated data use
seed 42.

| # | Test | Reference | Tolerance |
|---|---|---|---|
| T1 | Batch-adjusted partial Spearman with a single batch (k = 0) | `scipy.stats.spearmanr` (rho and p), simulated data | abs. difference <= 1e-10 |
| T2 | Batch-adjusted partial Spearman with two batches | `pingouin.partial_corr(method="spearman")` (rho and p, df = n - 2 - k), simulated data; version recorded | abs. difference <= 1e-8 |
| T2b | Partial Spearman with a continuous covariate (k = 2: batch and a simulated right-skewed, batch-related score), which enters as ranks (Section 5, step 2) | `pingouin.partial_corr(method="spearman")`, simulated data; the change in rho if the covariate entered raw is reported | abs. difference in rho and p <= 1e-8 |
| T3 | BH adjustment | `statsmodels.stats.multitest.multipletests(method="fdr_bh")` | abs. difference <= 1e-12 |
| T4 | Fisher-z and Bonett-Wright CI | Hand calculation: rho = 0.30, n = 100, k = 1: z = 0.309520, SE = 0.104333, 95% CI 0.1046 to 0.4731; k = 3: SE = 0.105437, CI 0.1025 to 0.4747 | 4 decimals |
| T5 | Schoenfeld residual tests (per covariate and global) | R `survival::cox.zph` (transform "km") on R `survival::lung` (age, sex) | abs. difference in p <= 1e-3. If they differ, R cox.zph is reported as primary and the difference is logged in DEVIATIONS.md |
| T6 | Known-answer log-rank; event-column inversion; survival sample identity | Section 6.5 | as stated there |
| T7 | Primary screen re-run in R | R: rank, `lm` residuals on batch, `cor`, on the same eligible genes | max abs. difference in rho <= 1e-8; identical rank order |
| T8 | Primary Cox model (Model 3) re-run in R | R `survival::coxph` (Efron ties) | relative difference in HR, CI and p <= 1e-4 |
| T9 | Guard checks (source-code scan, `tests/test_guards.py`). Exact patterns: (1) the substring `mannwhitneyu` in `pipeline/a2*.py` or `pipeline/a31_*.py`; (2) the substring `pearsonr` in `pipeline/a11_screen.py`; (3) a keyword argument `penalizer` whose value is not the literal 0, in any `pipeline/**/*.py` (syntax-tree scan); (4) in `pipeline/f1*.py` and `pipeline/plotstyle.py`, a string literal (including the literal parts of f-strings and concatenations) that is the first argument of `set_title`, `set_xlabel`, `set_ylabel`, `suptitle`, `text` or `annotate`, or a `label=` argument, and matches any of the regular expressions `\d+\.\d+`, `[=<>]\s*-?\d` or `\b[np]\s*[=<>]` (a decimal number; a comparison or equals sign followed by a number; n or p followed by a comparison or equals sign) | none | no match |
| T9b | T9 detects planted violations (`tests/test_guards_planted.py`): the pipeline directory is copied to a temporary directory, one pattern is planted in the copy (`mannwhitneyu` in `a21_km.py`; `pearsonr` in `a11_screen.py`; `penalizer=0.1` in `a20_cox.py`; `set_title("p = 0.002")` in `f10_main_figures.py`), and T9's scan is run on the copy, once per pattern | none | the unmodified copy scans clean, and each planted pattern is reported |
| T10 | Survival-data audit | Section 6.7 | n and events exactly equal |
| T11 | Null calibration of the batch-adjusted screen: simulated data (n = 252, 146 BI and 106 GDC-2024; 37,058 genes from a latent-factor model with batch shifts and ties), MGMT shuffled within batch and the full screen (`a11_screen.screen`) rerun 200 times | uniform distribution | maximum deviation of the empirical distribution of the pooled p-values from uniform <= 0.01; fraction of runs with any gene at q < 0.05 <= 0.09 |
| T12 | Planted effect, on the T11 simulated data: (a) one gene planted at batch-adjusted Spearman rho of about 0.45 with MGMT (n = 252); (b) a batch-driven pseudo-correlation (MGMT and one gene each shifted 2.5 SD between batches, independent within batch); (c) detection rate (q < 0.05) and rank-1 rate at rho 0.30 and 0.45 across the 200 null data sets, reported as a power check, not pass/fail | planted values | (a) rank 1 and q < 0.05; (b) pooled q < 0.05 and batch-adjusted q >= 0.05 |
| T13 | Cox recovery: 200 simulated data sets (n = 252; age, sex, GPR153 with a true HR of 1.25 per SD; Weibull baseline; about 17% censored; tied times), Model 3 fitted unpenalized with the a20 fitting function | true value | 95% CI coverage of the true log HR in 0.90-0.99 |
| T14 | ESTIMATE implementation (`a14_composition.py`) | tidyestimate 1.1.1 `estimate_score()` R code, sourced from the package tarball; simulated data with ties | abs. difference in stromal and immune scores <= 1e-8 |
| T15 | Vectorized univariate Cox (`common.cox_univariate`, used by a22) | R `survival::coxph` (Efron ties); simulated data with tied times | abs. difference in coefficient and SE <= 1e-6 |

### 13.2 Scope rule

- Only analyses with a row in `revision/methods_concordance.csv` enter the
  manuscript (main text or supplement).
- The MGMT-ranked enrichment, multiple imputation, CAMERA, the survival forest
  and DE are supplement-only.
- Each of these is included only if finished and verified (its 13.1 tests pass
  and it has a methods_concordance row) by 2026-10-29. An
  unfinished one is dropped and logged in DEVIATIONS.md, and the supplement is
  renumbered. The PLOS revision deadline is 2026-11-01.
- **Independent reproduction.** No number from the screen (Section 5,
  including the permutation null and the bootstrap), Model 3, the KM / log-rank
  analyses or the external meta-analyses (Sections 7 and 8) enters the
  manuscript unless an independent implementation reproduces it
  (`pipeline/a91_independent_check.py`, run in `run_all.sh`, which stops on a
  failure):
  - deterministic statistics: relative difference <= 1e-6; ranks: the
    pipeline's rank lies within the independent implementation's tie group
    (rho equal to 1e-12), since tied genes are ordered arbitrarily
  - bootstrap and permutation outputs: rerun with a separate implementation
    and another seed (43); absolute difference <= 3 Monte Carlo SE
    (proportions: binomial SE of the difference; quantiles: asymptotic SE with
    a kernel density estimate)
  - screen: batch-adjusted rho and p for GPR153, the top 100 genes and 500
    random eligible genes (seed 42) by `pingouin.partial_corr(method="spearman")`;
    rho, p, q and rank for all eligible genes by a separate pandas / numpy
    least-squares / statsmodels implementation; pooled Spearman rho and p for
    all non-constant genes by `scipy.stats.spearmanr`, with q and rank via
    statsmodels BH
  - Model 3: R `survival::coxph`. KM / log-rank: R `survival::survdiff` and
    `survfit` (log-log CI, as lifelines), numbers at risk, reverse-KM median
    follow-up. Meta-analyses: R `metafor::rma` (DerSimonian-Laird,
    Hartung-Knapp, Q-profile CI for tau², prediction interval with t on K - 2
    df)
  - both values of every comparison are logged in
    `revision/methods_concordance.csv`

## 14. Planned figures and tables

Fixed list, numbered consecutively in citation order. If an analysis is
dropped, its item is removed, the rest are renumbered, and the change is logged
in DEVIATIONS.md. Figures are TIFF or EPS, tables CSV/XLSX or in the text. No
zip files. Each item has a source-data CSV in `revision/source_data/` and a row
in `revision/figure_script_map.csv`.

**Main text**

| Item | Content | Section | Replaces (submitted) |
|---|---|---|---|
| Fig 1 | TCGA-GBM cohort flow (391 aliquots → 252 IDH-wildtype primary cohort) and external cohort selection | 3, 7 | new |
| Fig 2 | MGMT co-expression screen: (A) batch-adjusted rho of all eligible genes by rank, GPR153 marked; (B) the same for pooled rho; (C) GPR153 vs MGMT, log2(TPM + 1), colored by batch, both rho values in the legend; (D) bootstrap rank distribution of GPR153 and permutation null of the maximum rho | 5 | Fig 1A |
| Fig 3 | Kaplan-Meier OS, pooled median split, with at-risk tables: (A) GPR153; (B) MGMT | 6.4 | Fig 2 |
| Fig 4 | GPR153 HR per SD across Cox Models 1-7, Model 3 marked as primary | 6.1, 6.2 | Fig 3 |
| Fig 5 | External cohorts: (A) GPR153-MGMT rho per cohort with random-effects meta-analysis; (B) GPR153 survival HR per cohort with meta-analysis | 7, 8 | new |
| Fig 6 | Enrichment: (A) GSEA preranked, Hallmark and Reactome, by FDR, with the two targeted sets marked; (B) ORA by adjusted p | 10 | Fig 5 |
| Table 1 | Primary-cohort characteristics by batch | 3 | new |
| Table 2 | GPR153-MGMT: primary screen (batch-adjusted and pooled), decision conditions 1-4 and outcome branch, sensitivity analyses 5.1-5.3 | 5 | new |
| Table 3 | Cox Models 1-7: GPR153 and every covariate | 6.2 | Table 1 |
| Table 4 | External co-expression cohorts: n, Spearman and Pearson, 95% CI, power, GSE13041 runs (a)-(e), meta-analysis; sensitivity analyses (GSE13041 per platform with the cross-platform overlap check and pooling, GSE4290 all histologies); earlier same-session values, labeled "reported in same-session outputs; generating code unavailable; not reproduced" | 7 | Table 2 |
| Table 5 | Survival replication: audit counts, HR per cohort, meta-analysis | 6.7, 8 | new |
| Table 6 | Top enrichment results (ORA and GSEA) | 10 | Table 3 |

The submitted Fig 1B (differential expression) and Fig 4 (classifier) move to
the supplement (S7 Fig and S8 Fig).

**Supporting information**

| Item | Content | Section | Replaces (submitted) |
|---|---|---|---|
| S1 Fig | Distribution of mean log2(TPM + 1) across the 252 primary-cohort samples for all non-constant genes, with the eligibility filter marked; gene counts from the gene-universe table | 4 | S1 Fig (n = 391 aliquots, 59,569 genes) |
| S2 Fig | PCA of the primary cohort (top 2,000 most variable eligible genes), k-means Cluster 1-3, batch shown; no subtype names | 13 | S2 Fig (n = 115, subtype names) |
| S3 Fig | Kaplan-Meier, within-batch median split | 6.4 | new |
| S4 Fig | Cross-batch reproducibility in the 79 overlapping patients (all genes except GPR153 and MGMT); separate panel, labeled post hoc: GPR153 and MGMT (5.7) | 5.1, 5.7 | new |
| S5 Fig | Schoenfeld residuals for Model 3 | 6.2 | new |
| S6 Fig | CAMERA results for the targeted and top sets | 10 | new |
| S7 Fig | Differential expression volcano (DESeq2, shrunken log2FC) | 12 | Fig 1B |
| S8 Fig | Survival forest: C-index by feature set; permutation importance | 11 | Fig 4 |
| S1 Table | TCGA sample manifest | 3 | new |
| S2 Table | External-cohort sample manifest | 7 | new |
| S3 Table | Gene-level screen results, all non-constant genes | 5 | new |
| S4 Table | Per-analysis n and events | 13 | new |
| S5 Table | Survival-data audit (strings, mapping, units, counts) | 6.7 | new |
| S6 Table | Complete enrichment outputs (ORA, GSEA, CAMERA; MGMT-ranked enrichment if kept) | 10 | new |
| S7 Table | Multiple-imputation Cox results | 6.3 | new |
| S8 Table | Complete differential expression table | 12 | new |
| S9 Table | Claims ledger and replacement table | 13 | new |
| S1 File | This analysis plan, with its SHA-256 and registration record | 1 | new |
| S2 File | REMARK checklist | 13 | new |
| S3 File | STROBE checklist | 13 | new |
| S4 File | TRIPOD checklist (if the forest is kept) | 11 | new |

## 15. Tier assignment

One row per analysis.

| Analysis | Section | Tier | Placement |
|---|---|---|---|
| Cohort, manifests, batch comparison, annotation tables | 3 | 1 | main |
| Genome-wide screen, pooled and batch-adjusted, sensitivity runs | 5, 5.1 | 1 | main |
| Composition adjustment (condition 2) | 5.2 | 2, not droppable | main |
| Methylation adjustment (condition 4) and descriptive methylation analyses | 3.2, 5.3 | 2, condition 4 not droppable | main |
| Permutation null | 5.4 | 2 | main |
| Post hoc cross-batch diagnostic | 5.7 | 2 | supplement |
| Cox Models 1-7, KM, log-rank tests | 6.1, 6.2, 6.4, 6.5 | 1 | main |
| Multiple imputation | 6.3 | 2 | supplement |
| Genome-wide Cox percentile, incremental value | 6.6 | 2 | main |
| Survival-data audit | 6.7 | 1 | supplement table |
| External co-expression, meta-analysis, power, discordance | 7 | 1 | main |
| Survival replication | 8 | 2 | main |
| ORA and GSEA | 10 | 1 | main |
| CAMERA | 10 | 2 | supplement |
| MGMT-ranked ORA and GSEA (optional) | 10 | 2, last item; first to drop | supplement |
| Random survival forest | 11 | 3 | supplement |
| Differential expression | 12 | 3 | supplement |
| Reproducibility deliverables, verification tests and independent reproduction | 13, 13.1, 13.2 | 1 | supplement / repository |
| Reporting checklists | 13 | 2 | supplement |

The genome-wide screen (Section 5) is the only primary hypothesis test.
Everything else is exploratory.

## 16. Deviations

Any deviation from this plan, including any dropped analysis, is recorded in
`revision/DEVIATIONS.md` with date and reason, and reported in the revised
manuscript.

---

## Appendix A. Items not confirmed from a source file

Each item is stated in the plan but was not confirmed from a source file before
registration. IDs are kept stable across versions. Items marked "stated
limitation" are reported as limitations in the manuscript.

| # | Item | Status | Why unconfirmed | Handling |
|---|---|---|---|---|
| U1 | The earlier GSE13041 n = 163 is the 163 "Primary GBM" arrays (7.1) | stated limitation | The output has no sample list and the code is not on the machine | Confirmed or "not reproducible" by the criterion in 7.1, computed after registration |
| U2 | Number of primary-cohort patients with both STP27 probes non-missing (3.2) | stated limitation | Counting requires reading MGMT CpG values | Counted after registration; it decides between the primary and fallback measure for condition 4 (n >= 100 rule) |
| U3 | STP27 applied to SeSAMe beta-derived M-values reproduces the published GenomeStudio-based classification | stated limitation | Different preprocessing; agreement not measured | Kappa vs the Ceccarelli call after registration; below 0.60, the Ceccarelli call becomes the condition-4 measure |
| U4 | How the Ceccarelli 2016 and Brennan 2013 MGMT calls were derived | unconfirmed | Ceccarelli Table S1 (GDC publication page) has only the column "MGMT promoter status", with no method note. The Brennan 2013 publication-page files checked (DNA methylation call categories) contain no MGMT method. The supplementary methods text of both papers (PMC3910500, PMC4754110) is not open access: Europe PMC returns no full-text XML, and the PMC and Europe PMC article pages returned bot challenges, which were not bypassed | Described as "as provided by the publication"; the authors may supply the supplementary methods |
| U5 | (dropped in v1.5: Section 9 removed) | — | — | — |
| U6 | Rembrandt deposited values are on a linear or log scale (8) | stated limitation | Would require reading the expression matrix | Fixed rule in Section 8 (log2 if the matrix maximum exceeds 100) |
| U7 | (dropped in v1.3: GSE68848 is excluded, so its overlap with GSE108474 is not needed) | — | — | — |
| U8 | TCGA tumors in Toil come from the same aliquot as the selected GDC aliquot (12) | stated limitation | Toil uses 15-character sample barcodes without aliquot | Stated as a limitation of the DE analysis |
| U9 | (dropped in v1.5: Section 9 removed) | — | — | — |
| U10 | Extent of resection | resolved in v1.3 | Checked in GDC, TCGA-CDR, PanCanAtlas clinical data and the BCR clinical XML (Section 3.1): no field | Stated as "not available in any TCGA-GBM clinical source" |
| U11 | CGGA OS is in days (6.7) | unconfirmed | The CGGA clinical files give the column "OS" without a unit. Days are inferred from the range (19 to 5,028) | Stated as inferred; analyses use the values as days |
| U12 | Rembrandt EVENT_OS: 1 = death, 0 = censored (6.7) | unconfirmed | The data dictionary defines EVENT_OS as "Indicates death OR no record of death" without stating which code is which. 1 = death is inferred from the name and the event fraction (187 of 210) | Stated as inferred |
| U13 | Tumor composition does not explain the GPR153-MGMT correlation | stated limitation | Composition confounding is addressed only by the ESTIMATE and neuronal-score adjustments (5.2); with Section 9 removed there is no orthogonal check (cell lines, single cell) | Stated as a limitation |
