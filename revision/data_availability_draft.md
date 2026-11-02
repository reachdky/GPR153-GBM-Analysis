# Data Availability statement: draft

Draft for the authors to edit; not manuscript text. Items in **[brackets]** are unresolved. Sources for every
statement are listed at the end.

## Draft statement

All data analysed in this study are publicly available without access restriction; no controlled-access data
were used.

**TCGA-GBM (discovery cohort).** RNA-seq gene expression (391 "STAR - Counts" files, workflow "STAR - Counts",
GENCODE v36; values used: `tpm_unstranded`, `fpkm_unstranded` and unstranded counts), open-access masked somatic
mutation files (MAF), DNA methylation beta values (Illumina HumanMethylation27 and HumanMethylation450), clinical
case records and clinical supplements (BCR XML) were obtained from the NCI Genomic Data Commons (GDC,
https://portal.gdc.cancer.gov/), Data Release 46.0 (August 10, 2026). GDC metadata and clinical records were
retrieved on 24-26 October 2026; the expression files were re-retrieved and verified against their GDC md5
checksums on 1 November 2026. Clinical and molecular annotations were also taken from cBioPortal
(https://www.cbioportal.org/; studies `gbm_tcga_pub2013`, `gbm_tcga_pan_can_atlas_2018`, `lgggbm_tcga_pub`;
retrieved 24-25 October 2026) and from the supplementary files of the TCGA Pan-Cancer Clinical Data Resource
(Liu et al. 2018), Ceccarelli et al. 2016 and Brennan et al. 2013 (GDC publication pages; retrieved 26 October
2026).

**Normal brain reference.** GTEx brain cortex samples (with TCGA-GBM tumours) as RSEM expected counts from the
UCSC Toil recompute (https://xenabrowser.net/, `TcgaTargetGtex_gene_expected_count`, GENCODE v23; retrieved
24 October 2026).

**External cohorts.** Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/): GSE13041 (platform GPL96; its GPL570 arrays in a platform sensitivity analysis)
and GSE4290 (GPL570) (family SOFT files; re-retrieved and verified on 1 November 2026), and GSE108474 (Rembrandt;
GPL570 expression, clinical data and sample mapping; retrieved 25 October 2026). Chinese Glioma Genome Atlas
(https://www.cgga.org.cn/): mRNAseq_693 and mRNAseq_325, release 2020-05-06 (RSEM FPKM and clinical tables;
retrieved 24 October 2026). CPTAC-3 glioblastoma: STAR count files and clinical records from the GDC, Data
Release 46.0 (retrieved 25 October 2026). GSE16011 and GSE68848 were assessed as survival-replication
candidates only (retrieved 25 October 2026).

**Annotation and gene sets.** GEO platform tables GPL570, GPL8490, GPL13534 and GPL8542; MSigDB 2026.2.Hs
Hallmark and Reactome gene sets; Enrichr libraries Reactome_2022, KEGG_2021_Human and
GO_Biological_Process_2023; ESTIMATE gene sets from the tidyestimate 1.1.1 package (CRAN) (all retrieved 24-25
October 2026).

**Analysis plan.** The analysis plan (version 1.5) was registered on OSF before any revised GPR153 or MGMT value
was computed: https://osf.io/qpd2f/ (registered 27 October 2026, 23:44:45 UTC).

**Code and results.** All code, sample manifests (with inclusion and exclusion reasons), download records with
retrieval dates and checksums, GDC query parameters, the complete result tables, the source data for every
figure and table, and a figure-to-script map are publicly available at
https://github.com/reachdky/GPR153-GBM-Analysis.
The revised analysis is commit `55e41b3c030b95c965a6085c28e11e830006ada9` (the code of the final analysis run);
the archived release is commit `{ARCHIVE_COMMIT}`, which adds licence files, documentation and an archive build
script and changes no analysis code or result. The release is archived at Zenodo, https://doi.org/10.5281/zenodo.23872434 **[DOI reserved, record not
yet published: confirm before submission that the record is published and the DOI resolves]**. Software is
released under the MIT License, and figures, result tables and documentation under CC BY 4.0; third-party
materials keep their original terms. The code posted with the original submission is tagged `v1-submitted`
(commit `e0d6bb5`, 6 July 2026); it is retained in `legacy/` for reference and was not used for the revised
analysis.

## Unresolved

| Item | Status |
|---|---|
| Zenodo DOI | Reserved, not published: 10.5281/zenodo.23872434 (draft record saved by the authors; licences MIT and CC BY 4.0). The DOI resolves only after the record is published. Not resolved until then. |
| Public availability of the revised code | Resolved: the revised commits are on GitHub (`main`, public repository, first pushed 2026-11-01; tag `v1-submitted` = `e0d6bb5`). |
| Analysis version and archive commit | Analysis version: `55e41b3c030b95c965a6085c28e11e830006ada9` (final run 7). Archive commit: filled in the submission-package copy of this file from `REVISED_COMMIT.txt`; the repository copy keeps the placeholder `{ARCHIVE_COMMIT}` because a commit cannot contain its own hash. |
| Original retrieval of the TCGA STAR files and the GEO SOFT files | The files read by the analysis were first downloaded for the original submission (local file timestamps 24-25 December 2025; GDC file creation dates 2021-12-13 to 2024-10-04). The statement gives the 1 November 2026 verified re-retrieval; authors may add the original date. |
| Resources retrieved but not used by any analysis in the final run | DepMap 24Q4, GDSC release 8.5, the CELLxGENE GBmap core dataset, GSE131928 cell annotation, and the Bady et al. 2012 article files (`revision/manifests/download_dates.md`). Not named in the draft statement. |
| Submission form | The form currently says data will be made available on acceptance (journal requirement 3); it must be updated to match. |

## Sources

| Statement | File |
|---|---|
| Files, URLs, retrieval dates, md5 (r00) | `revision/manifests/DOWNLOAD_LOG.csv`, `revision/manifests/download_dates.md` |
| TCGA STAR and GEO files: retrieval, expected and observed md5 (r10) | `revision/manifests/raw_inputs_manifest.csv`, `revision/manifests/raw_inputs_expected_md5.csv`, `revision/logs/r10_fresh_download_2026-11-01.log` |
| GDC query parameters and data release | `revision/manifests/gdc_queries/*.json` |
| Gene annotation version and expression units | `revision/methods_parameters.json` |
| Sample manifests | `revision/manifests/tcga_samples.csv`, `revision/manifests/external_cohort_samples.csv` |
| Registration | `revision/PREREGISTRATION_REGISTRATION.md`, `revision/PREREGISTRATION.sha256` |
| Commits | `git show e0d6bb5`; `git show 55e41b3`; `REVISED_COMMIT.txt` in each archive |
| Licences | `LICENSE`, `LICENSE-CC-BY-4.0.txt`, `LICENSING.md` |
