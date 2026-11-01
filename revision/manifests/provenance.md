# Data provenance

Every input of the revised analysis, with its source, version, retrieval date and checksum record. A
**documented** date was written by a download script at retrieval time; an **inferred** date is a file-system
timestamp and is reported only as such. Checks of 2026-11-01 are in `revision/verification/`.

## Inputs fetched by `pipeline/r00_download_external.py`

| Inputs | Source and version | Retrieved (documented, UTC) | Checksum record |
|---|---|---|---|
| 63 files: GDC metadata and case records, cBioPortal JSON, CGGA, UCSC Toil, GEO series and platform files, MSigDB, Enrichr, ESTIMATE, DepMap, GDSC, CELLxGENE, TCGA publication supplements, article files | URL per file in `DOWNLOAD_LOG.csv`; GDC Data Release 46.0; CGGA release 2020-05-06; Toil GENCODE v23; MSigDB 2026.2.Hs; Enrichr Reactome_2022, KEGG_2021_Human, GO_Biological_Process_2023; tidyestimate 1.1.1; DepMap 24Q4; GDSC 8.5 | 2026-10-24 to 2026-10-26 (`retrieved_utc` per file) | md5 per file in `DOWNLOAD_LOG.csv`; all 63 re-checked 2026-11-01: match (`revision/verification/download_log_md5_check_2026-11-01.csv`) |
| 2,020 files from GDC queries: TCGA-GBM masked MAFs (464), methylation beta values (450), clinical supplements BCR XML (599); CPTAC-3 STAR counts (238) and masked MAFs (269) | GDC Data Release 46.0; query, release and file list in `gdc_queries/*.json` | per query, 2026-10-24 to 2026-10-26 (`retrieved_utc` in each query record) | GDC md5sum per file, checked at download; all 2,020 re-checked 2026-11-01: match (`revision/verification/gdc_query_files_md5_check_2026-11-01.csv`) |

## Inputs fetched by `pipeline/r10_download_raw_inputs.py`

| Inputs | Source and version | Retrieved | Checksum record |
|---|---|---|---|
| 391 TCGA-GBM STAR - Counts files (GENCODE v36; columns `tpm_unstranded`, `fpkm_unstranded`, `unstranded`) | GDC `/data/<file_id>`; query `gdc_queries/tcga_star__files_meta.json` (GDC Data Release 46.0; GDC file creation dates 2021-12-13 to 2024-10-04 in `tcga_samples.csv`) | Documented: re-retrieved 2026-11-01 (`raw_inputs_manifest.csv`). The copies read by the analysis were downloaded for the original submission; inferred from file timestamps: 2025-12-24 | GDC md5 (`tcga_samples.csv`, `md5_gdc`): 391/391 match for the analysis copies and for the 2026-11-01 retrieval; `a10` checks them again on every run |
| GSE13041 and GSE4290 family SOFT files | GEO FTP `https://ftp.ncbi.nlm.nih.gov/geo/series/<GSEnnn>/<GSE>/soft/` | Documented: re-retrieved 2026-11-01. Analysis copies: inferred from file timestamps 2025-12-25 | md5 of the analysis copies (`raw_inputs_expected_md5.csv`); the 2026-11-01 downloads are byte-identical |
| `expression.csv`, `clinical.csv` for GSE13041 and GSE4290 | Built from the SOFT files with GEOparse 2.0.3 (`pivot_samples("VALUE")`, `phenotype_data`), as by `legacy/03_download_geo.py` | Rebuilt 2026-11-01. Analysis copies: inferred 2025-12-25 | The rebuilt tables are byte-identical to the analysis copies (`raw_inputs_manifest.csv`) |

## Annotation versions

| Data | Annotation |
|---|---|
| TCGA STAR counts, CPTAC-3 STAR counts | GENCODE v36 (GDC harmonized pipeline) |
| UCSC Toil (TCGA + GTEx) | GENCODE v23 (`toil/gencode.v23.annotation.gene.probemap`) |
| GSE13041 (GPL96), GSE4290 and Rembrandt (GPL570) | GEO platform tables (GPL96 from the GSE13041 SOFT file; `platforms/GPL570_table.txt`) |
| Methylation | GPL8490 (HM27), GPL13534 (HM450) platform tables |
| Gene sets | MSigDB 2026.2.Hs; Enrichr libraries as named |

## Reproduction

`./run_all.sh --download` runs r00 and r10, then the chain; `README.md` lists the environment set-up and commands.
Files already present are checked, not downloaded again. A file whose md5 differs is not kept and the step stops.
