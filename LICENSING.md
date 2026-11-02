# Licensing

This repository, its source archive and the submission package contain three kinds of material, under different
terms. Where a file mixes them (for example a manifest with author-defined columns and values taken from a data
source), each part keeps its own terms.

| Material | Terms | Licence text |
|---|---|---|
| Software written by the authors | MIT License | `LICENSE` |
| Figures, result tables, source data and non-software documentation written by the authors | Creative Commons Attribution 4.0 International (CC BY 4.0) | `LICENSE-CC-BY-4.0.txt`; https://creativecommons.org/licenses/by/4.0/ |
| Third-party materials (listed below) | The original terms of each provider; not licensed by the authors | Provider pages below |

Copyright holder for the authors' material: the GPR153-GBM-Analysis authors (2025-2026).

## 1. Software (MIT License)

- `pipeline/` (all `.py` and `.R` files, including `pipeline/sim/`)
- `tests/`
- `run_all.sh`
- `legacy/*.py` (scripts posted with the original submission; kept unchanged)
- `revision/verification/scripts/`
- `requirements.txt` and the requirement and pip-freeze files in `revision/environment/`
- `.gitignore`

## 2. Figures, results and documentation (CC BY 4.0)

- Figures: `revision/figures/`
- Tables and source data: `revision/tables/`, `revision/source_data/`
- Analysis outputs computed by the pipeline: `revision/results/`, `revision/methods_parameters.json`,
  `revision/methods_concordance.csv`, `revision/figure_script_map.csv`
- Documentation: `README.md`, `legacy/README.md`, this file, the analysis plan and its records
  (`revision/PREREGISTRATION*`), `revision/DEVIATIONS.md`, the manuscript-support documents in `revision/`,
  `revision/verification/`, `revision/logs/`, the remaining files in `revision/environment/`,
  `revision/manifests/provenance.md`, `revision/manifests/download_dates.md`, `revision/literature_search_protocol.md`,
  `revision/literature/literature_search_summary.md`, `revision/archive_release/README.md`
- The authors' selection, arrangement and derived columns (inclusion flags, reasons, computed values) in
  `revision/manifests/` and `revision/literature/`

Suggested attribution: the GPR153-GBM-Analysis authors, https://github.com/reachdky/GPR153-GBM-Analysis, with the
archive DOI once published (`revision/archive_release/README.md`).

## 3. Third-party materials (original terms apply)

The following content was obtained from the providers named, is included only as needed to document or reproduce
the analysis, and is not covered by the MIT License or CC BY 4.0. Users must follow each provider's terms and
citation requirements.

| Content | Where it appears | Provider; terms |
|---|---|---|
| TCGA-GBM and CPTAC-3 identifiers, clinical and molecular annotations, and file metadata (GDC Data Release 46.0, open access only) | `revision/manifests/tcga_samples.csv`, `revision/manifests/gdc_queries/*.json`, `revision/manifests/DOWNLOAD_LOG.csv`, `revision/manifests/raw_inputs_*.csv`, per-sample columns of `revision/results/` | NCI Genomic Data Commons, https://gdc.cancer.gov/access-data/data-access-processes-and-tools |
| Clinical and molecular annotations from cBioPortal studies (`gbm_tcga_pub2013`, `gbm_tcga_pan_can_atlas_2018`, `lgggbm_tcga_pub`) and from the TCGA publications they derive from | `revision/manifests/tcga_samples.csv` | cBioPortal, https://www.cbioportal.org/; terms of the source publications |
| GEO sample accessions and annotations (GSE13041, GSE4290, GSE108474, GSE16011, GSE68848) | `revision/manifests/external_cohort_samples.csv`, `revision/results/` | NCBI Gene Expression Omnibus, https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html (submitters may claim rights in their data) |
| CGGA sample identifiers and clinical annotations (mRNAseq_693, mRNAseq_325) | `revision/manifests/external_cohort_samples.csv`, `revision/results/` | Chinese Glioma Genome Atlas, https://www.cgga.org.cn/ |
| GTEx and TCGA sample identifiers from the UCSC Toil recompute | `revision/results/de/` | UCSC Xena, https://xenabrowser.net/datapages/; GTEx Portal, https://gtexportal.org/home/ |
| Gene-set names and gene memberships (MSigDB 2026.2.Hs Hallmark and Reactome; Enrichr libraries Reactome_2022, KEGG_2021_Human, GO_Biological_Process_2023) | `revision/results/enrichment/`, `revision/tables/Table_6.csv`, `revision/tables/S6_Table.csv`, Fig 6, S6 Fig | MSigDB, https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp; Enrichr, https://maayanlab.cloud/Enrichr/; Reactome, https://reactome.org/license; KEGG, https://www.kegg.jp/kegg/legal.html; Gene Ontology, https://geneontology.org/docs/go-citation-policy/ |
| ESTIMATE gene sets (tidyestimate 1.1.1, CRAN) | Not included: read from the downloaded package at run time by `pipeline/a14_composition.py` and `tests/test_estimate.py` | https://cran.r-project.org/package=tidyestimate |
| PubMed bibliographic records, including abstracts | `revision/literature/pubmed_records.csv`, `revision/literature/pubmed_raw/` | U.S. National Library of Medicine, https://www.nlm.nih.gov/databases/download/terms_and_conditions.html (abstracts may be protected by the publishers' copyright) |
| Europe PMC search results (titles and identifiers) | `revision/literature/europepmc_fulltext_searches.csv` | Europe PMC, https://europepmc.org/ |
| Crossref bibliographic metadata | `revision/literature/citation_audit.csv` (columns `reg_*`) | Crossref, https://www.crossref.org/documentation/retrieve-metadata/rest-api/ |

## 4. Not included

- Raw and downloaded input data (`data/`): open access, retrieved from the providers by
  `pipeline/r00_download_external.py` and `pipeline/r10_download_raw_inputs.py`.
- Software dependencies (Python and R packages): installed from PyPI, CRAN and Bioconductor under their own
  licences (versions in `revision/environment/`).
- Font files: the figures were rendered with locally installed Arial font files, which are not distributed.

No third-party licence or copyright notice was present in the repository before this file was added; none was
removed or changed.
