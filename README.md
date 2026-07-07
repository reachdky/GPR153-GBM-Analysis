# GPR153-GBM-Analysis

Analysis scripts for the manuscript:

**"GPR153 as a potential biomarker of glioblastoma multiforme and its association with MGMT co-expression and overall survival"**

Yoog Patel — Independent Researcher, Matawan, NJ, USA

---

## Overview

This repository contains all Python scripts used to perform the computational analyses described in the manuscript. The pipeline covers data download, preprocessing, genome-wide correlation screening, survival analysis, machine learning, pathway enrichment, cross-dataset validation, and figure generation.

---

## Repository Structure

```
GPR153-GBM-Analysis/
├── download/
│   ├── download_tcga_gbm.py        # Download TCGA-GBM expression + clinical data via GDC API
│   ├── 02_download_gtex.py         # Download GTEx normal brain expression data
│   ├── 03_download_geo.py          # Download GEO validation cohorts (GSE13041, GSE4290)
│   └── 04_download_cgga.py         # Download CGGA Grade IV subset
│
├── processing/
│   ├── 05_clean_data.py            # Filter to common genes, clean expression matrices
│   ├── 06_normalize_expression.py  # log2(TPM+1) transformation
│   ├── 07_quality_control.py       # QC filtering and sample-level checks
│   └── combined_gbm_normal_expression.py  # Merge GBM and normal brain for DEG analysis
│
├── analysis/
│   ├── 08_expression_analysis.py   # Background expression distribution (Supp Fig S1)
│   ├── 09_correlation_analysis.py  # Genome-wide Spearman correlation + BH-FDR correction
│   ├── 10_survival_analysis.py     # Kaplan-Meier and Cox proportional hazards regression
│   ├── 11_pathway_analysis.py      # Enrichr pathway enrichment + GSEA
│   ├── 12_differential_expression.py  # GBM vs normal brain DEG analysis
│   ├── 13_clustering.py            # PCA and K-means molecular subtyping
│   ├── 14_predictive_modeling.py   # Random forest classifier + Gini feature importance
│   └── 15_validation.py            # Cross-dataset Pearson correlation validation
│
├── visualization/
│   ├── 16_generate_figure1.py      # Fig 1: GPR153-MGMT scatterplot + volcano plot
│   ├── 17_generate_figure2.py      # Fig 2: Kaplan-Meier survival curves
│   ├── 18_generate_figure3.py      # Fig 3: Cox proportional hazards forest plot
│   ├── 19_generate_figure4.py      # Fig 4: ROC curve + feature importance
│   └── 20_generate_figure5.py      # Fig 5: Pathway enrichment (GSEA + Enrichr Reactome)
│
└── cox_regression_v2.py            # Standalone multivariate Cox regression script
```

---

## Data Sources

All expression data used in this study are publicly available:

| Dataset | Source | Accession |
|---|---|---|
| TCGA-GBM (discovery) | GDC Data Portal | https://portal.gdc.cancer.gov/ |
| GSE13041 | NCBI GEO | GSE13041 |
| GSE4290 | NCBI GEO | GSE4290 |
| CGGA Grade IV | CGGA | http://www.cgga.org.cn |

No data files are included in this repository. All datasets must be downloaded using the scripts in the `download/` directory or accessed directly from the sources above.

---

## Requirements

Python 3.11 with the following packages:

```
pandas==2.0
numpy==1.25
scipy==1.11
scikit-learn==1.3
lifelines==0.27
matplotlib==3.7
seaborn==0.12
statsmodels
requests
```

Install all dependencies with:

```bash
pip install pandas==2.0 numpy==1.25 scipy==1.11 scikit-learn==1.3 lifelines==0.27 matplotlib==3.7 seaborn==0.12 statsmodels requests
```

---

## Running the Pipeline

Run scripts in numbered order:

```bash
# 1. Download data
python download/download_tcga_gbm.py
python download/03_download_geo.py
python download/04_download_cgga.py

# 2. Process
python processing/05_clean_data.py
python processing/06_normalize_expression.py
python processing/07_quality_control.py

# 3. Analyse
python analysis/09_correlation_analysis.py
python analysis/10_survival_analysis.py
python analysis/14_predictive_modeling.py
python analysis/15_validation.py

# 4. Generate figures
python visualization/16_generate_figure1.py
python visualization/17_generate_figure2.py
python visualization/18_generate_figure3.py
python visualization/19_generate_figure4.py
python visualization/20_generate_figure5.py
```

---

## Citation

If you use these scripts in your work, please cite:

> Patel Y. GPR153 as a potential biomarker of glioblastoma multiforme and its association with MGMT co-expression and overall survival. *PLOS ONE* (under review).

## License

This code is released under the MIT License. See LICENSE file for details.

---
