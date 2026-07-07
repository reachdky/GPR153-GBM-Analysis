#!/usr/bin/env python3
"""
Download TCGA GBM data using TCGAbiolinks
This script downloads RNA-seq expression and clinical data from TCGA
"""

import os
import sys

# Add project root to path
sys.path.append(os.path.expanduser('~/gpr153_project'))

from scripts.utils import create_output_dir, log_message

# Note: This script uses R/Bioconductor
# We'll create a wrapper that calls R script

def main():
    """Download TCGA GBM data"""
    
    log_message("Starting TCGA GBM data download")
    
    # Create output directories
    output_dir = os.path.join(os.path.expanduser('~/gpr153_project'), 
                             'data', 'raw', 'tcga')
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'expression'), exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'clinical'), exist_ok=True)
    
    # Create R script for download
    r_script = f"""
    library(TCGAbiolinks)
    library(SummarizedExperiment)
    
    # Set output directory
    output_dir <- "{output_dir}"
    
    # Query GBM RNA-seq data
    query_exp <- GDCquery(
        project = "TCGA-GBM",
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts"
    )
    
    # Download data
    GDCdownload(query_exp, directory = file.path(output_dir, "expression"))
    
    # Prepare data
    data_exp <- GDCprepare(query_exp, directory = file.path(output_dir, "expression"))
    
    # Extract expression matrix (TPM values)
    expr_matrix <- assay(data_exp, "tpm_unstrand")
    
    # Save expression data
    write.csv(expr_matrix, 
              file.path(output_dir, "expression", "tcga_gbm_expression_tpm.csv"))
    
    # Query clinical data
    query_clin <- GDCquery(
        project = "TCGA-GBM",
        data.category = "Clinical",
        file.type = "xml"
    )
    
    # Download clinical
    GDCdownload(query_clin, directory = file.path(output_dir, "clinical"))
    
    # Prepare clinical
    clinical <- GDCprepare_clinic(query_clin, clinical.info = "patient")
    
    # Save clinical data
    write.csv(clinical, 
              file.path(output_dir, "clinical", "tcga_gbm_clinical.csv"),
              row.names = FALSE)
    
    print("TCGA GBM download complete!")
    """
    
    # Save R script
    r_script_path = os.path.join(output_dir, 'download_tcga.R')
    with open(r_script_path, 'w') as f:
        f.write(r_script)
    
    print(f"R script created: {r_script_path}")
    print(f"\nTo download TCGA data, run:")
    print(f"  Rscript {r_script_path}")
    print(f"\nNote: This will download ~5-10 GB of data and may take 2-4 hours")
    
    # Execute R script
    print("\nExecuting R script...")
    os.system(f"Rscript {r_script_path}")
    
    log_message("TCGA GBM download complete")

if __name__ == "__main__":
    main()

