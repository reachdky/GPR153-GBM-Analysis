#!/usr/bin/env python3
"""
Normalize TCGA and GTEx expression data and create combined dataset
"""

import os
import pandas as pd
import numpy as np

from scripts.utils import log_message

def normalize_expression(df):
    """Log2-transform and z-score normalize by gene"""
    return df.apply(lambda x: (np.log2(x + 1) - np.mean(np.log2(x + 1))) / np.std(np.log2(x + 1)), axis=0)

def main():
    log_message("Loading expression data...")

    base_dir = os.path.expanduser('~/gpr153_project')
    processed_dir = os.path.join(base_dir, 'data', 'processed')

    tcga_file = os.path.join(processed_dir, 'tcga_gbm_expression_clean.csv')
    gtex_file = os.path.join(processed_dir, 'gtex_brain_expression_clean.csv')

    tcga = pd.read_csv(tcga_file, index_col=0)
    gtex = pd.read_csv(gtex_file, index_col=0)

    print(f"✓ TCGA shape: {tcga.shape}")
    print(f"✓ GTEx shape: {gtex.shape}")

    # Normalize
    print("Normalizing TCGA and GTEx data...")
    tcga_norm = normalize_expression(tcga)
    gtex_norm = normalize_expression(gtex)

    # Add 'group' column
    tcga_norm['group'] = 'GBM'
    gtex_norm['group'] = 'Normal'

    # Align genes
    common_genes = tcga_norm.columns.intersection(gtex_norm.columns)
    tcga_norm = tcga_norm[common_genes.tolist() + ['group']]
    gtex_norm = gtex_norm[common_genes.tolist() + ['group']]

    # Combine datasets
    combined_df = pd.concat([tcga_norm, gtex_norm], axis=0)

    # Transpose so genes are columns, samples are rows
    combined_df = combined_df.reset_index(drop=True)

    # Save normalized and combined datasets
    tcga_norm.to_csv(os.path.join(processed_dir, 'tcga_gbm_expression_normalized.csv'))
    gtex_norm.to_csv(os.path.join(processed_dir, 'gtex_expression_normalized.csv'))
    combined_df.to_csv(os.path.join(processed_dir, 'combined_gbm_normal_expression.csv'))

    print("\n======================================================================")
    print("✓ Combined dataset written to:", os.path.join(processed_dir, 'combined_gbm_normal_expression.csv'))
    print(f"✓ Shape: {combined_df.shape}\n")
    print("Sample counts:")
    print(combined_df['group'].value_counts())
    print("======================================================================")
    print("READY FOR SCRIPT 8")
    print("======================================================================")

if __name__ == "__main__":
    main()
