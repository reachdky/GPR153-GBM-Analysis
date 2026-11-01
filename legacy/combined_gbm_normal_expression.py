#!/usr/bin/env python3
"""
Create combined GBM and normal brain expression dataset
"""

import os
import pandas as pd

# Paths
base_dir = os.path.expanduser('~/gpr153_project')
processed_dir = os.path.join(base_dir, 'data', 'processed')

# TCGA GBM expression
tcga_file = os.path.join(processed_dir, 'tcga_gbm_expression_normalized.csv')
print("Loading TCGA GBM expression data...")
tcga = pd.read_csv(tcga_file, index_col=0)

# GTEx normal brain expression
gtex_file = os.path.join(processed_dir, 'gtex_expression_normalized.csv')
print("Loading GTEx Normal brain expression data...")
gtex = pd.read_csv(gtex_file, index_col=0)

# Make sure genes match
common_genes = tcga.index.intersection(gtex.index)
tcga = tcga.loc[common_genes]
gtex = gtex.loc[common_genes]

# Add 'group' column
tcga = tcga.T.copy()
tcga['group'] = 'GBM'

gtex = gtex.T.copy()
gtex['group'] = 'Normal'

# Combine datasets
combined = pd.concat([tcga, gtex], axis=0)

# Save combined dataset
combined_file = os.path.join(processed_dir, 'combined_gbm_normal_expression.csv')
combined.to_csv(combined_file)

print(f"Combined dataset shape: {combined.shape}")
print(f"✓ Combined dataset saved to {combined_file}")
