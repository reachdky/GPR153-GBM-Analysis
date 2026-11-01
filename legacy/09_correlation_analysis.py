#!/usr/bin/env python3
"""
Correlation Analysis - GPR153 and MGMT
Analyzes correlation between GPR153 and MGMT expression in GBM samples
"""

import os
import sys
import pandas as pd
import numpy as np
from scipy import stats
from datetime import datetime

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

def log_message(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")

def load_data():
    """Load processed expression data"""
    data_path = 'data/processed/combined_gbm_normal_expression.csv'
    
    log_message("Loading combined expression data...")
    data = pd.read_csv(data_path, index_col=0)
    
    # Filter for GBM samples only
    gbm_data = data[data['group'] == 'GBM'].copy()
    gbm_data = gbm_data.drop('group', axis=1)
    
    log_message(f"Loaded GBM data: {gbm_data.shape}")
    return gbm_data

def calculate_correlation(data, gene1, gene2):
    """Calculate correlation between two genes"""
    if gene1 not in data.columns or gene2 not in data.columns:
        return None
    
    # Get expression values
    expr1 = data[gene1].values
    expr2 = data[gene2].values
    
    # Remove any NaN values
    mask = ~(np.isnan(expr1) | np.isnan(expr2))
    expr1 = expr1[mask]
    expr2 = expr2[mask]
    
    if len(expr1) < 3:
        return None
    
    # Pearson correlation
    pearson_r, pearson_p = stats.pearsonr(expr1, expr2)
    
    # Spearman correlation
    spearman_r, spearman_p = stats.spearmanr(expr1, expr2)
    
    return {
        'gene1': gene1,
        'gene2': gene2,
        'n_samples': len(expr1),
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'spearman_r': spearman_r,
        'spearman_p': spearman_p
    }

def find_correlated_genes(data, target_gene, threshold=0.3, max_genes=100):
    """Find genes correlated with target gene"""
    if target_gene not in data.columns:
        return pd.DataFrame()
    
    target_expr = data[target_gene].values
    
    correlations = []
    for gene in data.columns:
        if gene == target_gene:
            continue
        
        gene_expr = data[gene].values
        
        # Remove NaN
        mask = ~(np.isnan(target_expr) | np.isnan(gene_expr))
        if mask.sum() < 3:
            continue
        
        r, p = stats.pearsonr(target_expr[mask], gene_expr[mask])
        
        if abs(r) >= threshold:
            correlations.append({
                'gene': gene,
                'correlation': r,
                'p_value': p,
                'abs_correlation': abs(r)
            })
    
    # Convert to DataFrame and sort
    corr_df = pd.DataFrame(correlations)
    if len(corr_df) == 0:
        return corr_df
    
    corr_df = corr_df.sort_values('abs_correlation', ascending=False)
    
    # Keep top genes
    return corr_df.head(max_genes)

def main():
    """Main correlation analysis pipeline"""
    log_message("Starting correlation analysis")
    print("=" * 70)
    print("CORRELATION ANALYSIS PIPELINE")
    print("=" * 70)
    
    # Load data
    data = load_data()
    
    # Create output directory
    output_dir = 'data/results/correlation'
    os.makedirs(output_dir, exist_ok=True)
    
    # Primary correlation: GPR153 vs MGMT
    print("\n" + "=" * 70)
    print("PRIMARY CORRELATION: GPR153 vs MGMT")
    print("=" * 70)
    
    gpr153_mgmt_corr = calculate_correlation(data, 'GPR153', 'MGMT')
    
    if gpr153_mgmt_corr:
        print(f"  Samples analyzed: {gpr153_mgmt_corr['n_samples']}")
        print(f"  Pearson r: {gpr153_mgmt_corr['pearson_r']:.4f}")
        print(f"  Pearson p-value: {gpr153_mgmt_corr['pearson_p']:.2e}")
        print(f"  Spearman r: {gpr153_mgmt_corr['spearman_r']:.4f}")
        print(f"  Spearman p-value: {gpr153_mgmt_corr['spearman_p']:.2e}")
        
        # Save
        corr_df = pd.DataFrame([gpr153_mgmt_corr])
        corr_df.to_csv(f'{output_dir}/gpr153_mgmt_correlation.csv', index=False)
        print(f"  ✓ Saved: gpr153_mgmt_correlation.csv")
    else:
        print("  ✗ Could not calculate correlation")
    
    # Find genes correlated with GPR153
    print("\n" + "=" * 70)
    print("GENES CORRELATED WITH GPR153")
    print("=" * 70)
    
    gpr153_correlated = find_correlated_genes(data, 'GPR153', threshold=0.3, max_genes=100)
    
    if len(gpr153_correlated) > 0:
        print(f"  Found {len(gpr153_correlated)} genes with |r| >= 0.3")
        print(f"  Top 5 correlations:")
        for idx, row in gpr153_correlated.head(5).iterrows():
            print(f"    {row['gene']}: r = {row['correlation']:.4f}, p = {row['p_value']:.2e}")
        
        # Save
        gpr153_correlated.to_csv(f'{output_dir}/gpr153_correlated_genes.csv', index=False)
        print(f"  ✓ Saved: gpr153_correlated_genes.csv")
    else:
        print("  No significantly correlated genes found")
    
    # Find genes correlated with MGMT
    print("\n" + "=" * 70)
    print("GENES CORRELATED WITH MGMT")
    print("=" * 70)
    
    mgmt_correlated = find_correlated_genes(data, 'MGMT', threshold=0.3, max_genes=100)
    
    if len(mgmt_correlated) > 0:
        print(f"  Found {len(mgmt_correlated)} genes with |r| >= 0.3")
        print(f"  Top 5 correlations:")
        for idx, row in mgmt_correlated.head(5).iterrows():
            print(f"    {row['gene']}: r = {row['correlation']:.4f}, p = {row['p_value']:.2e}")
        
        # Save
        mgmt_correlated.to_csv(f'{output_dir}/mgmt_correlated_genes.csv', index=False)
        print(f"  ✓ Saved: mgmt_correlated_genes.csv")
    else:
        print("  No significantly correlated genes found")
    
    # Multi-gene correlation matrix
    print("\n" + "=" * 70)
    print("MULTI-GENE CORRELATION MATRIX")
    print("=" * 70)
    
    genes_of_interest = ['GPR153', 'MGMT', 'TP53', 'EGFR', 'IDH1', 'PTEN']
    available_genes = [g for g in genes_of_interest if g in data.columns]
    
    if len(available_genes) >= 2:
        corr_matrix = data[available_genes].corr(method='pearson')
        
        print(f"  Genes analyzed: {', '.join(available_genes)}")
        print("\n  Correlation matrix:")
        print(corr_matrix.to_string())
        
        # Save
        corr_matrix.to_csv(f'{output_dir}/multi_gene_correlation_matrix.csv')
        print(f"\n  ✓ Saved: multi_gene_correlation_matrix.csv")
    
    print("\n" + "=" * 70)
    log_message("Correlation analysis complete")
    print("=" * 70)

if __name__ == '__main__':
    main()
