#!/usr/bin/env python3
"""
Differential Expression Analysis
Identifies differentially expressed genes between GBM and normal brain
"""

import os
import sys
import pandas as pd
import numpy as np
from scipy import stats
from datetime import datetime

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

def log_message(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")

def load_data():
    """Load combined expression data"""
    data_path = 'data/processed/combined_gbm_normal_expression.csv'
    
    log_message("Loading combined expression data...")
    data = pd.read_csv(data_path, index_col=0)
    
    log_message(f"Loaded data: {data.shape}")
    return data

def differential_expression_analysis(data, fdr_threshold=0.05, fold_change_threshold=2.0):
    """Perform differential expression analysis"""
    
    # Separate GBM and Normal samples
    gbm_data = data[data['group'] == 'GBM'].drop('group', axis=1)
    normal_data = data[data['group'] == 'Normal'].drop('group', axis=1)
    
    log_message(f"GBM samples: {len(gbm_data)}")
    log_message(f"Normal samples: {len(normal_data)}")
    
    results = []
    
    genes = gbm_data.columns
    total_genes = len(genes)
    
    log_message(f"Testing {total_genes} genes...")
    
    for i, gene in enumerate(genes):
        if (i + 1) % 5000 == 0:
            log_message(f"  Processed {i+1}/{total_genes} genes")
        
        gbm_expr = gbm_data[gene].values
        normal_expr = normal_data[gene].values
        
        # Remove NaN
        gbm_expr = gbm_expr[~np.isnan(gbm_expr)]
        normal_expr = normal_expr[~np.isnan(normal_expr)]
        
        if len(gbm_expr) < 3 or len(normal_expr) < 3:
            continue
        
        # Calculate fold change (log2)
        gbm_mean = np.mean(gbm_expr)
        normal_mean = np.mean(normal_expr)
        
        if normal_mean == 0:
            continue
        
        fold_change = gbm_mean / normal_mean
        log2_fold_change = np.log2(fold_change) if fold_change > 0 else 0
        
        # T-test
        t_stat, p_value = stats.ttest_ind(gbm_expr, normal_expr)
        
        results.append({
            'gene': gene,
            'gbm_mean': gbm_mean,
            'normal_mean': normal_mean,
            'fold_change': fold_change,
            'log2_fold_change': log2_fold_change,
            't_statistic': t_stat,
            'p_value': p_value,
            'n_gbm': len(gbm_expr),
            'n_normal': len(normal_expr)
        })
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # FDR correction (Benjamini-Hochberg)
    from scipy.stats import false_discovery_control
    results_df['fdr'] = false_discovery_control(results_df['p_value'].values)
    
    # Classify genes
    results_df['significant'] = (results_df['fdr'] < fdr_threshold) & \
                                 (np.abs(results_df['log2_fold_change']) >= np.log2(fold_change_threshold))
    
    results_df['regulation'] = 'unchanged'
    results_df.loc[(results_df['significant']) & (results_df['log2_fold_change'] > 0), 'regulation'] = 'upregulated'
    results_df.loc[(results_df['significant']) & (results_df['log2_fold_change'] < 0), 'regulation'] = 'downregulated'
    
    # Sort by p-value
    results_df = results_df.sort_values('p_value')
    
    return results_df

def main():
    """Main differential expression pipeline"""
    log_message("Starting differential expression analysis")
    print("=" * 70)
    print("DIFFERENTIAL EXPRESSION ANALYSIS")
    print("=" * 70)
    
    # Load data
    data = load_data()
    
    # Create output directory
    output_dir = 'data/results/differential'
    os.makedirs(output_dir, exist_ok=True)
    
    # Run analysis
    print("\n" + "=" * 70)
    print("PERFORMING DIFFERENTIAL EXPRESSION ANALYSIS")
    print("=" * 70)
    
    de_results = differential_expression_analysis(data, fdr_threshold=0.05, fold_change_threshold=2.0)
    
    # Summary statistics
    n_upregulated = (de_results['regulation'] == 'upregulated').sum()
    n_downregulated = (de_results['regulation'] == 'downregulated').sum()
    n_unchanged = (de_results['regulation'] == 'unchanged').sum()
    
    print(f"\nResults:")
    print(f"  Total genes tested: {len(de_results)}")
    print(f"  Upregulated in GBM: {n_upregulated}")
    print(f"  Downregulated in GBM: {n_downregulated}")
    print(f"  Unchanged: {n_unchanged}")
    
    # Top upregulated genes
    print(f"\nTop 10 upregulated genes:")
    top_up = de_results[de_results['regulation'] == 'upregulated'].head(10)
    for idx, row in top_up.iterrows():
        print(f"  {row['gene']}: FC = {row['fold_change']:.2f}, FDR = {row['fdr']:.2e}")
    
    # Top downregulated genes
    print(f"\nTop 10 downregulated genes:")
    top_down = de_results[de_results['regulation'] == 'downregulated'].head(10)
    for idx, row in top_down.iterrows():
        print(f"  {row['gene']}: FC = {row['fold_change']:.2f}, FDR = {row['fdr']:.2e}")
    
    # Check GPR153 and MGMT
    print(f"\nKey genes:")
    for gene in ['GPR153', 'MGMT']:
        gene_data = de_results[de_results['gene'] == gene]
        if len(gene_data) > 0:
            row = gene_data.iloc[0]
            print(f"  {gene}: FC = {row['fold_change']:.2f}, log2FC = {row['log2_fold_change']:.2f}, " +
                  f"FDR = {row['fdr']:.2e}, {row['regulation']}")
    
    # Save results
    de_results.to_csv(f'{output_dir}/differential_expression_all.csv', index=False)
    
    # Save significant genes only
    sig_results = de_results[de_results['significant']]
    sig_results.to_csv(f'{output_dir}/differential_expression_significant.csv', index=False)
    
    print(f"\n✓ Saved: differential_expression_all.csv")
    print(f"✓ Saved: differential_expression_significant.csv")
    
    print("\n" + "=" * 70)
    log_message("Differential expression analysis complete")
    print("=" * 70)

if __name__ == '__main__':
    main()
