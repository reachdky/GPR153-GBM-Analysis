#!/usr/bin/env python3
"""
Analyze GPR153 and MGMT expression in GBM vs Normal
Perform differential expression analysis
"""

import os
import sys
import pandas as pd
import numpy as np

sys.path.append(os.path.expanduser('~/gpr153_project'))
from scripts.utils import (
    log_message, 
    save_results, 
    perform_ttest,
    calculate_fold_change,
    load_processed_data
)


def analyze_gene_expression(data: pd.DataFrame, gene: str) -> dict:
    """
    Analyze single gene expression in GBM vs Normal
    
    Parameters
    ----------
    data : pd.DataFrame
        Combined expression data with 'group' column
    gene : str
        Gene name to analyze
    
    Returns
    -------
    dict
        Analysis results
    """
    
    if gene not in data.columns:
        print(f"  ✗ {gene} not found in data")
        return None
    
    print(f"\nAnalyzing {gene}...")
    
    # Split by group
    gbm_expr = data.loc[data['group'] == 'GBM', gene].values
    normal_expr = data.loc[data['group'] == 'Normal', gene].values
    
    # Calculate statistics
    gbm_mean = np.mean(gbm_expr)
    gbm_std = np.std(gbm_expr)
    normal_mean = np.mean(normal_expr)
    normal_std = np.std(normal_expr)
    
    # Perform t-test
    t_stat, p_val, cohens_d = perform_ttest(gbm_expr, normal_expr)
    
    # Calculate fold change
    fc, log2fc = calculate_fold_change(gbm_expr, normal_expr)
    
    # Print results
    print(f"  GBM (n={len(gbm_expr)}): {gbm_mean:.3f} ± {gbm_std:.3f}")
    print(f"  Normal (n={len(normal_expr)}): {normal_mean:.3f} ± {normal_std:.3f}")
    print(f"  Fold Change: {fc:.3f}")
    print(f"  Log2 FC: {log2fc:.3f}")
    print(f"  t-statistic: {t_stat:.3f}")
    print(f"  p-value: {p_val:.2e}")
    print(f"  Cohen's d: {cohens_d:.3f}")
    
    # Interpret effect size
    if abs(cohens_d) < 0.2:
        effect = "negligible"
    elif abs(cohens_d) < 0.5:
        effect = "small"
    elif abs(cohens_d) < 0.8:
        effect = "medium"
    else:
        effect = "large"
    
    print(f"  Effect size: {effect}")
    
    # Interpret significance
    if p_val < 0.001:
        sig = "highly significant (p < 0.001)"
    elif p_val < 0.01:
        sig = "very significant (p < 0.01)"
    elif p_val < 0.05:
        sig = "significant (p < 0.05)"
    else:
        sig = "not significant"
    
    print(f"  Significance: {sig}")
    
    results = {
        'gene': gene,
        'gbm_n': len(gbm_expr),
        'gbm_mean': gbm_mean,
        'gbm_std': gbm_std,
        'gbm_median': np.median(gbm_expr),
        'normal_n': len(normal_expr),
        'normal_mean': normal_mean,
        'normal_std': normal_std,
        'normal_median': np.median(normal_expr),
        'fold_change': fc,
        'log2_fold_change': log2fc,
        't_statistic': t_stat,
        'p_value': p_val,
        'cohens_d': cohens_d,
        'effect_size': effect,
        'significance': sig
    }
    
    return results


def compare_multiple_genes(data: pd.DataFrame, genes: list) -> pd.DataFrame:
    """
    Compare expression of multiple genes
    
    Parameters
    ----------
    data : pd.DataFrame
        Combined expression data
    genes : list
        List of genes to analyze
    
    Returns
    -------
    pd.DataFrame
        Comparison results
    """
    
    print("\n" + "="*70)
    print("MULTI-GENE EXPRESSION ANALYSIS")
    print("="*70)
    
    results_list = []
    
    for gene in genes:
        result = analyze_gene_expression(data, gene)
        if result is not None:
            results_list.append(result)
    
    results_df = pd.DataFrame(results_list)
    
    # Sort by p-value
    results_df = results_df.sort_values('p_value')
    
    return results_df


def analyze_expression_by_subgroup(data: pd.DataFrame, 
                                   clinical: pd.DataFrame,
                                   gene: str,
                                   subgroup_col: str) -> pd.DataFrame:
    """
    Analyze gene expression by clinical subgroups
    
    Parameters
    ----------
    data : pd.DataFrame
        Expression data
    clinical : pd.DataFrame
        Clinical data
    gene : str
        Gene to analyze
    subgroup_col : str
        Clinical variable for subgrouping (e.g., 'sex', 'age_group')
    
    Returns
    -------
    pd.DataFrame
        Subgroup analysis results
    """
    
    print(f"\nAnalyzing {gene} by {subgroup_col}...")
    
    # Merge expression and clinical
    merged = data.merge(clinical[[subgroup_col]], left_index=True, right_index=True, how='inner')
    
    # Get unique subgroups
    subgroups = merged[subgroup_col].unique()
    
    results = []
    
    for subgroup in subgroups:
        subgroup_data = merged.loc[merged[subgroup_col] == subgroup, gene]
        
        result = {
            'gene': gene,
            'subgroup_variable': subgroup_col,
            'subgroup': subgroup,
            'n': len(subgroup_data),
            'mean': subgroup_data.mean(),
            'std': subgroup_data.std(),
            'median': subgroup_data.median()
        }
        
        results.append(result)
        
        print(f"  {subgroup} (n={result['n']}): {result['mean']:.3f} ± {result['std']:.3f}")
    
    return pd.DataFrame(results)


def main():
    """Main expression analysis pipeline"""
    
    log_message("Starting expression analysis")
    
    base_dir = os.path.expanduser('~/gpr153_project')
    processed_dir = os.path.join(base_dir, 'data', 'processed')
    
    print("="*70)
    print("EXPRESSION ANALYSIS PIPELINE")
    print("="*70)
    
    # Load combined dataset
    print("\nLoading combined GBM vs Normal dataset...")
    combined_file = os.path.join(processed_dir, 'combined_gbm_normal_expression.csv')
    
    if not os.path.exists(combined_file):
        print(f"✗ Combined dataset not found: {combined_file}")
        print("Run 06_normalize_expression.py first")
        return
    
    data = pd.read_csv(combined_file, index_col=0)
    print(f"✓ Loaded: {data.shape}")
    
    # Analyze primary genes of interest
    print("\n" + "="*70)
    print("PRIMARY GENES OF INTEREST")
    print("="*70)
    
    gpr153_results = analyze_gene_expression(data, 'GPR153')
    mgmt_results = analyze_gene_expression(data, 'MGMT')
    
    # Save individual results
    if gpr153_results:
        save_results(
            pd.DataFrame([gpr153_results]),
            'gpr153_expression_analysis.csv',
            'expression'
        )
    
    if mgmt_results:
        save_results(
            pd.DataFrame([mgmt_results]),
            'mgmt_expression_analysis.csv',
            'expression'
        )
    
    # Analyze additional genes for context
    additional_genes = ['TP53', 'EGFR', 'IDH1', 'PTEN', 'RB1']
    
    all_genes = ['GPR153', 'MGMT'] + additional_genes
    
    multi_gene_results = compare_multiple_genes(data, all_genes)
    
    # Save multi-gene results
    save_results(
        multi_gene_results,
        'multi_gene_expression_analysis.csv',
        'expression'
    )
    
    # Print summary table
    print("\n" + "="*70)
    print("EXPRESSION ANALYSIS SUMMARY")
    print("="*70)
    
    summary = multi_gene_results[['gene', 'fold_change', 'log2_fold_change', 
                                  'p_value', 'cohens_d', 'significance']]
    
    print("\n" + summary.to_string(index=False))
    
    # Subgroup analysis (if clinical data available)
    clinical_file = os.path.join(processed_dir, 'tcga_gbm_clinical_clean.csv')
    
    if os.path.exists(clinical_file):
        print("\n" + "="*70)
        print("SUBGROUP ANALYSIS")
        print("="*70)
        
        clinical = pd.read_csv(clinical_file, index_col=0)
        
        # Filter to GBM samples only
        gbm_data = data[data['group'] == 'GBM'].copy()
        
        # Create age groups
        if 'age' in clinical.columns:
            clinical['age_group'] = pd.cut(
                clinical['age'],
                bins=[0, 50, 65, 100],
                labels=['<50', '50-65', '>65']
            )
        
        # Analyze by sex
        if 'sex' in clinical.columns:
            sex_results = analyze_expression_by_subgroup(
                gbm_data, clinical, 'GPR153', 'sex'
            )
            save_results(sex_results, 'gpr153_by_sex.csv', 'expression')
        
        # Analyze by age group
        if 'age_group' in clinical.columns:
            age_results = analyze_expression_by_subgroup(
                gbm_data, clinical, 'GPR153', 'age_group'
            )
            save_results(age_results, 'gpr153_by_age.csv', 'expression')
    
    print("\n" + "="*70)
    print("EXPRESSION ANALYSIS COMPLETE")
    print("="*70)
    
    # Key findings summary
    print("\nKEY FINDINGS:")
    
    if gpr153_results:
        print(f"\n1. GPR153 Expression:")
        print(f"   - Fold change in GBM: {gpr153_results['fold_change']:.2f}x")
        print(f"   - p-value: {gpr153_results['p_value']:.2e}")
        print(f"   - {gpr153_results['significance']}")
    
    if mgmt_results:
        print(f"\n2. MGMT Expression:")
        print(f"   - Fold change in GBM: {mgmt_results['fold_change']:.2f}x")
        print(f"   - p-value: {mgmt_results['p_value']:.2e}")
        print(f"   - {mgmt_results['significance']}")
    
    log_message("Expression analysis complete")


if __name__ == "__main__":
    main()
