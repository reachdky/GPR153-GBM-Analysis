#!/usr/bin/env python3
"""
Survival Analysis - GPR153 and MGMT
Analyzes patient survival based on gene expression levels
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
    """Load expression and clinical data"""
    expr_path = 'data/processed/combined_gbm_normal_expression.csv'
    clinical_path = 'data/processed/tcga_gbm_clinical_clean.csv'
    
    log_message("Loading expression data...")
    expr_data = pd.read_csv(expr_path, index_col=0)
    
    log_message("Loading clinical data...")
    clinical_data = pd.read_csv(clinical_path, index_col=0)
    
    # Filter for GBM samples only
    gbm_expr = expr_data[expr_data['group'] == 'GBM'].copy()
    gbm_expr = gbm_expr.drop('group', axis=1)
    
    log_message(f"Expression data: {gbm_expr.shape}")
    log_message(f"Clinical data: {clinical_data.shape}")
    
    # Create survival days column (use days_to_death or days_to_last_follow_up)
    clinical_data['survival_days'] = clinical_data['days_to_death'].fillna(
        clinical_data['days_to_last_follow_up']
    )
    
    # Create event column (1 = death, 0 = alive/censored)
    clinical_data['event'] = (clinical_data['vital_status'] == 'Dead').astype(int)
    
    # Use bcr_patient_barcode to match with expression data
    # Expression sample IDs are like: TCGA-06-6390-01A-11R-A96S-41
    # Clinical IDs are like: TCGA-06-6390
    # Extract patient ID from expression sample IDs
    
    def extract_patient_id(sample_id):
        """Extract patient ID from TCGA sample ID"""
        parts = sample_id.split('-')
        if len(parts) >= 3:
            return f"{parts[0]}-{parts[1]}-{parts[2]}"
        return sample_id
    
    # Map expression samples to patient IDs
    expr_patient_ids = [extract_patient_id(sid) for sid in gbm_expr.index]
    gbm_expr.index = expr_patient_ids
    
    # Use bcr_patient_barcode as index for clinical
    if 'bcr_patient_barcode' in clinical_data.columns:
        clinical_data = clinical_data.set_index('bcr_patient_barcode')
    
    # Find common samples
    common_samples = gbm_expr.index.intersection(clinical_data.index)
    
    # Handle duplicates by keeping first occurrence
    gbm_expr = gbm_expr[~gbm_expr.index.duplicated(keep='first')]
    clinical_data = clinical_data[~clinical_data.index.duplicated(keep='first')]
    
    # Re-find common samples after deduplication
    common_samples = gbm_expr.index.intersection(clinical_data.index)
    
    gbm_expr = gbm_expr.loc[common_samples]
    clinical_data = clinical_data.loc[common_samples]
    
    log_message(f"Aligned data: {len(common_samples)} samples")
    
    return gbm_expr, clinical_data

def categorize_expression(expression, method='median'):
    """Categorize expression into high/low groups"""
    if method == 'median':
        threshold = np.median(expression)
    elif method == 'tertile':
        threshold = np.percentile(expression, 33.33)
    else:
        threshold = np.mean(expression)
    
    categories = ['High' if x >= threshold else 'Low' for x in expression]
    return categories, threshold

def calculate_survival_statistics(clinical_data, gene_expr, gene_name):
    """Calculate survival statistics for a gene"""
    # Remove samples with missing survival data
    valid_mask = clinical_data['survival_days'].notna() & \
                 clinical_data['event'].notna() & \
                 gene_expr.notna()
    
    if valid_mask.sum() < 10:
        return None
    
    survival_days = clinical_data.loc[valid_mask, 'survival_days'].values
    events = clinical_data.loc[valid_mask, 'event'].values
    expression = gene_expr.loc[valid_mask].values
    
    # Categorize expression
    expr_categories, threshold = categorize_expression(expression, method='median')
    
    # Calculate median survival for each group
    high_expr_mask = np.array(expr_categories) == 'High'
    low_expr_mask = np.array(expr_categories) == 'Low'
    
    high_survival = survival_days[high_expr_mask]
    low_survival = survival_days[low_expr_mask]
    
    # Log-rank test approximation using Mann-Whitney U test
    statistic, p_value = stats.mannwhitneyu(high_survival, low_survival, alternative='two-sided')
    
    # Calculate medians
    high_median = np.median(high_survival)
    low_median = np.median(low_survival)
    
    # Count events (deaths)
    high_deaths = events[high_expr_mask].sum()
    low_deaths = events[low_expr_mask].sum()
    
    return {
        'gene': gene_name,
        'n_samples': valid_mask.sum(),
        'n_high': high_expr_mask.sum(),
        'n_low': low_expr_mask.sum(),
        'expression_threshold': threshold,
        'high_median_survival': high_median,
        'low_median_survival': low_median,
        'high_deaths': int(high_deaths),
        'low_deaths': int(low_deaths),
        'test_statistic': statistic,
        'p_value': p_value
    }

def create_survival_data(clinical_data, gene_expr, gene_name):
    """Create detailed survival data for plotting"""
    valid_mask = clinical_data['survival_days'].notna() & \
                 clinical_data['event'].notna() & \
                 gene_expr.notna()
    
    if valid_mask.sum() < 10:
        return None
    
    survival_days = clinical_data.loc[valid_mask, 'survival_days'].values
    events = clinical_data.loc[valid_mask, 'event'].values
    expression = gene_expr.loc[valid_mask].values
    sample_ids = clinical_data.loc[valid_mask].index.values
    
    # Categorize
    expr_categories, _ = categorize_expression(expression, method='median')
    
    # Create DataFrame
    surv_df = pd.DataFrame({
        'sample_id': sample_ids,
        'survival_days': survival_days,
        'event': events,
        'expression': expression,
        'expression_group': expr_categories
    })
    
    return surv_df

def main():
    """Main survival analysis pipeline"""
    log_message("Starting survival analysis")
    print("=" * 70)
    print("SURVIVAL ANALYSIS PIPELINE")
    print("=" * 70)
    
    # Load data
    expr_data, clinical_data = load_data()
    
    # Create output directory
    output_dir = 'data/results/survival'
    os.makedirs(output_dir, exist_ok=True)
    
    # Analyze GPR153
    print("\n" + "=" * 70)
    print("GPR153 SURVIVAL ANALYSIS")
    print("=" * 70)
    
    if 'GPR153' in expr_data.columns:
        gpr153_stats = calculate_survival_statistics(
            clinical_data, expr_data['GPR153'], 'GPR153'
        )
        
        if gpr153_stats:
            print(f"  Samples analyzed: {gpr153_stats['n_samples']}")
            print(f"  High expression (n={gpr153_stats['n_high']}): {gpr153_stats['high_median_survival']:.1f} days median")
            print(f"  Low expression (n={gpr153_stats['n_low']}): {gpr153_stats['low_median_survival']:.1f} days median")
            print(f"  p-value: {gpr153_stats['p_value']:.4f}")
            
            # Save statistics
            stats_df = pd.DataFrame([gpr153_stats])
            stats_df.to_csv(f'{output_dir}/gpr153_survival_stats.csv', index=False)
            
            # Save detailed data
            surv_data = create_survival_data(clinical_data, expr_data['GPR153'], 'GPR153')
            if surv_data is not None:
                surv_data.to_csv(f'{output_dir}/gpr153_survival_data.csv', index=False)
            
            print(f"  ✓ Saved results")
        else:
            print("  ✗ Insufficient data")
    else:
        print("  ✗ GPR153 not found")
    
    # Analyze MGMT
    print("\n" + "=" * 70)
    print("MGMT SURVIVAL ANALYSIS")
    print("=" * 70)
    
    if 'MGMT' in expr_data.columns:
        mgmt_stats = calculate_survival_statistics(
            clinical_data, expr_data['MGMT'], 'MGMT'
        )
        
        if mgmt_stats:
            print(f"  Samples analyzed: {mgmt_stats['n_samples']}")
            print(f"  High expression (n={mgmt_stats['n_high']}): {mgmt_stats['high_median_survival']:.1f} days median")
            print(f"  Low expression (n={mgmt_stats['n_low']}): {mgmt_stats['low_median_survival']:.1f} days median")
            print(f"  p-value: {mgmt_stats['p_value']:.4f}")
            
            # Save
            stats_df = pd.DataFrame([mgmt_stats])
            stats_df.to_csv(f'{output_dir}/mgmt_survival_stats.csv', index=False)
            
            surv_data = create_survival_data(clinical_data, expr_data['MGMT'], 'MGMT')
            if surv_data is not None:
                surv_data.to_csv(f'{output_dir}/mgmt_survival_data.csv', index=False)
            
            print(f"  ✓ Saved results")
        else:
            print("  ✗ Insufficient data")
    else:
        print("  ✗ MGMT not found")
    
    # Multi-gene survival analysis
    print("\n" + "=" * 70)
    print("MULTI-GENE SURVIVAL ANALYSIS")
    print("=" * 70)
    
    genes_of_interest = ['GPR153', 'MGMT', 'TP53', 'EGFR', 'IDH1', 'PTEN']
    available_genes = [g for g in genes_of_interest if g in expr_data.columns]
    
    all_stats = []
    for gene in available_genes:
        stats_result = calculate_survival_statistics(clinical_data, expr_data[gene], gene)
        if stats_result:
            all_stats.append(stats_result)
            print(f"  {gene}: p = {stats_result['p_value']:.4f}")
    
    if all_stats:
        multi_stats_df = pd.DataFrame(all_stats)
        multi_stats_df = multi_stats_df.sort_values('p_value')
        multi_stats_df.to_csv(f'{output_dir}/multi_gene_survival_stats.csv', index=False)
        print(f"\n  ✓ Saved: multi_gene_survival_stats.csv")
    
    print("\n" + "=" * 70)
    log_message("Survival analysis complete")
    print("=" * 70)

if __name__ == '__main__':
    main()
