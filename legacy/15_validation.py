#!/usr/bin/env python3
"""
Validation Analysis
Validates findings using independent datasets (CGGA, GEO)
"""

import os
import sys
import pandas as pd
import numpy as np
from scipy import stats
from datetime import datetime
import gzip

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

def log_message(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")

def validate_correlation(data, gene1, gene2, dataset_name):
    """Validate correlation in independent dataset"""
    if gene1 not in data.columns or gene2 not in data.columns:
        return None
    
    expr1 = data[gene1].values
    expr2 = data[gene2].values
    
    # Remove NaN
    mask = ~(np.isnan(expr1) | np.isnan(expr2))
    expr1 = expr1[mask]
    expr2 = expr2[mask]
    
    if len(expr1) < 3:
        return None
    
    # Pearson correlation
    r, p = stats.pearsonr(expr1, expr2)
    
    return {
        'dataset': dataset_name,
        'gene1': gene1,
        'gene2': gene2,
        'n_samples': len(expr1),
        'correlation': r,
        'p_value': p
    }

def load_cgga_data():
    """Load CGGA dataset"""
    cgga_dir = 'data/raw/cgga/CGGA_IDH_A_RNAseq_RSEM_20250915'
    
    if not os.path.exists(cgga_dir):
        return None
    
    # Look for expression file
    files = os.listdir(cgga_dir)
    expr_files = [f for f in files if 'exp' in f.lower() or 'rnaseq' in f.lower() or 'rsem' in f.lower()]
    
    if not expr_files:
        return None
    
    try:
        data = pd.read_csv(os.path.join(cgga_dir, expr_files[0]), sep='\t', index_col=0)
        log_message(f"Loaded CGGA: {data.shape}")
        return data
    except:
        try:
            data = pd.read_csv(os.path.join(cgga_dir, expr_files[0]), index_col=0)
            log_message(f"Loaded CGGA: {data.shape}")
            return data
        except Exception as e:
            log_message(f"Failed to load CGGA: {e}")
            return None

def load_geo_datasets():
    """Load GEO datasets"""
    datasets = {}
    
    geo_datasets = ['GSE107850', 'GSE13041', 'GSE4290']
    
    for geo_id in geo_datasets:
        geo_dir = f'data/raw/geo/{geo_id}'
        
        if not os.path.exists(geo_dir):
            continue
        
        # Look for series matrix file
        files = os.listdir(geo_dir)
        matrix_files = [f for f in files if 'series_matrix' in f.lower() or 'expression' in f.lower()]
        
        if matrix_files:
            try:
                file_path = os.path.join(geo_dir, matrix_files[0])
                
                # Handle gzipped files
                if file_path.endswith('.gz'):
                    with gzip.open(file_path, 'rt') as f:
                        # Skip header lines starting with !
                        lines = [line for line in f if not line.startswith('!')]
                        # Create DataFrame from remaining lines
                        from io import StringIO
                        data = pd.read_csv(StringIO(''.join(lines)), sep='\t', index_col=0)
                else:
                    data = pd.read_csv(file_path, sep='\t', index_col=0, comment='!')
                
                datasets[geo_id] = data
                log_message(f"Loaded {geo_id}: {data.shape}")
            except Exception as e:
                log_message(f"Failed to load {geo_id}: {e}")
                continue
    
    return datasets

def main():
    """Main validation pipeline"""
    log_message("Starting validation analysis")
    print("=" * 70)
    print("VALIDATION ANALYSIS")
    print("=" * 70)
    
    # Create output directory
    output_dir = 'data/results/validation'
    os.makedirs(output_dir, exist_ok=True)
    
    # Load validation datasets
    print("\n" + "=" * 70)
    print("LOADING VALIDATION DATASETS")
    print("=" * 70)
    
    datasets = {}
    
    # Try CGGA
    cgga_data = load_cgga_data()
    if cgga_data is not None:
        datasets['CGGA'] = cgga_data
    
    # Try GEO
    geo_datasets = load_geo_datasets()
    datasets.update(geo_datasets)
    
    if len(datasets) == 0:
        print("\n⚠️  Could not load validation datasets")
        print("Data exists but may need processing first.")
        print("\nNote: Validation is optional. Primary TCGA analysis is complete.")
        return
    
    print(f"\nLoaded {len(datasets)} validation dataset(s)")
    
    # Validate GPR153-MGMT correlation
    print("\n" + "=" * 70)
    print("VALIDATING GPR153-MGMT CORRELATION")
    print("=" * 70)
    
    # Load original TCGA result
    tcga_corr_path = 'data/results/correlation/gpr153_mgmt_correlation.csv'
    if os.path.exists(tcga_corr_path):
        tcga_corr = pd.read_csv(tcga_corr_path)
        print(f"\nTCGA (discovery):")
        print(f"  Pearson r = {tcga_corr['pearson_r'].values[0]:.4f}")
        print(f"  p-value = {tcga_corr['pearson_p'].values[0]:.2e}")
    
    # Validate in other datasets
    validation_results = []
    
    for dataset_name, data in datasets.items():
        result = validate_correlation(data, 'GPR153', 'MGMT', dataset_name)
        
        if result:
            validation_results.append(result)
            print(f"\n{dataset_name} (validation):")
            print(f"  n = {result['n_samples']}")
            print(f"  Pearson r = {result['correlation']:.4f}")
            print(f"  p-value = {result['p_value']:.2e}")
            
            if result['p_value'] < 0.05:
                print(f"  ✓ Significant correlation replicated")
            else:
                print(f"  ✗ Correlation not significant")
        else:
            print(f"\n{dataset_name}: GPR153 or MGMT not found in dataset")
    
    # Save validation results
    if validation_results:
        validation_df = pd.DataFrame(validation_results)
        validation_df.to_csv(f'{output_dir}/correlation_validation.csv', index=False)
        print(f"\n✓ Saved: correlation_validation.csv")
        
        # Meta-analysis summary
        print("\n" + "=" * 70)
        print("META-ANALYSIS SUMMARY")
        print("=" * 70)
        
        all_correlations = [r['correlation'] for r in validation_results]
        all_pvalues = [r['p_value'] for r in validation_results]
        
        print(f"\nAcross {len(validation_results)} validation dataset(s):")
        print(f"  Mean correlation: {np.mean(all_correlations):.4f}")
        print(f"  Median correlation: {np.median(all_correlations):.4f}")
        print(f"  Range: [{np.min(all_correlations):.4f}, {np.max(all_correlations):.4f}]")
        
        significant_count = sum(1 for p in all_pvalues if p < 0.05)
        print(f"\n  Significant in {significant_count}/{len(validation_results)} datasets")
        
        # Save summary
        summary = pd.DataFrame({
            'metric': ['mean_correlation', 'median_correlation', 'min_correlation', 
                      'max_correlation', 'n_significant', 'n_total'],
            'value': [np.mean(all_correlations), np.median(all_correlations),
                     np.min(all_correlations), np.max(all_correlations),
                     significant_count, len(validation_results)]
        })
        summary.to_csv(f'{output_dir}/validation_summary.csv', index=False)
        print(f"\n✓ Saved: validation_summary.csv")
    else:
        print("\n⚠️  Could not validate - genes not found in datasets")
        print("Note: Datasets may use different gene identifiers (Ensembl IDs, etc.)")
    
    print("\n" + "=" * 70)
    log_message("Validation analysis complete")
    print("=" * 70)

if __name__ == '__main__':
    main()
