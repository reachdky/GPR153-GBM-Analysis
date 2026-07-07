#!/usr/bin/env python3
"""
Quality control checks for processed data
Generate QC report and visualizations
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.path.expanduser('~/gpr153_project'))
from scripts.utils import log_message, create_output_dir, save_figure


def check_missing_values(data: pd.DataFrame, dataset_name: str):
    print(f"\n{dataset_name} - Missing Values:")
    missing = data.isnull().sum()
    missing_pct = 100 * missing / len(data)
    if missing.sum() == 0:
        print("  ✓ No missing values")
    else:
        print(f"  Columns with missing values:")
        for col in missing[missing > 0].index:
            print(f"    {col}: {missing[col]} ({missing_pct[col]:.1f}%)")
    return missing


def check_expression_distribution(data: pd.DataFrame, dataset_name: str, output_dir: str):
    print(f"\n{dataset_name} - Expression Distribution:")
    mean_expr = data.mean(axis=1)
    print(f"  Mean expression across samples: Min {mean_expr.min():.2f}, Max {mean_expr.max():.2f}, Median {mean_expr.median():.2f}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].hist(mean_expr, bins=50, edgecolor='black', alpha=0.7)
    axes[0].set_xlabel('Mean Expression (log2 TPM)', fontweight='bold')
    axes[0].set_ylabel('Number of Samples', fontweight='bold')
    axes[0].set_title(f'{dataset_name} - Mean Expression Distribution', fontweight='bold')
    axes[0].grid(alpha=0.3)

    sample_subset = data.iloc[:, :min(50, data.shape[1])].T
    axes[1].boxplot([sample_subset.iloc[i] for i in range(len(sample_subset))], showfliers=False)
    axes[1].set_xlabel('Samples (first 50)', fontweight='bold')
    axes[1].set_ylabel('Expression (log2 TPM)', fontweight='bold')
    axes[1].set_title(f'{dataset_name} - Expression Range', fontweight='bold')
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    save_figure(fig, f'qc_{dataset_name.lower()}_distribution', output_dir=output_dir)
    plt.close()


def check_genes_of_interest(data: pd.DataFrame, dataset_name: str):
    print(f"\n{dataset_name} - Genes of Interest:")
    genes = ['GPR153', 'MGMT', 'TP53', 'EGFR', 'IDH1']
    for gene in genes:
        if gene in data.columns:
            expr = data[gene]
            print(f"  ✓ {gene}: Mean {expr.mean():.2f}, SD {expr.std():.2f}, Range [{expr.min():.2f}, {expr.max():.2f}]")
        else:
            print(f"  ✗ {gene}: NOT FOUND")


def check_sample_correlation(data: pd.DataFrame, dataset_name: str, output_dir: str):
    print(f"\n{dataset_name} - Sample Correlation:")
    n_samples = min(50, len(data))
    sample_subset = data.iloc[:n_samples]
    corr_matrix = sample_subset.T.corr()
    avg_corr = corr_matrix.values[np.triu_indices_from(corr_matrix.values, k=1)].mean()
    print(f"  Average pairwise correlation: {avg_corr:.3f}")

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(corr_matrix, cmap='coolwarm', center=0, square=True, ax=ax, cbar_kws={'label': 'Correlation'})
    ax.set_title(f'{dataset_name} - Sample Correlation (n={n_samples})', fontweight='bold', pad=20)
    plt.tight_layout()
    save_figure(fig, f'qc_{dataset_name.lower()}_correlation', output_dir=output_dir)
    plt.close()


def check_clinical_data(clinical_file: str, dataset_name: str):
    if not os.path.exists(clinical_file):
        print(f"\n{dataset_name} - Clinical data not found")
        return

    print(f"\n{dataset_name} - Clinical Data:")
    data = pd.read_csv(clinical_file, index_col=0)
    print(f"  Patients: {len(data)}, Variables: {len(data.columns)}")

    key_vars = ['age', 'sex', 'survival_months', 'deceased']
    for var in key_vars:
        if var in data.columns:
            if var == 'age':
                print(f"  ✓ Age: {data[var].mean():.1f} ± {data[var].std():.1f} years")
            elif var == 'sex':
                print(f"  ✓ Sex: {data[var].value_counts().to_dict()}")
            elif var == 'survival_months':
                print(f"  ✓ Survival: median {data[var].median():.1f} months")
            elif var == 'deceased':
                n_dead = data[var].sum()
                pct = 100 * n_dead / len(data)
                print(f"  ✓ Deceased: {n_dead} ({pct:.1f}%)")
        else:
            print(f"  ✗ {var}: NOT FOUND")


def generate_qc_report(processed_dir: str, output_dir: str):
    print("\n" + "="*70)
    print("GENERATING QC REPORT")
    print("="*70)

    tcga_expr_file = os.path.join(processed_dir, 'tcga_gbm_expression_normalized.csv')
    if os.path.exists(tcga_expr_file):
        tcga_expr = pd.read_csv(tcga_expr_file, index_col=0)
        check_missing_values(tcga_expr, "TCGA GBM")
        check_expression_distribution(tcga_expr, "TCGA GBM", output_dir)
        check_genes_of_interest(tcga_expr, "TCGA GBM")
        check_sample_correlation(tcga_expr, "TCGA GBM", output_dir)

    tcga_clin_file = os.path.join(processed_dir, 'tcga_gbm_clinical_clean.csv')
    check_clinical_data(tcga_clin_file, "TCGA GBM")

    gtex_expr_file = os.path.join(processed_dir, 'gtex_brain_expression_normalized.csv')
    if os.path.exists(gtex_expr_file):
        gtex_expr = pd.read_csv(gtex_expr_file, index_col=0)
        check_missing_values(gtex_expr, "GTEx Brain")
        check_expression_distribution(gtex_expr, "GTEx Brain", output_dir)
        check_genes_of_interest(gtex_expr, "GTEx Brain")

    combined_file = os.path.join(processed_dir, 'combined_gbm_normal_expression.csv')
    if os.path.exists(combined_file):
        combined = pd.read_csv(combined_file, index_col=0)
        print(f"\nCombined Dataset Shape: {combined.shape}")
        print(f"  GBM samples: {(combined['group'] == 'GBM').sum()}")
        print(f"  Normal samples: {(combined['group'] == 'Normal').sum()}")


def main():
    log_message("Starting quality control checks")
    base_dir = os.path.expanduser('~/gpr153_project')
    processed_dir = os.path.join(base_dir, 'data', 'processed')
    qc_output_dir = os.path.join(base_dir, 'figures', 'supplementary')
    os.makedirs(qc_output_dir, exist_ok=True)

    print("="*70)
    print("QUALITY CONTROL PIPELINE")
    print("="*70)
    generate_qc_report(processed_dir, qc_output_dir)

    print("\n" + "="*70)
    print("QUALITY CONTROL COMPLETE")
    print("="*70)
    print(f"\nQC figures saved to: {qc_output_dir}/")
    log_message("Quality control complete")


if __name__ == "__main__":
    main()
