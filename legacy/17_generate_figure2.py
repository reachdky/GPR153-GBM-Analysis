#!/usr/bin/env python3
"""
Generate Figure 2: Correlation Analysis
Scatter plots showing GPR153-MGMT correlation
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from datetime import datetime

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

def log_message(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")

def load_data():
    """Load expression data"""
    data_path = 'data/processed/combined_gbm_normal_expression.csv'
    
    log_message("Loading expression data...")
    data = pd.read_csv(data_path, index_col=0)
    
    # Filter GBM only
    gbm_data = data[data['group'] == 'GBM'].drop('group', axis=1)
    
    return gbm_data

def create_correlation_scatter(data, gene1, gene2, output_path):
    """Create scatter plot for correlation"""
    
    if gene1 not in data.columns or gene2 not in data.columns:
        log_message(f"Genes not found: {gene1}, {gene2}")
        return False
    
    x = data[gene1].values
    y = data[gene2].values
    
    # Remove NaN
    mask = ~(np.isnan(x) | np.isnan(y))
    x = x[mask]
    y = y[mask]
    
    # Calculate correlation
    r, p = stats.pearsonr(x, y)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Scatter plot
    ax.scatter(x, y, alpha=0.5, s=50, color='#3498db', edgecolors='black', linewidth=0.5)
    
    # Regression line
    z = np.polyfit(x, y, 1)
    p_func = np.poly1d(z)
    x_line = np.linspace(x.min(), x.max(), 100)
    ax.plot(x_line, p_func(x_line), 'r-', linewidth=2, label=f'Linear fit')
    
    # Confidence interval
    from scipy.stats import t as t_dist
    n = len(x)
    predict = p_func(x)
    residuals = y - predict
    s_err = np.sqrt(np.sum(residuals**2) / (n - 2))
    
    # 95% CI
    ci = 1.96 * s_err
    ax.fill_between(x_line, p_func(x_line) - ci, p_func(x_line) + ci,
                     alpha=0.2, color='red', label='95% CI')
    
    # Correlation text
    text_str = f'Pearson r = {r:.3f}\np < {p:.2e}\nn = {len(x)}'
    ax.text(0.05, 0.95, text_str, transform=ax.transAxes,
            fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Labels
    ax.set_xlabel(f'{gene1} Expression (log2)', fontsize=14, fontweight='bold')
    ax.set_ylabel(f'{gene2} Expression (log2)', fontsize=14, fontweight='bold')
    ax.set_title(f'{gene1} vs {gene2} Correlation in GBM', fontsize=16, fontweight='bold', pad=20)
    
    # Legend
    ax.legend(loc='lower right', fontsize=10)
    
    # Styling
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    return True

def create_correlation_heatmap(data, genes, output_path):
    """Create correlation heatmap"""
    
    available_genes = [g for g in genes if g in data.columns]
    
    if len(available_genes) < 2:
        log_message("Not enough genes for heatmap")
        return False
    
    # Calculate correlation matrix
    corr_matrix = data[available_genes].corr()
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Heatmap
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
    sns.heatmap(corr_matrix, mask=mask, annot=True, fmt='.3f',
                cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                square=True, linewidths=1, cbar_kws={"shrink": 0.8},
                ax=ax)
    
    # Title
    ax.set_title('Gene Expression Correlation Matrix', fontsize=16, fontweight='bold', pad=20)
    
    plt.tight_layout()
    
    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    return True

def main():
    """Main figure generation"""
    log_message("Starting Figure 2 generation")
    print("=" * 70)
    print("FIGURE 2: CORRELATION ANALYSIS")
    print("=" * 70)
    
    # Load data
    data = load_data()
    
    # Create output directory
    output_dir = 'figures/publication'
    os.makedirs(output_dir, exist_ok=True)
    
    # Figure 2A: GPR153 vs MGMT scatter
    print("\nGenerating Figure 2A: GPR153-MGMT correlation...")
    success = create_correlation_scatter(
        data, 'GPR153', 'MGMT',
        f'{output_dir}/figure2a_gpr153_mgmt_correlation.pdf'
    )
    if success:
        print("  ✓ Saved: figure2a_gpr153_mgmt_correlation.pdf")
        print("  ✓ Saved: figure2a_gpr153_mgmt_correlation.png")
    
    # Figure 2B: Correlation heatmap
    print("\nGenerating Figure 2B: Correlation heatmap...")
    genes = ['GPR153', 'MGMT', 'TP53', 'EGFR', 'IDH1', 'PTEN']
    success = create_correlation_heatmap(
        data, genes,
        f'{output_dir}/figure2b_correlation_heatmap.pdf'
    )
    if success:
        print("  ✓ Saved: figure2b_correlation_heatmap.pdf")
        print("  ✓ Saved: figure2b_correlation_heatmap.png")
    
    print("\n" + "=" * 70)
    log_message("Figure 2 generation complete")
    print("=" * 70)

if __name__ == '__main__':
    main()
