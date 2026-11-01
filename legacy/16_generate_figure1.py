#!/usr/bin/env python3
"""
Generate Figure 1: Expression Analysis
Boxplots showing GPR153 and MGMT expression in GBM vs Normal brain
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
    
    return data

def create_expression_boxplot(data, gene, output_path):
    """Create boxplot for gene expression"""
    
    # Filter for samples with the gene
    if gene not in data.columns:
        log_message(f"Gene {gene} not found")
        return False
    
    # Prepare data - but we only have GBM samples now
    # Check if we have both groups
    groups = data['group'].unique()
    
    if len(groups) < 2:
        log_message(f"Warning: Only {groups} group(s) found. Expected GBM and Normal.")
        # Create plot anyway with available data
        plot_data = data[[gene, 'group']].copy()
        plot_data.columns = ['expression', 'group']
    else:
        plot_data = data[[gene, 'group']].copy()
        plot_data.columns = ['expression', 'group']
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Color palette
    colors = {'GBM': '#e74c3c', 'Normal': '#3498db'}
    available_colors = [colors.get(g, '#95a5a6') for g in plot_data['group'].unique()]
    
    # Boxplot
    sns.boxplot(data=plot_data, x='group', y='expression', 
                palette=available_colors, ax=ax, width=0.5)
    
    # Add individual points
    sns.stripplot(data=plot_data, x='group', y='expression', 
                  color='black', alpha=0.3, size=3, ax=ax)
    
    # Statistical test if we have both groups
    if len(groups) >= 2:
        gbm_vals = plot_data[plot_data['group'] == 'GBM']['expression'].values
        normal_vals = plot_data[plot_data['group'] == 'Normal']['expression'].values
        
        if len(normal_vals) > 0:
            t_stat, p_val = stats.ttest_ind(gbm_vals, normal_vals)
            
            # Add significance annotation
            y_max = plot_data['expression'].max()
            y_min = plot_data['expression'].min()
            y_range = y_max - y_min
            
            # Significance stars
            if p_val < 0.001:
                sig_text = '***'
            elif p_val < 0.01:
                sig_text = '**'
            elif p_val < 0.05:
                sig_text = '*'
            else:
                sig_text = 'n.s.'
            
            # Draw significance bracket
            x1, x2 = 0, 1
            y = y_max + y_range * 0.05
            h = y_range * 0.02
            
            ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], 'k-', linewidth=1.5)
            ax.text((x1+x2)/2, y+h, sig_text, ha='center', va='bottom', fontsize=14, fontweight='bold')
            ax.text((x1+x2)/2, y+h*2.5, f'p = {p_val:.2e}', ha='center', va='bottom', fontsize=10)
    
    # Labels and title
    ax.set_xlabel('Group', fontsize=14, fontweight='bold')
    ax.set_ylabel(f'{gene} Expression (log2)', fontsize=14, fontweight='bold')
    ax.set_title(f'{gene} Expression in Glioblastoma', fontsize=16, fontweight='bold', pad=20)
    
    # Sample sizes
    for i, group in enumerate(plot_data['group'].unique()):
        n = (plot_data['group'] == group).sum()
        ax.text(i, ax.get_ylim()[0], f'n={n}', ha='center', va='top', fontsize=10)
    
    # Styling
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    return True

def create_multi_gene_panel(data, genes, output_path):
    """Create multi-panel figure with multiple genes"""
    
    available_genes = [g for g in genes if g in data.columns]
    
    if len(available_genes) == 0:
        log_message("No genes available for multi-panel")
        return False
    
    n_genes = len(available_genes)
    n_cols = 3
    n_rows = (n_genes + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
    axes = axes.flatten() if n_genes > 1 else [axes]
    
    colors = {'GBM': '#e74c3c', 'Normal': '#3498db'}
    
    for idx, gene in enumerate(available_genes):
        ax = axes[idx]
        
        plot_data = data[[gene, 'group']].copy()
        plot_data.columns = ['expression', 'group']
        
        # Get colors for available groups
        available_colors = [colors.get(g, '#95a5a6') for g in plot_data['group'].unique()]
        
        # Boxplot
        sns.boxplot(data=plot_data, x='group', y='expression',
                   palette=available_colors, ax=ax, width=0.5)
        
        sns.stripplot(data=plot_data, x='group', y='expression',
                     color='black', alpha=0.2, size=2, ax=ax)
        
        # Title
        ax.set_title(gene, fontsize=14, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('Expression (log2)', fontsize=11)
        
        # Sample size
        for i, group in enumerate(plot_data['group'].unique()):
            n = (plot_data['group'] == group).sum()
            ax.text(i, ax.get_ylim()[0], f'n={n}', ha='center', va='top', fontsize=9)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Remove extra subplots
    for idx in range(len(available_genes), len(axes)):
        fig.delaxes(axes[idx])
    
    plt.suptitle('Gene Expression in Glioblastoma', fontsize=18, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    return True

def main():
    """Main figure generation"""
    log_message("Starting Figure 1 generation")
    print("=" * 70)
    print("FIGURE 1: EXPRESSION ANALYSIS")
    print("=" * 70)
    
    # Load data
    data = load_data()
    
    # Create output directory
    output_dir = 'figures/publication'
    os.makedirs(output_dir, exist_ok=True)
    
    # Figure 1A: GPR153 expression
    print("\nGenerating Figure 1A: GPR153 expression...")
    success = create_expression_boxplot(
        data, 'GPR153',
        f'{output_dir}/figure1a_gpr153_expression.pdf'
    )
    if success:
        print("  ✓ Saved: figure1a_gpr153_expression.pdf")
        print("  ✓ Saved: figure1a_gpr153_expression.png")
    
    # Figure 1B: MGMT expression
    print("\nGenerating Figure 1B: MGMT expression...")
    success = create_expression_boxplot(
        data, 'MGMT',
        f'{output_dir}/figure1b_mgmt_expression.pdf'
    )
    if success:
        print("  ✓ Saved: figure1b_mgmt_expression.pdf")
        print("  ✓ Saved: figure1b_mgmt_expression.png")
    
    # Figure 1C: Multi-gene panel
    print("\nGenerating Figure 1C: Multi-gene panel...")
    genes = ['GPR153', 'MGMT', 'TP53', 'EGFR', 'IDH1', 'PTEN']
    success = create_multi_gene_panel(
        data, genes,
        f'{output_dir}/figure1c_multi_gene_expression.pdf'
    )
    if success:
        print("  ✓ Saved: figure1c_multi_gene_expression.pdf")
        print("  ✓ Saved: figure1c_multi_gene_expression.png")
    
    print("\n" + "=" * 70)
    log_message("Figure 1 generation complete")
    print("=" * 70)

if __name__ == '__main__':
    main()
