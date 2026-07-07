#!/usr/bin/env python3
"""
Generate Figure 5: Machine Learning Results
Visualize clustering, PCA, and predictive model performance
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

def log_message(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")

def load_pca_data():
    """Load PCA results"""
    pca_path = 'data/results/ml/pca_results.csv'
    var_path = 'data/results/ml/pca_variance.csv'
    
    if not os.path.exists(pca_path):
        return None, None
    
    pca_data = pd.read_csv(pca_path, index_col=0)
    
    if os.path.exists(var_path):
        var_data = pd.read_csv(var_path)
    else:
        var_data = None
    
    return pca_data, var_data

def load_cluster_data():
    """Load cluster assignments"""
    cluster_path = 'data/results/ml/cluster_assignments.csv'
    char_path = 'data/results/ml/cluster_characteristics.csv'
    
    if not os.path.exists(cluster_path):
        return None, None
    
    cluster_data = pd.read_csv(cluster_path, index_col=0)
    
    if os.path.exists(char_path):
        char_data = pd.read_csv(char_path)
    else:
        char_data = None
    
    return cluster_data, char_data

def load_model_data():
    """Load predictive model results"""
    comparison_path = 'data/results/ml/model_comparison.csv'
    importance_path = 'data/results/ml/rf_feature_importance.csv'
    
    comparison = None
    importance = None
    
    if os.path.exists(comparison_path):
        comparison = pd.read_csv(comparison_path)
    
    if os.path.exists(importance_path):
        importance = pd.read_csv(importance_path)
    
    return comparison, importance

def create_pca_plot(pca_data, cluster_data, output_path):
    """Create PCA scatter plot colored by clusters"""
    
    if pca_data is None:
        return False
    
    # Merge with cluster assignments if available
    if cluster_data is not None:
        plot_data = pca_data.copy()
        plot_data['cluster'] = cluster_data['cluster'].values
    else:
        plot_data = pca_data.copy()
        plot_data['cluster'] = 0
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Color palette
    n_clusters = plot_data['cluster'].nunique()
    colors = sns.color_palette('Set2', n_clusters)
    
    # Scatter plot
    for i, cluster in enumerate(sorted(plot_data['cluster'].unique())):
        cluster_data = plot_data[plot_data['cluster'] == cluster]
        ax.scatter(cluster_data['PC1'], cluster_data['PC2'],
                  c=[colors[i]], label=f'Cluster {cluster}',
                  s=80, alpha=0.6, edgecolors='black', linewidth=0.5)
    
    # Labels
    ax.set_xlabel('PC1 (70.8% variance)', fontsize=13, fontweight='bold')
    ax.set_ylabel('PC2 (5.3% variance)', fontsize=13, fontweight='bold')
    ax.set_title('PCA of GBM Samples', fontsize=16, fontweight='bold', pad=20)
    
    # Legend
    ax.legend(loc='best', fontsize=11, framealpha=0.9)
    
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

def create_cluster_heatmap(char_data, output_path):
    """Create heatmap of cluster characteristics"""
    
    if char_data is None or len(char_data) == 0:
        return False
    
    # Select gene columns
    gene_cols = [col for col in char_data.columns if '_mean' in col]
    
    if len(gene_cols) == 0:
        return False
    
    # Prepare data for heatmap
    heatmap_data = char_data[gene_cols].T
    heatmap_data.columns = [f"Cluster {c}" for c in char_data['cluster']]
    
    # Clean row names
    heatmap_data.index = [idx.replace('_mean', '') for idx in heatmap_data.index]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Heatmap
    sns.heatmap(heatmap_data, annot=True, fmt='.2f', cmap='RdYlBu_r',
                center=heatmap_data.values.mean(), cbar_kws={'label': 'Mean Expression'},
                linewidths=1, linecolor='white', ax=ax)
    
    # Title
    ax.set_title('Gene Expression by Cluster', fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Cluster', fontsize=13, fontweight='bold')
    ax.set_ylabel('Gene', fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    
    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    return True

def create_feature_importance_plot(importance_data, output_path):
    """Create feature importance bar plot"""
    
    if importance_data is None or len(importance_data) == 0:
        return False
    
    # Sort by importance
    plot_data = importance_data.sort_values('importance', ascending=True)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Horizontal bar plot
    colors = ['#e74c3c' if f in ['GPR153', 'MGMT'] else '#3498db' 
              for f in plot_data['feature']]
    
    bars = ax.barh(range(len(plot_data)), plot_data['importance'],
                   color=colors, edgecolor='black', linewidth=0.5)
    
    # Y-axis
    ax.set_yticks(range(len(plot_data)))
    ax.set_yticklabels(plot_data['feature'], fontsize=12, fontweight='bold')
    
    # Labels
    ax.set_xlabel('Feature Importance', fontsize=13, fontweight='bold')
    ax.set_title('Random Forest Feature Importance\n(Survival Prediction)', 
                 fontsize=16, fontweight='bold', pad=20)
    
    # Styling
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    return True

def create_model_comparison(comparison_data, output_path):
    """Create model comparison bar plot"""
    
    if comparison_data is None or len(comparison_data) == 0:
        return False
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Accuracy
    ax = axes[0]
    bars = ax.bar(range(len(comparison_data)), comparison_data['accuracy'],
                  color=['#3498db', '#e74c3c'], edgecolor='black', linewidth=1.5, alpha=0.8)
    ax.set_xticks(range(len(comparison_data)))
    ax.set_xticklabels(comparison_data['model'], fontsize=12, fontweight='bold')
    ax.set_ylabel('Accuracy', fontsize=13, fontweight='bold')
    ax.set_title('Model Accuracy', fontsize=15, fontweight='bold')
    ax.set_ylim(0, 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add values on bars
    for i, (bar, val) in enumerate(zip(bars, comparison_data['accuracy'])):
        ax.text(i, val + 0.02, f'{val:.3f}', ha='center', va='bottom', 
                fontsize=11, fontweight='bold')
    
    # AUC
    ax = axes[1]
    bars = ax.bar(range(len(comparison_data)), comparison_data['auc'],
                  color=['#3498db', '#e74c3c'], edgecolor='black', linewidth=1.5, alpha=0.8)
    ax.set_xticks(range(len(comparison_data)))
    ax.set_xticklabels(comparison_data['model'], fontsize=12, fontweight='bold')
    ax.set_ylabel('AUC', fontsize=13, fontweight='bold')
    ax.set_title('Model AUC', fontsize=15, fontweight='bold')
    ax.set_ylim(0, 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add values on bars
    for i, (bar, val) in enumerate(zip(bars, comparison_data['auc'])):
        ax.text(i, val + 0.02, f'{val:.3f}', ha='center', va='bottom',
                fontsize=11, fontweight='bold')
    
    plt.suptitle('Predictive Model Performance\n(1-Year Survival Classification)',
                 fontsize=17, fontweight='bold', y=1.05)
    plt.tight_layout()
    
    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    return True

def main():
    """Main figure generation"""
    log_message("Starting Figure 5 generation")
    print("=" * 70)
    print("FIGURE 5: MACHINE LEARNING RESULTS")
    print("=" * 70)
    
    # Create output directory
    output_dir = 'figures/publication'
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    pca_data, var_data = load_pca_data()
    cluster_data, char_data = load_cluster_data()
    comparison_data, importance_data = load_model_data()
    
    # Figure 5A: PCA plot
    print("\nGenerating Figure 5A: PCA visualization...")
    success = create_pca_plot(
        pca_data, cluster_data,
        f'{output_dir}/figure5a_pca_clusters.pdf'
    )
    if success:
        print("  ✓ Saved: figure5a_pca_clusters.pdf")
        print("  ✓ Saved: figure5a_pca_clusters.png")
    else:
        print("  ✗ Could not generate PCA plot")
    
    # Figure 5B: Cluster characteristics heatmap
    print("\nGenerating Figure 5B: Cluster characteristics...")
    success = create_cluster_heatmap(
        char_data,
        f'{output_dir}/figure5b_cluster_heatmap.pdf'
    )
    if success:
        print("  ✓ Saved: figure5b_cluster_heatmap.pdf")
        print("  ✓ Saved: figure5b_cluster_heatmap.png")
    else:
        print("  ✗ Could not generate cluster heatmap")
    
    # Figure 5C: Feature importance
    print("\nGenerating Figure 5C: Feature importance...")
    success = create_feature_importance_plot(
        importance_data,
        f'{output_dir}/figure5c_feature_importance.pdf'
    )
    if success:
        print("  ✓ Saved: figure5c_feature_importance.pdf")
        print("  ✓ Saved: figure5c_feature_importance.png")
    else:
        print("  ✗ Could not generate feature importance plot")
    
    # Figure 5D: Model comparison
    print("\nGenerating Figure 5D: Model comparison...")
    success = create_model_comparison(
        comparison_data,
        f'{output_dir}/figure5d_model_comparison.pdf'
    )
    if success:
        print("  ✓ Saved: figure5d_model_comparison.pdf")
        print("  ✓ Saved: figure5d_model_comparison.png")
    else:
        print("  ✗ Could not generate model comparison")
    
    print("\n" + "=" * 70)
    log_message("Figure 5 generation complete")
    print("=" * 70)

if __name__ == '__main__':
    main()
