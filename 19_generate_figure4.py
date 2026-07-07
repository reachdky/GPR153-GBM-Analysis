#!/usr/bin/env python3
"""
Generate Figure 4: Pathway Analysis
Visualize enriched pathways from GPR153-correlated genes
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

def load_pathway_data(pathway_type):
    """Load pathway enrichment data"""
    
    pathway_files = {
        'kegg': 'data/results/pathways/kegg_pathways.csv',
        'go_bp': 'data/results/pathways/go_biological_process.csv',
        'reactome': 'data/results/pathways/reactome_pathways.csv'
    }
    
    if pathway_type not in pathway_files:
        return None
    
    path = pathway_files[pathway_type]
    
    if not os.path.exists(path):
        log_message(f"Pathway data not found: {path}")
        return None
    
    data = pd.read_csv(path)
    return data

def create_pathway_barplot(data, pathway_type, output_path, top_n=10):
    """Create horizontal bar plot for pathway enrichment"""
    
    if data is None or len(data) == 0:
        return False
    
    # Take top pathways
    plot_data = data.head(top_n).copy()
    
    # Clean pathway names (remove IDs, truncate if too long)
    plot_data['term_clean'] = plot_data['term'].str.replace(r'\(GO:\d+\)', '', regex=True)
    plot_data['term_clean'] = plot_data['term_clean'].str.replace(r'R-HSA-\d+', '', regex=True)
    plot_data['term_clean'] = plot_data['term_clean'].str.strip()
    
    # Truncate long names
    max_len = 50
    plot_data['term_clean'] = plot_data['term_clean'].apply(
        lambda x: x[:max_len] + '...' if len(x) > max_len else x
    )
    
    # -log10(p-value) for visualization
    plot_data['neg_log_p'] = -np.log10(plot_data['p_value'])
    
    # Sort by p-value (ascending)
    plot_data = plot_data.sort_values('p_value', ascending=True)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Color by significance
    colors = ['#e74c3c' if p < 0.01 else '#3498db' for p in plot_data['p_value']]
    
    # Horizontal bar plot
    bars = ax.barh(range(len(plot_data)), plot_data['neg_log_p'], color=colors, edgecolor='black', linewidth=0.5)
    
    # Y-axis labels
    ax.set_yticks(range(len(plot_data)))
    ax.set_yticklabels(plot_data['term_clean'], fontsize=11)
    
    # X-axis
    ax.set_xlabel('-log10(p-value)', fontsize=13, fontweight='bold')
    
    # Significance line
    ax.axvline(-np.log10(0.05), color='gray', linestyle='--', linewidth=1.5, alpha=0.7, label='p = 0.05')
    
    # Title
    titles = {
        'kegg': 'KEGG Pathway Enrichment',
        'go_bp': 'GO Biological Process Enrichment',
        'reactome': 'Reactome Pathway Enrichment'
    }
    ax.set_title(titles.get(pathway_type, 'Pathway Enrichment'), fontsize=16, fontweight='bold', pad=20)
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#e74c3c', edgecolor='black', label='p < 0.01'),
        Patch(facecolor='#3498db', edgecolor='black', label='p ≥ 0.01')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=10)
    
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

def create_pathway_dotplot(output_path):
    """Create dot plot combining multiple pathway databases"""
    
    # Load all pathway types
    kegg = load_pathway_data('kegg')
    go_bp = load_pathway_data('go_bp')
    reactome = load_pathway_data('reactome')
    
    # Combine top pathways from each
    combined = []
    
    if kegg is not None and len(kegg) > 0:
        top_kegg = kegg.head(5).copy()
        top_kegg['database'] = 'KEGG'
        combined.append(top_kegg)
    
    if go_bp is not None and len(go_bp) > 0:
        top_go = go_bp.head(5).copy()
        top_go['database'] = 'GO:BP'
        combined.append(top_go)
    
    if reactome is not None and len(reactome) > 0:
        top_reactome = reactome.head(5).copy()
        top_reactome['database'] = 'Reactome'
        combined.append(top_reactome)
    
    if len(combined) == 0:
        return False
    
    plot_data = pd.concat(combined, ignore_index=True)
    
    # Clean term names
    plot_data['term_clean'] = plot_data['term'].str.replace(r'\(GO:\d+\)', '', regex=True)
    plot_data['term_clean'] = plot_data['term_clean'].str.replace(r'R-HSA-\d+', '', regex=True)
    plot_data['term_clean'] = plot_data['term_clean'].str.strip()
    max_len = 40
    plot_data['term_clean'] = plot_data['term_clean'].apply(
        lambda x: x[:max_len] + '...' if len(x) > max_len else x
    )
    
    # -log10(p-value)
    plot_data['neg_log_p'] = -np.log10(plot_data['p_value'])
    
    # Sort
    plot_data = plot_data.sort_values(['database', 'p_value'])
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Create dot plot
    databases = plot_data['database'].unique()
    colors = {'KEGG': '#e74c3c', 'GO:BP': '#3498db', 'Reactome': '#2ecc71'}
    
    y_pos = 0
    yticks = []
    yticklabels = []
    
    for db in databases:
        db_data = plot_data[plot_data['database'] == db]
        
        for idx, row in db_data.iterrows():
            # Dot size based on number of genes
            size = row['n_genes'] * 20
            
            ax.scatter(row['neg_log_p'], y_pos, s=size, c=colors[db],
                      alpha=0.7, edgecolors='black', linewidth=1)
            
            yticks.append(y_pos)
            yticklabels.append(row['term_clean'])
            
            y_pos += 1
        
        y_pos += 0.5  # Gap between databases
    
    # Set ticks
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=10)
    
    # Labels
    ax.set_xlabel('-log10(p-value)', fontsize=13, fontweight='bold')
    ax.set_title('Enriched Pathways in GPR153-Correlated Genes', fontsize=16, fontweight='bold', pad=20)
    
    # Significance line
    ax.axvline(-np.log10(0.05), color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
    
    # Legend for databases
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=colors[db], edgecolor='black', label=db, alpha=0.7) 
                      for db in databases]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=11, title='Database')
    
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

def main():
    """Main figure generation"""
    log_message("Starting Figure 4 generation")
    print("=" * 70)
    print("FIGURE 4: PATHWAY ANALYSIS")
    print("=" * 70)
    
    # Create output directory
    output_dir = 'figures/publication'
    os.makedirs(output_dir, exist_ok=True)
    
    # Figure 4A: KEGG pathways
    print("\nGenerating Figure 4A: KEGG pathways...")
    kegg_data = load_pathway_data('kegg')
    success = create_pathway_barplot(
        kegg_data, 'kegg',
        f'{output_dir}/figure4a_kegg_pathways.pdf'
    )
    if success:
        print("  ✓ Saved: figure4a_kegg_pathways.pdf")
        print("  ✓ Saved: figure4a_kegg_pathways.png")
    
    # Figure 4B: GO Biological Process
    print("\nGenerating Figure 4B: GO Biological Process...")
    go_data = load_pathway_data('go_bp')
    success = create_pathway_barplot(
        go_data, 'go_bp',
        f'{output_dir}/figure4b_go_biological_process.pdf'
    )
    if success:
        print("  ✓ Saved: figure4b_go_biological_process.pdf")
        print("  ✓ Saved: figure4b_go_biological_process.png")
    
    # Figure 4C: Combined pathway dot plot
    print("\nGenerating Figure 4C: Combined pathway analysis...")
    success = create_pathway_dotplot(
        f'{output_dir}/figure4c_combined_pathways.pdf'
    )
    if success:
        print("  ✓ Saved: figure4c_combined_pathways.pdf")
        print("  ✓ Saved: figure4c_combined_pathways.png")
    
    print("\n" + "=" * 70)
    log_message("Figure 4 generation complete")
    print("=" * 70)

if __name__ == '__main__':
    main()
