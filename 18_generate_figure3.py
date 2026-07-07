#!/usr/bin/env python3
"""
Generate Figure 3: Survival Analysis
Kaplan-Meier curves for GPR153 and MGMT
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

def log_message(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")

def load_survival_data(gene):
    """Load survival data for a gene"""
    data_path = f'data/results/survival/{gene.lower()}_survival_data.csv'
    
    if not os.path.exists(data_path):
        log_message(f"Survival data not found: {data_path}")
        return None
    
    data = pd.read_csv(data_path, index_col=0)
    return data

def calculate_survival_function(times, events):
    """Calculate Kaplan-Meier survival function"""
    # Sort by time
    sorted_indices = np.argsort(times)
    times = times[sorted_indices]
    events = events[sorted_indices]
    
    # Unique event times
    unique_times = np.unique(times)
    
    n_at_risk = len(times)
    survival_prob = []
    survival_times = [0]
    survival_prob.append(1.0)
    
    current_survival = 1.0
    
    for t in unique_times:
        # Number of events at time t
        n_events = np.sum((times == t) & (events == 1))
        
        # Number at risk at time t
        n_at_risk_t = np.sum(times >= t)
        
        if n_events > 0 and n_at_risk_t > 0:
            # Update survival probability
            current_survival *= (1 - n_events / n_at_risk_t)
        
        survival_times.append(t)
        survival_prob.append(current_survival)
    
    return np.array(survival_times), np.array(survival_prob)

def create_km_curve(data, gene, output_path):
    """Create Kaplan-Meier survival curve"""
    
    if data is None:
        return False
    
    # Separate by expression group
    high_data = data[data['expression_group'] == 'High']
    low_data = data[data['expression_group'] == 'Low']
    
    if len(high_data) == 0 or len(low_data) == 0:
        log_message(f"Insufficient data for {gene}")
        return False
    
    # Calculate KM curves
    times_high, surv_high = calculate_survival_function(
        high_data['survival_days'].values,
        high_data['event'].values
    )
    
    times_low, surv_low = calculate_survival_function(
        low_data['survival_days'].values,
        low_data['event'].values
    )
    
    # Load statistics
    stats_path = f'data/results/survival/{gene.lower()}_survival_stats.csv'
    if os.path.exists(stats_path):
        stats_df = pd.read_csv(stats_path)
        p_value = stats_df['p_value'].values[0]
    else:
        p_value = None
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Plot KM curves
    ax.step(times_high, surv_high, where='post', linewidth=2.5, 
            color='#e74c3c', label=f'High {gene} (n={len(high_data)})')
    ax.step(times_low, surv_low, where='post', linewidth=2.5,
            color='#3498db', label=f'Low {gene} (n={len(low_data)})')
    
    # Confidence intervals (simplified - using standard error)
    # This is approximate
    
    # Median survival lines
    high_median = high_data['survival_days'].median()
    low_median = low_data['survival_days'].median()
    
    ax.axvline(high_median, color='#e74c3c', linestyle='--', alpha=0.5, linewidth=1.5)
    ax.axvline(low_median, color='#3498db', linestyle='--', alpha=0.5, linewidth=1.5)
    
    # P-value text
    if p_value is not None:
        if p_value < 0.001:
            p_text = f'p < 0.001'
        else:
            p_text = f'p = {p_value:.4f}'
        
        ax.text(0.98, 0.98, p_text, transform=ax.transAxes,
                fontsize=14, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                fontweight='bold')
    
    # Labels
    ax.set_xlabel('Time (days)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Survival Probability', fontsize=14, fontweight='bold')
    ax.set_title(f'{gene} Expression and Overall Survival in GBM', 
                 fontsize=16, fontweight='bold', pad=20)
    
    # Legend
    ax.legend(loc='upper right', fontsize=12, framealpha=0.9)
    
    # Styling
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, max(times_high.max(), times_low.max()) * 1.05)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Add median survival text
    text_y = 0.15
    ax.text(0.02, text_y, f'Median survival:\nHigh: {high_median:.0f} days\nLow: {low_median:.0f} days',
            transform=ax.transAxes, fontsize=11, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    return True

def create_combined_survival(output_path):
    """Create combined survival figure with both genes"""
    
    # Load data
    gpr153_data = load_survival_data('GPR153')
    mgmt_data = load_survival_data('MGMT')
    
    if gpr153_data is None or mgmt_data is None:
        log_message("Cannot create combined survival - data missing")
        return False
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # GPR153
    ax = axes[0]
    high_data = gpr153_data[gpr153_data['expression_group'] == 'High']
    low_data = gpr153_data[gpr153_data['expression_group'] == 'Low']
    
    times_high, surv_high = calculate_survival_function(
        high_data['survival_days'].values, high_data['event'].values
    )
    times_low, surv_low = calculate_survival_function(
        low_data['survival_days'].values, low_data['event'].values
    )
    
    ax.step(times_high, surv_high, where='post', linewidth=2.5,
            color='#e74c3c', label=f'High GPR153 (n={len(high_data)})')
    ax.step(times_low, surv_low, where='post', linewidth=2.5,
            color='#3498db', label=f'Low GPR153 (n={len(low_data)})')
    
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Survival Probability', fontsize=13, fontweight='bold')
    ax.set_title('GPR153 Expression', fontsize=15, fontweight='bold', pad=15)
    ax.legend(loc='upper right', fontsize=11)
    ax.set_ylim(0, 1.05)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # MGMT
    ax = axes[1]
    high_data = mgmt_data[mgmt_data['expression_group'] == 'High']
    low_data = mgmt_data[mgmt_data['expression_group'] == 'Low']
    
    times_high, surv_high = calculate_survival_function(
        high_data['survival_days'].values, high_data['event'].values
    )
    times_low, surv_low = calculate_survival_function(
        low_data['survival_days'].values, low_data['event'].values
    )
    
    ax.step(times_high, surv_high, where='post', linewidth=2.5,
            color='#e74c3c', label=f'High MGMT (n={len(high_data)})')
    ax.step(times_low, surv_low, where='post', linewidth=2.5,
            color='#3498db', label=f'Low MGMT (n={len(low_data)})')
    
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Survival Probability', fontsize=13, fontweight='bold')
    ax.set_title('MGMT Expression', fontsize=15, fontweight='bold', pad=15)
    ax.legend(loc='upper right', fontsize=11)
    ax.set_ylim(0, 1.05)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.suptitle('Overall Survival by Gene Expression in GBM', 
                 fontsize=17, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    return True

def main():
    """Main figure generation"""
    log_message("Starting Figure 3 generation")
    print("=" * 70)
    print("FIGURE 3: SURVIVAL ANALYSIS")
    print("=" * 70)
    
    # Create output directory
    output_dir = 'figures/publication'
    os.makedirs(output_dir, exist_ok=True)
    
    # Figure 3A: GPR153 survival
    print("\nGenerating Figure 3A: GPR153 survival curve...")
    gpr153_data = load_survival_data('GPR153')
    success = create_km_curve(
        gpr153_data, 'GPR153',
        f'{output_dir}/figure3a_gpr153_survival.pdf'
    )
    if success:
        print("  ✓ Saved: figure3a_gpr153_survival.pdf")
        print("  ✓ Saved: figure3a_gpr153_survival.png")
    
    # Figure 3B: MGMT survival
    print("\nGenerating Figure 3B: MGMT survival curve...")
    mgmt_data = load_survival_data('MGMT')
    success = create_km_curve(
        mgmt_data, 'MGMT',
        f'{output_dir}/figure3b_mgmt_survival.pdf'
    )
    if success:
        print("  ✓ Saved: figure3b_mgmt_survival.pdf")
        print("  ✓ Saved: figure3b_mgmt_survival.png")
    
    # Figure 3C: Combined survival
    print("\nGenerating Figure 3C: Combined survival curves...")
    success = create_combined_survival(
        f'{output_dir}/figure3c_combined_survival.pdf'
    )
    if success:
        print("  ✓ Saved: figure3c_combined_survival.pdf")
        print("  ✓ Saved: figure3c_combined_survival.png")
    
    print("\n" + "=" * 70)
    log_message("Figure 3 generation complete")
    print("=" * 70)

if __name__ == '__main__':
    main()
