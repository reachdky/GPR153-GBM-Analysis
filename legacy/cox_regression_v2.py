"""
Multivariate Cox Proportional Hazards Regression — Version 2
GPR153 Project — TCGA-GBM Cohort

Covariate rationale:
  - GPR153 expression: primary variable of interest
  - Age at diagnosis:  universally accepted GBM prognostic covariate (included for validity)
  - MGMT expression:  EXCLUDED — collinear with GPR153 by hypothesis (GPR153 regulates MGMT);
                       including it would absorb GPR153's signal and misrepresent the model
  - MGMT methylation, MMR, BER/GSC: ABSENT from TCGA clinical data — not available

Run from: /home/ypatel/gpr153_project
Usage:    python scripts/cox_regression_v2.py
Outputs:  data/results/survival/cox_v2_results.csv
          figures/survival/cox_v2_forest_plot.png
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

# ── 1. LOAD DATA ─────────────────────────────────────────────────────────────

print("Loading data...")

scatter = pd.read_csv('data/results/correlation/gpr153_mgmt_scatter.csv', index_col=0)
scatter = scatter.reset_index()
scatter.columns = ['full_barcode', 'GPR153', 'MGMT']
scatter['bcr_patient_barcode'] = scatter['full_barcode'].str[:12]
scatter = scatter.groupby('bcr_patient_barcode')[['GPR153']].mean().reset_index()

clin = pd.read_csv('data/processed/tcga_gbm_clinical_clean.csv')

# ── 2. MERGE ──────────────────────────────────────────────────────────────────

df = pd.merge(scatter, clin, on='bcr_patient_barcode', how='inner')
print(f"Merged: {df.shape[0]} patients")

# ── 3. SURVIVAL ENDPOINT ──────────────────────────────────────────────────────

df['os_days'] = np.where(
    df['vital_status'] == 'Dead',
    df['days_to_death'],
    df['days_to_last_follow_up']
)
df['event'] = (df['vital_status'] == 'Dead').astype(int)
df = df.dropna(subset=['os_days', 'event'])
df = df[df['os_days'] > 0]
print(f"After survival filtering: {df.shape[0]} patients, {int(df['event'].sum())} events")

# ── 4. COVARIATES ─────────────────────────────────────────────────────────────

df['age'] = df['age_at_diagnosis'].fillna(df['age_at_index'])

# Z-score both continuous covariates so HRs are comparable (per-SD change)
for col in ['GPR153', 'age']:
    df[col + '_z'] = (df[col] - df[col].mean()) / df[col].std()

# ── 5. MODEL DEFINITIONS ──────────────────────────────────────────────────────
#
#  Model A: GPR153 alone (univariate Cox — directly comparable to log-rank)
#  Model B: GPR153 + Age (primary multivariate model — age is only valid non-redundant covariate)
#
base = ['os_days', 'event']
dfA = df[base + ['GPR153_z']].dropna()
dfB = df[base + ['GPR153_z', 'age_z']].dropna()

models = {
    'A_GPR153_univariate':  dfA,
    'B_GPR153_age':         dfB,
}

# ── 6. FIT ────────────────────────────────────────────────────────────────────

print("\n" + "="*60)
print("COX PROPORTIONAL HAZARDS RESULTS")
print("="*60)

results_rows = []
fitted_models = {}

for name, data in models.items():
    covars = [c for c in data.columns if c not in base]
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(data, duration_col='os_days', event_col='event')
    fitted_models[name] = cph

    print(f"\n── Model {name} (n={len(data)}) ──")
    cph.print_summary(decimals=4)

    summary = cph.summary
    for covar in covars:
        row = summary.loc[covar]
        results_rows.append({
            'model':     name,
            'covariate': covar,
            'coef':      row['coef'],
            'HR':        row['exp(coef)'],
            'HR_lower':  row['exp(coef) lower 95%'],
            'HR_upper':  row['exp(coef) upper 95%'],
            'p':         row['p'],
            'n':         len(data),
        })

# ── 7. PH ASSUMPTION TEST (Model B) ──────────────────────────────────────────

print("\n── Schoenfeld Residuals Test (Model B) ──")
try:
    ph = proportional_hazard_test(fitted_models['B_GPR153_age'], dfB, time_transform='rank')
    ph.print_summary(decimals=4)
except Exception as e:
    print(f"PH test skipped: {e}")

# ── 8. SAVE CSV ───────────────────────────────────────────────────────────────

results_df = pd.DataFrame(results_rows)
results_df['sig'] = results_df['p'].apply(
    lambda p: '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'ns'))
)
out_csv = 'data/results/survival/cox_v2_results.csv'
results_df.to_csv(out_csv, index=False)
print(f"\nResults saved to {out_csv}")
print(results_df.to_string(index=False))

# ── 9. FOREST PLOT ────────────────────────────────────────────────────────────

print("\nGenerating forest plot...")

plot_rows = results_df[results_df['model'] == 'B_GPR153_age'].copy()

label_map = {
    'GPR153_z': 'GPR153 expression\n(z-scored)',
    'age_z':    'Age at diagnosis\n(z-scored)',
}
plot_rows['label'] = plot_rows['covariate'].map(label_map)
colors = {'GPR153_z': '#C0392B', 'age_z': '#27AE60'}

fig, ax = plt.subplots(figsize=(9, 3.5))

for i, (_, row) in enumerate(plot_rows.iterrows()):
    color = colors.get(row['covariate'], '#555555')
    ax.plot([row['HR_lower'], row['HR_upper']], [i, i],
            color=color, linewidth=2.5, solid_capstyle='round')
    ax.scatter(row['HR'], i, color=color, s=130, zorder=5)
    ax.text(row['HR_upper'] + 0.04, i,
            f"HR={row['HR']:.3f}  95% CI [{row['HR_lower']:.3f}–{row['HR_upper']:.3f}]  "
            f"p={row['p']:.4f}  {row['sig']}",
            va='center', fontsize=9.5)

ax.axvline(x=1.0, color='black', linestyle='--', linewidth=1, alpha=0.6)
ax.set_yticks(range(len(plot_rows)))
ax.set_yticklabels(plot_rows['label'].tolist(), fontsize=11)
ax.set_xlabel('Hazard Ratio (95% CI)', fontsize=11)
ax.set_title(
    f'Cox Proportional Hazards — GPR153 Adjusted for Age\n'
    f'TCGA-GBM  |  n={plot_rows["n"].iloc[0]}  |  '
    f'Covariates: GPR153 expression, Age at diagnosis',
    fontsize=11, fontweight='bold'
)
ax.set_xlim(0.5, 2.8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

import os
os.makedirs('figures/survival', exist_ok=True)
out_fig = 'figures/survival/cox_v2_forest_plot.png'
plt.savefig(out_fig, dpi=150, bbox_inches='tight')
plt.close()
print(f"Forest plot saved to {out_fig}")

# ── 10. SUMMARY ───────────────────────────────────────────────────────────────

gpr = results_df[(results_df['model'] == 'B_GPR153_age') &
                 (results_df['covariate'] == 'GPR153_z')].iloc[0]

print("\n" + "="*60)
print("SUMMARY FOR REPORTING")
print("="*60)
print(f"""
Model: GPR153 + Age at diagnosis  |  n={gpr['n']}

Note: MGMT expression excluded (collinear with GPR153 by hypothesis).
      MGMT methylation, MMR, BER/GSC absent from TCGA clinical data.

GPR153 (z-scored):
  HR  = {gpr['HR']:.3f}
  95% CI [{gpr['HR_lower']:.3f} – {gpr['HR_upper']:.3f}]
  p   = {gpr['p']:.4f}  {gpr['sig']}
""")

if gpr['p'] < 0.05:
    print("  GPR153 expression is an INDEPENDENT predictor of overall survival")
    print("  after adjustment for patient age (p < 0.05).")
    print("  Each 1-SD increase in GPR153 expression is associated with a")
    print(f"  {(gpr['HR']-1)*100:.1f}% increase in instantaneous hazard of death.")
else:
    print("  GPR153 does not reach conventional significance as an independent")
    print("  predictor after age adjustment. Univariate log-rank p=0.002 remains")
    print("  the primary survival finding. Cox result is reported transparently.")

print("\nDone.")
