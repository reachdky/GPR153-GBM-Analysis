#!/usr/bin/env python3
"""
r03_batch_comparison.py
Compare the two TCGA-GBM RNA-seq batches (BI vs GDC2024) on clinical and
molecular variables, one aliquot per patient. No expression values are read.

Input
  revision/manifests/tcga_samples.csv   (r01)
Outputs
  revision/results/cohort/tcga_batch_comparison.csv
  revision/results/cohort/tcga_batch_os.csv

Run from the repository root:
    python pipeline/r03_batch_comparison.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import median_survival_times
from scipy import stats

MAN = Path("revision/manifests/tcga_samples.csv")
RES = Path("revision/results/cohort")


def continuous(df, col, label):
    a, b = df.loc[df.batch == "BI", col].dropna(), df.loc[df.batch == "GDC2024", col].dropna()
    p = stats.mannwhitneyu(a, b).pvalue if len(a) and len(b) else np.nan
    fmt = lambda x: f"{x.median():.1f} [{x.quantile(.25):.1f}-{x.quantile(.75):.1f}] (n={len(x)})"
    return {"variable": label, "level": "median [IQR]", "BI": fmt(a), "GDC2024": fmt(b),
            "test": "Mann-Whitney U", "p": p}


def categorical(df, col, label):
    t = pd.crosstab(df[col].fillna("missing"), df.batch)
    obs = t.drop(index="missing", errors="ignore")
    if obs.shape[0] >= 2:
        p = stats.fisher_exact(obs.values).pvalue if obs.shape == (2, 2) else stats.chi2_contingency(obs.values)[1]
        test = "Fisher exact (missing excluded)" if obs.shape == (2, 2) else "chi-square (missing excluded)"
    else:
        p, test = np.nan, ""
    rows = []
    for lvl in t.index:
        rows.append({"variable": label, "level": lvl,
                     "BI": f"{t.loc[lvl, 'BI']} ({100 * t.loc[lvl, 'BI'] / t['BI'].sum():.0f}%)",
                     "GDC2024": f"{t.loc[lvl, 'GDC2024']} ({100 * t.loc[lvl, 'GDC2024'] / t['GDC2024'].sum():.0f}%)",
                     "test": test, "p": p if lvl == t.index[0] else np.nan})
    return rows


def survival(df, cohort):
    d = df[df.os_days.notna() & (df.os_days > 0) & df.os_event.notna()]
    out = []
    for b in ["BI", "GDC2024"]:
        x = d[d.batch == b]
        km = KaplanMeierFitter().fit(x.os_days, x.os_event)
        mci = median_survival_times(km.confidence_interval_)
        fu = KaplanMeierFitter().fit(x.os_days, 1 - x.os_event)  # reverse KM for follow-up
        out.append({"cohort": cohort, "batch": b, "n": len(x), "deaths": int(x.os_event.sum()),
                    "censored": int((1 - x.os_event).sum()),
                    "median_os_days": km.median_survival_time_,
                    "median_os_ci_low": float(mci.iloc[0, 0]), "median_os_ci_high": float(mci.iloc[0, 1]),
                    "median_followup_days_reverse_km": fu.median_survival_time_ if np.isfinite(fu.median_survival_time_) else "not reached"})
    a, b = d[d.batch == "BI"], d[d.batch == "GDC2024"]
    p = logrank_test(a.os_days, b.os_days, a.os_event, b.os_event).p_value
    for o in out:
        o["logrank_p_between_batches"] = p
    return out


def main():
    m = pd.read_csv(MAN)
    rows, surv = [], []
    for cohort, flag in [("sensitivity (all primary)", "in_screen_sensitivity"), ("primary (IDH-wildtype)", "in_screen_primary")]:
        df = m[m[flag]].copy()
        part = [{"variable": "n patients", "level": "", "BI": str((df.batch == "BI").sum()),
                 "GDC2024": str((df.batch == "GDC2024").sum()), "test": "", "p": np.nan},
                continuous(df, "age_at_diagnosis_years", "age at diagnosis (years)"),
                continuous(df, "kps", "Karnofsky performance score")]
        for col, label in [("sex", "sex"), ("idh_status", "IDH status"), ("mgmt_promoter", "MGMT promoter"),
                           ("tmz", "temozolomide"), ("radiation", "radiotherapy"), ("gcimp", "G-CIMP"),
                           ("expression_subtype_brennan2013", "expression subtype (Brennan 2013)"),
                           ("vital_status", "vital status")]:
            part += categorical(df, col, label)
        for r in part:
            r["cohort"] = cohort
        rows += part
        surv += survival(df, cohort)
    out = pd.DataFrame(rows)[["cohort", "variable", "level", "BI", "GDC2024", "test", "p"]]
    out.to_csv(RES / "tcga_batch_comparison.csv", index=False)
    s = pd.DataFrame(surv)
    s.to_csv(RES / "tcga_batch_os.csv", index=False)
    pd.set_option("display.width", 200)
    print(out[out.cohort.str.startswith("primary")].to_string(index=False))
    print(s.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
