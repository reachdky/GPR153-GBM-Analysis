#!/usr/bin/env python3
"""
a21_km.py
Kaplan-Meier and log-rank tests, TCGA primary cohort (PREREGISTRATION.md, Sections 6.4 and 6.5).
Median split of log2(TPM + 1): "high" = above the median, "low" = at or below. Two versions: pooled median
(main text, Fig 3) and within-batch median (each patient vs the median of their own batch; S3 Fig).
Per group: n, deaths, censored, median OS with 95% CI (lifelines), numbers at risk at 0, 6, 12, 18, 24 and 36
months (1 month = 30.4375 days), median follow-up by reverse Kaplan-Meier; log-rank test (lifelines
logrank_test with durations and event indicators). Curves (survival and 95% CI at every time) are written for the
figure scripts.

Inputs
  revision/work/survival_table_primary.csv (a20)
Outputs (revision/results/survival/)
  km_groups.csv, km_logrank.csv, km_at_risk.csv, km_curves.csv

Run from the repository root:
    python pipeline/a21_km.py
"""

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import median_survival_times

import common


def reverse_km_median(time, event):
    k = KaplanMeierFitter().fit(time, 1 - event)
    return k.median_survival_time_


def main():
    common.require_registration()
    out = common.outdir("survival")
    d = pd.read_csv(common.work_path("survival_table_primary.csv"))
    groups, tests, risk, curves = [], [], [], []
    for gene in ("gpr153", "mgmt"):
        x = d[f"{gene}_log2tpm"]
        splits = {"pooled median": x > x.median(),
                  "within-batch median": x > d.groupby("batch")[f"{gene}_log2tpm"].transform("median")}
        for split, high in splits.items():
            for lab, g in (("high", d[high]), ("low", d[~high])):
                k = KaplanMeierFitter().fit(g["time"], g["event"], label=lab)
                ci = median_survival_times(k.confidence_interval_)
                med = k.median_survival_time_
                groups.append({"gene": gene.upper(), "split": split, "group": lab, "n": len(g),
                               "deaths": int(g["event"].sum()), "censored": int((g["event"] == 0).sum()),
                               "median_os_days": med, "median_os_ci_low": ci.iloc[0, 0], "median_os_ci_high": ci.iloc[0, 1],
                               "median_os_months": med / common.DAYS_PER_MONTH,
                               "median_followup_days_reverse_km": reverse_km_median(g["time"], g["event"])})
                for mo in common.PARAMS["km_at_risk_months"]:
                    risk.append({"gene": gene.upper(), "split": split, "group": lab, "month": mo,
                                 "at_risk": int((g["time"] >= mo * common.DAYS_PER_MONTH).sum())})
                sf = k.survival_function_.join(k.confidence_interval_)
                sf.columns = ["survival", "ci_low", "ci_high"]
                sf = sf.reset_index().rename(columns={"timeline": "time_days"})
                sf.insert(0, "group", lab)
                sf.insert(0, "split", split)
                sf.insert(0, "gene", gene.upper())
                curves.append(sf)
            h, l_ = d[high], d[~high]
            r = logrank_test(h["time"], l_["time"], h["event"], l_["event"])
            tests.append({"gene": gene.upper(), "split": split, "n": len(d), "chisq": r.test_statistic,
                          "p": r.p_value, "median_followup_days_reverse_km": reverse_km_median(d["time"], d["event"])})
    pd.DataFrame(groups).to_csv(out / "km_groups.csv", index=False)
    pd.DataFrame(tests).to_csv(out / "km_logrank.csv", index=False)
    pd.DataFrame(risk).to_csv(out / "km_at_risk.csv", index=False)
    pd.concat(curves).to_csv(out / "km_curves.csv", index=False)
    print(pd.DataFrame(tests).to_string(index=False))


if __name__ == "__main__":
    main()
