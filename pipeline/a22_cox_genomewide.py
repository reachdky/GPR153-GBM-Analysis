#!/usr/bin/env python3
"""
a22_cox_genomewide.py
Selection context and incremental value (PREREGISTRATION.md, Section 6.6).
Univariate Cox (per SD of log2(TPM + 1), Efron ties) for every eligible gene in the survival-evaluable primary
cohort, by a vectorized Newton-Raphson fit (common.cox_univariate; checked against lifelines in
tests/test_cox_vectorized.py). GPR153's percentile by p and by |log HR|.
Incremental value of GPR153 over the clinical baseline (age, sex, KPS, MGMT methylation; Model 4 complete cases):
LR test of Model 4 vs baseline, Harrell's C of both and the change, and bootstrap optimism-corrected C for both
(1,000 resamples, seed 42; Harrell's method: refit in each resample, optimism = C in the resample minus C of that
fit in the original data).

Inputs
  revision/work/tcga_expr.npz (a10), revision/work/survival_table_primary.csv (a20),
  revision/results/screen/screen_gene_level.csv (a11)
Outputs (revision/results/survival/)
  cox_genomewide.csv, cox_genomewide_gpr153.csv, incremental_value.csv

Run from the repository root:
    python pipeline/a22_cox_genomewide.py
"""

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy import stats

import common
from common import GPR153, MGMT

BASE = ["age", "male", "kps", "meth"]


def cindex(c, x):
    return concordance_index(x["time"], -c.predict_partial_hazard(x), x["event"])


def incremental(d, base, label):
    x = d[["time", "event", "gpr153"] + base].dropna()
    full = CoxPHFitter(penalizer=0.0).fit(x, "time", "event", fit_options=common.COX_FIT_OPTIONS)
    red = CoxPHFitter(penalizer=0.0).fit(x.drop(columns="gpr153"), "time", "event", fit_options=common.COX_FIT_OPTIONS)
    stat = 2 * (full.log_likelihood_ - red.log_likelihood_)
    app_f, app_r = cindex(full, x), cindex(red, x.drop(columns="gpr153"))
    rng = np.random.default_rng(common.SEED)
    opt_f, opt_r = [], []
    for _ in range(common.PARAMS["n_bootstrap_optimism"]):
        b = x.iloc[rng.choice(len(x), len(x), replace=True)].reset_index(drop=True)
        try:
            bf = CoxPHFitter(penalizer=0.0).fit(b, "time", "event", fit_options=common.COX_FIT_OPTIONS)
            br = CoxPHFitter(penalizer=0.0).fit(b.drop(columns="gpr153"), "time", "event", fit_options=common.COX_FIT_OPTIONS)
        except Exception:
            continue
        opt_f.append(cindex(bf, b) - cindex(bf, x))
        opt_r.append(cindex(br, b.drop(columns="gpr153")) - cindex(br, x.drop(columns="gpr153")))
    return {"analysis": label, "n": len(x), "events": int(x["event"].sum()), "baseline": "+".join(base),
            "lr_stat": stat, "lr_p": stats.chi2.sf(stat, 1), "c_baseline": app_r, "c_with_gpr153": app_f,
            "delta_c": app_f - app_r, "c_baseline_optimism_corrected": app_r - np.mean(opt_r),
            "c_with_gpr153_optimism_corrected": app_f - np.mean(opt_f),
            "delta_c_optimism_corrected": (app_f - np.mean(opt_f)) - (app_r - np.mean(opt_r)),
            "bootstrap_resamples_used": len(opt_f)}


def main():
    common.set_seeds()
    common.require_registration()
    out = common.outdir("survival")
    t = common.TCGA()
    d = pd.read_csv(common.work_path("survival_table_primary.csv"))
    gl = pd.read_csv(common.RES / "screen/screen_gene_level.csv")
    genes = gl.loc[gl["eligible_20"], "gene_id"].to_numpy()
    X = t.log2tpm(d["aliquot"].values)[[t.row[g] for g in genes]]
    X = (X - X.mean(1, keepdims=True)) / X.std(1, ddof=1, keepdims=True)
    b, se, z, p = common.cox_univariate(X, d["time"].values, d["event"].values)
    gw = pd.DataFrame({"gene_id": genes, "symbol": gl.set_index("gene_id").loc[genes, "symbol"].values,
                       "log_hr": b, "se": se, "hr": np.exp(b), "p": p, "n": len(d), "events": int(d["event"].sum())})
    gw.to_csv(out / "cox_genomewide.csv", index=False)
    g = gw.set_index("gene_id").loc[GPR153]
    pd.DataFrame([{"gene": "GPR153", "log_hr": g.log_hr, "hr": g.hr, "p": g.p, "genes": len(gw),
                   "percentile_by_p": float((gw["p"] < g.p).mean() * 100),
                   "percentile_by_abs_log_hr": float((gw["log_hr"].abs() > abs(g.log_hr)).mean() * 100),
                   "mgmt_p": float(gw.set_index("gene_id").loc[MGMT, "p"]) if MGMT in set(gw.gene_id) else np.nan}]).to_csv(
        out / "cox_genomewide_gpr153.csv", index=False)
    inc = incremental(d, BASE, "TCGA primary, Model 4 vs baseline")
    pd.DataFrame([inc]).to_csv(out / "incremental_value.csv", index=False)
    print(pd.read_csv(out / "cox_genomewide_gpr153.csv").to_string(index=False))
    print({k: inc[k] for k in ("n", "lr_p", "delta_c", "delta_c_optimism_corrected")})


if __name__ == "__main__":
    main()
