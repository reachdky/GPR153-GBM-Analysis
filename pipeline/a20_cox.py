#!/usr/bin/env python3
"""
a20_cox.py
Cox models, TCGA (PREREGISTRATION.md, Sections 6.1-6.3). All unpenalized (lifelines CoxPHFitter,
penalizer = 0, Efron ties). GPR153 and MGMT expression per SD of log2(TPM + 1), standardized within the
survival-evaluable patients of each analysis cohort; age per year; sex as male = 1; KPS as recorded; MGMT
methylation = the measure selected for condition 4 (a15).

  Model 1  GPR153                                  Model 5  Model 4 + MGMT expression
  Model 2  GPR153 + age                            Model 6  GPR153 + age + sex + MGMT expression
  Model 3  GPR153 + age + sex  (primary)           Model 7  GPR153 + age + sex + MGMT methylation
  Model 4  GPR153 + age + sex + KPS + MGMT methylation

Each model is fitted on its complete cases; Models 1-3 are also refitted on the Model 4 complete-case set.
Likelihood-ratio tests: nested pairs on the larger model's sample, and GPR153 removed from each model.
Sensitivity: stratified by batch; batch as covariate; Models 4 and 5 + temozolomide + radiotherapy (immortal-
time bias; not for inference); IDH-not-mutant and sensitivity cohorts (sensitivity cohort with IDH-mutant and
IDH-unknown indicators). Proportional hazards: R survival::cox.zph (transform "km", per term and global) for every
fit, because T5 found lifelines and cox.zph differ (tests/test_schoenfeld.py). Age linearity: Model 3 with a
restricted cubic spline for age (3 knots at the 10th, 50th, 90th percentiles) vs linear, LR test.

Inputs
  revision/work/tcga_expr.npz (a10), revision/work/methylation_measure.csv (a15),
  revision/manifests/tcga_samples.csv (r01)
Outputs (revision/results/survival/)
  cox_models.csv, cox_lr_tests.csv, cox_ph_tests.csv, cox_age_spline.csv, cox_model3_schoenfeld.csv
  revision/work/survival_table_<cohort>.csv, revision/work/mi_input.csv

Run from the repository root:
    python pipeline/a20_cox.py
"""

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, NelsonAalenFitter
from scipy import stats

import common
from common import GPR153, MGMT

MODELS = {
    "Model 1": ["gpr153"],
    "Model 2": ["gpr153", "age"],
    "Model 3": ["gpr153", "age", "male"],
    "Model 4": ["gpr153", "age", "male", "kps", "meth"],
    "Model 5": ["gpr153", "age", "male", "kps", "meth", "mgmt"],
    "Model 6": ["gpr153", "age", "male", "mgmt"],
    "Model 7": ["gpr153", "age", "male", "meth"],
}
NESTED = [("Model 1", "Model 2"), ("Model 2", "Model 3"), ("Model 3", "Model 4"), ("Model 4", "Model 5"),
          ("Model 3", "Model 6"), ("Model 6", "Model 5"), ("Model 3", "Model 7"), ("Model 7", "Model 4")]
COHORTS = {"primary": "in_screen_primary", "idh_not_mutant": "in_screen_idh_not_mutant",
           "sensitivity": "in_screen_sensitivity"}


def survival_data(t, flag):
    m = t.m
    s = common.survival_table(m, flag)
    cov = m.set_index("aliquot").loc[s["aliquot"]]
    s["age"] = cov["age_at_diagnosis_years"].values
    s["male"] = (cov["sex"].values == "male").astype(float)
    s.loc[cov["sex"].isna().values, "male"] = np.nan
    s["kps"] = cov["kps"].values
    s["tmz"] = cov["tmz"].map({"yes": 1.0, "no": 0.0}).values
    s["radiation"] = cov["radiation"].map({"yes": 1.0, "no": 0.0}).values
    s["idh_mutant"] = (cov["idh_status"].values == "mutant").astype(float)
    s["idh_unknown"] = (cov["idh_status"].values == "unknown").astype(float)
    s["batch_gdc"] = (s["batch"] == "GDC2024").astype(float)
    s["gpr153_log2tpm"] = t.gene(GPR153, s["aliquot"].values)
    s["mgmt_log2tpm"] = t.gene(MGMT, s["aliquot"].values)
    for sym in common.FOREST_GENES[1:]:
        rows = np.where(t.symbol == sym)[0]
        if len(rows):
            r = rows[np.argmax(t.log2tpm(s["aliquot"].values)[rows].mean(1))]
            s[f"{sym.lower()}_log2tpm"] = np.log2(t.tpm[r, t.cols(s["aliquot"].values)] + 1)
    s["gpr153"] = common.standardize(s["gpr153_log2tpm"])
    s["mgmt"] = common.standardize(s["mgmt_log2tpm"])
    meth = pd.read_csv(common.work_path("methylation_measure.csv")).set_index("patient")
    s["meth"] = meth["measure"].reindex(s["patient"]).values
    s["meth_status"] = meth["status_desc"].reindex(s["patient"]).values
    s["meth_measure_name"] = meth["measure_name"].iloc[0]
    return s.rename(columns={"os_days": "time", "os_event": "event"})


def fit(d, covs, strata=None):
    """Complete cases on the model columns; a covariate constant in that sample is dropped (listed in
    cox_models.csv, column dropped_constant)."""
    cols = ["time", "event"] + covs + ([strata] if strata else [])
    x = d[cols].dropna()
    const = [c for c in covs if x[c].nunique() < 2]
    x = x.drop(columns=const)
    c = CoxPHFitter(penalizer=0.0).fit(x, "time", "event", strata=[strata] if strata else None,
                                   fit_options=common.COX_FIT_OPTIONS)
    c.dropped_constant = ";".join(const)
    return c, x


def rows_from(c, x, cohort, model, variant, sample_set):
    s = c.summary
    return [{"cohort": cohort, "model": model, "variant": variant, "sample_set": sample_set, "covariate": k,
             "coef": s.loc[k, "coef"], "se": s.loc[k, "se(coef)"], "hr": s.loc[k, "exp(coef)"],
             "hr_low": s.loc[k, "exp(coef) lower 95%"], "hr_high": s.loc[k, "exp(coef) upper 95%"],
             "p": s.loc[k, "p"], "n": len(x), "events": int(x["event"].sum()), "loglik": c.log_likelihood_,
             "dropped_constant": c.dropped_constant}
            for k in s.index]


def lr(big, small, df):
    stat = 2 * (big.log_likelihood_ - small.log_likelihood_)
    return stat, stats.chi2.sf(stat, df)


def main():
    common.set_seeds()
    common.require_registration()
    out = common.outdir("survival")
    t = common.TCGA()
    res, lrs, jobs = [], [], []
    tmp = Path(tempfile.mkdtemp(dir=common.work_path("")))
    for cohort, flag in COHORTS.items():
        d = survival_data(t, flag)
        d.to_csv(common.work_path(f"survival_table_{cohort}.csv"), index=False)
        extra = ["idh_mutant", "idh_unknown"] if cohort == "sensitivity" else []
        fits = {}
        for name, covs in MODELS.items():
            covs = covs + [e for e in extra if d[e].nunique() > 1]
            variants = [("main", covs, None)]
            variants += [("stratified by batch", covs, "batch_gdc"), ("batch covariate", covs + ["batch_gdc"], None)]
            if name in ("Model 4", "Model 5"):
                variants.append(("+ temozolomide + radiotherapy (immortal-time bias)", covs + ["tmz", "radiation"], None))
            for variant, cv, st in variants:
                c, x = fit(d, cv, st)
                res += rows_from(c, x, cohort, name, variant, "complete cases")
                cv = [v for v in cv if v in x.columns]
                if variant == "main":
                    fits[name] = (c, x, cv)
                tag = f"{cohort}_{name}_{variant}".replace(" ", "_").replace("+", "plus").replace("(", "").replace(")", "")
                x.to_csv(tmp / f"{tag}.csv", index=False)
                jobs.append({"input": tmp / f"{tag}.csv", "prefix": tmp / tag, "covariates": ",".join(cv),
                             "strata": st or "", "cohort": cohort, "model": name, "variant": variant})
        m4x = fits["Model 4"][1]
        for name in ("Model 1", "Model 2", "Model 3"):
            c, x = fit(m4x, [v for v in fits[name][2] if v in m4x.columns])
            res += rows_from(c, x, cohort, name, "main", "Model 4 complete cases")
        for small, big in NESTED:
            cb, xb, cvb = fits[big]
            sm = [v for v in fits[small][2] if v in cvb]
            cs, _ = fit(xb, sm)
            stat, p = lr(cb, cs, len(cvb) - len(sm))
            lrs.append({"cohort": cohort, "test": f"{small} vs {big}", "n": len(xb), "lr_stat": stat,
                        "df": len(cvb) - len(sm), "p": p})
        for name, (c, x, cv) in fits.items():
            if len(cv) > 1:
                cs, _ = fit(x, [v for v in cv if v != "gpr153"])
                stat, p = lr(c, cs, 1)
            else:
                stat = c.log_likelihood_ratio_test().test_statistic
                p = c.log_likelihood_ratio_test().p_value
            lrs.append({"cohort": cohort, "test": f"GPR153 in {name}", "n": len(x), "lr_stat": stat, "df": 1, "p": p})
        if cohort == "primary":
            c3, x3, cv3 = fits["Model 3"]
            knots = np.percentile(x3["age"], [10, 50, 90])
            xs = x3.assign(age_rcs1=common.rcs_basis(x3["age"], knots)[:, 0])
            cs, _ = fit(xs, cv3 + ["age_rcs1"])
            stat, p = lr(cs, c3, 1)
            pd.DataFrame([{"model": "Model 3", "knots": ";".join(f"{k:.2f}" for k in knots), "n": len(x3),
                           "lr_stat": stat, "df": 1, "p": p}]).to_csv(out / "cox_age_spline.csv", index=False)
            r = c3.compute_residuals(x3, "scaled_schoenfeld")
            r.insert(0, "time", x3.loc[r.index, "time"].values)
            r.to_csv(out / "cox_model3_schoenfeld.csv", index=False)
            mi = d[["time", "event", "gpr153", "age", "male", "kps", "meth_status", "mgmt", "batch_gdc"]].copy()
            mi["cumhaz"] = NelsonAalenFitter().fit(mi["time"], mi["event"]).cumulative_hazard_at_times(mi["time"]).values
            mi.to_csv(common.work_path("mi_input.csv"), index=False)
    pd.DataFrame(res).to_csv(out / "cox_models.csv", index=False)
    pd.DataFrame(lrs).to_csv(out / "cox_lr_tests.csv", index=False)
    jb = pd.DataFrame(jobs)
    jb[["input", "prefix", "covariates", "strata"]].to_csv(tmp / "jobs.csv", index=False)
    common.run_r("pipeline/a20_coxzph.R", tmp / "jobs.csv")
    ph = []
    for j in jb.itertuples():
        z = pd.read_csv(f"{j.prefix}_zph.csv")
        z.insert(0, "variant", j.variant)
        z.insert(0, "model", j.model)
        z.insert(0, "cohort", j.cohort)
        ph.append(z)
        rc = pd.read_csv(f"{j.prefix}_coxph.csv")
        py = pd.DataFrame(res)
        py = py[(py.cohort == j.cohort) & (py.model == j.model) & (py.variant == j.variant)
                & (py.sample_set == "complete cases")].set_index("covariate")
        assert np.allclose(py.loc[rc.term, "coef"].values, rc["coef"].values, rtol=1e-4, atol=1e-6), j.prefix
    pd.concat(ph).assign(source="R survival::cox.zph, transform km").to_csv(out / "cox_ph_tests.csv", index=False)
    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()
    r = pd.DataFrame(res)
    print(r[(r.cohort == "primary") & (r.variant == "main") & (r.covariate == "gpr153")
            & (r.sample_set == "complete cases")][["model", "hr", "hr_low", "hr_high", "p", "n", "events"]].to_string(index=False))


if __name__ == "__main__":
    main()
