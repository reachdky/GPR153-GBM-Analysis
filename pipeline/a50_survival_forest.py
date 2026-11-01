#!/usr/bin/env python3
"""
a50_survival_forest.py
Random survival forest (PREREGISTRATION.md, Section 11; Tier 3, supplement only). Runs in .venv_rsf
(scikit-survival 0.28.0 requires numpy >= 2; see revision/environment/).

RandomSurvivalForest(n_estimators = 500, min_samples_leaf = 15, max_features = "sqrt", random_state = 42).
Feature sets (labeled A-C, distinct from Cox Models 1-7): (A) age + sex; (B) A + MGMT, EGFR, PTEN, IDH1, TERT
expression (log2(TPM + 1)); (C) B + GPR153. Analyses: primary cohort; primary-cohort complete cases with KPS and
MGMT methylation added to every set; sensitivity cohort with IDH-mutant and IDH-unknown indicators added.
Validation: 10 repeats of 5-fold CV stratified by event (identical folds for every set); Harrell's C on held-out
folds. Paired difference C(C) - C(B) with the Nadeau-Bengio corrected interval
(variance x (1/J + n_test/n_train), t with J - 1 df, J = 50 folds). Permutation importance (10 repeats) of set C on
held-out folds.

Inputs
  revision/work/survival_table_{primary,sensitivity}.csv (a20)
Outputs (revision/results/forest/)
  forest_cv_folds.csv, forest_c_summary.csv, forest_paired_difference.csv, forest_permutation_importance.csv

Run from the repository root:
    .venv_rsf/bin/python pipeline/a50_survival_forest.py
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.inspection import permutation_importance
from sklearn.model_selection import RepeatedStratifiedKFold
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv

import common

GENES = ["mgmt_log2tpm", "egfr_log2tpm", "pten_log2tpm", "idh1_log2tpm", "tert_log2tpm"]


def sets(extra):
    a = ["age", "male"] + extra
    b = a + GENES
    return {"A": a, "B": b, "C": b + ["gpr153_log2tpm"]}


def run(d, fs, label, F):
    cols = sorted(set().union(*fs.values()))
    d = d[["time", "event"] + cols].dropna().reset_index(drop=True)
    y = Surv.from_arrays(d["event"].astype(bool), d["time"])
    cv = RepeatedStratifiedKFold(n_splits=F["cv_folds"], n_repeats=F["cv_repeats"], random_state=common.SEED)
    folds, imp = [], []
    for j, (tr, te) in enumerate(cv.split(d, d["event"])):
        for name, f in fs.items():
            m = RandomSurvivalForest(n_estimators=F["n_estimators"], min_samples_leaf=F["min_samples_leaf"],
                                     max_features=F["max_features"], random_state=F["random_state"], n_jobs=-1)
            m.fit(d.loc[tr, f], y[tr])
            c = concordance_index_censored(y[te]["event"], y[te]["time"], m.predict(d.loc[te, f]))[0]
            folds.append({"analysis": label, "fold": j + 1, "set": name, "c": c, "n_train": len(tr), "n_test": len(te)})
            if name == "C":
                pi = permutation_importance(m, d.loc[te, f], y[te], n_repeats=F["permutation_importance_repeats"],
                                            random_state=common.SEED)
                imp += [{"analysis": label, "fold": j + 1, "feature": f[k], "importance": pi.importances_mean[k]}
                        for k in range(len(f))]
    return pd.DataFrame(folds), pd.DataFrame(imp), len(d), int(d["event"].sum())


def main():
    common.set_seeds()
    common.require_registration()
    out = common.outdir("forest")
    F = common.PARAMS["forest"]
    prim = pd.read_csv(common.work_path("survival_table_primary.csv"))
    sens = pd.read_csv(common.work_path("survival_table_sensitivity.csv"))
    analyses = [("primary cohort", prim, sets([])), ("primary cohort, complete cases", prim, sets(["kps", "meth"])),
                ("sensitivity cohort", sens, sets(["idh_mutant", "idh_unknown"]))]
    all_f, all_i, summ, diff = [], [], [], []
    for label, d, fs in analyses:
        f, i, n, ev = run(d, fs, label, F)
        all_f.append(f)
        all_i.append(i)
        for s, g in f.groupby("set"):
            summ.append({"analysis": label, "set": s, "features": "+".join(fs[s]), "n": n, "events": ev,
                         "c_mean": g.c.mean(), "c_sd": g.c.std(ddof=1), "c_p025": g.c.quantile(0.025),
                         "c_p975": g.c.quantile(0.975)})
        w = f.pivot(index="fold", columns="set", values="c")
        dd = (w["C"] - w["B"]).values
        J = len(dd)
        nt, ntr = f["n_test"].mean(), f["n_train"].mean()
        se = np.sqrt((1 / J + nt / ntr) * dd.var(ddof=1))
        tq = stats.t.ppf(0.975, J - 1)
        diff.append({"analysis": label, "folds": J, "mean_diff_C_minus_B": dd.mean(), "nb_se": se,
                     "nb_ci_low": dd.mean() - tq * se, "nb_ci_high": dd.mean() + tq * se,
                     "nb_p": 2 * stats.t.sf(abs(dd.mean() / se), J - 1)})
    pd.concat(all_f).to_csv(out / "forest_cv_folds.csv", index=False)
    pd.DataFrame(summ).to_csv(out / "forest_c_summary.csv", index=False)
    pd.DataFrame(diff).to_csv(out / "forest_paired_difference.csv", index=False)
    pi = pd.concat(all_i)
    pi.groupby(["analysis", "feature"], as_index=False)["importance"].agg(["mean", "std"]).reset_index().to_csv(
        out / "forest_permutation_importance.csv", index=False)
    print(pd.DataFrame(diff).to_string(index=False))


if __name__ == "__main__":
    main()
