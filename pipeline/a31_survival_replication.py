#!/usr/bin/env python3
"""
a31_survival_replication.py
Survival replication (PREREGISTRATION.md, Section 8). In each usable cohort: Cox model (lifelines, penalizer 0,
Efron ties) of GPR153 per SD (standardized within the cohort's usable patients), (a) alone and (b) with age and
sex where available (GSE13041: characteristics "gender", M/F). Usable: survival time > 0, mapped vital status,
expression present; usable n must equal the survival-data audit (r07).

  CGGA mRNAseq_693, mRNAseq_325   log2(FPKM + 1), gene GPR153
  GSE13041 GPL96 primary (162)    as in a30 (highest-mean GPR153 probe set)
  Rembrandt (GSE108474)           dChip values; log2 if the matrix maximum exceeds 100, else as deposited;
                                  highest-mean GPR153 probe set; age = midpoint of AGE_RANGE
  CPTAC-3 GBM                     log2(TPM + 1); primary analysis IDH-wildtype, sensitivity all

Random-effects meta-analysis (DerSimonian-Laird on log HR) with tau^2 and I^2, with and without TCGA (Model 1 for
(a), Model 3 for (b), from a20); Hartung-Knapp sensitivity, Q-profile CI for tau^2, prediction interval.
Incremental value in CGGA (each release): baseline age + sex + MGMT promoter methylation, LR test, change in C,
optimism-corrected C (as a22).

Inputs
  revision/work/ext_*.npz and _clinical.csv (a29), revision/results/survival/cox_models.csv (a20),
  revision/results/cohort/survival_audit.csv (r07)
Outputs (revision/results/external/)
  survival_replication_cox.csv, survival_replication_meta.csv, survival_replication_incremental.csv

Run from the repository root:
    python pipeline/a31_survival_replication.py
"""

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

import common
from a22_cox_genomewide import incremental

AUDIT = {"CGGA_693": "CGGA mRNAseq_693", "CGGA_325": "CGGA mRNAseq_325", "GSE13041": "GSE13041 GPL96",
         "Rembrandt": "Rembrandt (GSE108474)", "CPTAC": "CPTAC-3 GBM"}


def cohort_table(name):
    z = common.load_matrix(f"ext_{name}.npz")
    X, feat, samp = z["X"], list(z["feature"]), list(z["sample"])
    c = pd.read_csv(common.work_path(f"ext_{name}_clinical.csv"))
    if "included" in c:
        c = c[c["included"].astype(bool)]
    if name.startswith("CGGA"):
        T, f = np.log2(X + 1), "GPR153"
    elif name == "CPTAC":
        T, f = np.log2(X + 1), common.GPR153
    else:
        T = np.log2(np.where(X > 0, X, np.nan)) if np.nanmax(X) > 100 else X
        idx = [samp.index(s) for s in c["sample_id"]]
        means = {p: np.nanmean(T[feat.index(p), idx]) for p in common.GPR153_PROBES if p in feat}
        f = max(means, key=means.get)
    col = {s: i for i, s in enumerate(samp)}
    c = c[c["sample_id"].isin(col)].copy()
    c["gpr153_raw"] = T[feat.index(f), [col[s] for s in c["sample_id"]]]
    c["male"] = c["sex"].map(lambda s: np.nan if pd.isna(s) else float(str(s).strip().lower() == "male"))
    c = c.rename(columns={"os_days": "time", "os_event": "event"})
    ok = c["time"].notna() & (c["time"] > 0) & c["event"].notna() & c["gpr153_raw"].notna()
    c = c[ok].copy()
    c["gpr153"] = common.standardize(c["gpr153_raw"])
    c["feature"] = f
    return c


def cox(d, covs):
    x = d[["time", "event"] + covs].dropna()
    s = CoxPHFitter(penalizer=0.0).fit(x, "time", "event", fit_options=common.COX_FIT_OPTIONS).summary.loc["gpr153"]
    return {"n": len(x), "events": int(x["event"].sum()), "log_hr": s["coef"], "se": s["se(coef)"],
            "hr": s["exp(coef)"], "hr_low": s["exp(coef) lower 95%"], "hr_high": s["exp(coef) upper 95%"], "p": s["p"]}


def main():
    common.require_registration()
    out = common.outdir("external")
    audit = pd.read_csv(common.RES / "cohort/survival_audit.csv").set_index("cohort")
    rows, tables = [], {}
    for name in AUDIT:
        d = cohort_table(name)
        if name != "CPTAC":
            assert len(d) == int(audit.loc[AUDIT[name], "n_usable"]), (name, len(d), audit.loc[AUDIT[name], "n_usable"])
        subsets = [("all", d)] if name != "CPTAC" else [("IDH-wildtype", d[d["idh_status"] == "wildtype"]), ("all", d)]
        for sub, g in subsets:
            g = g.copy()
            g["gpr153"] = common.standardize(g["gpr153_raw"])
            adj = [v for v in ("age", "male") if v in g and g[v].notna().any()]
            rows.append({"cohort": AUDIT[name], "subset": sub, "model": "(a) GPR153", "covariates": "",
                         "feature": g["feature"].iloc[0], **cox(g, ["gpr153"])})
            rows.append({"cohort": AUDIT[name], "subset": sub, "model": "(b) GPR153 + age + sex", "covariates": "+".join(adj),
                         "feature": g["feature"].iloc[0], **cox(g, ["gpr153"] + adj)})
            tables[(name, sub)] = g
    cx = pd.DataFrame(rows)
    cx.to_csv(out / "survival_replication_cox.csv", index=False)

    tcga = pd.read_csv(common.RES / "survival/cox_models.csv")
    tcga = tcga[(tcga.cohort == "primary") & (tcga.variant == "main") & (tcga.sample_set == "complete cases")
                & (tcga.covariate == "gpr153")].set_index("model")
    meta = []
    primary = cx[~((cx.cohort == "CPTAC-3 GBM") & (cx.subset == "all"))]
    for model, tm in (("(a) GPR153", "Model 1"), ("(b) GPR153 + age + sex", "Model 3")):
        e = primary[primary.model == model]
        meta.append({"model": model, "with_tcga": False,
                     **common.meta_random(e["log_hr"].values, e["se"].values ** 2, list(e["cohort"]))})
        meta.append({"model": model, "with_tcga": True,
                     **common.meta_random(np.r_[e["log_hr"].values, tcga.loc[tm, "coef"]],
                                          np.r_[e["se"].values ** 2, tcga.loc[tm, "se"] ** 2], list(e["cohort"]) + ["TCGA-GBM"])})
    md = pd.DataFrame(meta)
    for c in ("mu", "ci_low", "ci_high", "hk_ci_low", "hk_ci_high", "pi_low", "pi_high"):
        if c in md:
            md[f"{c}_hr"] = np.exp(md[c])
    md["note"] = "tau^2 and I^2 are imprecise with this few cohorts"
    md.to_csv(out / "survival_replication_meta.csv", index=False)

    inc = []
    for name in ("CGGA_693", "CGGA_325"):
        g = tables[(name, "all")].copy()
        g["meth"] = g["mgmt_promoter"].map({"methylated": 1.0, "un-methylated": 0.0})
        inc.append(incremental(g, ["age", "male", "meth"], f"{AUDIT[name]}: GPR153 + baseline vs baseline"))
    pd.DataFrame(inc).to_csv(out / "survival_replication_incremental.csv", index=False)
    print(cx[["cohort", "subset", "model", "n", "events", "hr", "p"]].to_string(index=False))


if __name__ == "__main__":
    main()
