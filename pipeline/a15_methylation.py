#!/usr/bin/env python3
"""
a15_methylation.py
MGMT promoter methylation (PREREGISTRATION.md, Sections 3.2 and 5.3; decision condition 4).

MGMT-STP27 (Bady 2012): M = log2(beta / (1 - beta)), beta clipped to [0.001, 0.999];
logit(y) = 4.3215 + 0.5271 M(cg12434587) + 0.9265 M(cg12981137); methylated if y >= 0.358.
STP27 population: primary-cohort patients with both probes non-missing.
Kappa: Cohen's kappa of STP27 status vs the Ceccarelli 2016 call, in patients with both.

Condition 4: partial Spearman of GPR153 with MGMT expression adjusted for batch, array platform (450K vs 27K)
and the methylation measure (k = 3; continuous measures enter as ranks), Fisher z with the Bonett-Wright SE.
Primary measure: STP27 linear predictor in the STP27 population if n >= 100 and kappa >= 0.60; otherwise the
Ceccarelli binary call in primary-cohort patients with a call; if that population has n < 100, condition 4 is
"not evaluable" (not met). Met if rho > 0 and the 95% CI excludes 0. A covariate that is constant in the
population is dropped and k is reduced; this is reported.
Sensitivity (do not change the decision): bootstrap percentile 95% CI (1,000 resamples, seed 42, re-ranked and
re-residualized in each); without platform (k = 2); binary STP27 status; the other measure.

Descriptive: binary status in fixed order (STP27 status, else Ceccarelli 2016, else Brennan 2013); GPR153-MGMT
within methylated and within unmethylated patients (batch-adjusted); GPR153 by status (Mann-Whitney U);
batch-adjusted Spearman of GPR153 with the STP27 linear predictor; mean beta of a single available STP27 CpG.

Inputs
  revision/work/tcga_expr.npz (a10), revision/work/methylation_inputs.csv (a19),
  revision/manifests/tcga_samples.csv (r01)
Outputs (revision/results/methylation/)
  stp27_calls.csv, methylation_concordance.csv, methylation_kappa_2x2.csv (kappa population: primary-cohort
  patients with both an STP27 call and a Ceccarelli call), condition4.csv, condition4_sensitivity.csv,
  methylation_descriptive.csv
  revision/work/methylation_measure.csv  (per patient of the primary, IDH-not-mutant and sensitivity cohorts:
  measure used for condition 4 and Cox Models 4, 5, 7)

Run from the repository root:
    python pipeline/a15_methylation.py
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import cohen_kappa_score

import common
from common import GPR153, MGMT

P = common.PARAMS["stp27"]
C4 = common.PARAMS["condition4"]


def stp27(b1, b2):
    lo, hi = P["beta_clip"]
    m1 = np.log2(np.clip(b1, lo, hi) / (1 - np.clip(b1, lo, hi)))
    m2 = np.log2(np.clip(b2, lo, hi) / (1 - np.clip(b2, lo, hi)))
    lp = P["intercept"] + P["b_cg12434587"] * m1 + P["b_cg12981137"] * m2
    return lp, (1 / (1 + np.exp(-lp)) >= P["cutoff"]).astype(int)


def calls(rows, mi):
    """STP27 linear predictor and status, Ceccarelli and Brennan calls, platform, per patient of `rows`."""
    d = rows[["patient", "aliquot", "batch", "mgmt_promoter", "mgmt_promoter_source"]].merge(
        mi.drop(columns="aliquot"), on="patient", how="left")
    d["batch_gdc"] = (d["batch"] == "GDC2024").astype(float)
    d["hm450"] = (d["platform_used"] == "HM450").astype(float)
    both = d["beta_cg12434587"].notna() & d["beta_cg12981137"].notna()
    d["stp27_lp"], d["stp27_status"] = np.nan, np.nan
    lp, st = stp27(d.loc[both, "beta_cg12434587"].values, d.loc[both, "beta_cg12981137"].values)
    d.loc[both, "stp27_lp"], d.loc[both, "stp27_status"] = lp, st
    d["ceccarelli"] = np.where(d["mgmt_promoter_source"] == "Ceccarelli 2016",
                               d["mgmt_promoter"].map({"methylated": 1, "unmethylated": 0}), np.nan)
    d["brennan"] = np.where(d["mgmt_promoter_source"] == "Brennan 2013",
                            d["mgmt_promoter"].map({"methylated": 1, "unmethylated": 0}), np.nan)
    d["status_desc"] = d["stp27_status"].fillna(d["ceccarelli"]).fillna(d["brennan"])
    return d, both


def partial(d, measure, use_platform=True):
    covs = [d["batch_gdc"].values]
    names = ["batch"]
    if use_platform:
        covs.append(d["hm450"].values)
        names.append("platform")
    covs.append(d[measure].values)
    names.append(measure)
    cov = np.column_stack(covs)
    keep = [i for i in range(cov.shape[1]) if np.unique(cov[:, i]).size > 1]
    dropped = [names[i] for i in range(len(names)) if i not in keep]
    rho, p, df = common.partial_spearman(d["gpr153"].values, d["mgmt"].values, cov[:, keep])
    k = len(keep)
    lo, hi, _, _ = common.bw_ci(rho[0], len(d), k)
    return {"n": len(d), "k": k, "df": df, "rho": rho[0], "p": p[0], "ci_low": lo, "ci_high": hi,
            "covariates": "+".join(names[i] for i in keep), "dropped_constant": "+".join(dropped)}, cov[:, keep]


def bootstrap_ci(d, cov):
    rng = np.random.default_rng(common.SEED)
    vals = []
    for _ in range(C4["n_bootstrap"]):
        i = rng.choice(len(d), len(d), replace=True)
        c = cov[i]
        c = c[:, [j for j in range(c.shape[1]) if np.unique(c[:, j]).size > 1]]
        vals.append(common.partial_spearman(d["gpr153"].values[i], d["mgmt"].values[i], c)[0][0])
    return np.percentile(vals, [2.5, 97.5])


def main():
    common.set_seeds()
    common.require_registration()
    out = common.outdir("methylation")
    t = common.TCGA()
    prim = common.cohort_aliquots(t.m, "in_screen_primary")
    mi = pd.read_csv(common.work_path("methylation_inputs.csv"))
    d, both = calls(prim, mi)
    d = d.drop(columns="status_desc")
    d["gpr153"] = t.gene(GPR153, d["aliquot"].values)
    d["mgmt"] = t.gene(MGMT, d["aliquot"].values)
    one = d["beta_cg12434587"].notna() ^ d["beta_cg12981137"].notna()
    d["single_cpg_beta"] = np.where(one, d[["beta_cg12434587", "beta_cg12981137"]].mean(1), np.nan)
    d.drop(columns=["gpr153", "mgmt"]).to_csv(out / "stp27_calls.csv", index=False)

    ov = d["stp27_status"].notna() & d["ceccarelli"].notna()
    kappa = cohen_kappa_score(d.loc[ov, "stp27_status"], d.loc[ov, "ceccarelli"]) if ov.sum() > 1 else np.nan
    conc = pd.crosstab(d.loc[ov, "stp27_status"], d.loc[ov, "ceccarelli"])
    lev = {1.0: "methylated", 0.0: "unmethylated"}
    pd.DataFrame([{"stp27": lev[a], "ceccarelli": lev[b], "n": int((d.loc[ov, "stp27_status"].eq(a)
                                                                    & d.loc[ov, "ceccarelli"].eq(b)).sum()),
                   "population_n": int(ov.sum()), "kappa": kappa} for a in lev for b in lev]).to_csv(
        out / "methylation_kappa_2x2.csv", index=False)
    n_stp = int(both.sum())
    use_stp = n_stp >= C4["min_n"] and not np.isnan(kappa) and kappa >= C4["min_kappa"]
    pd.DataFrame([{"stp27_population_n": n_stp, "overlap_with_ceccarelli": int(ov.sum()), "kappa": kappa,
                   "table": conc.to_dict(), "single_cpg_patients": int(one.sum()),
                   "single_cpg_mean_beta": float(d["single_cpg_beta"].mean()) if one.any() else np.nan,
                   "primary_measure": "STP27 linear predictor" if use_stp else "Ceccarelli 2016 call",
                   "reason": (f"STP27 population n = {n_stp} (>= {C4['min_n']}) and kappa = {kappa:.3f} (>= {C4['min_kappa']})"
                              if use_stp else f"STP27 population n = {n_stp}, kappa = {kappa}: fallback rule")}]).to_csv(
        out / "methylation_concordance.csv", index=False)

    pop_stp = d[both]
    pop_cec = d[d["ceccarelli"].notna() & d["platform_used"].notna()]
    pop, measure = (pop_stp, "stp27_lp") if use_stp else (pop_cec, "ceccarelli")
    if len(pop) < C4["min_n"]:
        res = {"measure": measure, "n": len(pop), "status": "not evaluable", "condition4_pass": False}
        cov = None
    else:
        res, cov = partial(pop, measure)
        res.update(measure=measure, status="evaluated",
                   condition4_pass=bool(res["rho"] > 0 and res["ci_low"] > 0))
    pd.DataFrame([res]).to_csv(out / "condition4.csv", index=False)

    sens = []
    if cov is not None:
        lo, hi = bootstrap_ci(pop, cov)
        sens.append({"analysis": "bootstrap percentile 95% CI", "measure": measure, "n": len(pop),
                     "rho": res["rho"], "ci_low": lo, "ci_high": hi})
        r, _ = partial(pop, measure, use_platform=False)
        sens.append({"analysis": "without platform (k = 2)", "measure": measure, **r})
    if n_stp >= 3:
        r, _ = partial(pop_stp, "stp27_status")
        sens.append({"analysis": "binary STP27 status", "measure": "stp27_status", **r})
    other = (pop_cec, "ceccarelli") if use_stp else (pop_stp, "stp27_lp")
    if len(other[0]) >= 5:
        r, _ = partial(*other)
        sens.append({"analysis": "other measure", "measure": other[1], **r})
    pd.DataFrame(sens).to_csv(out / "condition4_sensitivity.csv", index=False)

    d["status_desc"] = d["stp27_status"].fillna(d["ceccarelli"]).fillna(d["brennan"])
    d["status_source"] = np.select([d["stp27_status"].notna(), d["ceccarelli"].notna(), d["brennan"].notna()],
                                   ["STP27", "Ceccarelli 2016", "Brennan 2013"], "none")
    desc = [{"analysis": f"patients by status source: {s}", "n": int(n)} for s, n in d["status_source"].value_counts().items()]
    for lab, v in (("methylated", 1), ("unmethylated", 0)):
        g = d[d["status_desc"] == v]
        if len(g) >= 5:
            rho, p, df = common.partial_spearman(g["gpr153"].values, g["mgmt"].values, g[["batch_gdc"]].values
                                                 if g["batch_gdc"].nunique() > 1 else None)
            desc.append({"analysis": f"GPR153-MGMT within {lab}", "n": len(g), "rho": rho[0], "p": p[0]})
    g = d[d["status_desc"].notna()]
    if g["status_desc"].nunique() == 2:
        u = stats.mannwhitneyu(g.loc[g.status_desc == 1, "gpr153"], g.loc[g.status_desc == 0, "gpr153"])
        desc.append({"analysis": "GPR153 by binary status (Mann-Whitney U)", "n": len(g), "statistic": u.statistic,
                     "p": u.pvalue, "median_methylated": g.loc[g.status_desc == 1, "gpr153"].median(),
                     "median_unmethylated": g.loc[g.status_desc == 0, "gpr153"].median()})
    if n_stp >= 5:
        rho, p, _ = common.partial_spearman(pop_stp["gpr153"].values, pop_stp["stp27_lp"].values, pop_stp[["batch_gdc"]].values)
        desc.append({"analysis": "GPR153 vs STP27 linear predictor (batch-adjusted)", "n": n_stp, "rho": rho[0], "p": p[0]})
    pd.DataFrame(desc).to_csv(out / "methylation_descriptive.csv", index=False)
    # Measure for Cox Models 4, 5 and 7 in every analysis cohort (Section 6.2): selected in the primary cohort
    # (above), computed for all patients of the primary, IDH-not-mutant and sensitivity cohorts
    sel = t.m["selected_aliquot"] & (t.m["in_screen_primary"] | t.m["in_screen_idh_not_mutant"]
                                     | t.m["in_screen_sensitivity"])
    a, _ = calls(t.m.loc[sel].reset_index(drop=True), mi)
    a[["patient", "stp27_lp", "stp27_status", "ceccarelli", "status_desc", "hm450"]].assign(
        measure_name=measure, measure=a[measure]).to_csv(common.work_path("methylation_measure.csv"), index=False)
    print(pd.DataFrame([res]).to_string(index=False))


if __name__ == "__main__":
    main()
