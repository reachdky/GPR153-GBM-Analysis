#!/usr/bin/env python3
"""
a91_independent_check.py
Independent reproduction of the numbers that may enter the manuscript (PREREGISTRATION.md, Section 13.2): the
screen, Model 3, Kaplan-Meier / log-rank and the external meta-analyses. Each value is recomputed by a separate
implementation and compared with the pipeline's output. Deterministic statistics must agree to a relative
difference <= 1e-6; bootstrap and permutation outputs (re-run with another seed) must agree within 3 Monte Carlo
standard errors. The script stops (exit 1) if any comparison fails.

  Screen, batch-adjusted   pingouin.partial_corr(method="spearman", covar batch) for GPR153, the top 100 genes and
                           500 random eligible genes (seed 42); all eligible genes by pandas ranks, numpy lstsq
                           residuals and statsmodels BH (rho, p, q, rank)
  Screen, pooled           scipy.stats.spearmanr for every non-constant gene; statsmodels BH (rho, p, q, rank)
  Permutation null         10,000 within-batch permutations (seed 43), one permutation at a time (lstsq residuals)
  Bootstrap rank           1,000 resamples (seed 43), pandas ranks and lstsq residuals
  Model 3                  R survival::coxph (Efron ties, eps 1e-12)
  KM / log-rank            R survival::survdiff and survfit (log-log CI, as lifelines); at-risk numbers;
                           reverse-KM median follow-up
  Meta-analyses            R metafor::rma (DerSimonian-Laird; test "knha"; Q-profile CI, tol 1e-12;
                           prediction interval pi.type "Riley", t with K - 2 df)

Inputs
  revision/work/tcga_expr.npz, survival_table_primary.csv; revision/results/screen/*, survival/*, external/*
Outputs (revision/results/independent/)
  independent_check.csv          every comparison: pipeline value, independent value, criterion, pass
  independent_check_summary.csv  per analysis: values compared, largest difference, pass

Run from the repository root:
    python pipeline/a91_independent_check.py
"""

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pingouin
from joblib import Parallel, delayed
from scipy import stats
from scipy.stats import gaussian_kde
from statsmodels.stats.multitest import multipletests

import common
from common import GPR153, MGMT

P = common.PARAMS["independent_check"]
ROWS = []


def rel(a, b):
    a, b = float(a), float(b)
    if (np.isnan(a) or np.isinf(a)) and (np.isnan(b) or np.isinf(b)):
        return 0.0
    if np.isnan(a) or np.isnan(b) or np.isinf(a) or np.isinf(b):
        return np.inf
    m = max(abs(a), abs(b))
    return 0.0 if m == 0 else abs(a - b) / m


def det(analysis, quantity, pipe, ind, impl):
    d = rel(pipe, ind)
    ROWS.append({"analysis": analysis, "quantity": quantity, "pipeline_value": pipe, "independent_value": ind,
                 "independent_implementation": impl, "criterion": "relative difference", "tolerance": P["rel_tol"],
                 "difference": d, "pass": d <= P["rel_tol"]})


def mc(analysis, quantity, pipe, ind, se, impl):
    d = abs(float(pipe) - float(ind))
    tol = P["mc_z"] * se
    ROWS.append({"analysis": analysis, "quantity": quantity, "pipeline_value": pipe, "independent_value": ind,
                 "independent_implementation": impl, "criterion": f"abs. difference <= {P['mc_z']:g} Monte Carlo SE",
                 "tolerance": tol, "difference": d, "pass": d <= tol if tol > 0 else d == 0})


def tie_ranks(rho, tol=1e-12):
    """Smallest and largest rank (1 = largest rho) of each value's tie group; values within tol are tied."""
    order = np.argsort(-rho, kind="stable")
    v = -rho[order]
    new = np.r_[True, np.diff(v) > tol]
    start = np.maximum.accumulate(np.where(new, np.arange(len(v)), 0))
    end = np.minimum.accumulate(np.where(np.r_[new[1:], True], np.arange(len(v)), len(v))[::-1])[::-1]
    lo, hi = np.empty(len(v)), np.empty(len(v))
    lo[order], hi[order] = start + 1, end + 1
    return lo, hi


def rank_check(analysis, quantity, pipe, lo, hi, impl):
    """Ranks: pass if the pipeline's rank lies within the independent tie group (ties are ordered arbitrarily)."""
    pipe, lo, hi = np.atleast_1d(pipe).astype(float), np.atleast_1d(lo), np.atleast_1d(hi)
    d = np.where(pipe < lo, lo - pipe, np.where(pipe > hi, pipe - hi, 0.0))
    j = int(np.argmax(d))
    ROWS.append({"analysis": analysis, "quantity": quantity, "pipeline_value": pipe[j],
                 "independent_value": lo[j] if lo[j] == hi[j] else f"{lo[j]:.0f}-{hi[j]:.0f} (tied)",
                 "independent_implementation": impl, "criterion": "rank within the independent tie group (rho equal to 1e-12)",
                 "tolerance": 0, "difference": d[j], "pass": d[j] == 0})


def prop_se(a, b, na, nb):
    pbar = (a * na + b * nb) / (na + nb)
    return np.sqrt(pbar * (1 - pbar) * (1 / na + 1 / nb))


def quantile_se(draws, q):
    """Asymptotic Monte Carlo SE of a sample quantile: sqrt(q (1 - q) / B) / f(x_q), f by Gaussian KDE."""
    x = np.asarray(draws, float)
    return np.sqrt(q * (1 - q) / len(x)) / gaussian_kde(x)(np.quantile(x, q))[0]


def lstsq_resid(Y, Z):
    """Residuals of each column of Y (n x G) on Z by least squares."""
    return Y - Z @ np.linalg.lstsq(Z, Y, rcond=None)[0]


def adjusted_all(L, mg, batch):
    """Batch-adjusted Spearman of MGMT with every row of L: pandas ranks, lstsq residuals, explicit correlation."""
    n = L.shape[1]
    Z = np.column_stack([np.ones(n), batch])
    RY = lstsq_resid(pd.DataFrame(L.T).rank().to_numpy(), Z)
    rx = lstsq_resid(pd.Series(mg).rank().to_numpy()[:, None], Z)[:, 0]
    rho = (RY * rx[:, None]).sum(0) / np.sqrt((RY ** 2).sum(0) * (rx ** 2).sum())
    df = n - 3
    p = 2 * stats.t.sf(np.abs(rho * np.sqrt(df / (1 - rho ** 2))), df)
    return rho, p


def screen(t, prim, gl):
    al = prim["aliquot"].values
    L = t.log2tpm(al)
    mg = L[t.row[MGMT]]
    batch = (prim["batch"].values == "GDC2024").astype(float)
    g = gl.set_index("gene_id")
    el = g.index[g["eligible_20"] & (g.index != MGMT)].to_numpy()
    top = g.loc[el].sort_values("rank_adj").index[:P["top_genes"]].to_numpy()
    rest = np.setdiff1d(el, np.r_[top, [GPR153]])
    rnd = np.random.default_rng(common.SEED).choice(rest, P["random_genes"], replace=False)
    for label, genes in (("GPR153", [GPR153]), ("top 100", top), ("500 random", rnd)):
        for gid in genes:
            d = pd.DataFrame({"x": mg, "y": L[t.row[gid]], "batch": batch})
            r = pingouin.partial_corr(d, x="x", y="y", covar=["batch"], method="spearman")
            impl = f"pingouin {pingouin.__version__} partial_corr(method='spearman')"
            det(f"screen, batch-adjusted: {label}", f"{gid} rho", g.loc[gid, "rho_adj"], r["r"].iloc[0], impl)
            det(f"screen, batch-adjusted: {label}", f"{gid} p", g.loc[gid, "p_adj"], r["p-val"].iloc[0], impl)

    rows = np.array([t.row[x] for x in el])
    rho, p = adjusted_all(L[rows], mg, batch)
    q = multipletests(p, method="fdr_bh")[1]
    lo, hi = tie_ranks(rho)
    impl = "pandas ranks, numpy lstsq residuals, statsmodels fdr_bh"
    i = int(np.where(el == GPR153)[0][0])
    for k, v in (("rho", rho), ("p", p), ("q", q)):
        det("screen, batch-adjusted: GPR153 (all-gene implementation)", f"GPR153 {k}", g.loc[GPR153, f"{k}_adj"], v[i], impl)
    rank_check("screen, batch-adjusted: GPR153 (all-gene implementation)", "GPR153 rank", g.loc[GPR153, "rank_adj"],
               lo[i], hi[i], impl)
    gpl = g.loc[el]
    rank_check("screen, batch-adjusted: all eligible genes", "rank (largest difference)", gpl["rank_adj"].values, lo, hi, impl)
    for k, v in (("rho", rho), ("p", p), ("q", q)):
        ds = [rel(a, b) for a, b in zip(gpl[f"{k}_adj"].values, v)]
        j = int(np.argmax(ds))
        det("screen, batch-adjusted: all eligible genes", f"largest difference in {k} ({el[j]})",
            gpl[f"{k}_adj"].values[j], v[j], impl)

    nc = g.index[g["rho_pooled"].notna()].to_numpy()
    pr = np.array([stats.spearmanr(mg, L[t.row[x]]) for x in nc])
    for k, col in (("rho", 0), ("p", 1)):
        ds = [rel(a, b) for a, b in zip(g.loc[nc, f"{k}_pooled"].values, pr[:, col])]
        j = int(np.argmax(ds))
        det("screen, pooled: all non-constant genes", f"largest difference in {k} ({nc[j]})",
            g.loc[nc[j], f"{k}_pooled"], pr[j, col], f"scipy {__import__('scipy').__version__} spearmanr")
        det("screen, pooled: GPR153", f"GPR153 {k}", g.loc[GPR153, f"{k}_pooled"], pr[list(nc).index(GPR153), col],
            "scipy spearmanr")
    pe = pd.Series(pr[:, 1], index=nc).loc[el]
    qp = multipletests(pe.values, method="fdr_bh")[1]
    plo, phi = tie_ranks(pd.Series(pr[:, 0], index=nc).loc[el].values)
    ds = [rel(a, b) for a, b in zip(gpl["q_pooled"].values, qp)]
    j = int(np.argmax(ds))
    det("screen, pooled: all eligible genes", f"largest difference in q ({el[j]})", gpl["q_pooled"].values[j], qp[j],
        "statsmodels fdr_bh")
    det("screen, pooled: GPR153", "GPR153 q", g.loc[GPR153, "q_pooled"], qp[i], "statsmodels fdr_bh")
    rank_check("screen, pooled: all eligible genes", "rank (largest difference)", gpl["rank_pooled"].values, plo, phi,
               "scipy spearmanr rho, tie groups")
    rank_check("screen, pooled: GPR153", "GPR153 rank", g.loc[GPR153, "rank_pooled"], plo[i], phi[i],
               "scipy spearmanr rho, tie groups")
    return L, mg, batch, el


def permutation(L, mg, batch, el, t):
    s = pd.read_csv(common.RES / "screen/permutation_null_summary.csv").iloc[0]
    dr = pd.read_csv(common.RES / "screen/permutation_null_draws.csv")
    n = L.shape[1]
    Z = np.column_stack([np.ones(n), batch])
    RY = lstsq_resid(pd.DataFrame(L[[t.row[x] for x in el]].T).rank().to_numpy(), Z)
    RY /= np.sqrt((RY ** 2).sum(0))
    rng = np.random.default_rng(P["seed"])
    groups = [np.where(batch == b)[0] for b in (0.0, 1.0)]
    B, k = int(s["permutations"]), int(s["observed_gpr153_rank"])
    mx, at = np.empty(B), np.empty(B)
    for b in range(B):
        v = mg.copy()
        for g in groups:
            v[g] = v[rng.permutation(g)]
        rx = lstsq_resid(stats.rankdata(v)[:, None], Z)[:, 0]
        r = RY.T @ (rx / np.sqrt((rx ** 2).sum()))
        mx[b] = r.max()
        at[b] = np.sort(r)[-k]
    impl = f"one permutation at a time, lstsq residuals, seed {P['seed']}"
    for q, pv, iv in (("prop_max_ge_observed_max", s["prop_max_ge_observed_max"], (mx >= s["observed_max_rho"]).mean()),
                      ("prop_rank_value_ge_gpr153", s["prop_rank_value_ge_gpr153"], (at >= s["observed_gpr153_rho"]).mean())):
        mc("permutation null", q, pv, iv, prop_se(pv, iv, B, B), impl)
    for q, col in ((0.5, "max_rho_p50"), (0.95, "max_rho_p95"), (0.99, "max_rho_p99")):
        se = np.hypot(quantile_se(dr["max_rho"], q), quantile_se(mx, q))
        mc("permutation null", col, s[col], np.quantile(mx, q), se, impl)


def boot_draw(seed, X, mg, batch, gpos):
    idx = np.random.default_rng(seed).choice(X.shape[1], X.shape[1], replace=True)
    Xb, mb, bb = X[:, idx], mg[idx], batch[idx]
    Z = np.column_stack([np.ones(len(idx)), bb])
    RY = lstsq_resid(pd.DataFrame(Xb.T).rank().to_numpy(), Z)
    rx = lstsq_resid(pd.Series(mb).rank().to_numpy()[:, None], Z)[:, 0]
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = (RY * rx[:, None]).sum(0) / np.sqrt((RY ** 2).sum(0) * (rx ** 2).sum())
    rho = np.where(np.isnan(rho), -np.inf, rho)
    return rho[gpos], int((rho > rho[gpos]).sum() + 1)


def bootstrap(L, mg, batch, el, t):
    s = pd.read_csv(common.RES / "screen/bootstrap_rank_summary.csv").iloc[0]
    dr = pd.read_csv(common.RES / "screen/bootstrap_rank_draws.csv")
    X = L[[t.row[x] for x in el]]
    gpos = int(np.where(el == GPR153)[0][0])
    B = int(s["resamples"])
    seeds = np.random.SeedSequence(P["seed"]).generate_state(B)
    res = np.array(Parallel(n_jobs=-1)(delayed(boot_draw)(int(x), X, mg, batch, gpos) for x in seeds))
    impl = f"pandas ranks, lstsq residuals, seed {P['seed']}"
    for q, pv, iv in (("prop_rank1", s["prop_rank1"], (res[:, 1] == 1).mean()),
                      ("prop_top10", s["prop_top10"], (res[:, 1] <= 10).mean())):
        mc("bootstrap rank", q, pv, iv, prop_se(pv, iv, B, B), impl)
    for q, col in ((0.025, "rho_ci_low"), (0.975, "rho_ci_high")):
        se = np.hypot(quantile_se(dr["gpr153_rho"], q), quantile_se(res[:, 0], q))
        mc("bootstrap rank", col, s[col], np.percentile(res[:, 0], q * 100), se, impl)


def r_checks():
    d = pd.read_csv(common.work_path("survival_table_primary.csv"))
    cx = pd.read_csv(common.RES / "survival/cox_models.csv")
    m3 = cx[(cx.cohort == "primary") & (cx.model == "Model 3") & (cx.variant == "main")
            & (cx.sample_set == "complete cases")].set_index("covariate")
    meta_in, meta_pipe = [], []
    ce = pd.read_csv(common.RES / "external/external_coexpression.csv")
    tc = pd.read_csv(common.RES / "screen/screen_sensitivity.csv").set_index("run").loc["primary: batch-adjusted"]
    rz = lambda rho, n, k=0: (np.arctanh(rho), (1 + rho ** 2 / 2) / (n - 3 - k))
    em = pd.read_csv(common.RES / "external/external_meta.csv").set_index("analysis")
    for a, extra in (("external cohorts", False), ("external cohorts + TCGA", True)):
        y, v = rz(ce["spearman_rho"].values, ce["n"].values)
        if extra:
            yt, vt = rz(tc.gpr153_rho, tc.n, 1)
            y, v = np.r_[y, yt], np.r_[v, vt]
        meta_in += [{"analysis": f"co-expression: {a}", "yi": a_, "vi": b_} for a_, b_ in zip(y, v)]
        meta_pipe.append((f"co-expression: {a}", em.loc[a]))
    gp = pd.read_csv(common.RES / "external/gse13041_platforms.csv")
    for arr in ("all arrays", "Primary GBM"):
        e = gp[(gp.arrays == arr) & gp.platform.isin(common.PARAMS["gse13041_platforms"])]
        pooled = gp[(gp.arrays == arr) & ~gp.platform.isin(common.PARAMS["gse13041_platforms"])].iloc[0]
        if pd.notna(pooled.get("tau2")):
            y, v = rz(e["spearman_rho"].values, e["n"].values)
            meta_in += [{"analysis": f"GSE13041 platforms: {arr}", "yi": a_, "vi": b_} for a_, b_ in zip(y, v)]
            meta_pipe.append((f"GSE13041 platforms: {arr}", pd.Series({"mu": np.arctanh(pooled.spearman_rho),
                                                                       "p": pooled.spearman_p, "tau2": pooled.tau2,
                                                                       "I2": pooled.I2})))
    sx = pd.read_csv(common.RES / "external/survival_replication_cox.csv")
    sx = sx[~((sx.cohort == "CPTAC-3 GBM") & (sx.subset == "all"))]
    sm = pd.read_csv(common.RES / "external/survival_replication_meta.csv")
    tcx = cx[(cx.cohort == "primary") & (cx.variant == "main") & (cx.sample_set == "complete cases")
             & (cx.covariate == "gpr153")].set_index("model")
    for model, tm in (("(a) GPR153", "Model 1"), ("(b) GPR153 + age + sex", "Model 3")):
        e = sx[sx.model == model]
        for wt in (False, True):
            y, v = e["log_hr"].values, e["se"].values ** 2
            if wt:
                y, v = np.r_[y, tcx.loc[tm, "coef"]], np.r_[v, tcx.loc[tm, "se"] ** 2]
            a = f"survival: {model}, {'with' if wt else 'without'} TCGA"
            meta_in += [{"analysis": a, "yi": a_, "vi": b_} for a_, b_ in zip(y, v)]
            meta_pipe.append((a, sm[(sm.model == model) & (sm.with_tcga == wt)].iloc[0]))

    with tempfile.TemporaryDirectory(dir=common.work_path("")) as td:
        td = Path(td)
        d[["time", "event", "gpr153", "age", "male"]].dropna().to_csv(td / "m3.csv", index=False)
        d[["time", "event", "batch", "gpr153_log2tpm", "mgmt_log2tpm"]].to_csv(td / "km.csv", index=False)
        (td / "km_months.txt").write_text(",".join(map(str, common.PARAMS["km_at_risk_months"])))
        pd.DataFrame(meta_in).to_csv(td / "meta.csv", index=False)
        common.run_r("pipeline/a91_independent.R", td)
        r3 = pd.read_csv(td / "m3_out.csv").set_index("term")
        klr, kgr, kar = (pd.read_csv(td / f"km_{x}_out.csv") for x in ("logrank", "groups", "at_risk"))
        mo = pd.read_csv(td / "meta_out.csv").set_index("analysis")

    impl = "R survival::coxph (Efron, eps 1e-12)"
    for v in r3.index:
        for c in ("hr", "hr_low", "hr_high", "p"):
            det("Model 3", f"{v} {c}", m3.loc[v, c], r3.loc[v, c], impl)
    det("Model 3", "n", m3["n"].iloc[0], r3["n"].iloc[0], impl)
    det("Model 3", "events", m3["events"].iloc[0], r3["events"].iloc[0], impl)

    lr = pd.read_csv(common.RES / "survival/km_logrank.csv").set_index(["gene", "split"])
    for r in klr.itertuples():
        for c in ("n", "chisq", "p", "median_followup_days_reverse_km"):
            det("log-rank", f"{r.gene}, {r.split}: {c}", lr.loc[(r.gene, r.split), c], getattr(r, c), "R survival::survdiff")
    kg = pd.read_csv(common.RES / "survival/km_groups.csv").set_index(["gene", "split", "group"])
    for r in kgr.itertuples():
        for c in ("n", "deaths", "median_os_days", "median_os_ci_low", "median_os_ci_high", "median_followup_days_reverse_km"):
            det("Kaplan-Meier", f"{r.gene}, {r.split}, {r.group}: {c}", kg.loc[(r.gene, r.split, r.group), c],
                getattr(r, c), "R survival::survfit (conf.type log-log)")
    ka = pd.read_csv(common.RES / "survival/km_at_risk.csv").set_index(["gene", "split", "group", "month"])
    for r in kar.itertuples():
        det("Kaplan-Meier", f"{r.gene}, {r.split}, {r.group}: at risk, month {r.month}",
            ka.loc[(r.gene, r.split, r.group, r.month), "at_risk"], r.at_risk, "R summary.survfit(times)")

    for a, pipe in meta_pipe:
        for c in ("mu", "se", "ci_low", "ci_high", "p", "Q", "Q_p", "tau2", "I2", "hk_ci_low", "hk_ci_high", "hk_p",
                  "tau2_ci_low", "tau2_ci_high", "pi_low", "pi_high"):
            if c in pipe.index and pd.notna(pipe[c]):
                det(f"meta-analysis, {a}", c, pipe[c], mo.loc[a, c], f"R metafor {r_version('metafor')}")


def r_version(pkg):
    import subprocess
    return subprocess.run(["Rscript", "-e", f"cat(as.character(packageVersion('{pkg}')))"], capture_output=True, text=True,
                          env=dict(__import__("os").environ, R_LIBS_USER=str(Path(".Rlib").resolve()))).stdout


def main():
    common.set_seeds()
    common.require_registration()
    out = common.outdir("independent")
    t = common.TCGA()
    prim = common.cohort_aliquots(t.m, "in_screen_primary")
    gl = pd.read_csv(common.RES / "screen/screen_gene_level.csv")
    L, mg, batch, el = screen(t, prim, gl)
    permutation(L, mg, batch, el, t)
    bootstrap(L, mg, batch, el, t)
    r_checks()
    d = pd.DataFrame(ROWS)
    d.to_csv(out / "independent_check.csv", index=False)
    s = d.groupby("analysis", sort=False).agg(values_compared=("pass", "size"), failed=("pass", lambda x: int((~x).sum())),
                                              criterion=("criterion", "first"), largest_difference=("difference", "max"))
    s["pass"] = s["failed"] == 0
    s.reset_index().to_csv(out / "independent_check_summary.csv", index=False)
    print(s.to_string())
    if not s["pass"].all():
        print(d[~d["pass"]].to_string())
        sys.exit("Independent reproduction failed (Section 13.2): see independent_check.csv")


if __name__ == "__main__":
    main()
