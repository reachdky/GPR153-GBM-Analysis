"""
common.py
Shared paths, run mode, registration lock and statistical functions for the
analysis scripts (PREREGISTRATION.md, Sections 5-13).

Run mode (environment variable GBM_MODE):
  real (default)  inputs from data/ and revision/manifests, outputs under revision/.
                  Every script that reads expression or methylation values calls
                  require_registration() and stops unless the analysis plan has a
                  registration URL and its SHA-256 matches the recorded hash.
  sim             simulated inputs and outputs under revision/sim/ (pipeline/sim/).
"""

import hashlib
import json
import os
import random
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SEED = 42
MODE = os.environ.get("GBM_MODE", "real")
if MODE not in ("real", "sim"):
    sys.exit(f"GBM_MODE must be 'real' or 'sim', not {MODE!r}")

DL = Path("data/downloads")
ORIG = Path.home() / "gpr153_project_orig"  # original project, read-only
RAW = Path("data/raw")
BASE = Path("revision") if MODE == "real" else Path("revision/sim")
MAN = BASE / "manifests"
WORK = BASE / "work"
RES = BASE / "results"
SRC = BASE / "source_data"
FIG = BASE / "figures"
FIGMAP = BASE / "figure_script_map.csv"
PLAN = Path("revision/PREREGISTRATION.md")
PLAN_SHA = Path("revision/PREREGISTRATION.sha256")
PLAN_REG = Path("revision/PREREGISTRATION_REGISTRATION.md")

GPR153 = "ENSG00000158292"
MGMT = "ENSG00000170430"
GPR153_PROBES = ["221902_at", "64942_at"]
MGMT_PROBE = "204880_at"
DAYS_PER_MONTH = 365.25 / 12
Z975 = stats.norm.ppf(0.975)

NEURONAL = ["RBFOX3", "SNAP25", "SYT1", "STMN2", "GABRA1", "SYN1", "NEFL", "GAD1", "SLC17A7", "CAMK2A"]
FOREST_GENES = ["MGMT", "EGFR", "PTEN", "IDH1", "TERT"]

PARAMS = {
    "seed": SEED,
    "expression_unit": "log2(TPM + 1), tpm_unstranded, GENCODE v36",
    "eligibility_tpm": 1.0,
    "eligibility_fraction_primary": 0.20,
    "eligibility_fraction_sensitivity": [0.10, 0.50, 0.0],
    "outlier_abs_z": 3.0,
    "n_bootstrap_rank": 1000,
    "n_permutations": 10000,
    "permutation_block": 500,
    "stp27": {"intercept": 4.3215, "b_cg12434587": 0.5271, "b_cg12981137": 0.9265, "cutoff": 0.358,
              "beta_clip": [0.001, 0.999]},
    "condition4": {"min_n": 100, "min_kappa": 0.60, "k": 3, "n_bootstrap": 1000},
    "cox": {"penalizer": 0.0, "ties": "efron", "precision": 1e-12},
    "km_at_risk_months": [0, 6, 12, 18, 24, 36],
    "n_bootstrap_optimism": 1000,
    "mi": {"m": 20, "maxit": 20, "seed": SEED},
    "ora": {"max_genes": 200, "q": 0.05, "min_size": 10, "max_size": 500,
            "libraries": ["Reactome_2022", "KEGG_2021_Human", "GO_Biological_Process_2023"]},
    "gsea": {"min_size": 15, "max_size": 500, "permutations": 1000,
             "libraries": ["h.all.v2026.2.Hs.symbols", "c2.cp.reactome.v2026.2.Hs.symbols"]},
    "targeted_sets": ["HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_DNA_REPAIR"],
    "forest": {"n_estimators": 500, "min_samples_leaf": 15, "max_features": "sqrt", "random_state": SEED,
               "cv_folds": 5, "cv_repeats": 10, "permutation_importance_repeats": 10},
    "de": {"padj": 0.05, "abs_lfc": 1.0},
    "pca": {"n_genes": 2000, "k": 3, "n_init": 50},
    "power_rhos": [0.12, 0.15, 0.20],
    "gse13041_platforms": ["GPL96", "GPL570"],
    "mgmt_enrichment_max_genes": 200,
    "independent_check": {"rel_tol": 1e-6, "random_genes": 500, "top_genes": 100, "mc_z": 3.0, "seed": SEED + 1},
}

# lifelines Newton-Raphson convergence: the default stops at a relative difference of about 1e-5 from R coxph;
# 1e-12 converges to the same optimum as R (Section 13.2 requires agreement to 1e-6)
COX_FIT_OPTIONS = {"precision": PARAMS["cox"]["precision"]}


def set_seeds():
    random.seed(SEED)
    np.random.seed(SEED)
    os.environ.setdefault("PYTHONHASHSEED", str(SEED))


def require_registration():
    """Stop unless the analysis plan is registered (URL recorded, hash matches). No-op in sim mode."""
    if MODE == "sim":
        return
    reg = PLAN_REG.read_text() if PLAN_REG.exists() else ""
    url = re.search(r"\|\s*OSF project / registration URL\s*\|\s*(\S[^|]*?)\s*\|", reg)
    recorded = PLAN_SHA.read_text().split()[0] if PLAN_SHA.exists() else ""
    current = hashlib.sha256(PLAN.read_bytes()).hexdigest()
    if not url or "osf.io" not in url.group(1) or recorded != current:
        sys.exit("Locked: the analysis plan is not registered (no OSF URL in "
                 f"{PLAN_REG}, or {PLAN_SHA} does not match {PLAN}). "
                 "Real expression and methylation values are not read. Use GBM_MODE=sim to test.")


def outdir(*parts):
    d = RES.joinpath(*parts)
    d.mkdir(parents=True, exist_ok=True)
    return d


def work_path(name):
    WORK.mkdir(parents=True, exist_ok=True)
    return WORK / name


def save_matrix(name, **arrays):
    np.savez(work_path(name), **arrays)


def load_matrix(name):
    with np.load(work_path(name), allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def tcga_manifest():
    return pd.read_csv(MAN / "tcga_samples.csv")


def external_manifest():
    return pd.read_csv(MAN / "external_cohort_samples.csv")


def cohort_aliquots(m, flag):
    """Selected aliquots of an analysis cohort, in manifest order."""
    return m.loc[m[flag] & m["selected_aliquot"]].reset_index(drop=True)


def survival_table(m, flag="in_screen_primary"):
    """Patient-level survival table for an analysis cohort, from the manifest only (Section 6.7 rules)."""
    p = cohort_aliquots(m, flag)
    ok = p["os_event"].notna() & p["os_days"].notna() & (p["os_days"] > 0)
    s = p.loc[ok, ["patient", "aliquot", "batch", "os_days", "os_event"]].copy()
    s["os_event"] = s["os_event"].astype(int)
    return s.reset_index(drop=True)


# ── Correlation ───────────────────────────────────────────────────────────────

def rank_rows(x):
    """Average ranks along the last axis (ties averaged); NaN propagates as NaN rows."""
    return stats.rankdata(x, axis=-1)


def design(covariates, n):
    """Intercept plus covariate columns (n x (1 + k))."""
    cols = [np.ones(n)]
    if covariates is not None:
        c = np.asarray(covariates, dtype=float)
        cols += list(c.reshape(n, -1).T)
    return np.column_stack(cols)


def residualize(y, z):
    """Residuals of each row of y (... x n) after OLS on the columns of z (n x p)."""
    q, _ = np.linalg.qr(z)
    return y - (y @ q) @ q.T


def rank_covariates(covariates, n):
    """Continuous covariates enter partial Spearman as ranks; binary indicators are unchanged by ranking up to a
    linear transform, so ranking all columns gives the same residuals."""
    if covariates is None:
        return None
    c = np.asarray(covariates, dtype=float).reshape(n, -1)
    return np.column_stack([stats.rankdata(col) for col in c.T])


def partial_spearman(y, x, covariates=None):
    """Batch- (or covariate-) adjusted Spearman of each row of y (G x n) with x (n).
    Ranks y and x, regresses each rank vector on an intercept and the ranked covariates (k columns), and
    correlates the residuals. t = rho * sqrt(df / (1 - rho^2)), df = n - 2 - k. Returns rho, p, df."""
    y = np.atleast_2d(np.asarray(y, dtype=float))
    n = y.shape[1]
    zc = rank_covariates(covariates, n)
    z = design(zc, n)
    k = z.shape[1] - 1
    ry = residualize(rank_rows(y), z)
    rx = residualize(stats.rankdata(x), z)
    rho = (ry @ rx) / (np.linalg.norm(ry, axis=1) * np.linalg.norm(rx))
    df = n - 2 - k
    with np.errstate(divide="ignore", invalid="ignore"):
        t = rho * np.sqrt(df / (1 - rho ** 2))
    p = 2 * stats.t.sf(np.abs(t), df)
    return rho, p, df


def bh(p):
    """Benjamini-Hochberg adjusted p-values; NaN entries are ignored and returned as NaN."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    m = len(pv)
    if m == 0:
        return q
    order = np.argsort(pv)
    adj = pv[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


def bw_ci(rho, n, k=0, level=0.95):
    """Fisher-z CI with the Bonett-Wright SE sqrt((1 + rho^2 / 2) / (n - 3 - k))."""
    rho = np.asarray(rho, dtype=float)
    z = np.arctanh(rho)
    se = np.sqrt((1 + rho ** 2 / 2) / (n - 3 - k))
    zc = stats.norm.ppf(0.5 + level / 2)
    return np.tanh(z - zc * se), np.tanh(z + zc * se), z, se


def corr_power(rho, n, alpha=0.05):
    """Two-sided power of a test of rho = 0 via Fisher z."""
    zc = stats.norm.ppf(1 - alpha / 2)
    d = np.abs(np.arctanh(rho)) * np.sqrt(n - 3)
    return stats.norm.cdf(d - zc) + stats.norm.cdf(-d - zc)


def signed_rank(rho, eligible):
    """Rank 1 = largest positive rho among eligible entries; NaN elsewhere."""
    r = np.full(len(rho), np.nan)
    idx = np.where(eligible & ~np.isnan(rho))[0]
    order = idx[np.argsort(-rho[idx], kind="stable")]
    r[order] = np.arange(1, len(order) + 1)
    return r


def abs_rank(rho, eligible):
    return signed_rank(np.abs(rho), eligible)


# ── Meta-analysis ─────────────────────────────────────────────────────────────

def meta_random(y, v, labels=None):
    """DerSimonian-Laird random effects with Hartung-Knapp(-Sidik-Jonkman) CI, Q-profile CI for tau^2 and a
    95% prediction interval. y: effect estimates (Fisher z or log HR); v: their variances."""
    from scipy.optimize import brentq
    y, v = np.asarray(y, float), np.asarray(v, float)
    K = len(y)
    w = 1 / v
    mu_fe = np.sum(w * y) / np.sum(w)
    Q = float(np.sum(w * (y - mu_fe) ** 2))
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - (K - 1)) / c) if K > 1 else 0.0
    ws = 1 / (v + tau2)
    mu = float(np.sum(ws * y) / np.sum(ws))
    se = float(np.sqrt(1 / np.sum(ws)))
    out = {"k": K, "mu": mu, "se": se, "ci_low": mu - Z975 * se, "ci_high": mu + Z975 * se,
           "p": float(2 * stats.norm.sf(abs(mu / se))), "Q": Q, "Q_df": K - 1,
           "Q_p": float(stats.chi2.sf(Q, K - 1)) if K > 1 else np.nan,
           "tau2": tau2, "I2": max(0.0, (Q - (K - 1)) / Q) if Q > 0 else 0.0}
    if K > 1:
        q = float(np.sum(ws * (y - mu) ** 2) / (K - 1))
        se_hk = np.sqrt(q / np.sum(ws))
        tq = stats.t.ppf(0.975, K - 1)
        out.update(hk_se=float(se_hk), hk_ci_low=mu - tq * se_hk, hk_ci_high=mu + tq * se_hk,
                   hk_p=float(2 * stats.t.sf(abs(mu / se_hk), K - 1)))

        def qgen(t):
            wt = 1 / (v + t)
            m_t = np.sum(wt * y) / np.sum(wt)
            return np.sum(wt * (y - m_t) ** 2)

        hi_crit, lo_crit = stats.chi2.ppf(0.975, K - 1), stats.chi2.ppf(0.025, K - 1)
        upper_t = 1e4 * (np.var(y) + np.max(v)) + 1
        out["tau2_ci_low"] = 0.0 if qgen(0) <= hi_crit else brentq(lambda t: qgen(t) - hi_crit, 0, upper_t)
        out["tau2_ci_high"] = 0.0 if qgen(0) <= lo_crit else brentq(lambda t: qgen(t) - lo_crit, 0, upper_t)
    if K > 2:
        tp = stats.t.ppf(0.975, K - 2)
        half = tp * np.sqrt(tau2 + se ** 2)
        out.update(pi_low=mu - half, pi_high=mu + half)
    if labels is not None:
        out["cohorts"] = ";".join(labels)
    return out


# ── Survival ──────────────────────────────────────────────────────────────────

def cox_univariate(x, time, event, max_iter=50, tol=1e-9):
    """Vectorized single-covariate Cox model (Efron ties) for every row of x (G x n), by Newton-Raphson with
    step halving. Returns beta, se (from the observed information), Wald z and p."""
    x = np.atleast_2d(np.asarray(x, float))
    time, event = np.asarray(time, float), np.asarray(event, int)
    order = np.argsort(time, kind="stable")
    x, time, event = x[:, order], time[order], event[order]
    G, n = x.shape
    ev_times = np.unique(time[event == 1])
    J = len(ev_times)
    start = np.searchsorted(time, ev_times, side="left")
    E = np.zeros((n, J))
    E[np.where(event == 1)[0], np.searchsorted(ev_times, time[event == 1])] = 1
    d = E.sum(0)
    xd_sum = x @ E.sum(1)

    def pieces(b):
        eta = b[:, None] * x
        eta -= eta.max(1, keepdims=True)
        w = np.exp(eta)
        wx, wxx = w * x, w * x * x
        rc = lambda a: np.cumsum(a[:, ::-1], 1)[:, ::-1][:, start]
        S0, S1, S2 = rc(w), rc(wx), rc(wxx)
        T0, T1, T2 = w @ E, wx @ E, wxx @ E
        ll = eta @ E.sum(1)  # linear predictor over deaths, on the shifted scale used for phi
        U, I = xd_sum.copy(), np.zeros(G)
        logterm = np.zeros(G)
        for l in range(int(d.max())):
            m = d > l
            f = np.where(m, l / d, 0.0)
            phi = S0 - f * T0
            a1 = S1 - f * T1
            a2 = S2 - f * T2
            with np.errstate(divide="ignore", invalid="ignore"):
                U -= np.where(m, a1 / phi, 0).sum(1)
                I += np.where(m, a2 / phi - (a1 / phi) ** 2, 0).sum(1)
                logterm += np.where(m, np.log(phi), 0).sum(1)
        return ll - logterm, U, I

    b = np.zeros(G)
    ll, U, I = pieces(b)
    for _ in range(max_iter):
        step = U / I
        newb = b + step
        nll, nU, nI = pieces(newb)
        for _ in range(20):
            worse = nll < ll - 1e-12
            if not worse.any():
                break
            step = np.where(worse, step / 2, step)
            newb = b + step
            nll, nU, nI = pieces(newb)
        b, ll, U, I = newb, nll, nU, nI
        if np.nanmax(np.abs(step)) < tol:
            break
    se = 1 / np.sqrt(I)
    z = b / se
    return b, se, z, 2 * stats.norm.sf(np.abs(z))


def standardize(v):
    v = np.asarray(v, float)
    return (v - np.nanmean(v)) / np.nanstd(v, ddof=1)


def rcs_basis(x, knots):
    """Restricted cubic spline nonlinear term(s) (Harrell), for 3 knots one column."""
    x = np.asarray(x, float)
    t = np.asarray(knots, float)
    k = len(t)
    pos = lambda u: np.maximum(u, 0) ** 3
    cols = []
    for j in range(k - 2):
        cols.append((pos(x - t[j]) - pos(x - t[k - 2]) * (t[k - 1] - t[j]) / (t[k - 1] - t[k - 2])
                     + pos(x - t[k - 1]) * (t[k - 2] - t[j]) / (t[k - 1] - t[k - 2])) / (t[k - 1] - t[0]) ** 2)
    return np.column_stack(cols)


# ── Records ───────────────────────────────────────────────────────────────────

def record_figure(item, script, output, source_data):
    """Add or replace a row in figure_script_map.csv."""
    row = {"item": item, "script": script, "output": str(output), "source_data": ";".join(map(str, source_data))}
    if FIGMAP.exists():
        f = pd.read_csv(FIGMAP)
        f = f[f["item"] != item]
        f = pd.concat([f, pd.DataFrame([row])], ignore_index=True)
    else:
        f = pd.DataFrame([row])
    f.sort_values("item").to_csv(FIGMAP, index=False)


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, default=str))


def run_r(script, *args):
    env = dict(os.environ, R_LIBS_USER=str(Path(".Rlib").resolve()))
    subprocess.run(["Rscript", str(script), *map(str, args)], check=True, env=env)


# ── TCGA expression helpers ───────────────────────────────────────────────────

class TCGA:
    """TCGA expression (from a10_tcga_matrix.py) aligned to the manifest."""

    def __init__(self):
        z = load_matrix("tcga_expr.npz")
        self.gene_id, self.symbol, self.gene_type = z["gene_id"], z["symbol"], z["gene_type"]
        self.aliquot = z["aliquot"]
        self.tpm, self.fpkm, self.counts = z["tpm"], z["fpkm"], z["counts"]
        self.col = {a: i for i, a in enumerate(self.aliquot)}
        self.row = {g: i for i, g in enumerate(self.gene_id)}
        self.m = tcga_manifest()

    def cols(self, aliquots):
        return np.array([self.col[a] for a in aliquots])

    def log2tpm(self, aliquots):
        return np.log2(self.tpm[:, self.cols(aliquots)] + 1)

    def log2fpkm(self, aliquots):
        return np.log2(self.fpkm[:, self.cols(aliquots)] + 1)

    def gene(self, gid, aliquots, unit="tpm"):
        x = self.tpm if unit == "tpm" else self.fpkm
        return np.log2(x[self.row[gid], self.cols(aliquots)] + 1)

    def eligible(self, aliquots, fraction):
        """TPM >= 1 in at least `fraction` of the samples, and non-constant log2(TPM + 1)."""
        t = self.tpm[:, self.cols(aliquots)]
        nonconst = np.log2(t + 1).std(1) > 0
        if fraction == 0:
            return nonconst
        return nonconst & ((t >= PARAMS["eligibility_tpm"]).mean(1) >= fraction)


def second_batch_aliquots(m, patients):
    """For patients with primary-tumor aliquots in both batches: the GDC-2024 aliquot chosen by the r01 rule
    (lowest vial, portion, plate, then full barcode)."""
    g = m[(m["sample_type"] == "Primary Tumor") & (m["batch"] == "GDC2024") & m["patient"].isin(patients)].copy()
    g["vial"], g["portion"], g["plate"] = g["aliquot"].str[15], g["aliquot"].str[17:19].astype(int), g["aliquot"].str[21:25]
    g = g.sort_values(["patient", "vial", "portion", "plate", "aliquot"]).groupby("patient").head(1)
    return g.set_index("patient")["aliquot"]
