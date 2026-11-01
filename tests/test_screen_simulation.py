#!/usr/bin/env python3
"""
T11 and T12 (PREREGISTRATION.md, Section 13.1), on simulated data (seed 42), through the screen function of
pipeline/a11_screen.py (batch-adjusted partial Spearman, BH over eligible genes, signed rank).

Simulated data: n = 252 (146 BI, 106 GDC-2024), 37,058 genes plus MGMT, log2-scale expression from a latent-factor
model (6 factors; 30% of loadings non-zero), batch shifts in 10% of genes and in MGMT, and 5% of genes with half
their values at zero (ties).

  T11 null calibration: MGMT is shuffled within batch and the full screen is rerun, 200 times. Pass if the maximum
      deviation of the empirical distribution of the pooled p-values (200 runs x tested genes) from uniform is
      <= 0.01 and the fraction of runs with any gene at q < 0.05 is <= 0.09.
  T12 planted effect, in the first null data set: one gene is replaced by a planted gene with batch-adjusted
      Spearman rho of about 0.45 with MGMT (normal-scores construction, Pearson r = 2 sin(pi rho / 6)). Pass if
      it is rank 1 with q < 0.05. A batch-driven pseudo-correlation is planted separately (MGMT and one gene each
      shifted by 2.5 SD between batches, independent within batch): pass if the pooled screen gives q < 0.05 and
      the batch-adjusted screen gives q >= 0.05 for that gene. Reported, not pass/fail: detection rates (q < 0.05;
      rank 1) at rho 0.30 and 0.45 across the 200 null data sets, by replacing one gene's p-value and rho with the
      planted gene's (each gene's statistic is computed independently, so this equals rerunning the screen).
Run from the repository root:  python tests/test_screen_simulation.py
"""

import sys
import time

import numpy as np

from _report import record, common

sys.path.insert(0, "pipeline")
from a11_screen import screen  # noqa: E402

rng = np.random.default_rng(common.SEED)
N_BI, N_GDC, G, RUNS = 146, 106, 37058, 200
n = N_BI + N_GDC
batch = np.r_[np.zeros(N_BI), np.ones(N_GDC)]


class SimTCGA:
    """The attributes of common.TCGA that a11_screen.screen uses."""

    def __init__(self, L):
        self.L = L
        self.row = {g: i for i, g in enumerate([f"g{i}" for i in range(G)] + [common.MGMT])}

    def log2tpm(self, aliquots):
        return self.L[:, aliquots]


def simulate():
    f = rng.normal(size=(6, n))
    load = rng.normal(0, 0.35, (G + 1, 6)) * (rng.random((G + 1, 6)) < 0.3)
    load[G] = [0.7, 0, 0, 0, 0, 0]
    L = rng.normal(4, 1.5, (G + 1, 1)) + load @ f + rng.normal(0, 0.5, (G + 1, n))
    shift = rng.normal(0, 0.3, G + 1) * (rng.random(G + 1) < 0.1)
    shift[G] = 0.3
    L += shift[:, None] * batch[None, :]
    low = np.where(rng.random(G) < 0.05)[0]
    L[low] *= rng.random((len(low), n)) < 0.5
    return np.clip(L, 0, None)


def shuffle_within_batch(v):
    out = v.copy()
    for b in (0.0, 1.0):
        i = np.where(batch == b)[0]
        out[i] = v[rng.permutation(i)]
    return out


def planted(mgmt, rho_s):
    """Gene with batch-adjusted Spearman of about rho_s with mgmt (normal scores, Pearson r = 2 sin(pi rho / 6))."""
    m = mgmt - np.array([mgmt[batch == b].mean() for b in batch])
    m = m / m.std()
    r = 2 * np.sin(np.pi * rho_s / 6)
    return 4 + r * m + np.sqrt(1 - r ** 2) * rng.normal(size=n)


t0 = time.time()
L0 = simulate()
t = SimTCGA(L0.copy())
al = np.arange(n)
eligible = np.ones(G + 1, bool)
cov = batch[:, None]
mg0 = L0[G].copy()
P, RHO = np.empty((RUNS, G)), np.empty((RUNS, G))
any_q, mgmts = 0, []
for k in range(RUNS):
    t.L[G] = shuffle_within_batch(mg0)
    mgmts.append(t.L[G].copy())
    r = screen(t, al, eligible, cov)
    P[k], RHO[k] = r["p"][:G], r["rho"][:G]
    any_q += bool(np.nanmin(r["q"]) < 0.05)
n_untested = int(np.isnan(P[0]).sum())  # genes constant in the simulated data (all zero) are not tested
p_all = np.sort(P[~np.isnan(P)])
N = len(p_all)
ecdf_dev = max(np.max(np.arange(1, N + 1) / N - p_all), np.max(p_all - np.arange(N) / N))
frac_any = any_q / RUNS
ok11 = record("T11", "null calibration of the batch-adjusted screen (MGMT shuffled within batch, 200 runs)",
              ecdf_dev <= 0.01 and frac_any <= 0.09, "simulated",
              f"max ECDF deviation of {N} pooled p-values {ecdf_dev:.4f} (tolerance 0.01); fraction of runs with any "
              f"q < 0.05 {frac_any:.3f} (tolerance 0.09); fraction p < 0.05 {np.mean(P < 0.05):.4f}; "
              f"{G - n_untested} genes tested ({n_untested} constant), n = {n}, {time.time() - t0:.0f} s")

# T12a: planted rho ~ 0.45 in the first null data set, through the screen function
t.L[G] = mgmts[0]
t.L[0] = planted(mgmts[0], 0.45)
r = screen(t, al, eligible, cov)
rank_a, q_a, rho_a = r["rank"][0], r["q"][0], r["rho"][0]
pass_a = rank_a == 1 and q_a < 0.05

# T12b: batch-driven pseudo-correlation
t.L[0] = L0[0]
t.L[G] = shuffle_within_batch(mg0) + 2.5 * mg0.std() * batch
t.L[1] = 4 + rng.normal(size=n) + 2.5 * batch
rp = screen(t, al, eligible, None)
ra = screen(t, al, eligible, cov)
pass_b = rp["q"][1] < 0.05 and ra["q"][1] >= 0.05

# Detection rates across the 200 null data sets (gene 0 replaced by the planted gene)
rates = {}
for rho_s in (0.30, 0.45):
    det = top = 0
    for k in range(RUNS):
        rr, pp, _ = common.partial_spearman(planted(mgmts[k], rho_s), mgmts[k], cov)
        p = P[k].copy()
        p[0] = pp[0]
        det += common.bh(p)[0] < 0.05
        top += rr[0] > np.nanmax(RHO[k][1:])
    rates[rho_s] = (det / RUNS, top / RUNS)
ok12 = record("T12", "planted effect: rank 1 and q < 0.05 at rho 0.45; batch pseudo-correlation removed by adjustment",
              pass_a and pass_b, "simulated",
              f"planted rho {rho_a:.3f}: rank {rank_a:.0f}, q {q_a:.2e}; batch-driven gene: pooled rho {rp['rho'][1]:.3f} "
              f"q {rp['q'][1]:.2e}, adjusted rho {ra['rho'][1]:.3f} q {ra['q'][1]:.2f}; detection (q < 0.05) / rank 1 "
              f"over {RUNS} data sets: rho 0.30 {rates[0.30][0]:.2f} / {rates[0.30][1]:.2f}, rho 0.45 "
              f"{rates[0.45][0]:.2f} / {rates[0.45][1]:.2f} (power check, not pass/fail)")
sys.exit(0 if ok11 and ok12 else 1)
