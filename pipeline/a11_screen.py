#!/usr/bin/env python3
"""
a11_screen.py
Genome-wide MGMT co-expression screen (PREREGISTRATION.md, Section 5) and its sensitivity runs (5.1), except
the bootstrap (a12), permutation null (a13), composition (a14) and methylation (a15) analyses.

Primary run: IDH-wildtype primary cohort (n = 252), log2(TPM + 1), batch-adjusted partial Spearman (k = 1) of
MGMT with every eligible gene (TPM >= 1 in >= 20%), MGMT excluded, BH over eligible genes. The pooled
(unadjusted) Spearman is computed alongside for every gene.

Sensitivity runs (one row each in screen_sensitivity.csv): pooled; BI only and GDC-2024 only (sign check,
Bonett-Wright CI with k = 0, within-batch rank); FPKM; MGMT outlier removal (screen) and GPR153/MGMT outlier
removal (single value); filters 10%, 50%, none; IDH-not-mutant cohort; sensitivity cohort; primary cohort
excluding IDH source disagreement. Subsets of the primary cohort keep the primary eligibility set; the
IDH-not-mutant and sensitivity cohorts use eligibility computed within the cohort.

Cross-batch reproducibility: per-gene Spearman between the BI and the GDC-2024 aliquot of the 79 primary-cohort
patients with both (85 across all primary tumors). GPR153 and MGMT are removed before this computation; their
own values are the post hoc diagnostic in a17_posthoc_crossbatch.py.

Inputs
  revision/work/tcga_expr.npz (a10), revision/manifests/tcga_samples.csv (r01)
Outputs (revision/results/screen/)
  screen_gene_level.csv          every non-constant gene: eligibility flags, mean log2(TPM + 1), fraction TPM >= 1,
                                 batch-adjusted and pooled rho, p, q, signed and absolute rank
  screen_sensitivity.csv         GPR153 rho, p, q, rank (and CI where defined) for every run
  screen_gene_universe.csv       gene counts at each step
  screen_gpr153_mgmt_values.csv  per-sample GPR153 and MGMT log2(TPM + 1) (Fig 2C)
  crossbatch_reproducibility.csv per-gene cross-batch Spearman (GPR153 and MGMT excluded)
  crossbatch_summary.csv

Run from the repository root:
    python pipeline/a11_screen.py
"""

import numpy as np
import pandas as pd

import common
from common import GPR153, MGMT


def screen(t, aliquots, eligible, covariate=None, unit="tpm", label=""):
    """Partial (or pooled) Spearman of MGMT with every gene non-constant in these samples. q and ranks over
    eligible genes, MGMT excluded."""
    x = t.log2tpm(aliquots) if unit == "tpm" else t.log2fpkm(aliquots)
    mg = x[t.row[MGMT]]
    tested = (x.std(1) > 0)
    rho = np.full(len(x), np.nan)
    p = np.full(len(x), np.nan)
    r, pv, df = common.partial_spearman(x[tested], mg, covariate)
    rho[tested], p[tested] = r, pv
    rho[t.row[MGMT]] = p[t.row[MGMT]] = np.nan
    el = eligible & tested
    el[t.row[MGMT]] = False
    q = np.full(len(x), np.nan)
    q[el] = common.bh(p[el])
    return {"rho": rho, "p": p, "q": q, "rank": common.signed_rank(rho, el), "abs_rank": common.abs_rank(rho, el),
            "n": len(aliquots), "df": df, "n_tested": int(el.sum()), "eligible": el, "label": label}


def gpr_row(res, run, cohort, k, unit="log2(TPM + 1)", filt="20%", note=""):
    i = res["_gpr"]
    row = {"run": run, "cohort": cohort, "unit": unit, "filter": filt, "k": k, "n": res["n"], "df": res["df"],
           "genes_tested": res["n_tested"], "gpr153_rho": res["rho"][i], "gpr153_p": res["p"][i],
           "gpr153_q": res["q"][i], "gpr153_rank": res["rank"][i], "gpr153_abs_rank": res["abs_rank"][i],
           "rank1_gene": res["_top"], "rank1_rho": np.nanmax(np.where(res["eligible"], res["rho"], np.nan)),
           "note": note}
    return row


def batch_cov(m_rows):
    return (m_rows["batch"].values == "GDC2024").astype(float)[:, None]


def main():
    common.set_seeds()
    common.require_registration()
    out = common.outdir("screen")
    t = common.TCGA()
    m = t.m
    prim = common.cohort_aliquots(m, "in_screen_primary")
    al = prim["aliquot"].values
    gi = t.row[GPR153]
    elig = {f: t.eligible(al, f) for f in (0.20, 0.10, 0.50, 0.0)}

    def run(rows, eligible, cov, unit="tpm"):
        r = screen(t, rows["aliquot"].values, eligible, cov, unit)
        r["_gpr"] = gi
        top = np.nanargmax(np.where(r["eligible"], r["rho"], -np.inf))
        r["_top"] = t.symbol[top] or t.gene_id[top]
        return r

    primary = run(prim, elig[0.20], batch_cov(prim))
    pooled = run(prim, elig[0.20], None)

    L = t.log2tpm(al)
    gl = pd.DataFrame({"gene_id": t.gene_id, "symbol": t.symbol, "gene_type": t.gene_type,
                       "mean_log2tpm": L.mean(1), "frac_tpm_ge_1": (t.tpm[:, t.cols(al)] >= 1).mean(1),
                       "eligible_20": elig[0.20], "eligible_10": elig[0.10], "eligible_50": elig[0.50],
                       "eligible_none": elig[0.0], "n": len(al)})
    for tag, r in (("adj", primary), ("pooled", pooled)):
        for k in ("rho", "p", "q", "rank", "abs_rank"):
            gl[f"{k}_{tag}"] = r[k]
    gl = gl[elig[0.0]]
    gl.to_csv(out / "screen_gene_level.csv", index=False)
    pd.DataFrame({"aliquot": al, "patient": prim["patient"].values, "batch": prim["batch"].values,
                  "gpr153_log2tpm": L[gi], "mgmt_log2tpm": L[t.row[MGMT]]}).to_csv(
        out / "screen_gpr153_mgmt_values.csv", index=False)

    rows = [gpr_row(primary, "primary: batch-adjusted", "primary", 1),
            gpr_row(pooled, "pooled (unadjusted)", "primary", 0)]
    for b in ("BI", "GDC2024"):
        sub = prim[prim["batch"] == b]
        r = run(sub, elig[0.20], None)
        lo, hi, _, _ = common.bw_ci(r["rho"][gi], r["n"], 0)
        row = gpr_row(r, f"{b} only (sign check)", "primary", 0)
        row.update(ci_low=lo, ci_high=hi, note="condition 3; Bonett-Wright 95% CI, k = 0; rank within batch")
        rows.append(row)
    rows.append(gpr_row(run(prim, elig[0.20], batch_cov(prim), "fpkm"), "FPKM", "primary", 1, "log2(FPKM + 1)"))

    mg = L[t.row[MGMT]]
    zm = (mg - mg.mean()) / mg.std(ddof=1)
    keep = np.abs(zm) <= common.PARAMS["outlier_abs_z"]
    rows.append(gpr_row(run(prim[keep], elig[0.20], batch_cov(prim[keep])), "MGMT outliers removed (|z| > 3)",
                        "primary", 1, note=f"{int((~keep).sum())} samples removed"))
    gp = L[gi]
    zg = (gp - gp.mean()) / gp.std(ddof=1)
    keep2 = keep & (np.abs(zg) <= common.PARAMS["outlier_abs_z"])
    r2, p2, df2 = common.partial_spearman(gp[keep2], mg[keep2], batch_cov(prim[keep2]))
    rows.append({"run": "GPR153 or MGMT outliers removed (single value)", "cohort": "primary",
                 "unit": "log2(TPM + 1)", "filter": "n/a", "k": 1, "n": int(keep2.sum()), "df": df2,
                 "gpr153_rho": r2[0], "gpr153_p": p2[0], "note": f"{int((~keep2).sum())} samples removed"})
    for f, lab in ((0.10, "10%"), (0.50, "50%"), (0.0, "none")):
        rows.append(gpr_row(run(prim, elig[f], batch_cov(prim)), f"filter {lab}", "primary", 1, filt=lab))
    for flag, lab in (("in_screen_idh_not_mutant", "IDH-not-mutant"), ("in_screen_sensitivity", "sensitivity")):
        c = common.cohort_aliquots(m, flag)
        rows.append(gpr_row(run(c, t.eligible(c["aliquot"].values, 0.20), batch_cov(c)), f"{lab} cohort",
                            lab, 1, note="eligibility computed within the cohort"))
    nd = prim[~prim["idh_source_disagreement"].astype(bool)]
    rows.append(gpr_row(run(nd, elig[0.20], batch_cov(nd)), "excluding IDH source disagreement", "primary", 1,
                        note=f"{len(prim) - len(nd)} patients excluded"))
    pd.DataFrame(rows).to_csv(out / "screen_sensitivity.csv", index=False)

    pd.DataFrame([
        {"step": "GDC v36 genes (after PAR_Y removal)", "n": len(t.gene_id)},
        {"step": "non-constant in primary cohort", "n": int(elig[0.0].sum())},
        {"step": "eligible, 10% filter", "n": int(elig[0.10].sum())},
        {"step": "eligible, 20% filter (primary)", "n": int(elig[0.20].sum())},
        {"step": "eligible, 50% filter", "n": int(elig[0.50].sum())},
        {"step": "tested in primary screen (eligible, MGMT excluded)", "n": primary["n_tested"]},
        {"step": "GPR153 eligible at 10% / 20% / 50%", "n": f"{elig[0.10][gi]} / {elig[0.20][gi]} / {elig[0.50][gi]}"},
        {"step": "MGMT eligible at 10% / 20% / 50%", "n": "{} / {} / {}".format(*(elig[f][t.row[MGMT]] for f in (0.10, 0.20, 0.50)))},
    ]).to_csv(out / "screen_gene_universe.csv", index=False)

    crossbatch(t, m, prim, out)
    print(pd.DataFrame(rows)[["run", "n", "gpr153_rho", "gpr153_q", "gpr153_rank"]].to_string(index=False))


def crossbatch(t, m, prim, out):
    drop = np.isin(t.gene_id, [GPR153, MGMT])
    summ = []
    for lab, pts in (("primary cohort", prim.loc[prim["patient_in_both_batches"], "patient"]),
                     ("all primary tumors", m.loc[m["selected_aliquot"] & m["in_screen_sensitivity"]
                                                  & m["patient_in_both_batches"], "patient"])):
        sel = m[m["selected_aliquot"] & m["patient"].isin(pts)].set_index("patient")
        sel = sel[sel["batch"] == "BI"]
        other = common.second_batch_aliquots(m, sel.index)
        sel = sel.loc[other.index]
        a = t.log2tpm(sel["aliquot"].values)[~drop]
        b = t.log2tpm(other.values)[~drop]
        ok = (a.std(1) > 0) & (b.std(1) > 0)
        ra, rb = common.rank_rows(a[ok]), common.rank_rows(b[ok])
        ra -= ra.mean(1, keepdims=True)
        rb -= rb.mean(1, keepdims=True)
        rho = (ra * rb).sum(1) / np.sqrt((ra ** 2).sum(1) * (rb ** 2).sum(1))
        d = pd.DataFrame({"gene_id": t.gene_id[~drop][ok], "symbol": t.symbol[~drop][ok], "set": lab,
                          "n_patients": len(sel), "crossbatch_spearman": rho})
        d.to_csv(out / f"crossbatch_reproducibility{'' if lab == 'primary cohort' else '_all_primary'}.csv",
                 index=False)
        q = np.quantile(rho, [0.05, 0.25, 0.5, 0.75, 0.95])
        summ.append({"set": lab, "n_patients": len(sel), "genes": int(ok.sum()), "p05": q[0], "p25": q[1],
                     "median": q[2], "p75": q[3], "p95": q[4]})
    pd.DataFrame(summ).to_csv(out / "crossbatch_summary.csv", index=False)


if __name__ == "__main__":
    main()
