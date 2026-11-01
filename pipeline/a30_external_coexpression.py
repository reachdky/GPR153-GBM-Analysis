#!/usr/bin/env python3
"""
a30_external_coexpression.py
External co-expression cohorts (PREREGISTRATION.md, Section 7).

Expression: microarrays (GSE13041, GSE4290) are log2-transformed if the matrix maximum exceeds 100, otherwise
used as deposited (values <= 0 become missing before log2); CGGA log2(FPKM + 1). GPR153 probe set: the one with
the highest mean expression in the cohort (221902_at, 64942_at); MGMT 204880_at. Every probe combination is
reported as a sensitivity analysis.

Per cohort: Spearman rho (primary) and Pearson r, 95% CI by Fisher z with the Bonett-Wright SE
sqrt((1 + rho^2 / 2) / (n - 3)); BH across the four cohorts. GSE13041 runs (a) 162 primary arrays, (b) duplicate
retained (163), (c) all 191 GPL96 arrays, (d) all GPL96 arrays after the duplicate rule (Spearman r >= 0.98 across
all probes except the GPR153 and MGMT probes; the later GSM is excluded), (e) without arrays from patients treated
before resection. U1 criterion (7.1): Pearson r for each GPR153 probe set x MGMT 204880_at on the deposited scale
and on log2, on the 163 "Primary GBM" arrays and all 191; confirmed only if one combination rounds to +0.116 and
+0.119 (3 decimals). CGGA pooled (sensitivity): Spearman pooled, and partial Spearman adjusted for release.
GSE13041 by platform (sensitivity; platforms with both GPR153 and MGMT probe sets: GPL96, GPL570): all arrays and
"Primary GBM" arrays per platform, each with its own highest-mean GPR153 probe set. Before pooling, a cross-platform
overlap check: shared tumor ID (title), identical clinical record (age, sex, survival time, vital status), and
whole-profile Spearman r >= 0.98 on the probe sets common to both platforms (GPR153 and MGMT probe sets excluded).
The platforms are pooled by random effects (as below) only if all three find no overlap. GSE4290 with all
histologies (sensitivity, all 180 arrays). Earlier same-session values (original project outputs; generating code
unavailable) are listed beside them, labeled "not reproduced".
Meta-analysis on Fisher z with variance (1 + rho^2 / 2) / (n - 3 - k), DerSimonian-Laird, with and without TCGA
(batch-adjusted primary estimate, k = 1), Hartung-Knapp sensitivity, Q-profile CI for tau^2, prediction interval.
Power at rho = 0.12, 0.15, 0.20 and the TCGA estimate. Discordance diagnostics: per-sample percentile of the GPR153
and MGMT features within each array or sample, and composition tables from the manifest.

Inputs
  revision/work/ext_{GSE13041,GSE4290,CGGA_693,CGGA_325}.npz and _clinical.csv (a29),
  revision/results/screen/screen_sensitivity.csv (a11)
Outputs (revision/results/external/)
  external_coexpression.csv, external_probe_combinations.csv, gse13041_runs.csv, gse13041_u1_criterion.csv,
  cgga_pooled.csv, external_meta.csv, external_power.csv, external_discordance.csv, external_composition.csv,
  gse13041_platforms.csv, gse13041_platform_overlap.csv, gse4290_all_histologies.csv, earlier_same_session_values.csv

Run from the repository root:
    python pipeline/a30_external_coexpression.py
"""

import itertools

import numpy as np
import pandas as pd
from scipy import stats

import common

ARRAY = {"GSE13041": "GSE13041 GPL96 primary", "GSE4290": "GSE4290 GBM"}
EARLIER = [  # (file under the original project, line, label); line 1 is the header
    ("data/results/validation/cross_dataset_validation.csv", 4, "GSE13041, all platforms pooled"),
    ("data/results/validation/cross_dataset_validation.csv", 3, "GSE4290, all histologies"),
    ("data/results/validation/cross_dataset_validation_filtered.csv", 4, "GSE13041 GPL570"),
    ("data/results/validation/cross_dataset_validation_filtered.csv", 5, "GSE13041 'Primary' (n matches GPL570 Primary GBM)"),
]
RNA = {"CGGA_693": "CGGA mRNAseq_693", "CGGA_325": "CGGA mRNAseq_325"}
EXCLUDE = set(common.GPR153_PROBES) | {common.MGMT_PROBE}


class Cohort:
    def __init__(self, name, array=None):
        z = common.load_matrix(f"ext_{name}.npz")
        self.name, self.X, self.feature, self.sample = name, z["X"], z["feature"], z["sample"]
        self.clin = pd.read_csv(common.work_path(f"ext_{name}_clinical.csv"))
        self.row = {f: i for i, f in enumerate(self.feature)}
        self.array = name in ARRAY if array is None else array
        if self.array:
            self.log2 = bool(np.nanmax(self.X) > 100)
            self.T = np.log2(np.where(self.X > 0, self.X, np.nan)) if self.log2 else self.X
        else:
            self.log2 = True
            self.T = np.log2(self.X + 1)

    def cols(self, samples):
        idx = {s: i for i, s in enumerate(self.sample)}
        return np.array([idx[s] for s in samples])

    def values(self, feature, samples, scale="analysis"):
        x = self.T if scale == "analysis" else self.X
        return x[self.row[feature], self.cols(samples)]

    def gpr153_feature(self, samples):
        if not self.array:
            return "GPR153"
        means = {p: np.nanmean(self.values(p, samples)) for p in common.GPR153_PROBES if p in self.row}
        return max(means, key=means.get)

    def mgmt_feature(self):
        return common.MGMT_PROBE if self.array else "MGMT"


def corr(x, y, k=0):
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    n = len(x)
    s, sp = stats.spearmanr(x, y)
    p, pp = stats.pearsonr(x, y)
    lo, hi, _, _ = common.bw_ci(s, n, k)
    plo, phi, _, _ = common.bw_ci(p, n, k)
    return {"n": n, "spearman_rho": s, "spearman_p": sp, "spearman_ci_low": lo, "spearman_ci_high": hi,
            "pearson_r": p, "pearson_p": pp, "pearson_ci_low": plo, "pearson_ci_high": phi}


def duplicates(c, samples):
    keep = np.array([f not in EXCLUDE for f in c.feature])
    x = c.X[keep][:, c.cols(samples)]
    x = x[~np.isnan(x).any(1)]
    r = np.corrcoef(common.rank_rows(x.T))
    iu = np.triu_indices_from(r, 1)
    hit = r[iu] >= 0.98
    pairs = [(samples[i], samples[j], r[i, j]) for i, j in zip(iu[0][hit], iu[1][hit])]
    later = {max(a, b) for a, b, _ in pairs}
    return pairs, later


def platform_overlap(a, b):
    """Cross-platform overlap check between GSE13041 GPL96 (a) and GPL570 (b) arrays."""
    ta, tb = (set(c.clin["title"].str.replace("Tumor_fresh-frozen_", "", regex=False)) for c in (a, b))
    key = ["age", "sex", "os_days", "os_event"]
    ka = set(map(tuple, a.clin[key].astype(str).values))
    kb = set(map(tuple, b.clin[key].astype(str).values))
    common_p = [f for f in a.feature if f in b.row and f not in EXCLUDE]
    xa = a.X[[a.row[f] for f in common_p]]
    xb = b.X[[b.row[f] for f in common_p]]
    ok = ~(np.isnan(xa).any(1) | np.isnan(xb).any(1))
    ra, rb = common.rank_rows(xa[ok].T), common.rank_rows(xb[ok].T)
    ra = (ra - ra.mean(1, keepdims=True)) / np.linalg.norm(ra - ra.mean(1, keepdims=True), axis=1, keepdims=True)
    rb = (rb - rb.mean(1, keepdims=True)) / np.linalg.norm(rb - rb.mean(1, keepdims=True), axis=1, keepdims=True)
    r = ra @ rb.T
    rows = [{"check": "shared tumor ID (title)", "pairs_found": len(ta & tb), "detail": "; ".join(sorted(ta & tb))},
            {"check": "identical clinical record (age, sex, survival time, vital status)", "pairs_found": len(ka & kb),
             "detail": ""},
            {"check": f"whole-profile Spearman r >= 0.98 on {int(ok.sum())} common probe sets", "pairs_found": int((r >= 0.98).sum()),
             "detail": f"maximum cross-platform r {r.max():.3f}"}]
    return pd.DataFrame(rows)


def earlier_values():
    rows = []
    for f, line, lab in EARLIER:
        d = pd.read_csv(common.ORIG / f)
        r = d.iloc[line - 2]
        rows.append({"analysis": lab, "dataset_label": r["Dataset"], "n": int(r["N"]), "pearson_r": r["r"], "p": r["p"],
                     "source": f"{f} line {line}",
                     "status": "reported in same-session outputs; generating code unavailable; not reproduced"})
    return pd.DataFrame(rows)


def main():
    common.require_registration()
    out = common.outdir("external")
    co = {n: Cohort(n) for n in (*ARRAY, *RNA)}
    rows, combos = [], []
    for n, c in co.items():
        inc = c.clin.loc[c.clin["included"].astype(bool), "sample_id"].to_numpy() if "included" in c.clin else c.sample
        gf, mf = c.gpr153_feature(inc), c.mgmt_feature()
        rows.append({"cohort": ARRAY.get(n, RNA.get(n)), "gpr153_feature": gf, "mgmt_feature": mf,
                     "scale": "log2" if (c.array and c.log2) else ("as deposited" if c.array else "log2(FPKM + 1)"),
                     **corr(c.values(gf, inc), c.values(mf, inc))})
        if c.array:
            for p in common.GPR153_PROBES:
                if p in c.row:
                    combos.append({"cohort": ARRAY[n], "gpr153_feature": p, "mgmt_feature": mf,
                                   "primary": p == gf, **corr(c.values(p, inc), c.values(mf, inc))})
    ce = pd.DataFrame(rows)
    ce["spearman_q_bh"] = common.bh(ce["spearman_p"].values)
    ce.to_csv(out / "external_coexpression.csv", index=False)
    pd.DataFrame(combos).to_csv(out / "external_probe_combinations.csv", index=False)

    g = co["GSE13041"]
    cl = g.clin
    prim = cl.loc[cl["included"].astype(bool), "sample_id"].to_numpy()
    desc163 = cl.loc[cl["histology"] == "Primary GBM", "sample_id"].to_numpy()
    all191 = cl["sample_id"].to_numpy()
    pairs, later = duplicates(g, list(all191))
    after = np.array([s for s in all191 if s not in later])
    untreated = cl.loc[cl["included"].astype(bool) & ~cl["prior_treatment"].astype(bool), "sample_id"].to_numpy()
    gf = g.gpr153_feature(prim)
    runs = []
    for lab, s in (("(a) primary arrays", prim), ("(b) duplicate retained", desc163), ("(c) all GPL96 arrays", all191),
                   ("(d) all GPL96 arrays after duplicate exclusion", after),
                   ("(e) primary arrays without prior treatment", untreated)):
        runs.append({"run": lab, "gpr153_feature": gf, **corr(g.values(gf, s), g.values(common.MGMT_PROBE, s))})
    rd = pd.DataFrame(runs)
    rd["duplicate_pairs_in_191"] = "; ".join(f"{a}/{b} r={r:.3f}" for a, b, r in pairs)
    rd.to_csv(out / "gse13041_runs.csv", index=False)

    u1 = []
    for p, scale in itertools.product(common.GPR153_PROBES, ("deposited", "log2")):
        rec = {"gpr153_probe": p, "mgmt_probe": common.MGMT_PROBE, "scale": scale}
        for lab, s, target in (("163", desc163, 0.116), ("191", all191, 0.119)):
            x = g.values(p, s, "deposited")
            y = g.values(common.MGMT_PROBE, s, "deposited")
            if scale == "log2":
                x, y = np.log2(np.where(x > 0, x, np.nan)), np.log2(np.where(y > 0, y, np.nan))
            ok = ~(np.isnan(x) | np.isnan(y))
            r = stats.pearsonr(x[ok], y[ok])[0]
            rec.update({f"n_{lab}": int(ok.sum()), f"pearson_{lab}": r, f"reproduces_{lab}": round(r, 3) == target})
        rec["reproduces_both"] = rec["reproduces_163"] and rec["reproduces_191"]
        u1.append(rec)
    u1 = pd.DataFrame(u1)
    u1["ledger"] = "confirmed" if u1["reproduces_both"].any() else "not reproducible"
    u1.to_csv(out / "gse13041_u1_criterion.csv", index=False)

    g5 = Cohort("GSE13041_GPL570", array=True)
    ov = platform_overlap(g, g5)
    ov.to_csv(out / "gse13041_platform_overlap.csv", index=False)
    plat = []
    for lab, c in (("GPL96", g), ("GPL570", g5)):
        for sub_, s in (("all arrays", c.clin["sample_id"].to_numpy()),
                        ("Primary GBM", c.clin.loc[c.clin["histology"] == "Primary GBM", "sample_id"].to_numpy())):
            if lab == "GPL96" and sub_ == "Primary GBM":
                s = prim  # the 162 primary arrays after the duplicate rule
            f_ = c.gpr153_feature(s)
            plat.append({"platform": lab, "arrays": sub_, "gpr153_feature": f_,
                         "scale": "log2" if c.log2 else "as deposited", **corr(c.values(f_, s), c.values(common.MGMT_PROBE, s))})
    pl = pd.DataFrame(plat)
    for sub_ in ("all arrays", "Primary GBM"):
        e = pl[pl.arrays == sub_]
        rec = {"platform": "GPL96 + GPL570, random effects", "arrays": sub_}
        if ov["pairs_found"].sum() == 0:
            rho_ = e["spearman_rho"].values
            mr = common.meta_random(np.arctanh(rho_), (1 + rho_ ** 2 / 2) / (e["n"].values - 3), list(e["platform"]))
            rec.update(n=int(e["n"].sum()), spearman_rho=np.tanh(mr["mu"]), spearman_ci_low=np.tanh(mr["ci_low"]),
                       spearman_ci_high=np.tanh(mr["ci_high"]), spearman_p=mr["p"], tau2=mr["tau2"], I2=mr["I2"],
                       note="pooled: no cross-platform overlap found")
        else:
            rec["note"] = "not pooled: cross-platform overlap found (gse13041_platform_overlap.csv)"
        plat.append(rec)
    pd.DataFrame(plat).to_csv(out / "gse13041_platforms.csv", index=False)

    h = co["GSE4290"]
    hall = h.clin["sample_id"].to_numpy()
    hf = h.gpr153_feature(hall)
    hd = pd.DataFrame([{"analysis": "GSE4290, all histologies", "gpr153_feature": hf,
                        **corr(h.values(hf, hall), h.values(common.MGMT_PROBE, hall))}])
    hd["histology_counts"] = "; ".join(f"{k}: {v}" for k, v in h.clin["histology"].fillna("missing").value_counts().items())
    hd.to_csv(out / "gse4290_all_histologies.csv", index=False)
    earlier_values().to_csv(out / "earlier_same_session_values.csv", index=False)

    a, b = co["CGGA_693"], co["CGGA_325"]
    ia = a.clin.loc[a.clin["included"].astype(bool), "sample_id"].to_numpy()
    ib = b.clin.loc[b.clin["included"].astype(bool), "sample_id"].to_numpy()
    gx = np.concatenate([a.values("GPR153", ia), b.values("GPR153", ib)])
    mx = np.concatenate([a.values("MGMT", ia), b.values("MGMT", ib)])
    rel = np.r_[np.zeros(len(ia)), np.ones(len(ib))]
    pooled = corr(gx, mx)
    adj, adjp, _ = common.partial_spearman(gx, mx, rel[:, None])
    lo, hi, _, _ = common.bw_ci(adj[0], len(gx), 1)
    pd.DataFrame([{"analysis": "CGGA 693 + 325 pooled", **pooled},
                  {"analysis": "CGGA 693 + 325, adjusted for release (k = 1)", "n": len(gx), "spearman_rho": adj[0],
                   "spearman_p": adjp[0], "spearman_ci_low": lo, "spearman_ci_high": hi}]).to_csv(out / "cgga_pooled.csv", index=False)

    tc = pd.read_csv(common.RES / "screen/screen_sensitivity.csv").set_index("run").loc["primary: batch-adjusted"]
    y = np.arctanh(ce["spearman_rho"].values)
    v = (1 + ce["spearman_rho"].values ** 2 / 2) / (ce["n"].values - 3)
    meta = [{"analysis": "external cohorts", **common.meta_random(y, v, list(ce["cohort"]))}]
    yt = np.r_[y, np.arctanh(tc.gpr153_rho)]
    vt = np.r_[v, (1 + tc.gpr153_rho ** 2 / 2) / (tc.n - 3 - 1)]
    meta.append({"analysis": "external cohorts + TCGA", **common.meta_random(yt, vt, list(ce["cohort"]) + ["TCGA-GBM"])})
    md = pd.DataFrame(meta)
    for c in [c for c in md.columns if c in ("mu", "ci_low", "ci_high", "hk_ci_low", "hk_ci_high", "pi_low", "pi_high")]:
        md[f"{c}_rho"] = np.tanh(md[c])
    md["note"] = "tau^2 and I^2 are imprecise with this few cohorts"
    md.to_csv(out / "external_meta.csv", index=False)

    pw = []
    for r in ce.itertuples():
        for rho in common.PARAMS["power_rhos"] + [tc.gpr153_rho]:
            pw.append({"cohort": r.cohort, "n": r.n, "rho": rho, "power": common.corr_power(rho, r.n),
                       "rho_source": "TCGA estimate" if rho == tc.gpr153_rho else "fixed"})
    pd.DataFrame(pw).to_csv(out / "external_power.csv", index=False)

    dis, comp = [], []
    for n, c in co.items():
        inc = c.clin.loc[c.clin["included"].astype(bool), "sample_id"].to_numpy()
        T = c.T[:, c.cols(inc)]
        pct = (common.rank_rows(np.nan_to_num(T, nan=-np.inf).T).T / T.shape[0]) * 100
        for f in (c.gpr153_feature(inc), c.mgmt_feature()):
            v_ = pct[c.row[f]]
            dis.append({"cohort": ARRAY.get(n, RNA.get(n)), "feature": f, "n": len(inc),
                        "percentile_median": np.median(v_), "percentile_q25": np.percentile(v_, 25),
                        "percentile_q75": np.percentile(v_, 75), "frac_below_sample_median": float((v_ < 50).mean()),
                        "frac_fpkm_ge_1": float((c.X[c.row[f], c.cols(inc)] >= 1).mean()) if not c.array else np.nan})
        cl_ = c.clin[c.clin["included"].astype(bool)]
        for col in ("histology", "grade", "subtype_deposited", "idh_status"):
            if col in cl_:
                for val, k in cl_[col].fillna("missing").value_counts().items():
                    comp.append({"cohort": ARRAY.get(n, RNA.get(n)), "variable": col, "value": val, "n": int(k)})
    pd.DataFrame(dis).to_csv(out / "external_discordance.csv", index=False)
    pd.DataFrame(comp).to_csv(out / "external_composition.csv", index=False)
    print(ce[["cohort", "n", "spearman_rho", "spearman_q_bh", "pearson_r"]].to_string(index=False))
    print(u1[["gpr153_probe", "scale", "pearson_163", "pearson_191", "reproduces_both"]].to_string(index=False))


if __name__ == "__main__":
    main()
