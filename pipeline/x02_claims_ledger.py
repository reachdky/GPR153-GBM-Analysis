#!/usr/bin/env python3
"""
x02_claims_ledger.py
Revised value for each claim of the submitted manuscript (revision/results/claims_trace_v1.csv), read from the
revised-pipeline outputs. Reads result files only; computes no statistic.

Writes
  revision/results/claims_trace_v2.csv    v1 columns + revised_value, revised_source, revised_status
  revision/methods_concordance.csv        rows "Revised value: <claim>" replaced or appended (pipeline_value, and
                                          independent_value where Section 13.2 reproduced it)
  revision/tables/S9_Table.csv            claims ledger for the supplement: value reported, value in the original
                                          code or outputs, revised value, with neutral comparison labels (the free-text
                                          notes stay in claims_trace_v2.csv); recorded in figure_script_map.csv

Status: "consistent" (same direction and significance as claimed), "not supported" (claim does not hold in the
revised analysis), "changed" (estimate differs in size or definition), "no revised analysis" (outside the plan).

Run from the repository root, after run_all.sh --analysis:
    python pipeline/x02_claims_ledger.py
"""

import sys

import pandas as pd

sys.path.insert(0, "pipeline")
import common  # noqa: E402

R = common.RES
V1 = R / "claims_trace_v1.csv"
CONC = common.BASE / "methods_concordance.csv"
SECTION = {**{c: "5" for c in ("C01", "C02", "C03")}, **{c: "6.4" for c in ("C04", "C05", "C06")},
           **{c: "6.1, 6.2" for c in ("C07", "C08", "C09", "C10", "C11")}, "C12": "11", "C13": "11",
           **{c: "7" for c in ("C14", "C15", "C16", "C17", "C18")}, **{c: "10" for c in ("C19", "C20", "C21")},
           "C22": "12", "C23": "", "C24": "13", "C25": "2"}


# Neutral labels for S9 Table
ORIGINAL_LABEL = {"traced": "located in original code outputs",
                  "traced (output only)": "located in original outputs; generating code not available",
                  "traced (sample list inferred)": "located in original outputs; sample list inferred",
                  "inferred": "inferred from original code",
                  "not reproducible": "not located in original code or outputs"}
REVISED_LABEL = {"consistent": "same direction and significance as reported",
                 "changed": "revised value differs from the reported value",
                 "not supported": "reported result not reproduced in the revised analysis",
                 "no revised analysis": "not re-analysed (outside the analysis plan)"}


def rd(p):
    return pd.read_csv(R / p)


def f(x, d=3):
    return f"{x:.{d}f}" if abs(x) >= 10 ** -d else f"{x:.2e}"


def main():
    common.require_registration()
    v1 = pd.read_csv(V1)
    dec = rd("decision/decision.csv").iloc[0]
    sens = rd("screen/screen_sensitivity.csv").set_index("run")
    boot = rd("screen/bootstrap_rank_summary.csv").iloc[0]
    km = rd("survival/km_logrank.csv").set_index(["gene", "split"])
    kmg = rd("survival/km_groups.csv").set_index(["gene", "split", "group"])
    cm = rd("survival/cox_models.csv")
    ph = rd("survival/cox_ph_tests.csv")
    ph = ph[ph.sample_set == "complete cases"]
    ind = rd("independent/independent_check.csv")
    ext = rd("external/external_coexpression.csv").set_index("cohort")
    runs = rd("external/gse13041_runs.csv").set_index("run")
    g4290 = rd("external/gse4290_all_histologies.csv").iloc[0]
    meta = rd("external/external_meta.csv").set_index("analysis")
    cg = rd("external/cgga_pooled.csv").iloc[0]
    fc = rd("forest/forest_c_summary.csv")
    fd = rd("forest/forest_paired_difference.csv").set_index("analysis")
    fi = rd("forest/forest_permutation_importance.csv")
    de = rd("de/de_summary.csv").set_index("comparison")
    pca = rd("pca/pca_variance.csv").set_index("component")
    gs = rd("enrichment/gsea_h.all.v2026.2.Hs.symbols.csv").set_index("Term")
    gr = rd("enrichment/gsea_c2.cp.reactome.v2026.2.Hs.symbols.csv").set_index("Term")
    cam = rd("enrichment/camera_targeted_sets.csv")
    ora = {lib: rd(f"enrichment/ora_{lib}_positive.csv") for lib in ("Reactome_2022", "KEGG_2021_Human",
                                                                     "GO_Biological_Process_2023")}

    def hr(model, cov="gpr153", sample="complete cases"):
        r = cm[(cm.cohort == "primary") & (cm.model == model) & (cm.variant == "main") & (cm.sample_set == sample)
               & (cm.covariate == cov)].iloc[0]
        return r, f"HR {r.hr:.3f} [{r.hr_low:.3f}, {r.hr_high:.3f}] per {'SD' if cov == 'gpr153' else 'unit'}, p {f(r.p, 4)}, n {r.n:.0f}, events {r.events:.0f}"

    def indep(analysis, quantity):
        r = ind[(ind.analysis == analysis) & (ind.quantity == quantity)]
        return float(r.independent_value.iloc[0]) if len(r) else None

    m1, s1 = hr("Model 1")
    m2, s2 = hr("Model 2")
    m3, s3 = hr("Model 3")
    a2, sa2 = hr("Model 2", "age")
    m6, s6 = hr("Model 6")
    php = ph[(ph.cohort == "primary") & (ph.model == "Model 2") & (ph.variant == "main")].set_index("term").p
    ph3 = ph[(ph.cohort == "primary") & (ph.model == "Model 3") & (ph.variant == "main")].set_index("term").p
    rs = ora["Reactome_2022"]
    wnt = rs[rs.Term == "Signaling By WNT R-HSA-195721"]
    imp = fi[fi.analysis == "primary cohort"].sort_values("mean", ascending=False).reset_index(drop=True)
    gi = imp.index[imp.feature == "gpr153_log2tpm"][0]
    fcc = fc[(fc.analysis == "primary cohort") & (fc.set == "C")].iloc[0]
    n_ext = int(sum(ext.loc[c, "n"] for c in ext.index))
    camw = cam[(cam.term == "HALLMARK_WNT_BETA_CATENIN_SIGNALING") & (cam.inter_gene_cor == "estimated")].iloc[0]

    rev = {
        "C01": (f"batch-adjusted partial Spearman rho {f(dec.adjusted_rho)}, p {f(sens.loc['primary: batch-adjusted', 'gpr153_p'], 4)}, "
                f"q {f(dec.adjusted_q, 4)}, n {dec.n:.0f} patients (one aliquot each); pooled (unadjusted) rho "
                f"{f(dec.pooled_rho)}, q {f(dec.pooled_q)}",
                "results/decision/decision.csv; results/screen/screen_sensitivity.csv", "changed",
                ("screen, batch-adjusted: GPR153 (all-gene implementation)", "GPR153 rho"), dec.adjusted_rho),
        "C02": (f"rank {dec.adjusted_rank:.0f} of {dec.genes_tested:.0f} eligible genes (batch-adjusted), q {f(dec.adjusted_q, 4)}; "
                f"rank 1 {dec.rank1_gene} (rho {f(dec.rank1_rho)}; gap {f(dec.rho_gap_to_rank1)}); pooled rank "
                f"{dec.pooled_rank:.0f}; bootstrap: rank 1 in {boot.prop_rank1:.3f} of 1,000 resamples, median rank "
                f"{boot.rank_median:.0f}; outcome branch {dec.branch}",
                "results/decision/decision.csv; results/screen/bootstrap_rank_summary.csv", "not supported",
                ("screen, batch-adjusted: GPR153 (all-gene implementation)", "GPR153 rank"), dec.adjusted_rank),
        "C03": (f"FPKM: rho {f(sens.loc['FPKM', 'gpr153_rho'])}, q {f(sens.loc['FPKM', 'gpr153_q'])}, rank "
                f"{sens.loc['FPKM', 'gpr153_rank']:.0f}; single GPR153/MGMT outlier removed: rho "
                f"{f(sens.loc['GPR153 or MGMT outliers removed (single value)', 'gpr153_rho'])} "
                f"(MGMT |z| > 3: 0 samples removed)",
                "results/screen/screen_sensitivity.csv", "changed", None, sens.loc["FPKM", "gpr153_rho"]),
        "C04": (f"log-rank, pooled median split: chi-square {f(km.loc[('GPR153', 'pooled median'), 'chisq'])}, p "
                f"{f(km.loc[('GPR153', 'pooled median'), 'p'], 4)}; within-batch median: p "
                f"{f(km.loc[('GPR153', 'within-batch median'), 'p'], 4)}; median OS high "
                f"{kmg.loc[('GPR153', 'pooled median', 'high'), 'median_os_months']:.1f} vs low "
                f"{kmg.loc[('GPR153', 'pooled median', 'low'), 'median_os_months']:.1f} months; n 249",
                "results/survival/km_logrank.csv; km_groups.csv", "changed",
                ("log-rank", "GPR153, pooled median: p"), km.loc[("GPR153", "pooled median"), "p"]),
        "C05": (f"MGMT expression, pooled median split: p {f(km.loc[('MGMT', 'pooled median'), 'p'], 4)}; within-batch median: "
                f"p {f(km.loc[('MGMT', 'within-batch median'), 'p'], 4)}; n 249",
                "results/survival/km_logrank.csv", "changed",
                ("log-rank", "MGMT, pooled median: p"), km.loc[("MGMT", "pooled median"), "p"]),
        "C06": ("n 249 patients (207 deaths), one primary-tumor aliquot per patient", "results/survival/km_groups.csv",
                "changed", ("log-rank", "GPR153, pooled median: n"), 249.0),
        "C07": (f"Model 1 (GPR153 alone): {s1}", "results/survival/cox_models.csv", "changed", None, m1.hr),
        "C08": (f"Model 2 (GPR153 + age): {s2}; primary Model 3 (GPR153 + age + sex): {s3}",
                "results/survival/cox_models.csv", "consistent", ("Model 3", "gpr153 hr"), m3.hr),
        "C09": (f"age in Model 2: {sa2} (per year)", "results/survival/cox_models.csv", "changed", None, a2.hr),
        "C10": (f"R cox.zph (km): Model 2 GPR153 p {f(php['gpr153'])}, age p {f(php['age'])}, global p {f(php['GLOBAL'])}; "
                f"Model 3 GPR153 p {f(ph3['gpr153'])}, global p {f(ph3['GLOBAL'])}",
                "results/survival/cox_ph_tests.csv", "consistent", None, php["gpr153"]),
        "C11": (f"Model 6 (GPR153 + age + sex + MGMT expression): {s6}", "results/survival/cox_models.csv",
                "changed", None, m6.hr),
        "C12": (f"random survival forest, set C (with GPR153): C-index {f(fcc.c_mean)} (50 folds, SD {f(fcc.c_sd)}); "
                f"C minus B (without GPR153) {f(fd.loc['primary cohort', 'mean_diff_C_minus_B'], 4)}, p "
                f"{f(fd.loc['primary cohort', 'nb_p'])}",
                "results/forest/forest_c_summary.csv; forest_paired_difference.csv", "not supported", None,
                fcc.c_mean),
        "C13": (f"permutation importance (set C, primary cohort): GPR153 rank {gi + 1} of {len(imp)} "
                f"({f(imp.loc[gi, 'mean'], 4)}); rank 1 {imp.loc[0, 'feature']} ({f(imp.loc[0, 'mean'], 4)})",
                "results/forest/forest_permutation_importance.csv", "not supported", None, imp.loc[gi, "mean"]),
        "C14": (f"GSE13041 all GPL96 arrays: Spearman {f(runs.loc['(c) all GPL96 arrays', 'spearman_rho'])} "
                f"(p {f(runs.loc['(c) all GPL96 arrays', 'spearman_p'])}), Pearson "
                f"{f(runs.loc['(c) all GPL96 arrays', 'pearson_r'])}, n 191 (64942_at, log2)",
                "results/external/gse13041_runs.csv", "changed", None, runs.loc["(c) all GPL96 arrays", "spearman_rho"]),
        "C15": (f"GSE13041 GPL96 primary arrays: Spearman {f(ext.loc['GSE13041 GPL96 primary', 'spearman_rho'])} "
                f"[{f(ext.loc['GSE13041 GPL96 primary', 'spearman_ci_low'])}, {f(ext.loc['GSE13041 GPL96 primary', 'spearman_ci_high'])}], "
                f"p {f(ext.loc['GSE13041 GPL96 primary', 'spearman_p'])}, n 162 (163 with the duplicate retained: "
                f"{f(runs.loc['(b) duplicate retained', 'spearman_rho'])})",
                "results/external/external_coexpression.csv; gse13041_runs.csv", "changed", None,
                ext.loc["GSE13041 GPL96 primary", "spearman_rho"]),
        "C16": (f"GSE4290 GBM: Spearman {f(ext.loc['GSE4290 GBM', 'spearman_rho'])} [{f(ext.loc['GSE4290 GBM', 'spearman_ci_low'])}, "
                f"{f(ext.loc['GSE4290 GBM', 'spearman_ci_high'])}], Pearson {f(ext.loc['GSE4290 GBM', 'pearson_r'])}, n 77; all "
                f"histologies: Spearman {f(g4290.spearman_rho)}, n {g4290.n:.0f}",
                "results/external/external_coexpression.csv; gse4290_all_histologies.csv", "changed", None,
                ext.loc["GSE4290 GBM", "spearman_rho"]),
        "C17": (f"CGGA mRNAseq_693 GBM IDH-wildtype: Spearman {f(ext.loc['CGGA mRNAseq_693', 'spearman_rho'])}, n 109; "
                f"mRNAseq_325: {f(ext.loc['CGGA mRNAseq_325', 'spearman_rho'])}, n 74; pooled {f(cg.spearman_rho)}, n {cg.n:.0f}",
                "results/external/external_coexpression.csv; cgga_pooled.csv", "not supported", None,
                ext.loc["CGGA mRNAseq_693", "spearman_rho"]),
        "C18": (f"4 external cohorts, n {n_ext} (no double counting); random-effects rho "
                f"{f(meta.loc['external cohorts', 'mu_rho'])} [{f(meta.loc['external cohorts', 'ci_low_rho'])}, "
                f"{f(meta.loc['external cohorts', 'ci_high_rho'])}], p {f(meta.loc['external cohorts', 'p'])}, HKSJ p "
                f"{f(meta.loc['external cohorts', 'hk_p'])}, I2 {meta.loc['external cohorts', 'I2'] * 100:.0f}%",
                "results/external/external_meta.csv", "not supported",
                ("meta-analysis, co-expression: external cohorts", "p"), meta.loc["external cohorts", "p"]),
        "C19": (f"ORA, Reactome_2022, GPR153-positive list (200 genes, eligible-gene background): Signaling By WNT "
                f"adjusted p (BH over all tested sets) {f(wnt['q_bh'].iloc[0])}; sets with adjusted p < 0.05: "
                f"{int((rs['q_bh'] < 0.05).sum())}",
                "results/enrichment/ora_Reactome_2022_positive.csv", "not supported", None,
                wnt["q_bh"].iloc[0]),
        "C20": ("ORA, positive list: adjusted p < 0.05 in " + "; ".join(
                    f"{lib} {int((d['q_bh'] < 0.05).sum())} sets" for lib, d in ora.items()),
                "results/enrichment/ora_*_positive.csv", "not supported", None,
                float(sum(int((d["q_bh"] < 0.05).sum()) for d in ora.values()))),
        "C21": (f"GSEA preranked (Hallmark): WNT_BETA_CATENIN NES {f(gs.loc['HALLMARK_WNT_BETA_CATENIN_SIGNALING', 'NES'])}, "
                f"FDR {f(gs.loc['HALLMARK_WNT_BETA_CATENIN_SIGNALING', 'FDR q-val'])}, rank "
                f"{int(gs.NES.rank(ascending=False)['HALLMARK_WNT_BETA_CATENIN_SIGNALING'])} of {len(gs)} by NES; DNA_REPAIR "
                f"FDR {f(gs.loc['HALLMARK_DNA_REPAIR', 'FDR q-val'])}; CAMERA (estimated inter-gene correlation) WNT p "
                f"{f(camw.PValue)}, FDR {f(camw.FDR)}",
                "results/enrichment/gsea_h.all.*.csv; camera_targeted_sets.csv", "changed", None,
                gs.loc["HALLMARK_WNT_BETA_CATENIN_SIGNALING", "NES"]),
        "C22": (f"DESeq2 (apeglm), TCGA primary tumors vs GTEx cortex (Toil): GPR153 log2FC "
                f"{f(de.loc['primary', 'gpr153_log2fc_apeglm'])}, padj {f(de.loc['primary', 'gpr153_padj'])}, rank by padj "
                f"{de.loc['primary', 'gpr153_rank_padj']:.0f} of {de.loc['primary', 'genes_tested']:.0f}; vs TCGA normals (n 5): "
                f"log2FC {f(de.loc['tcga_normals', 'gpr153_log2fc_apeglm'])}",
                "results/de/de_summary.csv", "not supported", None, de.loc["primary", "gpr153_log2fc_apeglm"]),
        "C23": ("no revised analysis (KM by MGMT methylation status is not in the plan; KM is by GPR153 and MGMT "
                "expression, C04/C05)", "", "no revised analysis", None, None),
        "C24": (f"PCA of the 2,000 most variable genes, n 252: PC1 {pca.loc['PC1', 'variance_explained'] * 100:.1f}%, PC2 "
                f"{pca.loc['PC2', 'variance_explained'] * 100:.1f}%; no subtype labels",
                "results/pca/pca_variance.csv", "changed", None, pca.loc["PC1", "variance_explained"]),
        "C25": ("downloads 2026-10-24 to 2026-10-26 UTC (data/downloads/DOWNLOAD_LOG.csv); GDC Data Release 46.0",
                "data/downloads/DOWNLOAD_LOG.csv", "changed", None, None),
    }
    v2 = v1.copy()
    v2["revised_value"] = v2.claim_id.map(lambda c: rev[c][0])
    v2["revised_source"] = v2.claim_id.map(lambda c: rev[c][1])
    v2["revised_status"] = v2.claim_id.map(lambda c: rev[c][2])
    v2.to_csv(R / "claims_trace_v2.csv", index=False)
    s9 = pd.DataFrame({"claim": v2.claim_id, "location in submitted manuscript": v2.claim_location,
                       "statement": v2.claim_text, "reported value": v2.claimed_value,
                       "original code or output": v2.traced_value, "original source": v2.source_file_or_script,
                       "original code comparison": v2.trace_status.map(ORIGINAL_LABEL),
                       "revised value": v2.revised_value, "revised source": v2.revised_source,
                       "revised comparison": v2.revised_status.map(REVISED_LABEL),
                       "analysis plan section": v2.claim_id.map(SECTION)})
    assert s9["original code comparison"].notna().all() and s9["revised comparison"].notna().all()
    tab = common.BASE / "tables"
    tab.mkdir(parents=True, exist_ok=True)
    s9.to_csv(tab / "S9_Table.csv", index=False)
    common.record_figure("S9 Table", "pipeline/x02_claims_ledger.py", tab / "S9_Table.csv",
                         [tab / "S9_Table.csv", V1, R / "claims_trace_v2.csv"])

    conc = pd.read_csv(CONC)
    conc = conc[~conc.statement.str.startswith("Revised value: ")]
    rows = []
    for _, r in v1.iterrows():
        val, src, status, chk, pv = rev[r.claim_id]
        iv = indep(*chk) if chk else None
        rows.append({"statement": f"Revised value: {r.claim_id} {r.claim_text}", "plan_section": SECTION[r.claim_id],
                     "script": "pipeline/x02_claims_ledger.py", "output": "revision/results/claims_trace_v2.csv",
                     "output_exists": True, "pipeline_value": pv, "independent_value": iv,
                     "independent_implementation": f"a91: {chk[0]} / {chk[1]}" if chk else "",
                     "criterion": "relative difference" if chk else "", "tolerance": 1e-6 if chk else None,
                     "difference": abs(pv - iv) / max(abs(iv), 1e-300) if chk else None,
                     "pass": (abs(pv - iv) <= 1e-6 * max(abs(iv), 1e-300)) if chk else None})
    pd.concat([conc, pd.DataFrame(rows)], ignore_index=True).to_csv(CONC, index=False)
    print(v2[["claim_id", "revised_status"]].value_counts("revised_status").to_string())


if __name__ == "__main__":
    main()
