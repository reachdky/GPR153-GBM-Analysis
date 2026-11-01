#!/usr/bin/env python3
"""
f12_tables.py
Main and supporting tables (PREREGISTRATION.md, Section 14) assembled from result CSVs. A table whose inputs are
absent (analysis dropped under Section 13.2) is skipped and listed.

  Table 1  primary-cohort characteristics by batch            Table 4  external co-expression, incl. platform and histology sensitivities
  Table 2  screen, decision conditions, branch, sensitivity   Table 5  survival audit and replication
  Table 3  Cox Models 1-7                                     Table 6  top enrichment results
  S1-S8 Table  manifests, gene-level screen, per-analysis n, survival audit, enrichment, MI Cox, DE
S9 Table (claims ledger) is written by pipeline/x02_claims_ledger.py.
Outputs: revision/tables/<item>.csv; a row per table in revision/figure_script_map.csv (the table CSV is its own
source data; the input result files are listed with it)

Run from the repository root:
    python pipeline/f12_tables.py
"""

import pandas as pd

import common

R = common.RES
OUT = common.BASE / "tables"
ME = "pipeline/f12_tables.py"


def read(*parts):
    p = R.joinpath(*parts)
    return pd.read_csv(p) if p.exists() else None


def per_analysis_n():
    rows = []
    s = read("screen", "screen_sensitivity.csv")
    if s is not None:
        rows += [{"analysis": f"screen: {r.run}", "n": r.n, "events": ""} for r in s.itertuples()]
    c = read("survival", "cox_models.csv")
    if c is not None:
        c = c[c.covariate == "gpr153"].drop_duplicates(["cohort", "model", "variant", "sample_set"])
        rows += [{"analysis": f"Cox {r.cohort} {r.model} {r.variant} ({r.sample_set})", "n": r.n, "events": r.events}
                 for r in c.itertuples()]
    for f, lab in (("external/external_coexpression.csv", "co-expression"), ("external/survival_replication_cox.csv", "survival")):
        e = read(f)
        if e is not None:
            rows += [{"analysis": f"external {lab}: {r.cohort}" + (f" {r.subset} {r.model}" if lab == "survival" else ""),
                      "n": r.n, "events": getattr(r, "events", "")} for r in e.itertuples()]
    c4 = read("methylation", "condition4.csv")
    if c4 is not None:
        rows.append({"analysis": f"condition 4 ({c4.measure.iloc[0]})", "n": c4.n.iloc[0], "events": ""})
    return pd.DataFrame(rows)


def main():
    common.require_registration()
    OUT.mkdir(parents=True, exist_ok=True)
    made, skipped = [], []

    def put(item, df, sources):
        if df is None or (isinstance(df, list) and any(x is None for x in df)):
            skipped.append(item)
            return
        if isinstance(df, list):
            df = pd.concat(df, ignore_index=True)
        path = OUT / f"{item.replace(' ', '_')}.csv"
        df.to_csv(path, index=False)
        common.record_figure(item, ME, path, [path, *sources])
        made.append(item)

    b = read("cohort", "tcga_batch_comparison.csv")
    put("Table 1", None if b is None else b[b.cohort.str.startswith("primary")], [R / "cohort/tcga_batch_comparison.csv"])
    parts = [read("decision", "decision.csv"), read("decision", "decision_conditions.csv"),
             read("screen", "screen_sensitivity.csv"), read("composition", "composition_adjusted.csv"),
             read("methylation", "condition4.csv"), read("methylation", "condition4_sensitivity.csv")]
    put("Table 2", None if any(p is None for p in parts) else
        [p.assign(section=s) for p, s in zip(parts, ["branch", "conditions", "screen and sensitivity (5.1)",
                                                   "composition (5.2)", "condition 4 (5.3)", "condition 4 sensitivity"])],
        [R / "decision/decision.csv", R / "decision/decision_conditions.csv", R / "screen/screen_sensitivity.csv",
         R / "composition/composition_adjusted.csv", R / "methylation/condition4.csv",
         R / "methylation/condition4_sensitivity.csv"])
    c = read("survival", "cox_models.csv")
    put("Table 3", None if c is None else c[(c.cohort == "primary") & (c.variant == "main")], [R / "survival/cox_models.csv"])
    ext_files = ("external_coexpression.csv", "gse13041_runs.csv", "gse13041_platforms.csv",
                                         "gse13041_platform_overlap.csv", "gse4290_all_histologies.csv",
                                         "cgga_pooled.csv", "external_meta.csv", "external_power.csv",
                 "earlier_same_session_values.csv")
    ext = [read("external", f) for f in ext_files]
    put("Table 4", None if any(e is None for e in ext) else
        [e.assign(part=p) for e, p in zip(ext, ["cohorts", "GSE13041 runs", "GSE13041 by platform (sensitivity)",
                                                "GSE13041 cross-platform overlap check",
                                                "GSE4290 all histologies (sensitivity)", "CGGA pooled", "meta-analysis",
                                                "power", "earlier same-session values (not reproduced)"])],
        [R / "external" / f for f in ext_files])
    sv = [read("cohort", "survival_audit.csv"), read("external", "survival_replication_cox.csv"),
          read("external", "survival_replication_meta.csv")]
    put("Table 5", None if any(e is None for e in sv) else
        [e.assign(part=p) for e, p in zip(sv, ["audit", "Cox per cohort", "meta-analysis"])],
        [R / "cohort/survival_audit.csv", R / "external/survival_replication_cox.csv",
         R / "external/survival_replication_meta.csv"])
    enr, enr_files = [], []
    for lib in common.PARAMS["ora"]["libraries"]:
        for d in ("positive", "negative"):
            e = read("enrichment", f"ora_{lib}_{d}.csv")
            if e is not None:
                enr.append(e.nsmallest(10, "p").assign(method="ORA", library=lib, direction=d))
                enr_files.append(R / f"enrichment/ora_{lib}_{d}.csv")
    for lib in common.PARAMS["gsea"]["libraries"]:
        e = read("enrichment", f"gsea_{lib}.csv")
        if e is not None:
            enr.append(e.assign(fdr=e["FDR q-val"].astype(float)).nsmallest(10, "fdr").assign(method="GSEA", library=lib))
            enr_files.append(R / f"enrichment/gsea_{lib}.csv")
    put("Table 6", enr or None, enr_files)

    for item, src in (("S1 Table", common.MAN / "tcga_samples.csv"), ("S2 Table", common.MAN / "external_cohort_samples.csv"),
                      ("S3 Table", R / "screen/screen_gene_level.csv")):
        put(item, pd.read_csv(src) if src.exists() else None, [src])
    put("S4 Table", per_analysis_n(), [R / f for f in ("screen/screen_sensitivity.csv", "survival/cox_models.csv",
                                                       "external/external_coexpression.csv",
                                                       "external/survival_replication_cox.csv", "methylation/condition4.csv")])
    a = [read("cohort", "survival_audit_strings.csv"), read("cohort", "survival_audit.csv")]
    put("S5 Table", None if any(x is None for x in a) else [a[0].assign(part="raw strings"), a[1].assign(part="counts")],
        [R / "cohort/survival_audit_strings.csv", R / "cohort/survival_audit.csv"])
    files = sorted((R / "enrichment").glob("ora_*_*tive.csv")) + sorted((R / "enrichment").glob("gsea_*.csv")) + \
        sorted((R / "enrichment").glob("camera_*_igc_*.csv")) + sorted((R / "enrichment").glob("mgmt_ora_*_positive.csv")) + \
        sorted((R / "enrichment").glob("mgmt_gsea_*.csv"))
    put("S6 Table", [pd.read_csv(f).assign(source_file=f.name) for f in files] or None, files)
    put("S7 Table", read("survival", "cox_multiple_imputation.csv"), [R / "survival/cox_multiple_imputation.csv"])
    put("S8 Table", read("de", "deseq2_primary.csv"), [R / "de/deseq2_primary.csv"])
    print("tables written:", ", ".join(made), "| skipped (no results):", ", ".join(skipped) or "none",
          "| S9 Table is written by x02_claims_ledger.py")


if __name__ == "__main__":
    main()
