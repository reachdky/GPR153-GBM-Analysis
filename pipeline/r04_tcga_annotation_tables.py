#!/usr/bin/env python3
"""
r04_tcga_annotation_tables.py
Tabulate diagnostic and molecular annotations of the TCGA-GBM cohorts by batch
(PREREGISTRATION.md, Section 3.1). No expression values are read.

Input
  revision/manifests/tcga_samples.csv   (r01)
Outputs
  revision/results/cohort/tcga_annotation_tables.csv
  revision/results/cohort/tcga_cohort_overlap.csv

Run from the repository root:
    python pipeline/r04_tcga_annotation_tables.py
"""

from pathlib import Path

import pandas as pd

MAN = Path("revision/manifests/tcga_samples.csv")
RES = Path("revision/results/cohort")

COHORTS = [("primary (IDH-wildtype)", "in_screen_primary"),
           ("IDH not mutant (wildtype + unknown)", "in_screen_idh_not_mutant"),
           ("sensitivity (all primary)", "in_screen_sensitivity")]
VARIABLES = [("gdc_primary_diagnosis", "GDC primary_diagnosis"), ("gdc_morphology", "GDC morphology (ICD-O-3)"),
             ("histology_ceccarelli2016", "histology (Ceccarelli 2016)"),
             ("gdc_residual_disease", "GDC residual_disease (extent of resection)"),
             ("idh_status", "IDH status"), ("idh_source", "IDH source"),
             ("idh_source_disagreement", "IDH source disagreement"),
             ("mgmt_promoter", "MGMT promoter"), ("mgmt_promoter_source", "MGMT promoter source"),
             ("methylation_platforms", "methylation array (open beta file)"),
             ("egfr_amplified", "EGFR amplification (GISTIC = 2)"), ("egfr_gistic_source", "EGFR GISTIC source"),
             ("chr7gain_chr10loss_ceccarelli2016", "chr7 gain / chr10 loss (Ceccarelli 2016)"),
             ("tert_promoter_ceccarelli2016", "TERT promoter (Ceccarelli 2016)")]


def main():
    m = pd.read_csv(MAN)
    rows = []
    for cohort, flag in COHORTS:
        df = m[m[flag]]
        for col, label in VARIABLES:
            t = pd.crosstab(df[col].fillna("missing").astype(str), df["batch"])
            for lvl in t.index:
                rows.append({"cohort": cohort, "variable": label, "level": lvl,
                             **{f"n_{b}": int(t.loc[lvl, b]) for b in t.columns}, "n_total": int(t.loc[lvl].sum())})
    out = pd.DataFrame(rows).fillna(0)
    out.to_csv(RES / "tcga_annotation_tables.csv", index=False)

    prim = m[m["is_primary"]]
    both = set(m.loc[m["patient_in_both_batches"], "patient"])
    ov = [{"cohort": "all primary-tumor aliquots", "patients_with_BI_and_GDC2024_aliquots": len(both)}]
    for cohort, flag in COHORTS + [("survival primary", "in_survival_primary")]:
        ov.append({"cohort": cohort, "patients_with_BI_and_GDC2024_aliquots": len(set(m.loc[m[flag], "patient"]) & both)})
    pd.DataFrame(ov).to_csv(RES / "tcga_cohort_overlap.csv", index=False)
    assert prim.loc[prim["selected_aliquot"], "patient"].is_unique
    pd.set_option("display.width", 200)
    print(out[out.cohort.str.startswith("primary")].to_string(index=False))
    print(pd.DataFrame(ov).to_string(index=False))


if __name__ == "__main__":
    main()
