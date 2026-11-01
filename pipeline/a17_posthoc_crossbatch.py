#!/usr/bin/env python3
"""
a17_posthoc_crossbatch.py
POST HOC diagnostic (PREREGISTRATION.md, Section 5.7; not part of the decision rule): cross-batch Spearman of
GPR153 and of MGMT between the BI and the GDC-2024 aliquot in the 79 overlapping primary-cohort patients, and
their percentile among all genes in a11's cross-batch reproducibility table.

Inputs
  revision/work/tcga_expr.npz (a10), revision/manifests/tcga_samples.csv (r01),
  revision/results/screen/crossbatch_reproducibility.csv (a11)
Outputs
  revision/results/screen/posthoc_crossbatch_gpr153_mgmt.csv

Run from the repository root:
    python pipeline/a17_posthoc_crossbatch.py
"""

import pandas as pd
from scipy import stats

import common
from common import GPR153, MGMT


def main():
    common.require_registration()
    t = common.TCGA()
    m = t.m
    prim = common.cohort_aliquots(m, "in_screen_primary")
    sel = prim[prim["patient_in_both_batches"] & (prim["batch"] == "BI")].set_index("patient")
    other = common.second_batch_aliquots(m, sel.index)
    sel = sel.loc[other.index]
    ref = pd.read_csv(common.RES / "screen/crossbatch_reproducibility.csv")["crossbatch_spearman"]
    rows = []
    for gid, name in ((GPR153, "GPR153"), (MGMT, "MGMT")):
        r = stats.spearmanr(t.gene(gid, sel["aliquot"].values), t.gene(gid, other.values)).statistic
        rows.append({"label": "post hoc", "gene": name, "gene_id": gid, "n_patients": len(sel),
                     "crossbatch_spearman": r, "percentile_among_genes": float((ref < r).mean() * 100),
                     "genes_in_reference": len(ref)})
    pd.DataFrame(rows).to_csv(common.outdir("screen") / "posthoc_crossbatch_gpr153_mgmt.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
