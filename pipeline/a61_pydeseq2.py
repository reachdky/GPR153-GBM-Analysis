#!/usr/bin/env python3
"""
a61_pydeseq2.py
Cross-check of the primary DESeq2 comparison with PyDESeq2 (PREREGISTRATION.md, Section 12): same counts and
design (~ condition), Wald test, BH, apeglm-type shrinkage (DeseqStats.lfc_shrink). Reports agreement with R:
Spearman correlation of log2FC and of shrunken log2FC, maximum absolute log2FC difference, agreement of the
significance call, and GPR153's values from both.

Inputs
  revision/work/de_primary_counts.csv.gz, de_primary_samples.csv (a60); revision/results/de/deseq2_primary.csv (a60)
Outputs (revision/results/de/)
  pydeseq2_primary.csv, pydeseq2_vs_deseq2.csv

Run from the repository root:
    python pipeline/a61_pydeseq2.py
"""

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from scipy import stats

import common
from common import GPR153

P = common.PARAMS["de"]


def main():
    common.require_registration()
    out = common.outdir("de")
    cts = pd.read_csv(common.work_path("de_primary_counts.csv.gz"), index_col=0)
    meta = pd.read_csv(common.work_path("de_primary_samples.csv")).set_index("sample")
    meta["condition"] = pd.Categorical(meta["condition"], categories=["normal", "tumor"])
    keep = cts.sum(1) > 0
    dds = DeseqDataSet(counts=cts[keep].T, metadata=meta, design="~condition", refit_cooks=True,
                       inference=DefaultInference(n_cpus=8), quiet=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=["condition", "tumor", "normal"], inference=DefaultInference(n_cpus=8), quiet=True)
    ds.summary()
    res = ds.results_df.copy()
    ds.lfc_shrink(coeff="condition[T.tumor]")
    res["log2FoldChange_shrunk"] = ds.results_df["log2FoldChange"]
    res = res.rename_axis("gene_id").reset_index()
    res.to_csv(out / "pydeseq2_primary.csv", index=False)
    r = pd.read_csv(out / "deseq2_primary.csv").merge(res, on="gene_id", suffixes=("_r", "_py"))
    ok = r[["log2FoldChange_r", "log2FoldChange_py"]].notna().all(1)
    sig_r = (r.padj_r < P["padj"]) & (r.log2FoldChange_apeglm.abs() >= P["abs_lfc"])
    sig_p = (r.padj_py < P["padj"]) & (r.log2FoldChange_shrunk.abs() >= P["abs_lfc"])
    g = r.set_index("gene_id").loc[GPR153] if GPR153 in set(r.gene_id) else None
    pd.DataFrame([{"genes_compared": int(ok.sum()),
                   "spearman_log2fc": stats.spearmanr(r.loc[ok, "log2FoldChange_r"], r.loc[ok, "log2FoldChange_py"])[0],
                   "spearman_shrunken_log2fc": stats.spearmanr(r.loc[ok, "log2FoldChange_apeglm"], r.loc[ok, "log2FoldChange_shrunk"], nan_policy="omit")[0],
                   "max_abs_log2fc_difference": float(np.nanmax(np.abs(r.log2FoldChange_r - r.log2FoldChange_py))),
                   "significance_agreement": float((sig_r == sig_p).mean()),
                   "gpr153_log2fc_r": None if g is None else g.log2FoldChange_r,
                   "gpr153_log2fc_py": None if g is None else g.log2FoldChange_py,
                   "gpr153_padj_r": None if g is None else g.padj_r, "gpr153_padj_py": None if g is None else g.padj_py}]).to_csv(
        out / "pydeseq2_vs_deseq2.csv", index=False)
    print(pd.read_csv(out / "pydeseq2_vs_deseq2.csv").T)


if __name__ == "__main__":
    main()
