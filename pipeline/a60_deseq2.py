#!/usr/bin/env python3
"""
a60_deseq2.py
Differential expression (PREREGISTRATION.md, Section 12; Tier 3, supplement only), R DESeq2 with apeglm
shrinkage (pipeline/a60_deseq2.R). No ComBat-seq: tissue source and batch are fully confounded.
  primary       TCGA primary-cohort tumors in Toil vs GTEx brain cortex, one sample per donor
  sensitivity   the same tumors vs all GTEx brain cortex and frontal cortex (BA9) samples
  tcga_normals  all primary-cohort tumors vs the 5 TCGA solid-tissue normals (GDC STAR unstranded counts)
Thresholds: BH-adjusted p < 0.05 and |apeglm log2FC| >= 1. Reports GPR153's rank by padj and by shrunken log2FC,
the numbers of tumors and normals, and gene-ID handling (v23 vs v36 overlap on version-stripped IDs).

Inputs
  revision/work/de_toil.npz, de_toil_samples.csv (a59); revision/work/tcga_expr.npz (a10)
Outputs (revision/results/de/)
  deseq2_<comparison>.csv, de_summary.csv, de_gene_id_handling.csv

Run from the repository root:
    python pipeline/a60_deseq2.py
"""

import numpy as np
import pandas as pd

import common
from common import GPR153

P = common.PARAMS["de"]


def write_inputs(name, counts, gene_id, samples, condition):
    pd.DataFrame(counts, index=gene_id, columns=samples).rename_axis("gene_id").to_csv(
        common.work_path(f"de_{name}_counts.csv.gz"))
    pd.DataFrame({"sample": samples, "condition": condition}).to_csv(common.work_path(f"de_{name}_samples.csv"), index=False)


def summarize(name, r, n_t, n_n, normals):
    r["significant"] = (r["padj"] < P["padj"]) & (r["log2FoldChange_apeglm"].abs() >= P["abs_lfc"])
    r["rank_padj"] = r["padj"].rank(method="min")
    r["rank_shrunken_lfc"] = (-r["log2FoldChange_apeglm"]).rank(method="min")
    g = r.set_index("gene_id").loc[GPR153] if GPR153 in set(r.gene_id) else None
    return {"comparison": name, "tumors": n_t, "normals": n_n, "normal_tissue": normals, "genes_tested": int(r.padj.notna().sum()),
            "genes": len(r), "significant": int(r.significant.sum()),
            "gpr153_log2fc_apeglm": None if g is None else g.log2FoldChange_apeglm,
            "gpr153_padj": None if g is None else g.padj, "gpr153_rank_padj": None if g is None else g.rank_padj,
            "gpr153_rank_shrunken_lfc": None if g is None else g.rank_shrunken_lfc,
            "most_significant_up_gene": r[r.log2FoldChange_apeglm > 0].nsmallest(1, "padj").gene_id.iloc[0]}


def main():
    common.require_registration()
    out = common.outdir("de")
    z = common.load_matrix("de_toil.npz")
    s = pd.read_csv(common.work_path("de_toil_samples.csv"))
    assert list(s["sample"]) == list(z["sample"])
    t = common.TCGA()
    comps = []
    for name, mask in (("primary", (s.group == "tumor") | s.normal_one_per_donor), ("sensitivity", pd.Series(True, index=s.index))):
        write_inputs(name, z["counts"][:, mask.values], z["gene_id"], s.loc[mask, "sample"].values, s.loc[mask, "group"].values)
        comps.append((name, int((s.group[mask] == "tumor").sum()), int((s.group[mask] == "normal").sum()),
                      "; ".join(f"{k} {v}" for k, v in s.loc[mask & (s.group == "normal"), "tissue"].value_counts().items())))
    prim = common.cohort_aliquots(t.m, "in_screen_primary")
    nor = t.m[t.m["sample_type"] == "Solid Tissue Normal"]
    al = np.r_[prim["aliquot"].values, nor["aliquot"].values]
    write_inputs("tcga_normals", t.counts[:, t.cols(al)], t.gene_id, al, ["tumor"] * len(prim) + ["normal"] * len(nor))
    comps.append(("tcga_normals", len(prim), len(nor), f"TCGA solid-tissue normal {len(nor)}"))
    summ = []
    for name, n_t, n_n, normals in comps:
        common.run_r("pipeline/a60_deseq2.R", common.work_path(f"de_{name}_counts.csv.gz"),
                     common.work_path(f"de_{name}_samples.csv"), out / f"deseq2_{name}.csv")
        summ.append(summarize(name, pd.read_csv(out / f"deseq2_{name}.csv"), n_t, n_n, normals))
    pd.DataFrame(summ).to_csv(out / "de_summary.csv", index=False)
    v23, v36 = set(z["gene_id"]), set(t.gene_id)
    pd.DataFrame([{"item": "Toil v23 genes (version-stripped, PAR_Y dropped)", "n": len(v23)},
                  {"item": "GDC v36 genes", "n": len(v36)}, {"item": "overlap (stripped Ensembl ID)", "n": len(v23 & v36)},
                  {"item": "GPR153 in v23", "n": int(GPR153 in v23)}]).to_csv(out / "de_gene_id_handling.csv", index=False)
    print(pd.DataFrame(summ).to_string(index=False))


if __name__ == "__main__":
    main()
