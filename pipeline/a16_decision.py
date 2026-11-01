#!/usr/bin/env python3
"""
a16_decision.py
Decision rule and outcome branch (PREREGISTRATION.md, Sections 5.5 and 5.6), from the outputs of a11, a14 and
a15. All four conditions are evaluated and reported. Branch, in order: (d) GPR153 not positive with q < 0.05 in
the primary analysis; else (b) condition 1 fails; else (c) any of conditions 2-4 fails; else (a).
Also reports whether GPR153's rank-1 status differs between the batch-adjusted and the pooled screen (if so, the
abstract states both).

Inputs (revision/results/)
  screen/screen_sensitivity.csv, composition/composition_adjusted.csv, methylation/condition4.csv
Outputs (revision/results/decision/)
  decision_conditions.csv, decision.csv

Run from the repository root:
    python pipeline/a16_decision.py
"""

import pandas as pd

import common


def main():
    common.require_registration()
    out = common.outdir("decision")
    s = pd.read_csv(common.RES / "screen/screen_sensitivity.csv").set_index("run")
    comp = pd.read_csv(common.RES / "composition/composition_adjusted.csv")
    c4 = pd.read_csv(common.RES / "methylation/condition4.csv").iloc[0]
    pr, po = s.loc["primary: batch-adjusted"], s.loc["pooled (unadjusted)"]
    bi, gdc = s.loc["BI only (sign check)"], s.loc["GDC2024 only (sign check)"]
    positive = bool(pr.gpr153_rho > 0 and pr.gpr153_q < 0.05)
    conds = [
        {"condition": 1, "description": "GPR153 rank 1 with q < 0.05 (batch-adjusted)",
         "value": f"rank {pr.gpr153_rank:.0f}, rho {pr.gpr153_rho:.4f}, q {pr.gpr153_q:.3g}",
         "met": bool(pr.gpr153_rank == 1 and pr.gpr153_q < 0.05)},
        {"condition": 2, "description": "positive with q < 0.05 after each composition adjustment",
         "value": "; ".join(f"{r.adjustment}: rho {r.gpr153_rho:.4f}, q {r.gpr153_q_adjusted_screen:.3g}"
                            for r in comp.itertuples()),
         "met": bool(comp["condition2_pass"].all())},
        {"condition": 3, "description": "positive in BI only and in GDC-2024 only (sign check)",
         "value": f"BI rho {bi.gpr153_rho:.4f} [{bi.ci_low:.3f}, {bi.ci_high:.3f}]; "
                  f"GDC-2024 rho {gdc.gpr153_rho:.4f} [{gdc.ci_low:.3f}, {gdc.ci_high:.3f}]",
         "met": bool(bi.gpr153_rho > 0 and gdc.gpr153_rho > 0)},
        {"condition": 4, "description": "methylation-adjusted partial rho > 0 with 95% CI excluding 0",
         "value": (f"{c4.measure}: n {c4.n}, rho {c4.rho:.4f} [{c4.ci_low:.3f}, {c4.ci_high:.3f}], k {c4.k}"
                   if c4.status == "evaluated" else f"not evaluable (n {c4.n})"),
         "met": bool(c4.condition4_pass)},
    ]
    cd = pd.DataFrame(conds)
    cd.to_csv(out / "decision_conditions.csv", index=False)
    met = dict(zip(cd.condition, cd.met))
    if not positive:
        branch = "d"
    elif not met[1]:
        branch = "b"
    elif not (met[2] and met[3] and met[4]):
        branch = "c"
    else:
        branch = "a"
    row = {"branch": branch, "gpr153_positive_q05": positive, "all_conditions_met": all(met.values()),
           "failing_conditions": ";".join(str(c) for c, v in met.items() if not v),
           "adjusted_rho": pr.gpr153_rho, "adjusted_q": pr.gpr153_q, "adjusted_rank": pr.gpr153_rank,
           "pooled_rho": po.gpr153_rho, "pooled_q": po.gpr153_q, "pooled_rank": po.gpr153_rank,
           "rank1_status_differs": bool((pr.gpr153_rank == 1) != (po.gpr153_rank == 1)),
           "genes_tested": pr.genes_tested, "n": pr.n,
           "rank1_gene": pr.rank1_gene, "rank1_rho": pr.rank1_rho, "rho_gap_to_rank1": pr.rank1_rho - pr.gpr153_rho}
    pd.DataFrame([row]).to_csv(out / "decision.csv", index=False)
    print(cd.to_string(index=False))
    print(f"branch ({branch})")


if __name__ == "__main__":
    main()
