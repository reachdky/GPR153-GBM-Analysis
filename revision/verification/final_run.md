# Final analysis run

## Identification

| Item | Value |
|---|---|
| Final run | Run 5, `./run_all.sh --analysis`, 2026-10-30 18:49:00 to 19:22:43 UTC (33 min 43 s), exit 0 |
| Log | `revision/logs/run_analysis_2026-10-30.log` (start time `revision/logs/run_analysis_start.txt`) |
| Code | Analysis scripts as in the revised commit (`revision/submission_package/REVISED_COMMIT.txt`); no analysis script changed after run 5 (`revision/DEVIATIONS.md`, entry 2026-11-01) |
| Inputs | Verified by md5 on 2026-11-01: `revision/manifests/raw_inputs_manifest.csv` (391 TCGA STAR files, 2 GEO SOFT files, 4 GEO tables: all match); r00 downloads: `revision/manifests/DOWNLOAD_LOG.csv` |
| Results | `revision/results/`, `revision/figures/`, `revision/tables/`, `revision/source_data/`, `revision/figure_script_map.csv`, `revision/methods_parameters.json`, `revision/methods_concordance.csv` |

## Tests (`revision/results/tests/test_results.csv`)

16 of 17 tests pass. T5 fails, as recorded before registration (`revision/DEVIATIONS.md`, 2026-10-26): lifelines
and R `survival::cox.zph` disagree on R `lung` (age p 0.7229 vs 0.6474; sex p 0.1180 vs 0.1063). Section 13.1
prescribes R `cox.zph` for this case, and every proportional-hazards test in the results uses it
(`source` column of `revision/results/survival/cox_ph_tests.csv`).

## Independent reproduction (`revision/results/independent/`)

`independent_check_summary.csv`: 21 analyses, 1,457 values compared with independent implementations (R
`survival`, `metafor`, separate Python code), 0 failures; per-value detail in `independent_check.csv`. T7 and T8
(screen and Model 3 re-run in R) pass.

## Deviations

All entries in `revision/DEVIATIONS.md` (2026-10-26 to 2026-11-01). Entries made after the first real run:
a15 and a29 implementation errors (2026-10-28, chain rerun), figure presentation (2026-10-29, rerun), x02 ORA
column (2026-10-29, rerun), T9 stop (2026-10-29, rerun), f12 map rows (2026-10-29, rerun), S9 Table
(2026-10-30, rerun = run 5), r10 download step and the two summary scripts (2026-11-01, no rerun).

## Run-to-run comparisons

| Comparison | File | Result |
|---|---|---|
| Run 1 vs run 2 (row level) | `revision/verification/run1_vs_run2_diff.csv` | Changes confined to the methylation-dependent Cox models (Models 4, 5, 7 and Models 1-3 on the Model 4 sample) in the IDH-not-mutant and sensitivity cohorts, CPTAC-3 model (b), the model (b) meta-analyses, and the tables, source data and figure built from them (S4 Table, Table 5, Fig 5). Run 1 tables, figures and source data were not archived; they were rebuilt from the archived run 1 results with the current figure and table code (`scripts/rebuild_outputs.py`) |
| Run 2 vs run 5, results (SHA-256) | `revision/verification/run_comparison_sha256.csv` | 108 of 110 files identical. Differ: `tests/test_results.csv` (run times and the T11 run time in seconds) and `claims_trace_v2.csv` (C19 revised value 0.648 to 1.000, from the ORA column correction logged 2026-10-29) |
| Run 2 vs run 5, figures, tables, source data | — | Not comparable: run 2 figures, tables and source data were not archived (run 2 figure files were overwritten before archiving). Equivalence of run 2 and run 5 is established for `results/` only |
| Run 2 vs run 5, `methods_concordance.csv` | `run_comparison_sha256.csv` | Differs: the C19 revised-value row (as above) and, in 7 rows, the last printed digit of the computed relative difference |
| Run 4 vs run 5 (SHA-256) | `run_comparison_sha256.csv` | Results 109 of 110 identical (`test_results.csv` times); figures 15 of 15, source data 15 of 15, tables 14 of 14 identical; S9 Table new in run 5; `figure_script_map.csv` adds the S9 Table row; `methods_concordance.csv` differs only in the last printed digit of the relative difference in 7 rows |

## Consistency of the final outputs (`revision/verification/final_outputs_regenerated_sha256.csv`)

Figures, tables, source data, the methods record and the claims ledger were regenerated on 2026-11-01 from the
final `revision/results/` in a scratch directory (`scripts/regen_check.py`: f10, f11, f12, a90, x02) and compared
by SHA-256: all 14 figures, 15 tables (incl. S9), 14 source-data files, `methods_parameters.json`,
`methods_concordance.csv` and `claims_trace_v2.csv` are identical; `figure_script_map.csv` is identical after
the scratch path is replaced by `revision/`. The tables, figures, claims ledger, methods record and concordance
therefore all reflect the same final results.

`revision/figures/FigS1_cohort_flow.png` and `revision/source_data/FigS1_cohort_flow.csv` are written by
`pipeline/r01_tcga_cohort.py` (cohort construction), not by the figure scripts; they are not a numbered
submission figure and are not in `figure_script_map.csv` (the cohort flow is Fig 1).

## Unresolved

| Item | Note |
|---|---|
| Last-digit differences in `methods_concordance.csv` | The relative-difference column is recomputed by a90 from identical inputs; 7 of 315 rows differed in the 16th-17th significant digit between runs 4 and 5 (not in the 2026-11-01 regeneration). Pass/fail is unaffected (largest relative difference 4.4e-11 against a tolerance of 1e-6) |
| Run 2 figures, tables, source data | Not archived; see above |
