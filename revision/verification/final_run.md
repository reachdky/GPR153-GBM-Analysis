# Final analysis run

## Identification

| Item | Value |
|---|---|
| Final run | Run 7, `./run_all.sh --analysis`, 2026-11-01 21:51:15 to 22:25:03 UTC (33 min 48 s), exit 0 |
| Log | `revision/logs/run_analysis_2026-11-01_run7.log` (start and end: `revision/logs/run_analysis_start_run7.txt`) |
| Code | The revised commit (`revision/submission_package/REVISED_COMMIT.txt`); changes since run 5 are in `revision/DEVIATIONS.md` (2026-11-01) |
| Inputs | md5-verified 2026-11-01: 63 logged downloads, 2,020 GDC query files, 391 TCGA STAR files, GEO SOFT files and tables (`revision/verification/download_log_md5_check_2026-11-01.csv`, `gdc_query_files_md5_check_2026-11-01.csv`, `revision/manifests/raw_inputs_manifest.csv`) |
| Previous runs | Run 6 (2026-11-01, `revision/archive/run6_2026-11-01/`, log `revision/logs/run_analysis_2026-11-01.log`); run 5 (2026-10-30, `revision/archive/run5_2026-10-30/`) |

## Tests (`revision/results/tests/test_results.csv`)

16 of 17 pass. **T5 fails**: it compares lifelines' Schoenfeld test with R `survival::cox.zph` on R `lung`
(age p 0.7229 vs 0.6474; sex p 0.1180 vs 0.1063) and the two differ beyond the tolerance. The plan (Section 13.1)
prescribes R `cox.zph` as the source of every PH test in that case (`DEVIATIONS.md`, 2026-10-26). Verified: all
439 rows of `results/survival/cox_ph_tests.csv` have `source` = "R survival::cox.zph, transform km", and no
pipeline script calls lifelines' PH test.

## Independent reproduction (`revision/results/independent/`)

21 analyses, 1,457 values compared with independent implementations (R `survival`, `metafor`, separate Python
code), 0 failures. T7 and T8 (screen and Model 3 in R) pass.

## Run 7 vs run 6 (`run6_vs_run7_sha256.csv`)

Run 7 differs from run 6 only in the figure font (Arial, read from the four Arial 7.07 font files; DEVIATIONS
2026-11-01): the 14 figure TIFFs differ; all 110 result files except `tests/test_results.csv` (run times), all 15
tables, all source-data files, `figure_script_map.csv`, `methods_parameters.json` and `methods_concordance.csv`
are identical by SHA-256, and `numbers_for_manuscript.md` and `claims_not_reproduced.md` regenerate identically.
Every text element of every figure resolves to an Arial font file (checked in `plotstyle.save()`; the files used
are recorded in each TIFF).

## Run 6 vs run 5

| Output | Result | File |
|---|---|---|
| Results (110 files) | 108 identical by SHA-256. `survival/cox_ph_tests.csv`: the 406 run-5 rows are unchanged; 33 rows added (Models 1-3 on the Model 4 sample) with a `sample_set` column. `tests/test_results.csv`: run times only | `run5_vs_run6_sha256.csv` |
| Tables (15), `methods_parameters.json`, `methods_concordance.csv`, `claims_trace_v2.csv` | Identical | — |
| Figure source data | Values identical; S5 Fig source data gains the `sample_set` column | `figure_source_data_run5_vs_run6.csv` |
| Figures | Redrawn to the PLOS specifications from the same values; main figures renamed `Fig1.tif` ... | `DEVIATIONS.md` 2026-11-01 |

New PH results (not previously computed): on the Model 4 complete-case sample, Models 2 and 3 have global
`cox.zph` p 0.003-0.008 in all three cohorts, driven by age (p 0.0008-0.003) and sex (p 0.018-0.026); the GPR153
term has p 0.34-0.54. The same pattern was already present for Models 4 and 5 on that sample (global p 0.0014
and 0.0031 in run 5 and run 6).

## Earlier comparisons

| Comparison | File | Result |
|---|---|---|
| Run 1 vs run 2 (row level) | `run1_vs_run2_diff.csv` | Changes confined to the methylation-dependent Cox models in the IDH-not-mutant and sensitivity cohorts, CPTAC-3 model (b), the model (b) meta-analyses, and S4 Table, Table 5, Fig 5 (run 1 tables and figures rebuilt from archived results with `scripts/rebuild_outputs.py`) |
| Run 2, run 4 vs run 5 (SHA-256) | `run_comparison_sha256.csv` | Run 2 vs run 5 results: 108 of 110 identical (test times; C19 ORA-column correction). Run 2 figures, tables and source data were not archived, so equivalence is shown for results only |

## Consistency of the final outputs (`final_outputs_regenerated_sha256.csv`)

Figures, tables, source data, methods record and claims ledger regenerated on 2026-11-01 from the final
`revision/results/` in a scratch directory (`scripts/regen_check.py`) and compared by SHA-256: 14 figures, 15
tables, 14 source-data files, `methods_parameters.json`, `methods_concordance.csv` and `claims_trace_v2.csv`
identical; `figure_script_map.csv` identical after the scratch path is replaced by `revision/`.
`FigS1_cohort_flow.png/.csv` are written by `r01` (cohort construction), not by the figure scripts, and are not a
numbered figure.

## Unresolved

| Item | Note |
|---|---|
| NAAS | Not run |
| Run 2 figures, tables, source data | Not archived |
