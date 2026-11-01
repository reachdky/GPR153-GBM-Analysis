#!/usr/bin/env python3
"""
a23_mice.py
Runs the multiple-imputation sensitivity analysis (PREREGISTRATION.md, Section 6.3) in R (pipeline/a23_mice.R).

Inputs   revision/work/mi_input.csv (a20)
Outputs  revision/results/survival/cox_multiple_imputation.csv

Run from the repository root:
    python pipeline/a23_mice.py
"""

import common


def main():
    common.require_registration()
    common.run_r("pipeline/a23_mice.R", common.work_path("mi_input.csv"),
                 common.outdir("survival") / "cox_multiple_imputation.csv")


if __name__ == "__main__":
    main()
