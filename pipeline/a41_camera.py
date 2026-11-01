#!/usr/bin/env python3
"""
a41_camera.py
Runs CAMERA (PREREGISTRATION.md, Section 10; supplement only) in R (pipeline/a41_camera.R).

Inputs   revision/work/camera_expression.csv.gz, camera_covariates.csv (a40); data/downloads/genesets/
Outputs  revision/results/enrichment/camera_<library>_igc_<estimated|0.01>.csv, camera_targeted_sets.csv

Run from the repository root:
    python pipeline/a41_camera.py
"""

import common


def main():
    common.require_registration()
    common.run_r("pipeline/a41_camera.R", common.work_path("camera_expression.csv.gz"),
                 common.work_path("camera_covariates.csv"), common.DL / "genesets", common.outdir("enrichment"))


if __name__ == "__main__":
    main()
