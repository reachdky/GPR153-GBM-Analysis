#!/usr/bin/env python3
"""
a19_methylation_inputs.py
MGMT-STP27 input values (PREREGISTRATION.md, Section 3.2): beta values of cg12434587 and cg12981137 from the
GDC SeSAMe level-3 beta files of TCGA-GBM primary tumors. Each file is md5-checked against GDC metadata. If a
patient has 450K data the 450K file is used, else the 27K file; with more than one file on that platform, the
lowest aliquot barcode. Locked until the analysis plan is registered.

Inputs
  data/downloads/gdc_methylation/files_meta.json and <file_id>/<file_name>       (r00)
Outputs
  revision/work/methylation_inputs.csv   patient, platform_used, aliquot, file_id, platforms_available,
                                         beta_cg12434587, beta_cg12981137

Run from the repository root:
    python pipeline/a19_methylation_inputs.py
"""

import hashlib
import json
from collections import defaultdict

import pandas as pd

import common

PROBES = ("cg12434587", "cg12981137")


def betas(path):
    v = {p: None for p in PROBES}
    with open(path) as f:
        for line in f:
            pid, _, val = line.rstrip("\n").partition("\t")
            if pid in v:
                v[pid] = float(val) if val not in ("", "NA", "NaN") else None
    return v


def main():
    common.require_registration()
    hits = json.load(open(common.DL / "gdc_methylation/files_meta.json"))["hits"]
    by_pt = defaultdict(list)
    for h in hits:
        s = h["cases"][0]["samples"][0]
        if s["sample_type"] != "Primary Tumor":
            continue
        al = s["portions"][0]["analytes"][0]["aliquots"][0]["submitter_id"]
        plat = h["platform"].replace("Illumina Human Methylation ", "HM")
        by_pt[h["cases"][0]["submitter_id"]].append((0 if plat == "HM450" else 1, al, plat, h))
    rows = []
    for pt, files in sorted(by_pt.items()):
        files.sort(key=lambda x: (x[0], x[1]))
        _, al, plat, h = files[0]
        path = common.DL / "gdc_methylation" / h["file_id"] / h["file_name"]
        if hashlib.md5(path.read_bytes()).hexdigest() != h["md5sum"]:
            raise SystemExit(f"md5 mismatch: {path}")
        b = betas(path)
        rows.append({"patient": pt, "platform_used": plat, "aliquot": al, "file_id": h["file_id"],
                     "platforms_available": "+".join(sorted({f[2] for f in files})),
                     "beta_cg12434587": b["cg12434587"], "beta_cg12981137": b["cg12981137"]})
    pd.DataFrame(rows).to_csv(common.work_path("methylation_inputs.csv"), index=False)
    print(f"methylation_inputs.csv: {len(rows)} patients")


if __name__ == "__main__":
    main()
