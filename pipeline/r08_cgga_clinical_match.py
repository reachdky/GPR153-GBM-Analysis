#!/usr/bin/env python3
"""
r08_cgga_clinical_match.py
Clinical-record match between CGGA mRNAseq_325 and mRNAseq_693 (PREREGISTRATION.md,
Section 7). A pair is a match if age, sex, grade, histology, OS time and censor
status are all equal (text fields stripped of surrounding whitespace). Partial
matches on age, sex and OS are listed for information. Only clinical tables are read.

Outputs
  revision/results/cohort/cgga_clinical_match.csv

Run from the repository root:
    python pipeline/r08_cgga_clinical_match.py
"""

import io
import zipfile
from pathlib import Path

import pandas as pd

DL = Path("data/downloads/cgga")
RES = Path("revision/results/cohort")
FULL = ["Age", "Gender", "Grade", "Histology", "OS", "Censor (alive=0; dead=1)"]
PARTIAL = ["Age", "Gender", "OS"]


def clinical(n):
    z = zipfile.ZipFile(DL / f"CGGA.mRNAseq_{n}_clinical.20200506.txt.zip")
    name = [x for x in z.namelist() if not x.startswith("__")][0]
    t = pd.read_csv(io.StringIO(z.read(name).decode().replace("\r\n", "\n").replace("\r", "\n")), sep="\t")
    t = t.apply(lambda s: s.str.strip() if s.dtype == object else s)
    t["Age"] = t["Age"].astype(float)
    t["selected"] = (t.PRS_type == "Primary") & (t.Grade == "WHO IV") & (t.IDH_mutation_status == "Wildtype")
    return t


def main():
    a, b = clinical("693"), clinical("325")
    rows = []
    for level, keys in [("all six fields", FULL), ("age, sex, OS only", PARTIAL)]:
        m = b.merge(a, on=keys, suffixes=("_325", "_693"))
        for _, r in m.iterrows():
            rows.append({"match_level": level, "id_325": r["CGGA_ID_325"], "id_693": r["CGGA_ID_693"],
                         "selected_325": r["selected_325"], "selected_693": r["selected_693"]})
    out = pd.DataFrame(rows, columns=["match_level", "id_325", "id_693", "selected_325", "selected_693"])
    out.to_csv(RES / "cgga_clinical_match.csv", index=False)
    print(out.to_string(index=False) if len(out) else "no matches")
    print("full matches:", int((out.match_level == "all six fields").sum()),
          "| among selected:", int(((out.match_level == "all six fields") & out.selected_325 & out.selected_693).sum()))


if __name__ == "__main__":
    main()
