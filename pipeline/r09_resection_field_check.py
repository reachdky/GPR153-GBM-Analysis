#!/usr/bin/env python3
"""
r09_resection_field_check.py
Look for an extent-of-resection field for the TCGA-GBM primary cohort in four
sources (PREREGISTRATION.md, Section 3.1): GDC diagnoses (residual_disease),
the TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR, Table S1), the PanCanAtlas
merged BCR clinical table, and the original BCR clinical XML. Field names are
matched on resect / residual / extent / surg / margin; for each matching field the
number of primary-cohort patients with a non-missing value is counted.

Outputs
  revision/results/cohort/resection_field_check.csv

Run from the repository root:
    python pipeline/r09_resection_field_check.py
"""

import re
from pathlib import Path

import pandas as pd

DL = Path("data/downloads")
RES = Path("revision/results/cohort")
PATTERN = re.compile(r"resect|residual|extent|surg|margin", re.I)
MISSING = {"", "[Not Available]", "[Not Applicable]", "[Unknown]", "[Not Evaluated]", "[Discrepancy]", "nan"}


def present(values):
    return int(sum(str(v).strip() not in MISSING for v in values))


def main():
    m = pd.read_csv("revision/manifests/tcga_samples.csv")
    prim = set(m.loc[m["in_screen_primary"], "patient"])
    rows = [{"source": "GDC diagnoses", "field": "residual_disease", "patients_in_source": len(prim),
             "patients_with_value": present(m.loc[m["in_screen_primary"], "gdc_residual_disease"].fillna(""))}]

    cdr = pd.read_excel(DL / "tcga_pub/TCGA-CDR-SupplementalTableS1.xlsx", sheet_name="TCGA-CDR")
    cdr = cdr[cdr["bcr_patient_barcode"].isin(prim)]
    for c in [c for c in cdr.columns if PATTERN.search(c)]:
        rows.append({"source": "TCGA-CDR Table S1", "field": c, "patients_in_source": len(cdr),
                     "patients_with_value": present(cdr[c].fillna(""))})

    pc = pd.read_csv(DL / "tcga_pub/clinical_PANCAN_patient_with_followup.tsv", sep="\t", low_memory=False,
                     encoding="latin1", dtype=str)
    pc = pc[pc["bcr_patient_barcode"].isin(prim)]
    for c in [c for c in pc.columns if PATTERN.search(c)]:
        rows.append({"source": "PanCanAtlas clinical_PANCAN_patient_with_followup", "field": c,
                     "patients_in_source": len(pc), "patients_with_value": present(pc[c].fillna(""))})

    xml_files = [p for p in (DL / "gdc_bcr_xml").glob("*/*clinical*.xml")]
    tag = re.compile(r"<([a-z_]+:[a-z_0-9]+)[ >]")
    hits, seen = {}, set()
    for p in xml_files:
        patient = re.search(r"(TCGA-\w\w-\w{4})", p.name).group(1)
        if patient not in prim:
            continue
        seen.add(patient)
        text = p.read_text(errors="replace")
        for t in set(tag.findall(text)):
            if PATTERN.search(t.split(":")[1]):
                filled = re.search(rf"<{t}[^>]*>\s*[^<\s][^<]*</{t}>", text)
                hits.setdefault(t, set())
                if filled:
                    hits[t].add(patient)
    for t in sorted(hits):
        rows.append({"source": "GDC BCR clinical XML", "field": t, "patients_in_source": len(seen),
                     "patients_with_value": len(hits[t])})
    if not hits:
        rows.append({"source": "GDC BCR clinical XML", "field": "(no matching element)", "patients_in_source": len(seen),
                     "patients_with_value": 0})
    out = pd.DataFrame(rows)
    out.to_csv(RES / "resection_field_check.csv", index=False)
    pd.set_option("display.width", 200)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
