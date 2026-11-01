#!/usr/bin/env python3
"""
x05_pubmed_search.py
PubMed searches requested in the decision letter (Section 10) and searches for GPR153 modulators, through the NCBI
E-utilities (esearch, esummary, efetch). Records the exact query, the run time (UTC), the count and every PMID, and
for each record the title, journal, year, DOI and abstract. Each GPR153 record is flagged for the terms it mentions
in its title or abstract. Not part of the analysis; writes under revision/literature/ only.

Outputs (revision/literature/)
  pubmed_searches.csv        query, purpose, run_utc, count, pmids
  pubmed_records.csv         pmid, title, journal, year, doi, flags, abstract (all records found)
  pubmed_raw/<n>.json        esearch responses

Run from the repository root:
    python pipeline/x05_pubmed_search.py
"""

import json
import re
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

OUT = Path("revision/literature")
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL = {"tool": "gpr153_revision_search"}
SEARCHES = [
    ("GPR153 AND glioblastoma", "decision letter Section 10"),
    ("GPR153 AND glioma", "decision letter Section 10"),
    ("GPR153 AND cancer", "decision letter Section 10"),
    ("GPR153 AND MGMT", "decision letter Section 10"),
    ("GPR153 AND (GLI OR GLI1)", "decision letter Section 10"),
    ("GPR153 AND WNT", "decision letter Section 10"),
    ("GPR153 AND temozolomide", "decision letter Section 10"),
    ("GPR153 AND (leukemia OR AML)", "decision letter Section 10"),
    ("orphan GPCR AND glioblastoma", "decision letter Section 10"),
    ("GPR153", "all GPR153 records (screening set)"),
    ("GPR153 AND (agonist OR antagonist OR ligand OR inhibitor OR modulator OR inverse agonist)", "GPR153 modulators"),
    ("GPR153 AND (small molecule OR compound OR drug OR screening)", "GPR153 modulators"),
    ("GPR153 AND (cAMP OR signaling OR signalling)", "GPR153 signaling"),
]
FLAGS = {"glioblastoma/glioma": r"glioblastoma|glioma|\bGBM\b", "cancer": r"cancer|tumou?r|carcinoma|leuk|AML",
         "MGMT": r"\bMGMT\b", "GLI/hedgehog": r"\bGLI1?\b|hedgehog", "WNT": r"\bWNT", "temozolomide": r"temozolomide",
         "modulator": r"agonist|antagonist|ligand|inhibitor|small.molecule", "cAMP": r"\bcAMP\b"}


def get(endpoint, **params):
    for attempt in range(4):
        r = requests.get(EU + endpoint, params={**TOOL, **params}, timeout=60)
        if r.status_code == 200:
            time.sleep(0.4)
            return r
        time.sleep(2 * (attempt + 1))
    r.raise_for_status()


def main():
    (OUT / "pubmed_raw").mkdir(parents=True, exist_ok=True)
    rows, pmids = [], set()
    for i, (q, purpose) in enumerate(SEARCHES, 1):
        run = datetime.now(timezone.utc).isoformat(timespec="seconds")
        j = get("esearch.fcgi", db="pubmed", term=q, retmax=10000, retmode="json").json()
        (OUT / "pubmed_raw" / f"{i:02d}.json").write_text(json.dumps({"query": q, "run_utc": run, "response": j}, indent=1))
        ids = j["esearchresult"]["idlist"]
        rows.append({"n": i, "query": q, "purpose": purpose, "run_utc": run, "count": int(j["esearchresult"]["count"]),
                     "query_translation": j["esearchresult"].get("querytranslation", ""), "pmids": ";".join(ids)})
        pmids |= set(ids)
    pd.DataFrame(rows).to_csv(OUT / "pubmed_searches.csv", index=False)
    recs = []
    ids = sorted(pmids, key=int)
    for k in range(0, len(ids), 200):
        x = ET.fromstring(get("efetch.fcgi", db="pubmed", id=",".join(ids[k:k + 200]), retmode="xml").content)
        for a in x.iter("PubmedArticle"):
            pmid = a.findtext(".//PMID")
            title = "".join(a.find(".//ArticleTitle").itertext()) if a.find(".//ArticleTitle") is not None else ""
            abst = " ".join("".join(t.itertext()) for t in a.findall(".//Abstract/AbstractText"))
            year = a.findtext(".//JournalIssue/PubDate/Year") or a.findtext(".//JournalIssue/PubDate/MedlineDate") or ""
            doi = next((e.text for e in a.findall(".//ArticleId") if e.get("IdType") == "doi"), "")
            text = f"{title} {abst}"
            recs.append({"pmid": pmid, "title": title, "journal": a.findtext(".//Journal/ISOAbbreviation") or "",
                         "year": year[:4], "doi": doi,
                         "mentions_gpr153": bool(re.search(r"GPR153", text, re.I)),
                         "flags": ";".join(f for f, pat in FLAGS.items() if re.search(pat, text, re.I)),
                         "found_by": ";".join(str(r["n"]) for r in rows if pmid in r["pmids"].split(";")),
                         "abstract": abst})
    pd.DataFrame(recs).to_csv(OUT / "pubmed_records.csv", index=False)
    print(pd.DataFrame(rows)[["n", "query", "count"]].to_string(index=False))
    print(f"{len(recs)} unique records")


if __name__ == "__main__":
    main()
