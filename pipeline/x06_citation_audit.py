#!/usr/bin/env python3
"""
x06_citation_audit.py
Check each reference of the submitted manuscript and each reference recommended in the decision letter against
Crossref (DOI metadata) and PubMed (DOI or title search): first author, title, journal, year, volume and pages as
cited vs as registered. Not part of the analysis; writes under revision/literature/ only.

Inputs
  revision/literature/references_to_audit.csv   references as printed (transcribed from the manuscript and letter)
Outputs
  revision/literature/citation_audit.csv          cited vs registered fields, per-field match, PMID, retrieval time

Run from the repository root:
    python pipeline/x06_citation_audit.py
"""

import difflib
import re
import time
import unicodedata
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

LIT = Path("revision/literature")
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CR = "https://api.crossref.org/works/"
HDR = {"User-Agent": "gpr153-revision-citation-audit"}


def norm(s):
    s = unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode().lower()
    return re.sub(r"[^a-z0-9]+", " ", s).strip()


def first_page(p):
    return str(p).split("-")[0].strip().lower() if pd.notna(p) else ""


def crossref(doi):
    r = requests.get(CR + doi, headers=HDR, timeout=60)
    time.sleep(0.3)
    if r.status_code != 200:
        return None
    m = r.json()["message"]
    a = (m.get("author") or [{}])[0]
    year = (m.get("published-print") or m.get("published-online") or m.get("issued") or {}).get("date-parts", [[None]])[0][0]
    return {"reg_first_author": a.get("family", ""), "reg_title": (m.get("title") or [""])[0],
            "reg_journal": (m.get("container-title") or [""])[0], "reg_year": year, "reg_volume": m.get("volume", ""),
            "reg_pages": m.get("page", "") or m.get("article-number", "")}


def pubmed(term):
    r = requests.get(EU + "esearch.fcgi", params={"db": "pubmed", "term": term, "retmode": "json",
                                                  "tool": "gpr153_revision_search"}, timeout=60)
    time.sleep(0.4)
    ids = r.json()["esearchresult"]["idlist"]
    return ids


def pubmed_record(pmid):
    x = ET.fromstring(requests.get(EU + "efetch.fcgi", params={"db": "pubmed", "id": pmid, "retmode": "xml"},
                                   timeout=60).content)
    time.sleep(0.4)
    a = x.find(".//PubmedArticle")
    return {"reg_first_author": a.findtext(".//AuthorList/Author/LastName") or "",
            "reg_title": "".join(a.find(".//ArticleTitle").itertext()),
            "reg_journal": a.findtext(".//Journal/ISOAbbreviation") or "",
            "reg_year": a.findtext(".//JournalIssue/PubDate/Year") or "", "reg_volume": a.findtext(".//JournalIssue/Volume") or "",
            "reg_pages": a.findtext(".//Pagination/MedlinePgn") or ""}


def main():
    refs = pd.read_csv(LIT / "references_to_audit.csv", dtype=str)
    out = []
    for r in refs.itertuples():
        row = {"ref": r.ref, "source": r.source, "cited": f"{r.first_author}, {r.title}, {r.journal} {r.year};{r.volume}:{r.pages}",
               "cited_doi": r.doi if pd.notna(r.doi) else "", "checked_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
        reg, how = None, ""
        if pd.notna(r.doi):
            reg, how = crossref(r.doi), "Crossref DOI"
            ids = pubmed(f"{r.doi}[doi]")
        else:
            ids = pubmed(f"\"{r.title}\"[ti]")
            if ids:
                reg, how = pubmed_record(ids[0]), "PubMed title search"
        row["pmid"] = ";".join(ids)
        if pd.notna(r.url):
            resp = requests.get(r.url, headers=HDR, timeout=60, allow_redirects=True)
            row.update({"registered_by": "URL request", "url_status": resp.status_code, "url_final": resp.url})
            out.append(row)
            continue
        if reg is None:
            row.update({"registered_by": how, "status": "DOI not found" if pd.notna(r.doi) else "not found by title"})
            out.append(row)
            continue
        row.update(reg)
        row["registered_by"] = how
        row["title_similarity"] = round(difflib.SequenceMatcher(None, norm(r.title), norm(reg["reg_title"])).ratio(), 3)
        row["first_author_match"] = norm(r.first_author) == norm(reg["reg_first_author"])
        row["year_match"] = str(r.year) == str(reg["reg_year"])
        row["volume_match"] = (pd.isna(r.volume) and not reg["reg_volume"]) or str(r.volume) == str(reg["reg_volume"])
        row["first_page_match"] = first_page(r.pages) == first_page(reg["reg_pages"])
        bad = [k.replace("_match", "") for k in ("first_author_match", "year_match", "volume_match", "first_page_match")
               if not row[k]] + (["title"] if row["title_similarity"] < 0.9 else [])
        row["status"] = "matches" if not bad else "differs: " + ", ".join(bad)
        out.append(row)
    pd.DataFrame(out).to_csv(LIT / "citation_audit.csv", index=False)
    print(pd.DataFrame(out)[["ref", "status", "pmid"]].to_string(index=False))


if __name__ == "__main__":
    main()
