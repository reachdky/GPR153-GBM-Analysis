#!/usr/bin/env python3
"""
r10_download_raw_inputs.py
Retrieve the analysis inputs under data/raw/ that r00_download_external.py does not fetch, and verify each file
against its expected md5:
  tcga  the 391 GDC STAR - Counts files listed in revision/manifests/tcga_samples.csv (read by a10). Expected md5:
        the GDC md5 recorded in the manifest (column md5_gdc), which must also equal the md5sum the GDC files
        endpoint returns today.
  geo   GSE13041 and GSE4290 family SOFT files from the GEO FTP site (read by r02), and the expression.csv and
        clinical.csv tables built from them with GEOparse exactly as at initial submission (read by r02, r07, a29):
        expression = pivot_samples("VALUE"), clinical = phenotype_data. Expected md5: the files the final analysis
        run read, listed in revision/manifests/raw_inputs_expected_md5.csv.

A file already present at its destination is checked, not downloaded again. A downloaded file is written to a
temporary name and moved into place only if its md5 matches. On any mismatch the step stops, reports the file and
leaves existing files unchanged.

Outputs
  <dest>/tcga/expression/TCGA-GBM/Transcriptome_Profiling/Gene_Expression_Quantification/<file_id>/<file_name>
  <dest>/geo/GSE13041_family.soft.gz, <dest>/geo/GSE4290_family.soft.gz
  <dest>/geo/{GSE13041,GSE4290}/{expression,clinical}.csv
  revision/manifests/gdc_queries/tcga_star__files_meta.json   GDC query, data release, retrieval time, file records
  revision/manifests/raw_inputs_manifest.csv                  path, source URL, retrieval time, bytes, md5 expected
                                                              and observed, status

Run from the repository root:
    python pipeline/r10_download_raw_inputs.py [--dest DIR] [tcga] [geo]
--dest defaults to data/raw (the paths the pipeline reads). Steps default to both.
"""

import argparse
import csv
import hashlib
import json
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

MAN = Path("revision/manifests")
EXPECTED = MAN / "raw_inputs_expected_md5.csv"
REPORT = MAN / "raw_inputs_manifest.csv"
GDC = "https://api.gdc.cancer.gov"
STAR_SUB = Path("tcga/expression/TCGA-GBM/Transcriptome_Profiling/Gene_Expression_Quantification")
GEO_SOFT = "https://ftp.ncbi.nlm.nih.gov/geo/series/{p}/{g}/soft/{g}_family.soft.gz"
GEO_SERIES = {"GSE13041": "GSE13nnn", "GSE4290": "GSE4nnn"}
FIELDS = ["path", "source_url", "retrieved_utc", "bytes", "md5_expected", "md5_expected_source", "md5_observed",
          "status"]


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def fetch(url, path, method="GET", payload=None, retries=3):
    for attempt in range(1, retries + 1):
        try:
            with requests.request(method, url, json=payload, stream=True, timeout=300) as r:
                r.raise_for_status()
                with open(path, "wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
            return
        except requests.RequestException as e:
            print(f"  attempt {attempt} failed: {e}")
            time.sleep(10 * attempt)
    sys.exit(f"Download failed: {url}")


class Report:
    def __init__(self):
        self.rows = {}
        if REPORT.exists():
            with open(REPORT) as f:
                self.rows = {r["path"]: r for r in csv.DictReader(f)}

    def add(self, path, url, expected, source, observed, status, retrieved=""):
        self.rows[str(path)] = {"path": str(path), "source_url": url, "retrieved_utc": retrieved,
                                "bytes": path.stat().st_size if path.exists() else "", "md5_expected": expected,
                                "md5_expected_source": source, "md5_observed": observed, "status": status}
        with open(REPORT, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS)
            w.writeheader()
            for k in sorted(self.rows):
                w.writerow(self.rows[k])


def place(rep, url, dest, expected, source, method="GET", payload=None):
    """Check dest if present; otherwise download to a temporary file and move it into place on an md5 match.
    Returns True on a match."""
    if dest.exists():
        obs = md5(dest)
        ok = obs == expected
        rep.add(dest, url, expected, source, obs, "present, md5 match" if ok else "present, MD5 MISMATCH",
                rep.rows.get(str(dest), {}).get("retrieved_utc", ""))
        return ok
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    retrieved = now()
    fetch(url, tmp, method, payload)
    obs = md5(tmp)
    if obs != expected:
        tmp.unlink()
        rep.add(dest, url, expected, source, obs, "downloaded, MD5 MISMATCH, not kept", retrieved)
        return False
    tmp.rename(dest)
    rep.add(dest, url, expected, source, obs, "downloaded, md5 match", retrieved)
    return True


def step_tcga(dest, rep):
    m = pd.read_csv(MAN / "tcga_samples.csv", usecols=["file_id", "file_name", "md5_gdc"])
    query = {"filters": {"op": "and", "content": [
                 {"op": "in", "content": {"field": "file_id", "value": sorted(m["file_id"])}},
                 {"op": "=", "content": {"field": "access", "value": "open"}}]},
             "fields": "file_id,file_name,md5sum,file_size,access,state,data_type,analysis.workflow_type,"
                       "cases.submitter_id",
             "format": "JSON", "size": len(m)}
    retrieved = now()
    hits = requests.post(GDC + "/files", json=query, timeout=300).json()["data"]["hits"]
    status = requests.get(GDC + "/status", timeout=60).json()
    (MAN / "gdc_queries").mkdir(parents=True, exist_ok=True)
    with open(MAN / "gdc_queries" / "tcga_star__files_meta.json", "w") as f:
        json.dump({"endpoint": GDC + "/files", "query": query, "retrieved_utc": retrieved, "gdc_status": status,
                   "hits": hits}, f, indent=1)
    gdc = {h["file_id"]: h for h in hits}
    bad = []
    for r in m.itertuples():
        h = gdc.get(r.file_id)
        if h is None or h["md5sum"] != r.md5_gdc or h["file_name"] != r.file_name:
            bad.append(f"{r.file_id}: not returned as an open file with the recorded name and md5")
    if bad:
        sys.exit("tcga: GDC records differ from the manifest, nothing downloaded\n  " + "\n  ".join(bad))
    print(f"tcga: {len(hits)} open files, GDC {status.get('data_release', '')}")
    for i, r in enumerate(m.itertuples(), 1):
        path = dest / STAR_SUB / r.file_id / r.file_name
        if not place(rep, f"{GDC}/data/{r.file_id}", path, r.md5_gdc, "revision/manifests/tcga_samples.csv md5_gdc"):
            sys.exit(f"tcga: md5 mismatch for {path}; see {REPORT}")
        if i % 50 == 0:
            print(f"  {i}/{len(m)}")
    print(f"tcga: {len(m)} files, all md5 match")


def step_geo(dest, rep):
    import GEOparse

    exp = pd.read_csv(EXPECTED).set_index("path")["md5"].to_dict()
    failed = []
    for g, p in GEO_SERIES.items():
        url = GEO_SOFT.format(p=p, g=g)
        soft = dest / "geo" / f"{g}_family.soft.gz"
        key = f"geo/{g}_family.soft.gz"
        if not place(rep, url, soft, exp[key], str(EXPECTED)):
            failed.append(str(soft))
            continue
        out = dest / "geo" / g
        targets = {"expression.csv": None, "clinical.csv": None}
        if all((out / n).exists() for n in targets):
            for n in targets:
                obs = md5(out / n)
                ok = obs == exp[f"geo/{g}/{n}"]
                rep.add(out / n, f"built from {soft}", exp[f"geo/{g}/{n}"], str(EXPECTED), obs,
                        "present, md5 match" if ok else "present, MD5 MISMATCH")
                if not ok:
                    failed.append(str(out / n))
            continue
        tmpdir = dest / "geo" / f".{g}.build"
        shutil.rmtree(tmpdir, ignore_errors=True)
        tmpdir.mkdir(parents=True)
        gse = GEOparse.get_GEO(filepath=str(soft), silent=True)
        gse.pivot_samples("VALUE").to_csv(tmpdir / "expression.csv")
        gse.phenotype_data.to_csv(tmpdir / "clinical.csv")
        built = now()
        obs = {n: md5(tmpdir / n) for n in targets}
        if any(obs[n] != exp[f"geo/{g}/{n}"] for n in targets):
            for n in targets:
                rep.add(out / n, f"built from {soft}", exp[f"geo/{g}/{n}"], str(EXPECTED), obs[n],
                        "built, md5 match" if obs[n] == exp[f"geo/{g}/{n}"] else "built, MD5 MISMATCH, not kept",
                        built)
            failed.append(str(out))
            shutil.rmtree(tmpdir)
            continue
        out.mkdir(parents=True, exist_ok=True)
        for n in targets:
            (tmpdir / n).rename(out / n)
            rep.add(out / n, f"built from {soft}", exp[f"geo/{g}/{n}"], str(EXPECTED), obs[n], "built, md5 match",
                    built)
        tmpdir.rmdir()
    if failed:
        sys.exit("geo: md5 mismatch, files not kept:\n  " + "\n  ".join(failed) + f"\nsee {REPORT}")
    print("geo: all md5 match")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", default="data/raw")
    ap.add_argument("steps", nargs="*", default=["tcga", "geo"])
    a = ap.parse_args()
    dest, rep = Path(a.dest), Report()
    for s in a.steps:
        {"tcga": step_tcga, "geo": step_geo}[s](dest, rep)


if __name__ == "__main__":
    main()
