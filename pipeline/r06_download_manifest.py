#!/usr/bin/env python3
"""
r06_download_manifest.py
Copy the download log to revision/manifests/ and write the per-file retrieval
table used in PREREGISTRATION.md Section 2. Files fetched through a GDC query
(MAFs, methylation, CPTAC) are listed by their query metadata file, with the
file count and the GDC data release.

Outputs
  revision/manifests/DOWNLOAD_LOG.csv
  revision/manifests/download_dates.md
  revision/manifests/gdc_queries/*.json   (query, release and retrieval time)

Run from the repository root:
    python pipeline/r06_download_manifest.py
"""

import json
import shutil
from pathlib import Path

import pandas as pd

DL = Path("data/downloads")
MAN = Path("revision/manifests")


def main():
    (MAN / "gdc_queries").mkdir(parents=True, exist_ok=True)
    shutil.copy(DL / "DOWNLOAD_LOG.csv", MAN / "DOWNLOAD_LOG.csv")
    log = pd.read_csv(DL / "DOWNLOAD_LOG.csv")
    rows = []
    for _, r in log.sort_values("path").iterrows():
        note = ""
        if r["path"].endswith(".json") and "gdc" in r["url"]:
            meta = json.load(open(r["path"]))
            rel = meta.get("gdc_status", {}).get("data_release", "")
            note = f"{len(meta.get('hits', []))} records; {rel}"
            q = {k: meta[k] for k in ("query", "retrieved_utc", "gdc_status") if k in meta}
            name = r["path"].replace(str(DL) + "/", "").replace("/", "__")
            with open(MAN / "gdc_queries" / name, "w") as f:
                json.dump(q, f, indent=1)
        rows.append(f"| `{r['path'].replace(str(DL) + '/', '')}` | {r['description']} | "
                    f"{r['retrieved_utc'][:10]} | {note} |")
    md = ["| File (under data/downloads/) | Content | Retrieved (UTC) | Note |", "|---|---|---|---|", *rows]
    (MAN / "download_dates.md").write_text("\n".join(md) + "\n")
    print(f"{len(rows)} files")


if __name__ == "__main__":
    main()
