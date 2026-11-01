"""Test result recorder: appends one row per test to revision/results/tests/test_results.csv."""

import csv
import datetime
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))
import common  # noqa: E402

OUT = common.outdir("tests") / "test_results.csv"
FIELDS = ["test", "name", "result", "data", "detail", "run_utc"]


def record(test, name, passed, data, detail):
    new = not OUT.exists()
    rows = []
    if not new:
        with open(OUT) as f:
            rows = [r for r in csv.DictReader(f) if r["test"] != test]
    rows.append({"test": test, "name": name, "result": "PASS" if passed else "FAIL", "data": data,
                 "detail": detail, "run_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")})
    rows.sort(key=lambda r: (len(r["test"]), r["test"]))
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"{test} {name}: {'PASS' if passed else 'FAIL'} ({detail})")
    return passed
