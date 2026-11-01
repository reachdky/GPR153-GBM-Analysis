"""Regenerate figures, tables, source data, methods record and claims ledger from the final results into a scratch
base and compare with revision/ by SHA-256. Usage (repo root): python regen_check.py <scratch base> <out csv>"""
import hashlib, shutil, sys
from pathlib import Path
sys.path.insert(0, "pipeline")
import common  # noqa: E402
out, csv = Path(sys.argv[1]), Path(sys.argv[2])
shutil.rmtree(out, ignore_errors=True)
shutil.copytree("revision/results", out / "results")
shutil.copytree("revision/manifests", out / "manifests")
common.BASE, common.RES, common.MAN = out, out / "results", out / "manifests"
common.SRC, common.FIG, common.FIGMAP = out / "source_data", out / "figures", out / "figure_script_map.csv"
import f10_main_figures, f11_supplement_figures, f12_tables  # noqa: E402,E401
f10_main_figures.R = f11_supplement_figures.R = common.RES
f12_tables.R, f12_tables.OUT = common.RES, out / "tables"
f10_main_figures.main(); f11_supplement_figures.main(); f12_tables.main()
import a90_methods_record  # noqa: E402
a90_methods_record.main()
import x02_claims_ledger  # noqa: E402
x02_claims_ledger.R, x02_claims_ledger.V1 = common.RES, common.RES / "claims_trace_v1.csv"
x02_claims_ledger.CONC = out / "methods_concordance.csv"
x02_claims_ledger.main()
h = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
rows = ["file,sha256_final,sha256_regenerated,status"]
for d in ["figures", "tables", "source_data"]:
    for p in sorted((Path("revision") / d).iterdir()):
        q = out / d / p.name
        rows.append(f"{d}/{p.name},{h(p)},{h(q) if q.exists() else ''},"
                    f"{'identical' if q.exists() and h(p) == h(q) else ('not regenerated' if not q.exists() else 'differs')}")
for f in ["figure_script_map.csv", "methods_parameters.json", "methods_concordance.csv", "results/claims_trace_v2.csv"]:
    p, q = Path("revision") / f, out / f
    rows.append(f"{f},{h(p)},{h(q)},{'identical' if h(p) == h(q) else 'differs'}")
csv.write_text("\n".join(rows) + "\n")
