"""SHA-256 comparison of two result trees. Usage: python sha_compare.py <label> <tree A> <tree B>; prints CSV rows."""
import hashlib, sys
from pathlib import Path
lab, a, b = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
h = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
fa = {p.relative_to(a) for p in a.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
fb = {p.relative_to(b) for p in b.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
for f in sorted(fa | fb):
    sa = h(a / f) if f in fa else ""
    sb = h(b / f) if f in fb else ""
    st = "identical" if sa == sb else ("only in A" if not sb else "only in B" if not sa else "differs")
    print(f"{lab},{f},{sa},{sb},{st}")
