"""Summarise runs.jsonl: per tag -> runs, total evals, best G / Gw / G2, face counts."""
import json, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
rows = [json.loads(l) for l in open(os.path.join(HERE, "runs.jsonl")) if l.strip()]
by = defaultdict(list)
for r in rows:
    by[r["tag"]].append(r)
tot = 0
for tag, rs in by.items():
    ev = sum(r["evals"] for r in rs)
    tot += ev
    key = "G2" if "G2" in rs[0] else "G"
    best = max(rs, key=lambda r: r.get(key, -1e9))
    npos = sum(1 for r in rs if r.get(key, -1) > 0)
    print(f"{tag:14s} runs={len(rs):4d} evals={ev / 1e6:8.1f}M  best {key}={best.get(key):+.6f} "
          f"(F={best['F']}, nfaces={best['nfaces']}, drop={best.get('drop')}, split={best.get('split')}, "
          f"aspect={best.get('aspect', 0):.1f})  #{key}>0: {npos}")
    if len(sys.argv) > 1 and sys.argv[1] == "-v":
        for r in sorted(rs, key=lambda r: -r.get(key, -1e9))[:5]:
            print("     ", {k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()
                           if k in ("seed", "drop", "sigma", "G", "Gw", "G2", "nfaces", "nstable", "aspect", "split", "pattern")})
print(f"total evaluations logged: {tot / 1e6:.1f}M")
