"""Structural statistics of the trees not certified by the plain joint SDP.

    python analyze_V8.py [V] [uncertified json]

Prints counts by (graph, root), the distribution of tree height and of the
depth profile among uncertified trees versus all trees, and tests simple
candidate characterisations."""
import json
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import certify_joint  # noqa: F401,E402  (sets sys.path for lowerbound)
import enumerate_cases as ec  # noqa: E402


def depths(parent, v):
    d = []
    for i in range(len(parent)):
        x, s = i, 0
        while x != v:
            x = parent[x]; s += 1
        d.append(s)
    return d


def main(V=8, path=None):
    path = path or os.path.join(HERE, f"uncertified_V{V}.json")
    unc = json.load(open(path))
    gds = ec.graph_data(V)
    U = {(c["graph"], c["root"], tuple(c["parent"])) for c in unc}
    print(f"V={V}: {len(U)} uncertified trees")
    byr = Counter((g, r) for g, r, _ in U)
    for gd in gds:
        for v in gd.roots:
            n = len(gd.trees[v])
            print(f"  graph {gd.gid} (degrees {gd.degrees}) root {v} (deg {gd.degrees[v]}): "
                  f"{byr[(gd.gid, v)]} / {n} uncertified")
    hall, hunc = Counter(), Counter()
    kall, kunc = Counter(), Counter()     # number of vertices at depth >= 4
    for gd in gds:
        for v in gd.roots:
            for p in gd.trees[v]:
                d = depths(p, v)
                h = max(d)
                k = sum(1 for t in d if t >= 4)
                hall[h] += 1; kall[k] += 1
                if (gd.gid, v, p) in U:
                    hunc[h] += 1; kunc[k] += 1
    print("  tree height: (#uncertified / #trees)")
    for h in sorted(hall):
        print(f"    height {h}: {hunc[h]} / {hall[h]}")
    print("  number of vertices at depth >= 4: (#uncertified / #trees)")
    for k in sorted(kall):
        print(f"    {k}: {kunc[k]} / {kall[k]}")
    # characterisation: some face F not containing v with all three vertices at depth >= t
    for t in (3, 4, 5):
        hit_all = hit_unc = 0
        for gd in gds:
            for v in gd.roots:
                tf = gd.decomposition[v]
                for p in gd.trees[v]:
                    d = depths(p, v)
                    if any(min(d[a] for a in f) >= t for f in tf):
                        hit_all += 1
                        if (gd.gid, v, p) in U:
                            hit_unc += 1
        print(f"  trees having a face (avoiding v) with all 3 vertices at depth >= {t}: "
              f"{hit_all}; uncertified among them: {hit_unc} (of {len(U)} uncertified)")


if __name__ == "__main__":
    V = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    main(V, sys.argv[2] if len(sys.argv) > 2 else None)
