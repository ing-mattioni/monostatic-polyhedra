"""
Step 4: re-test the systems that the baseline SDP (Algorithm 2 of
arXiv:2401.17906) does not certify, with the valid strengthening inequalities
of certify.py (level "b": |r_i|^2 <= 1, r_{i,1} <= 1; level "bc": in addition
|r_i|^2 <= |r_a|^2 and r_i.r_a <= |r_a|^2 for all tree ancestors a of i).

Reads  results_V<V>.csv  (written by `certify.py run V`), writes
uncertified_V<V>.json with, for every baseline-uncertified system, the
optimal z and the exact-verification flag under each level, and prints the
uncertified systems grouped by graph and root.

Usage:  python strengthen.py [V] [--workers 4]
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import certify as C  # noqa: E402
import enumerate_cases as ec  # noqa: E402


def load_results(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            rows.append(dict(graph=int(r["graph"]), root=int(r["root"]),
                             parent=[int(t) for t in r["parent"].split()],
                             face=[int(t) for t in r["tetrahedron"].split()][1:],
                             z=float(r["z"]), certified=bool(int(r["certified"])),
                             y=[float(t) for t in r["y"].split()]))
    return rows


def depth(parent, i):
    d = 0
    while parent[i] != -1:
        i = parent[i]
        d += 1
    return d


def describe(case):
    """Short combinatorial description: depths of the face vertices in the
    tree and the tree's depth profile."""
    p = case["parent"]
    V = len(p)
    depths = [depth(p, i) for i in range(V)]
    return dict(face_depths=sorted(depths[k] for k in case["face"]),
                tree_depth=max(depths),
                is_path=(max(depths) == V - 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("V", nargs="?", type=int, default=7)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    V = a.V
    rows = load_results(os.path.join(HERE, f"results_V{V}.csv"))
    gds = ec.graph_data(V, with_trees=False)
    unc = [r for r in rows if not r["certified"]]
    print(f"V={V}: {len(rows)} systems, {len(rows) - len(unc)} certified at baseline, "
          f"{len(unc)} not certified")
    cases = [(r["graph"], r["root"], r["parent"], r["face"]) for r in unc]
    res = {"base": unc}
    for lev in ("b", "bc"):
        res[lev] = C.run_all(V, lev, 2.0, workers=min(a.workers, 4), gds=gds,
                             cases=cases, verbose=False)
        n_unc = sum(not r["certified"] for r in res[lev])
        print(f"  level {lev}: {len(unc) - n_unc} of the {len(unc)} newly certified, "
              f"{n_unc} remain uncertified")
    out = []
    for k, r in enumerate(unc):
        d = dict(graph=r["graph"], root=r["root"], parent=r["parent"],
                 face=r["face"], tetrahedron=[r["root"]] + r["face"],
                 z_base=r["z"], cert_base=r["certified"],
                 z_b=res["b"][k]["z"], cert_b=res["b"][k]["certified"],
                 z_bc=res["bc"][k]["z"], cert_bc=res["bc"][k]["certified"],
                 y_bc=res["bc"][k]["y"])
        d.update(describe(d))
        out.append(d)
    with open(os.path.join(HERE, f"uncertified_V{V}.json"), "w") as f:
        json.dump(out, f, indent=1)

    # ---- grouped summary --------------------------------------------------
    print("\nBaseline-uncertified systems grouped by (graph, root):")
    groups = defaultdict(list)
    for d in out:
        groups[(d["graph"], d["root"])].append(d)
    for (g, v), L in sorted(groups.items()):
        gd = gds[g]
        print(f"  graph {g} (degrees {gd.degrees}), root {v} (deg {gd.degrees[v]}): "
              f"{len(L)} uncertified at baseline, "
              f"{sum(not d['cert_b'] for d in L)} after (b), "
              f"{sum(not d['cert_bc'] for d in L)} after (b)+(c)")
        fd = Counter((tuple(d["face_depths"]), d["tree_depth"]) for d in L)
        print("     (face depths, tree depth): counts ->",
              dict(sorted(fd.items())))
    print("\nFace depth profiles over all baseline-uncertified systems:",
          dict(sorted(Counter(tuple(d["face_depths"]) for d in out).items())))
    print("Tree is a Hamiltonian path from the root:",
          sum(d["is_path"] for d in out), "of", len(out))


if __name__ == "__main__":
    main()
