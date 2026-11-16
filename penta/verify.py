"""
verify.py -- TRUSTED verifier of the certificate files certs_Z.json.gz and certs_P.json.gz.

    python verify.py [--workers 14]

It uses only shapes.py (exact symbolic derivation of the chart polynomials, with internal identity
checks) and core.py (rigorous ball-arithmetic certificate check); it does not import the search code.

What is checked (exit status 0 iff everything passes):
 1. the two files use the same dyadic rho0 in (0,1);
 2. chart Z: the top cells are exactly the cells of the stated grid on [0,1]^3;
    chart P: exactly the cells of the stated grid on [0,rho0] x [0,1]^2;
 3. in every cell the stored bisection tree is traversed; the boxes of the children are COMPUTED here
    (bisection of the stated axis at the exact midpoint), so the leaves partition the cell;
 4. along every root-to-leaf path each of the 375 rooted spanning trees of the face-adjacency graph
    receives a certificate that passes core.check on the box of the node where it is stored -- unless
    the leaf is a "skip" leaf of chart Z, which is accepted only if lo_u*lo_v*lo_w > 1 - rho0 holds
    exactly (those shapes have t_1 = uvw > 1 - rho0 and are covered by chart P).
Consequently, for every shape of chart Z with t_1 <= 1 - rho0 and every shape of chart P, and every
rooted tree, a certificate holds at that shape; see README.md for why this proves the theorem.
"""
import argparse
import gzip
import json
import math
import os
import sys
import time
from fractions import Fraction

os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

NTREES = 375


class VerificationError(Exception):
    pass


def require(cond, msg):
    """explicit check (unlike `assert`, not removed by `python -O`)"""
    if not cond:
        raise VerificationError(msg)


def is_dyadic_float(x):
    return isinstance(x, float) and math.isfinite(x)


def verify_cell(args):
    chart, lo, hi, node, rho0 = args
    from core import BoxData, check, TREES
    require(len(TREES) == NTREES, "wrong number of rooted trees")
    stats = {"nodes": 0, "leaves": 0, "skip_leaves": 0, "certs": 0, "maxdepth": 0}

    def fail(msg):
        raise VerificationError(f"chart {chart} cell {lo}-{hi}: {msg}")

    def rec(lo, hi, node, pending, depth):
        stats["nodes"] += 1
        stats["maxdepth"] = max(stats["maxdepth"], depth)
        if node.get("skip"):
            if chart != "Z":
                fail("skip leaf outside chart Z")
            if not (Fraction(lo[0]) * Fraction(lo[1]) * Fraction(lo[2]) > 1 - Fraction(rho0)):
                fail(f"illegal skip leaf at lo = {lo}")
            stats["leaves"] += 1
            stats["skip_leaves"] += 1
            return
        if "FAILED" in node:
            fail("node marked FAILED by the search")
        pending = set(pending)
        certs = node.get("certs", [])
        if certs:
            bd = BoxData(chart, lo, hi)
            for tid, y0, Y1 in certs:
                if not (isinstance(tid, int) and 0 <= tid < NTREES):
                    fail("bad tree id")
                if not (len(y0) == 4 and all(is_dyadic_float(t) for t in y0)
                        and len(Y1) == 4 and all(len(r) == 3 and all(is_dyadic_float(t) for t in r) for r in Y1)):
                    fail("malformed multipliers")
                if not check(bd, tid, y0, Y1):
                    fail(f"certificate for tree {tid} does not verify on box {lo}-{hi}")
                pending.discard(tid)
                stats["certs"] += 1
        if "children" in node:
            axis = node["axis"]
            if axis not in (0, 1, 2) or len(node["children"]) != 2:
                fail("bad bisection record")
            mid = (lo[axis] + hi[axis]) / 2
            if Fraction(mid) * 2 != Fraction(lo[axis]) + Fraction(hi[axis]) or not (lo[axis] < mid < hi[axis]):
                fail("midpoint not exact")
            hi1 = list(hi); hi1[axis] = mid
            lo2 = list(lo); lo2[axis] = mid
            rec(list(lo), hi1, node["children"][0], pending, depth + 1)
            rec(lo2, list(hi), node["children"][1], pending, depth + 1)
        else:
            if pending:
                fail(f"leaf with {len(pending)} uncertified trees")
            stats["leaves"] += 1

    rec(list(lo), list(hi), node, range(NTREES), 0)
    return stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--dir", default=os.path.dirname(os.path.abspath(__file__)))
    a = ap.parse_args()
    from concurrent.futures import ProcessPoolExecutor
    t0 = time.time()
    data = {}
    for chart in ("Z", "P"):
        with gzip.open(os.path.join(a.dir, f"certs_{chart}.json.gz"), "rt") as fh:
            data[chart] = json.load(fh)
        require(isinstance(data[chart], dict) and data[chart].get("chart") == chart, f"file of chart {chart}: wrong header")
    rho0 = data["Z"]["rho0"]
    require(is_dyadic_float(rho0) and rho0 == data["P"]["rho0"] and 0 < rho0 < 1, "rho0 missing, different in the two files, or out of range")
    print(f"verify.py: rho0 = {rho0} (= {Fraction(rho0)}), trees = {NTREES}")
    import shapes  # noqa  (runs the exact symbolic derivation and its identity checks)
    for name in ("Z", "P"):
        shapes.get_chart(name)
    print("chart polynomials derived and identity-checked symbolically (shapes.py)")
    ok = True
    for chart in ("Z", "P"):
        d = data[chart]
        g = d["grid"]
        dlo = [0.0, 0.0, 0.0]
        dhi = [1.0, 1.0, 1.0] if chart == "Z" else [rho0, 1.0, 1.0]
        require(isinstance(g, list) and len(g) == 3 and all(isinstance(n, int) and n >= 1 for n in g), "bad grid")
        expected = {(i, j, k) for i in range(g[0]) for j in range(g[1]) for k in range(g[2])}
        got = [tuple(c["index"]) for c in d["cells"]]
        require(len(got) == len(set(got)) and set(got) == expected, f"chart {chart}: cells do not form the grid")
        tasks = []
        for c in d["cells"]:
            idx = c["index"]
            lo, hi = [], []
            for m in range(3):
                w = Fraction(dhi[m]) - Fraction(dlo[m])
                l = Fraction(dlo[m]) + w * idx[m] / g[m]
                h = Fraction(dlo[m]) + w * (idx[m] + 1) / g[m]
                lf, hf = float(l), float(h)
                require(Fraction(lf) == l and Fraction(hf) == h, "grid point not exactly representable")
                lo.append(lf)
                hi.append(hf)
            tasks.append((chart, lo, hi, c["node"], rho0))
        tot = {}
        try:
            with ProcessPoolExecutor(a.workers) as ex:
                for st in ex.map(verify_cell, tasks, chunksize=1):
                    for key, val in st.items():
                        tot[key] = max(tot.get(key, 0), val) if key == "maxdepth" else tot.get(key, 0) + val
        except Exception as e:           # any failure in a worker (VerificationError or malformed data)
            print("VERIFICATION FAILED:", e)
            ok = False
            break
        print(f"chart {chart}: grid {g}, {len(tasks)} cells: {tot}  [{time.time() - t0:.0f}s]")
    if ok:
        print("ALL CHECKS PASSED: for every shape of chart Z with t_1 <= 1 - rho0 and every shape of chart P, "
              "each of the 375 rooted trees is refuted.")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    try:
        main()
    except VerificationError as e:
        print("VERIFICATION FAILED:", e)
        sys.exit(1)
