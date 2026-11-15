"""
search.py -- UNTRUSTED branch-and-bound that produces the certificate file of a chart.

    python search.py Z|P [--rho0 0.5] [--workers 14] [--maxdepth 30]

The parameter domain of the chart is cut into a fixed grid of top cells; on each cell the 375 rooted
trees are certified recursively: at a node (box) every pending tree is tried with the multipliers
inherited from the parent box and, failing that, with multipliers proposed by propose.py; the trees
that pass core.check are recorded at the node, the box is bisected (longest edge) for the others.
Chart Z: boxes with lo_u * lo_v * lo_w > 1 - rho0 are skipped (that region is covered by chart P).
Output: certs_<chart>.json.gz   (verified independently by verify.py).
"""
import argparse
import gzip
import json
import os
import sys
import time
from fractions import Fraction

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

GRID = {"Z": (8, 8, 8), "P": (4, 8, 8)}


def domain(chart, rho0):
    return ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]) if chart == "Z" else ([0.0, 0.0, 0.0], [rho0, 1.0, 1.0])


def in_skip_region(chart, lo, rho0):
    """exact test of  lo_u * lo_v * lo_w > 1 - rho0  (chart Z only)"""
    if chart != "Z":
        return False
    return Fraction(lo[0]) * Fraction(lo[1]) * Fraction(lo[2]) > 1 - Fraction(rho0)


def work(args):
    chart, idx, lo, hi, rho0, maxdepth = args
    import warnings
    warnings.filterwarnings("ignore")
    import numpy as np
    from core import BoxData, check, TREES
    from propose import Proposer, edge_mats_float, SAMPLES
    pr = Proposer()
    stats = {"nodes": 0, "certs": 0, "sdp": 0, "maxdepth": 0, "failed": 0, "skipped": 0}

    def rec(lo, hi, pending, depth):
        stats["nodes"] += 1
        stats["maxdepth"] = max(stats["maxdepth"], depth)
        if in_skip_region(chart, lo, rho0):
            stats["skipped"] += 1
            return {"skip": True}
        bd = BoxData(chart, lo, hi)
        c = [(a + b) / 2 for a, b in zip(lo, hi)]
        h = [(b - a) / 2 for a, b in zip(lo, hi)]
        sm = None
        certs, rest = [], {}
        for tid, cand in pending.items():
            ok = False
            if cand is not None:
                y0, Y1, c_old = cand
                y0n = [y0[i] + sum(Y1[i][k] * (c[k] - c_old[k]) for k in range(3)) for i in range(4)]
                if check(bd, tid, y0n, Y1):
                    certs.append([tid, y0n, Y1])
                    ok = True
            if not ok:
                if sm is None:
                    sm = [edge_mats_float(chart, np.array(c) + np.array(h) * xi) for xi in SAMPLES]
                stats["sdp"] += 1
                r = pr.propose(chart, lo, hi, tid, sm)
                if r is not None:
                    lam, y0, Y1 = r
                    if check(bd, tid, y0, Y1):
                        certs.append([tid, y0, Y1])
                        ok = True
                    else:
                        rest[tid] = (y0, Y1, c)
                else:
                    rest[tid] = None
        stats["certs"] += len(certs)
        node = {"certs": certs}
        if rest:
            if depth >= maxdepth:
                node["FAILED"] = sorted(rest)
                stats["failed"] += len(rest)
                return node
            widths = [b - a for a, b in zip(lo, hi)]
            axis = max(range(3), key=lambda k: (widths[k], -k))
            mid = (lo[axis] + hi[axis]) / 2
            hi1 = list(hi); hi1[axis] = mid
            lo2 = list(lo); lo2[axis] = mid
            node["axis"] = axis
            node["children"] = [rec(list(lo), hi1, dict(rest), depth + 1), rec(lo2, list(hi), dict(rest), depth + 1)]
        return node

    t0 = time.time()
    node = rec(list(lo), list(hi), {tid: None for tid in range(len(TREES))}, 0)
    stats["time"] = time.time() - t0
    return idx, node, stats


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("chart", choices=["Z", "P"])
    ap.add_argument("--rho0", type=float, default=0.5)
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--maxdepth", type=int, default=30)
    ap.add_argument("--only", type=int, default=None, help="process only the first N cells (testing)")
    a = ap.parse_args()
    from concurrent.futures import ProcessPoolExecutor, as_completed
    dlo, dhi = domain(a.chart, a.rho0)
    g = GRID[a.chart]
    tasks = []
    for i in range(g[0]):
        for j in range(g[1]):
            for k in range(g[2]):
                idx = (i, j, k)
                lo = [dlo[m] + (dhi[m] - dlo[m]) * idx[m] / g[m] for m in range(3)]
                hi = [dlo[m] + (dhi[m] - dlo[m]) * (idx[m] + 1) / g[m] for m in range(3)]
                tasks.append((a.chart, idx, lo, hi, a.rho0, a.maxdepth))
    if a.only:
        tasks = tasks[:a.only]
    t0 = time.time()
    cells, tot = [], {}
    with ProcessPoolExecutor(a.workers) as ex:
        futs = [ex.submit(work, t) for t in tasks]
        for n, f in enumerate(as_completed(futs)):
            idx, node, st = f.result()
            cells.append({"index": list(idx), "node": node})
            for key, val in st.items():
                tot[key] = max(tot.get(key, 0), val) if key == "maxdepth" else tot.get(key, 0) + val
            if (n + 1) % 16 == 0 or n + 1 == len(tasks):
                print(f"[{a.chart}] {n + 1}/{len(tasks)} cells, {time.time() - t0:.0f}s, totals {tot}", flush=True)
    cells.sort(key=lambda c: c["index"])
    out = {"chart": a.chart, "rho0": a.rho0, "grid": list(g), "cells": cells}
    name = f"certs_{a.chart}.json.gz" if not a.only else f"certs_{a.chart}_partial.json.gz"
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), name)
    with gzip.open(path, "wt") as fh:
        json.dump(out, fh)
    print(f"wrote {path}; failed (box,tree) pairs: {tot.get('failed', 0)}")
