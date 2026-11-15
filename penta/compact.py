"""
compact.py -- UNTRUSTED post-processing: shorten the multipliers in a certificate file.

    python compact.py Z|P [--digits 6] [--workers 14]

Every multiplier (y0, Y1) is rounded to a few significant decimal digits; the rounded certificate is kept
only if it still passes core.check on its box (otherwise the original is kept).  The result replaces
certs_<chart>.json.gz (the original is kept as certs_<chart>.full.json.gz).  Run verify.py afterwards:
that is the actual proof check, compact.py is only a size optimisation.
"""
import argparse
import gzip
import json
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def rnd(x, digits):
    return float(f"{x:.{digits}g}")


def compact_cell(args):
    chart, lo, hi, node, digits = args
    from core import BoxData, check
    stats = {"rounded": 0, "kept": 0}

    def rec(lo, hi, node):
        if node.get("skip"):
            return
        certs = node.get("certs", [])
        if certs:
            bd = BoxData(chart, lo, hi)
            for cert in certs:
                tid, y0, Y1 = cert
                y0r = [rnd(t, digits) for t in y0]
                Y1r = [[rnd(t, digits) for t in row] for row in Y1]
                if check(bd, tid, y0r, Y1r):
                    cert[1], cert[2] = y0r, Y1r
                    stats["rounded"] += 1
                else:
                    stats["kept"] += 1
        if "children" in node:
            axis = node["axis"]
            mid = (lo[axis] + hi[axis]) / 2
            hi1 = list(hi); hi1[axis] = mid
            lo2 = list(lo); lo2[axis] = mid
            rec(list(lo), hi1, node["children"][0])
            rec(lo2, list(hi), node["children"][1])

    rec(list(lo), list(hi), node)
    return node, stats


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("chart", choices=["Z", "P"])
    ap.add_argument("--digits", type=int, default=6)
    ap.add_argument("--workers", type=int, default=14)
    a = ap.parse_args()
    from concurrent.futures import ProcessPoolExecutor
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, f"certs_{a.chart}.json.gz")
    with gzip.open(path, "rt") as fh:
        d = json.load(fh)
    g, rho0 = d["grid"], d["rho0"]
    dhi = [1.0, 1.0, 1.0] if a.chart == "Z" else [rho0, 1.0, 1.0]
    tasks = []
    for c in d["cells"]:
        idx = c["index"]
        lo = [dhi[m] * idx[m] / g[m] for m in range(3)]
        hi = [dhi[m] * (idx[m] + 1) / g[m] for m in range(3)]
        tasks.append((a.chart, lo, hi, c["node"], a.digits))
    tot = {"rounded": 0, "kept": 0}
    with ProcessPoolExecutor(a.workers) as ex:
        for c, (node, st) in zip(d["cells"], ex.map(compact_cell, tasks, chunksize=1)):
            c["node"] = node
            for k in tot:
                tot[k] += st[k]
    tmp = path + ".tmp"
    with gzip.open(tmp, "wt", compresslevel=9) as fh:
        json.dump(d, fh, separators=(",", ":"))
    os.replace(path, os.path.join(here, f"certs_{a.chart}.full.json.gz"))
    os.replace(tmp, path)                      # atomic: readers never see a partial file
    print(f"chart {a.chart}: {tot}; wrote {path} ({os.path.getsize(path) / 1e6:.1f} MB)")
