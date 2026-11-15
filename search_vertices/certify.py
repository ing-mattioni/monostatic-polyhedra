"""Certify a float solution exactly with verify/exact_verify.py.
usage: python certify.py sol.npy NAME [--denoms 1e6,1e7,1e8,1e9]
Normalises (centroid->0, max radius 1), rationalises, runs the exact verifier,
copies certified polyhedra to certified/V{n}_NAME_d{D}.json (+ report)."""
import os, sys, json, subprocess, argparse, re
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import core

PY = sys.executable
VER = os.path.join(ROOT, "verify", "exact_verify.py")
CERT = os.path.join(HERE, "certified")


def certify(npy, name, denoms=(10**6, 10**7, 10**8, 10**9), keep_all=False):
    X = core.normalize(np.load(npy))
    inf = core.info(X)
    os.makedirs(CERT, exist_ok=True)
    tmp = os.path.join(HERE, "tmp_cert"); os.makedirs(tmp, exist_ok=True)
    npn = os.path.join(tmp, name + ".npy"); np.save(npn, X)
    results = []
    for D in denoms:
        tag = f"V{len(X)}_{name}_d{int(round(np.log10(D)))}"
        pj = os.path.join(tmp, tag + ".json"); rj = os.path.join(tmp, tag + ".report.json")
        subprocess.run([PY, VER, "--from-npy", npn, "--denom", str(D), "--out", pj, "--name", tag],
                       check=True, capture_output=True, text=True)
        r = subprocess.run([PY, VER, pj, "--report", rj], capture_output=True, text=True)
        txt = r.stdout
        m = re.search(r"\(S, H, U\) = \((\d+), (\d+), (\d+)\)", txt)
        shu = tuple(int(v) for v in m.groups()) if m else None
        mono = "MONO-UNSTABLE: True" in txt
        hv = re.search(r"convex hull\s+: V = (\d+)", txt)
        nv = int(hv.group(1)) if hv else None
        ok = (r.returncode == 0) and mono
        results.append(dict(D=D, code=r.returncode, shu=shu, mono=mono, nvert_exact=nv, ok=ok))
        print(f"  denom 1e{int(round(np.log10(D)))}: exit={r.returncode} (S,H,U)={shu} mono-unstable={mono} V_exact={nv}", flush=True)
        if ok:
            for src, ext in ((pj, ".json"), (rj, ".report.json")):
                with open(src) as fh, open(os.path.join(CERT, tag + ext), "w") as gh:
                    gh.write(fh.read())
            with open(os.path.join(CERT, tag + ".report.txt"), "w") as gh:
                gh.write(txt)
            if not keep_all:
                break
    return inf, results


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npy"); ap.add_argument("name")
    ap.add_argument("--denoms", default="1e6,1e7,1e8,1e9")
    a = ap.parse_args()
    ds = [int(float(s)) for s in a.denoms.split(",")]
    inf, res = certify(a.npy, a.name, ds)
    print(f"float F={inf['F']:.4e}  nvert={inf['nvert']}  far={inf['far']}  parent={inf['parent'].tolist()}")
