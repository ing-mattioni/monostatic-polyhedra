"""Exact certification: integer plane data -> exact rational vertices -> verify/exact_verify.py.

Planes are given in Reshetov's convention (p_i integer, plane p_i . x = |p_i|^2) or generally as
integer (A, b) with A x <= b.  Vertices are computed EXACTLY by enumerating all plane triples
(Fraction arithmetic, exact feasibility test), so no float information enters the certificate.

usage: python certify.py in.json out_name      (in.json: {"p": [[..],..]} or {"A":..,"b":..})
"""
import os, sys, json, subprocess, itertools
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PY = sys.executable


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def exact_vertices(A, b):
    """A: list of integer 3-vectors, b: list of integers (or Fractions). Returns list of Fraction triples."""
    F = len(A)
    A = [tuple(int(x) for x in a) for a in A]
    b = [Fraction(x) for x in b]
    verts = {}
    for i, j, k in itertools.combinations(range(F), 3):
        D = det3(A[i], A[j], A[k])
        if D == 0:
            continue
        # Cramer
        bi, bj, bk = b[i], b[j], b[k]
        cols = []
        for col in range(3):
            M = [list(A[i]), list(A[j]), list(A[k])]
            M[0][col] = bi
            M[1][col] = bj
            M[2][col] = bk
            cols.append(Fraction(det3(M[0], M[1], M[2])) / D)
        x = tuple(cols)
        if x in verts:
            continue
        if all(A[l][0] * x[0] + A[l][1] * x[1] + A[l][2] * x[2] <= b[l] for l in range(F)):
            verts[x] = True
    return list(verts.keys())


def write_json(verts, name, path):
    with open(path, "w") as fh:
        json.dump({"name": name, "vertices": [[str(c) for c in v] for v in verts]}, fh, indent=0)


def run_verifier(path, report):
    r = subprocess.run([PY, os.path.join(ROOT, "verify", "exact_verify.py"), path, "--report", report],
                       capture_output=True, text=True)
    txt = r.stdout
    S = H = U = None
    mono = None
    for line in txt.splitlines():
        if "MONO-STABLE" in line:
            mono = "True" in line
    return r.returncode, mono, txt


def certify_planes(A, b, name, outdir=None, quiet=False):
    outdir = outdir or os.path.join(HERE, "certified")
    os.makedirs(outdir, exist_ok=True)
    verts = exact_vertices(A, b)
    path = os.path.join(outdir, name + ".json")
    rep = os.path.join(outdir, name + ".report.json")
    write_json(verts, name, path)
    code, mono, txt = run_verifier(path, rep)
    with open(os.path.join(outdir, name + ".report.txt"), "w") as fh:
        fh.write(txt)
    with open(os.path.join(outdir, name + ".planes.json"), "w") as fh:
        json.dump({"A": [[int(x) for x in a] for a in A], "b": [str(x) for x in b]}, fh)
    if not quiet:
        print("\n".join(l for l in txt.splitlines() if any(k in l for k in
                        ("convex hull", "face sizes", "S =", "MONO", "CERTIF", "certified", "(S", "equilibria"))))
    return code, mono, txt, len(verts)


def certify_p(P, name, **kw):
    P = [[int(x) for x in p] for p in P]
    b = [p[0] * p[0] + p[1] * p[1] + p[2] * p[2] for p in P]
    return certify_planes(P, b, name, **kw)


def round_p(P, digits):
    """Scale P so that max |entry| = 10**digits and round to integers."""
    P = np.asarray(P, float)
    s = 10 ** digits / np.abs(P).max()
    return np.rint(P * s).astype(np.int64)


def certify_float(P, name, digits_list=(3, 4, 5, 6, 7), quiet=True):
    """Round a float p-parametrised solution to integer p (smallest digits that keep G > 0 in float),
    then certify exactly.  Returns (digits, exit code, mono, float G of rounded body)."""
    from core import evaluate_p
    P = np.asarray(P, float)
    F0 = evaluate_p(P)["nfaces"]
    for dg in digits_list:
        Q = round_p(P, dg)
        r = evaluate_p(Q.astype(float))
        if r["status"] == 0 and r["G"] > 0 and r["nfaces"] == F0:
            code, mono, txt, nv = certify_p(Q.tolist(), name, quiet=quiet)
            return dict(digits=dg, code=code, mono=mono, G_rounded=float(r["G"]), nfaces=r["nfaces"],
                        p=Q.tolist(), txt=txt)
    return None


if __name__ == "__main__":
    d = json.load(open(sys.argv[1]))
    name = sys.argv[2]
    if "p" in d:
        code, mono, txt, nv = certify_p(d["p"], name)
    else:
        code, mono, txt, nv = certify_planes(d["A"], [Fraction(x) for x in d["b"]], name)
    print("exit", code, "mono", mono)

