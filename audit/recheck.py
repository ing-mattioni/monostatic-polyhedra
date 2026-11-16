"""
recheck.py -- INDEPENDENT exact re-check of certificates (audit).  Shares no code with penta/core.py or
penta/shapes.py.  Chart polynomials nu/eta are read as DATA from shapes_data.json; derive_check.py verifies
them exactly against covectors derived from the vertices.  Here everything is exact integer/rational
arithmetic (Python ints, fractions.Fraction):

  * the box [lo,hi] of every node is reconstructed from the stored bisection trees (Fractions);
  * on a box, with x = c + h*s, s in [-1,1]^3, every needed polynomial is Taylor-expanded exactly in s;
  * N(s) = sum_i y_i(s) sym(nu_i eta_e^T)(s) is formed exactly (common positive denominator dropped);
  * N = N0 + sum_k s_k N_k + R(s); |R_ab| <= E_ab := sum of |coefficients| of degree >= 2;
    ||R||_2 <= r := max_a sum_b E_ab  (symmetric matrix: spectral norm <= infinity norm);
  * -(N0 + sum_k sigma_k N_k) - r I is checked positive definite at the 8 corners sigma (exact Sylvester),
    and y_i >= 0 on the box exactly.  If this cruder test fails, the box is subdivided (8 children,
    same affine multiplier) up to a few levels before reporting failure.

Usage:  python recheck.py sample Z <certs.json.gz> <n_uniform> <seed>      (box certificates)
        python recheck.py points Z <certsZ> <certsP> <n_points> <seed>    (pointwise end-to-end)
"""
import gzip
import itertools
import json
import random
import sys
import time
from fractions import Fraction

DATA = json.load(open(r"C:\Projects\Claude4\audit\shapes_data.json"))
ADJ_PAIRS = [(i, j) for i in range(5) for j in range(i + 1, 5) if (i, j) != (3, 4)]
ADJ = {i: sorted(j for j in range(5) if j != i and (min(i, j), max(i, j)) in ADJ_PAIRS) for i in range(5)}


# ---------------------------------------------------------------- trees (own enumeration, same order convention)
def my_trees():
    out = []
    for s in range(5):
        others = [i for i in range(5) if i != s]
        for choice in itertools.product(*[ADJ[i] for i in others]):
            par = {i: j for i, j in zip(others, choice)}
            good = True
            for i in others:
                seen, x = set(), i
                while x != s:
                    if x in seen:
                        good = False
                        break
                    seen.add(x)
                    x = par[x]
                if not good:
                    break
            if good:
                out.append((s, par))
    return out


TREES = my_trees()
assert len(TREES) == 375 and len({(s, tuple(sorted(p.items()))) for s, p in TREES}) == 375


def kirchhoff_count():
    # number of spanning trees of K5 - {A,B} by the matrix-tree theorem (exact)
    Lm = [[Fraction(0)] * 5 for _ in range(5)]
    for i, j in ADJ_PAIRS:
        Lm[i][j] -= 1; Lm[j][i] -= 1; Lm[i][i] += 1; Lm[j][j] += 1
    M = [row[1:] for row in Lm[1:]]
    n, det = 4, Fraction(1)
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]; det = -det
        det *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return det


assert kirchhoff_count() == 75


def tree_edge_list(tid):
    s, par = TREES[tid]
    return [(i, par[i]) for i in range(5) if i != s]


# ---------------------------------------------------------------- polynomials
def load_polys(chart):
    D = DATA[chart]
    conv = lambda terms: {tuple(m): c for m, c in terms}
    nu = {int(f): [conv(t) for t in D["nu"][f]] for f in D["nu"]}
    eta = {tuple(int(a) for a in k.split(",")): [conv(t) for t in D["eta"][k]] for k in D["eta"]}
    return nu, eta


POLYS = {ch: load_polys(ch) for ch in ("Z", "P")}
DMAX = {ch: max(sum(m) for vecs in (POLYS[ch][0].values(), POLYS[ch][1].values()) for vec in vecs for comp in vec
                for m in comp) for ch in ("Z", "P")}


def edge_vectors(chart, i, j):
    """(nu_i, eta_{i->j}) as lists of 3 polynomial dicts; eta_{i->j} = +eta_(i,j) if i<j else -eta_(j,i)"""
    nu, eta = POLYS[chart]
    if i < j:
        e = eta[(i, j)]
    else:
        e = [{m: -c for m, c in comp.items()} for comp in eta[(j, i)]]
    return nu[i], e


def eval_poly(pd, x):
    tot = Fraction(0)
    for m, c in pd.items():
        tot += c * x[0] ** m[0] * x[1] ** m[1] * x[2] ** m[2]
    return tot


def dyadic_exp(fr):
    d = fr.denominator
    e = d.bit_length() - 1
    assert d == 1 << e, "not dyadic"
    return e


class Box:
    """exact Taylor data of the chart polynomials on the box: x = c + h s; integer polys scaled by 2^(E*DMAX)."""

    def __init__(self, chart, lo, hi):
        self.chart, self.lo, self.hi = chart, lo, hi
        self.c = [(a + b) / 2 for a, b in zip(lo, hi)]
        self.h = [(b - a) / 2 for a, b in zip(lo, hi)]
        assert all(hk >= 0 for hk in self.h)
        self.E = max(dyadic_exp(t) for t in self.c + self.h)
        S = 1 << self.E
        self.C = [int(t * S) for t in self.c]
        self.H = [int(t * S) for t in self.h]
        assert all(Fraction(Ck, S) == ck for Ck, ck in zip(self.C, self.c))
        self.dmax = DMAX[chart]
        self.pow = []
        for k in range(3):
            lst = [[1]]
            for n in range(1, self.dmax + 1):
                prev = lst[-1]
                nxt = [0] * (n + 1)
                for a, coef in enumerate(prev):     # multiply by (C + H s)
                    nxt[a] += coef * self.C[k]
                    nxt[a + 1] += coef * self.H[k]
                lst.append(nxt)
            self.pow.append(lst)
        self.cache = {}

    def tpoly(self, pd):
        """2^(E*dmax) * p(c + h s) as integer dict over exponents of s"""
        out = {}
        S = self.E
        for m, a in pd.items():
            sc = a << (S * (self.dmax - sum(m)))
            P0, P1, P2 = self.pow[0][m[0]], self.pow[1][m[1]], self.pow[2][m[2]]
            for i, c0 in enumerate(P0):
                if c0 == 0:
                    continue
                for j, c1 in enumerate(P1):
                    if c1 == 0:
                        continue
                    c01 = sc * c0 * c1
                    for k, c2 in enumerate(P2):
                        if c2:
                            key = (i, j, k)
                            out[key] = out.get(key, 0) + c01 * c2
        return out

    def edge_C(self, i, j):
        """6 integer polys = 2 * 2^(2 E dmax) * sym(nu eta^T) entries (order xx,xy,xz,yy,yz,zz)"""
        key = (i, j)
        if key not in self.cache:
            nu, eta = edge_vectors(self.chart, i, j)
            a = [self.tpoly(t) for t in nu]
            b = [self.tpoly(t) for t in eta]
            ent = []
            for (pp, qq) in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)):
                ent.append(padd(pmul(a[pp], b[qq]), pmul(a[qq], b[pp])))
            self.cache[key] = ent
        return self.cache[key]


def pmul(A, B):
    out = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            k = (ma[0] + mb[0], ma[1] + mb[1], ma[2] + mb[2])
            out[k] = out.get(k, 0) + ca * cb
    return out


def padd(A, B):
    out = dict(A)
    for k, c in B.items():
        out[k] = out.get(k, 0) + c
    return out


def posdef3(M):
    """exact Sylvester on symmetric 3x3 given as (xx,xy,xz,yy,yz,zz) (ints or Fractions)"""
    xx, xy, xz, yy, yz, zz = M
    if not xx > 0:
        return False
    if not xx * yy - xy * xy > 0:
        return False
    det = xx * (yy * zz - yz * yz) - xy * (xy * zz - yz * xz) + xz * (xy * yz - yy * xz)
    return det > 0


def check_on_box(bx, tid, yc, Yh):
    """yc[i] = y_i at the box centre, Yh[i][k] = coefficient of s_k (i.e. Y1[i][k]*h_k); exact Fractions.
    True iff y >= 0 on the box and the (sufficient) exact test proves N negative definite on the box."""
    for i in range(4):
        if yc[i] - sum(abs(t) for t in Yh[i]) < 0:
            return False
    # integerise multipliers with a common dyadic denominator
    den = 1
    for t in yc + [t for row in Yh for t in row]:
        den = max(den, t.denominator)
    Y0 = [int(t * den) for t in yc]
    YK = [[int(t * den) for t in row] for row in Yh]
    assert all(Fraction(Y0[i], den) == yc[i] for i in range(4))
    edges = tree_edge_list(tid)
    Ntot = [dict() for _ in range(6)]
    for i, (a, b) in enumerate(edges):
        Ce = bx.edge_C(a, b)
        ypoly = {(0, 0, 0): Y0[i], (1, 0, 0): YK[i][0], (0, 1, 0): YK[i][1], (0, 0, 1): YK[i][2]}
        ypoly = {m: c for m, c in ypoly.items() if c}
        for m in range(6):
            Ntot[m] = padd(Ntot[m], pmul(ypoly, Ce[m]))
    N0 = [Ntot[m].get((0, 0, 0), 0) for m in range(6)]
    NK = [[Ntot[m].get(e, 0) for m in range(6)] for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    Eb = [sum(abs(c) for mono, c in Ntot[m].items() if sum(mono) >= 2) for m in range(6)]
    # rows of the symmetric bound matrix: (xx,xy,xz),(xy,yy,yz),(xz,yz,zz)
    rows = [Eb[0] + Eb[1] + Eb[2], Eb[1] + Eb[3] + Eb[4], Eb[2] + Eb[4] + Eb[5]]
    r = max(rows)
    for sg in itertools.product((-1, 1), repeat=3):
        M = [-(N0[m] + sg[0] * NK[0][m] + sg[1] * NK[1][m] + sg[2] * NK[2][m]) for m in range(6)]
        M[0] -= r; M[3] -= r; M[5] -= r
        if not posdef3(M):
            return False
    return True


def check_cert(chart, lo, hi, tid, y0, Y1, maxsub=3, bx=None):
    """exact check of the affine certificate y(x) = y0 + Y1 (x - c_node) on [lo,hi]; subdivides on failure.
    Returns (ok, number of sub-boxes used)."""
    c0 = [(a + b) / 2 for a, b in zip(lo, hi)]
    y0f = [Fraction(t) for t in y0]
    Y1f = [[Fraction(t) for t in row] for row in Y1]

    def rec(blo, bhi, depth, bxx):
        bxx = bxx or Box(chart, blo, bhi)
        yc = [y0f[i] + sum(Y1f[i][k] * (bxx.c[k] - c0[k]) for k in range(3)) for i in range(4)]
        Yh = [[Y1f[i][k] * bxx.h[k] for k in range(3)] for i in range(4)]
        if check_on_box(bxx, tid, yc, Yh):
            return True, 1
        if depth >= maxsub:
            return False, 1
        # y >= 0 failure cannot be cured by subdivision only if it fails at a corner; subdivide anyway
        mid = bxx.c
        tot = 0
        for corner in itertools.product((0, 1), repeat=3):
            sl = [blo[k] if corner[k] == 0 else mid[k] for k in range(3)]
            sh = [mid[k] if corner[k] == 0 else bhi[k] for k in range(3)]
            ok, n = rec(sl, sh, depth + 1, None)
            tot += n
            if not ok:
                return False, tot
        return True, tot

    return rec(lo, hi, 0, bx)


# ---------------------------------------------------------------- traversal of certificate files
def domain(chart, rho0):
    return ([Fraction(0)] * 3, [Fraction(1)] * 3) if chart == "Z" else ([Fraction(0)] * 3, [Fraction(rho0), Fraction(1), Fraction(1)])


def traverse(d, collect=True):
    """own traversal: exact boxes, own coverage logic.  Yields records and checks coverage.
    Returns (records, stats); record = (lo, hi, depth, cert, cellindex)."""
    chart, rho0, g = d["chart"], Fraction(d["rho0"]), d["grid"]
    dlo, dhi = domain(chart, d["rho0"])
    idxs = sorted(tuple(c["index"]) for c in d["cells"])
    assert idxs == sorted(itertools.product(range(g[0]), range(g[1]), range(g[2]))), "grid mismatch"
    recs, st = [], {"nodes": 0, "leaves": 0, "skip": 0, "certs": 0, "maxdepth": 0, "dup_tid_on_path": 0}
    for cell in d["cells"]:
        idx = cell["index"]
        lo = [dlo[k] + (dhi[k] - dlo[k]) * Fraction(idx[k], g[k]) for k in range(3)]
        hi = [dlo[k] + (dhi[k] - dlo[k]) * Fraction(idx[k] + 1, g[k]) for k in range(3)]
        stack = [(lo, hi, cell["node"], frozenset(range(375)), 0)]
        while stack:
            lo, hi, node, pending, depth = stack.pop()
            st["nodes"] += 1
            st["maxdepth"] = max(st["maxdepth"], depth)
            if node.get("skip"):
                assert chart == "Z" and lo[0] * lo[1] * lo[2] > 1 - rho0, "bad skip"
                st["skip"] += 1; st["leaves"] += 1
                continue
            assert "FAILED" not in node
            p = set(pending)
            for cert in node.get("certs", []):
                tid = cert[0]
                if tid not in p:
                    st["dup_tid_on_path"] += 1
                p.discard(tid)
                st["certs"] += 1
                if collect:
                    recs.append((lo, hi, depth, cert, tuple(idx)))
            if "children" in node:
                ax = node["axis"]
                assert ax in (0, 1, 2) and len(node["children"]) == 2
                mid = (lo[ax] + hi[ax]) / 2
                h1 = list(hi); h1[ax] = mid
                l2 = list(lo); l2[ax] = mid
                stack.append((lo, h1, node["children"][0], frozenset(p), depth + 1))
                stack.append((l2, hi, node["children"][1], frozenset(p), depth + 1))
            else:
                assert not p, f"leaf with {len(p)} pending trees"
                st["leaves"] += 1
    return recs, st


def _work(args):
    chart, lo, hi, tid, y0, Y1, tag = args
    t0 = time.time()
    ok, nsub = check_cert(chart, lo, hi, tid, y0, Y1)
    return tag, ok, nsub, time.time() - t0


def run_sample(chart, path, n_uniform, seed, workers=int(__import__("os").environ.get("AUDIT_WORKERS", "3"))):
    t0 = time.time()
    with gzip.open(path, "rt") as fh:
        d = json.load(fh)
    assert d["chart"] == chart
    recs, st = traverse(d)
    print(f"[{chart}] traversal (own code): {st}; {len(recs)} certificate records; {time.time() - t0:.0f}s", flush=True)
    rng = random.Random(seed)
    lo_dom, hi_dom = domain(chart, d["rho0"])
    maxd = st["maxdepth"]
    deep = [r for r in recs if r[2] >= maxd - 1]
    bnd = {(k, side): [r for r in recs if (r[0][k] == lo_dom[k] if side == 0 else r[1][k] == hi_dom[k])]
           for k in range(3) for side in (0, 1)}
    chosen = {}
    for r in rng.sample(recs, min(n_uniform, len(recs))):
        chosen[id(r)] = (r, "uniform")
    for r in rng.sample(deep, min(200, len(deep))):
        chosen.setdefault(id(r), (r, f"depth>={maxd - 1}"))
    for key, lst in bnd.items():
        for r in rng.sample(lst, min(60, len(lst))):
            chosen.setdefault(id(r), (r, f"boundary x{key[0]}={'lo' if key[1] == 0 else 'hi'}"))
    print(f"[{chart}] sampled {len(chosen)} certificates (deep pool {len(deep)}, boundary pools "
          f"{ {k: len(v) for k, v in bnd.items()} })", flush=True)
    tasks = [(chart, r[0], r[1], r[3][0], r[3][1], r[3][2], (tag, r[2], [str(t) for t in r[0]], [str(t) for t in r[1]], r[3][0]))
             for r, tag in chosen.values()]
    del d, recs
    from concurrent.futures import ProcessPoolExecutor
    res = {"ok": 0, "fail": 0, "subdivided": 0}
    bytag = {}
    fails = []
    tmax = 0
    with ProcessPoolExecutor(workers) as ex:
        for n, (tag, ok, nsub, dt) in enumerate(ex.map(_work, tasks, chunksize=4)):
            tmax = max(tmax, dt)
            key = tag[0].split(" ")[0]
            b = bytag.setdefault(key, [0, 0])
            if ok:
                res["ok"] += 1; b[0] += 1
                if nsub > 1:
                    res["subdivided"] += 1
            else:
                res["fail"] += 1; b[1] += 1
                fails.append(tag)
            if (n + 1) % 250 == 0:
                print(f"  {n + 1}/{len(tasks)} {res} [{time.time() - t0:.0f}s]", flush=True)
    print(f"[{chart}] RESULT {res}; by category (ok, fail): {bytag}; slowest cert {tmax:.2f}s; total {time.time() - t0:.0f}s")
    for f in fails[:20]:
        print("  FAIL:", f)


# ---------------------------------------------------------------- pointwise end-to-end check
def frac_det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]))


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


ORIENTED = [("B3", "B2", "A2", "A3"), ("B1", "B3", "A3", "A1"), ("B2", "B1", "A1", "A2"), ("A1", "A3", "A2"), ("B1", "B2", "B3")]


def verts_exact(chart, x):
    if chart == "Z":
        uu, vv, ww = x
        t = [uu * vv * ww, vv * ww, ww]
        B = [[Fraction(1), 0, 0], [0, Fraction(1), 0], [0, 0, Fraction(1)]]
        A = [[t[i] * B[i][k] for k in range(3)] for i in range(3)]
    else:
        r, pp, qq = x
        sig = [r, r * pp, r * pp * qq]
        s = [Fraction(1), pp, pp * qq]
        Bxy = [(0, 0), (1, 0), (0, 1)]
        G = (Fraction(1, 3), Fraction(1, 3))
        B = [[Fraction(a), Fraction(b), Fraction(0)] for a, b in Bxy]
        A = [[(1 - sig[i]) * Bxy[i][0] + sig[i] * G[0], (1 - sig[i]) * Bxy[i][1] + sig[i] * G[1], s[i]] for i in range(3)]
    return {"A1": A[0], "A2": A[1], "A3": A[2], "B1": B[0], "B2": B[1], "B3": B[2]}


def covectors_exact(V):
    vol, mom = Fraction(0), [Fraction(0)] * 3
    for f in ORIENTED:
        for k in range(1, len(f) - 1):
            a, b, c = V[f[0]], V[f[k]], V[f[k + 1]]
            dd = frac_det3(a, b, c)
            vol += dd
            mom = [m + dd * (a[i] + b[i] + c[i]) for i, m in enumerate(mom)]
    g = [m / (4 * vol) for m in mom]
    L = []
    for f in ORIENTED:
        P = [V[nm] for nm in f]
        n = cross(sub(P[1], P[0]), sub(P[2], P[0])) if len(P) == 3 else cross(sub(P[2], P[0]), sub(P[3], P[1]))
        dd = dot(n, P[0])
        assert all(dot(n, X) == dd for X in P)
        e = dd - dot(n, g)
        assert e != 0
        L.append([t / e for t in n])
    return L, g, vol


def pos_multiple(a, b):
    """exact: a = lam * b with lam > 0 ?  (b != 0)"""
    k = next(i for i in range(3) if b[i] != 0)
    lam = a[k] / b[k]
    return lam > 0 and all(a[i] == lam * b[i] for i in range(3))


def locate_and_check(d, chart, x):
    """walk down the cell tree containing x; for every certificate on the path, check y(x) >= 0 and
    N(x) < 0 exactly at x; return set of trees certified at x"""
    rho0 = Fraction(d["rho0"])
    g = d["grid"]
    dlo, dhi = domain(chart, d["rho0"])
    idx = []
    for k in range(3):
        j = int((x[k] - dlo[k]) / (dhi[k] - dlo[k]) * g[k])
        idx.append(min(j, g[k] - 1))
    cell = next(c for c in d["cells"] if c["index"] == idx)
    lo = [dlo[k] + (dhi[k] - dlo[k]) * Fraction(idx[k], g[k]) for k in range(3)]
    hi = [dlo[k] + (dhi[k] - dlo[k]) * Fraction(idx[k] + 1, g[k]) for k in range(3)]
    assert all(lo[k] <= x[k] <= hi[k] for k in range(3))
    nu, eta = POLYS[chart]
    nuv = {f: [eval_poly(t, x) for t in nu[f]] for f in nu}
    etav = {e: [eval_poly(t, x) for t in eta[e]] for e in eta}
    node, certified, bad = cell["node"], set(), []
    while True:
        if node.get("skip"):
            return None, bad, "skip"
        c = [(a + b) / 2 for a, b in zip(lo, hi)]
        for tid, y0, Y1 in node.get("certs", []):
            yx = [Fraction(y0[i]) + sum(Fraction(Y1[i][k]) * (x[k] - c[k]) for k in range(3)) for i in range(4)]
            Nm = [Fraction(0)] * 6
            for i, (a, b) in enumerate(tree_edge_list(tid)):
                n_ = nuv[a]
                e_ = etav[(a, b)] if a < b else [-t for t in etav[(b, a)]]
                for m, (pp, qq) in enumerate(((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))):
                    Nm[m] += yx[i] * (n_[pp] * e_[qq] + n_[qq] * e_[pp]) / 2
            if all(t >= 0 for t in yx) and posdef3([-t for t in Nm]):
                certified.add(tid)
            else:
                bad.append(tid)
        if "children" not in node:
            break
        ax = node["axis"]
        mid = (lo[ax] + hi[ax]) / 2
        if x[ax] <= mid:
            hi = list(hi); hi[ax] = mid; node = node["children"][0]
        else:
            lo = list(lo); lo[ax] = mid; node = node["children"][1]
    return certified, bad, (nuv, etav)


def run_points(pathZ, pathP, npts, seed):
    t0 = time.time()
    dZ = json.load(gzip.open(pathZ, "rt"))
    dP = json.load(gzip.open(pathP, "rt")) if pathP != "-" else None
    rho0 = Fraction(dZ["rho0"])
    rng = random.Random(seed)
    stats = {"Z": 0, "P_from_Z": 0, "P_direct": 0, "fully_certified": 0, "bad_cert_at_point": 0,
             "positivity_ok": 0, "skipped_no_P": 0}

    def rnd01(kind):
        if kind == 0:
            return Fraction(rng.randrange(1, 2 ** 20), 2 ** 20)
        if kind == 1:   # near 1
            return 1 - Fraction(rng.randrange(1, 2 ** 12), 2 ** 30)
        return Fraction(rng.randrange(1, 2 ** 12), 2 ** 30)  # near 0

    def check_point(chart, x, d):
        V = verts_exact(chart, x)
        L, gc, vol = covectors_exact(V)
        cert, bad, data = locate_and_check(d, chart, x)
        if cert is None:
            raise AssertionError(f"point {x} in chart {chart} lands in a skip leaf")
        nuv, etav = data
        okpos = all(pos_multiple(nuv[f], L[f]) for f in range(5)) and \
            all(pos_multiple(etav[(i, j)], sub(L[j], L[i])) for (i, j) in ADJ_PAIRS)
        stats["positivity_ok"] += okpos
        if not okpos:
            print("  POSITIVITY FAILURE at", chart, [float(t) for t in x])
        stats["bad_cert_at_point"] += len(bad)
        if bad:
            print("  certificate fails pointwise:", chart, [float(t) for t in x], bad[:5])
        if len(cert) == 375:
            stats["fully_certified"] += 1
        else:
            print("  NOT all trees certified at", chart, [float(t) for t in x], 375 - len(cert))

    for it in range(npts):
        # chart Z point (mixture of interior / near-boundary / exact boundary values incl. pyramids w=1)
        kinds = [rng.choice((0, 0, 1, 2)) for _ in range(3)]
        x = [rnd01(k) for k in kinds]
        if rng.random() < 0.15:
            x[2] = Fraction(1)              # pyramid
        if x[1] == 1 and x[2] == 1:
            continue
        if x[0] * x[1] * x[2] <= 1 - rho0:
            stats["Z"] += 1
            check_point("Z", x, dZ)
        else:
            if dP is None:
                stats["skipped_no_P"] += 1
                continue
            t1, t2, t3 = x[0] * x[1] * x[2], x[1] * x[2], x[2]
            y = [1 - t1, (1 - t2) / (1 - t1), (1 - t3) / (1 - t2)]
            assert 0 < y[0] < rho0 and 0 <= y[1] <= 1 and 0 <= y[2] <= 1
            # affine equivalence of the two frames: O=0 -> apex (G, 1/rho), e_i -> B_i^P
            VZ, VP = verts_exact("Z", x), verts_exact("P", y)
            apex = [Fraction(1, 3), Fraction(1, 3), 1 / y[0]]
            Bs = [VP["B1"], VP["B2"], VP["B3"]]
            T = lambda z: [apex[k] + sum(z[i] * (Bs[i][k] - apex[k]) for i in range(3)) for k in range(3)]
            assert all(T(VZ[nm]) == VP[nm] for nm in VZ), "frames not affinely equivalent"
            stats["P_from_Z"] += 1
            check_point("P", y, dP)
        # chart P point directly (incl. prisms rho=0 and pyramids q=0)
        if dP is not None:
            y = [Fraction(rng.randrange(0, 2 ** 20), 2 ** 20) * rho0, rnd01(rng.choice((0, 0, 1))), rnd01(rng.choice((0, 0, 1, 2)))]
            if rng.random() < 0.2:
                y[0] = Fraction(0)
            if rng.random() < 0.15:
                y[2] = Fraction(0)
            stats["P_direct"] += 1
            check_point("P", y, dP)
    print(f"pointwise end-to-end: {stats}  [{time.time() - t0:.0f}s]")


def _work_node(args):
    chart, lo, hi, certs = args
    bx = Box(chart, lo, hi)
    out = []
    for tid, y0, Y1 in certs:
        ok, nsub = check_cert(chart, lo, hi, tid, y0, Y1, bx=bx)
        if not ok or nsub > 1:
            out.append((tid, ok, nsub))
    return len(certs), out, [str(t) for t in lo], [str(t) for t in hi]


def run_all(chart, path, workers):
    t0 = time.time()
    with gzip.open(path, "rt") as fh:
        d = json.load(fh)
    assert d["chart"] == chart
    recs, st = traverse(d)
    print(f"[{chart}] traversal (own code): {st} [{time.time() - t0:.0f}s]", flush=True)
    nodes = {}
    for lo, hi, depth, cert, idx in recs:
        key = (tuple(lo), tuple(hi))
        nodes.setdefault(key, []).append(cert)
    tasks = [(chart, list(k[0]), list(k[1]), v) for k, v in nodes.items()]
    del recs, d
    from concurrent.futures import ProcessPoolExecutor
    done, nfail, nsubd, fails = 0, 0, 0, []
    with ProcessPoolExecutor(workers) as ex:
        for n, (cnt, out, lo, hi) in enumerate(ex.map(_work_node, tasks, chunksize=8)):
            done += cnt
            for tid, ok, nsub in out:
                if ok:
                    nsubd += 1
                else:
                    nfail += 1
                    fails.append((lo, hi, tid))
            if (n + 1) % 2000 == 0:
                print(f"  nodes {n + 1}/{len(tasks)}, certs {done}, fail {nfail}, subdivided {nsubd} "
                      f"[{time.time() - t0:.0f}s]", flush=True)
    print(f"[{chart}] ALL: {done} certificates on {len(tasks)} nodes checked exactly: fail {nfail}, "
          f"needed subdivision {nsubd}; total {time.time() - t0:.0f}s")
    for f in fails[:20]:
        print("  FAIL:", f)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "all":
        run_all(sys.argv[2], sys.argv[3], int(sys.argv[4]))
        sys.exit(0)
    if mode == "sample":
        run_sample(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
    elif mode == "points":
        run_points(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
    elif mode == "selftest":
        # quick timing on one certificate
        d = json.load(gzip.open(sys.argv[2], "rt"))
        recs, st = traverse(d)
        print(st)
        rng = random.Random(1)
        for r in rng.sample(recs, 5):
            t0 = time.time()
            print(r[2], [float(t) for t in r[0]], [float(t) for t in r[1]], r[3][0],
                  check_cert(d["chart"], r[0], r[1], r[3][0], r[3][1], r[3][2]), f"{time.time() - t0:.2f}s")
