"""
Infeasibility certificates (Theorem 8 / SDP (3.3) of Bozoki-Domokos-Papp-Regos,
arXiv:2401.17906) for the quadratic systems of their Algorithm 2.

System for a case (graph G, candidate unstable vertex v, rooted spanning tree
with parent map j, face F not containing v); unknowns r_i in R^3, i != v,
with r_v = (1,0,0) substituted:

    (shadow)  |r_i|^2 - r_i . r_j(i)            <= 0     for every i != v
    (CoM)     sum_{k in F u {v}} r_{k,1}         <= 0     (= 1 + sum_{k in F} r_{k,1})

Optional valid (implied, hence redundant for feasibility, but possibly useful
for the relaxation) inequalities:

  level "b":  |r_i|^2 - |r_v|^2 <= 0  and  (r_i - r_v) . r_v <= 0   (i != v)
  level "bc": additionally, for every i and every proper ancestor a != v of i
              in the tree:  |r_i|^2 - |r_a|^2 <= 0  and  r_i . r_a - |r_a|^2 <= 0
              (consequences of Lemma 6 / Cauchy-Schwarz along tree paths;
              for a = v these are exactly the level-"b" inequalities).

Every inequality has the form  c + sum_{p,q} Gm[p,q] r_p.r_q + sum_p h[p] r_{p,1} <= 0
with a symmetric V x V "Gram coefficient" matrix Gm.  After substituting
r_v = e_1 its quadratic part is (A kron I_3) in the variables
(r_{i,1}, r_{i,2}, r_{i,3})_{i != v}, and only first coordinates occur linearly.

DECOUPLING.  Order the 3(V-1) variables as (all first coordinates, all second
coordinates, all third coordinates).  Then for every y >= 0

    sum_k y_k Qhat_k  =  blockdiag( M(y), A(y), A(y) ),
    M(y) = sum_k y_k [[a_k, bx_k^T/2], [bx_k/2, A_k]]   (size V),
    A(y) = sum_k y_k A_k                                  (size V-1),

and A(y) is the trailing principal submatrix of M(y).  By Cauchy interlacing
lambda_min(A(y)) >= lambda_min(M(y)), so lambda_min(full) = lambda_min(M(y)):
the constraint  sum y_k Qhat_k >= z I  of (3.3) is *equivalent* to M(y) >= z I.
Hence the SDP can be solved with V x V matrices with the same optimal z and y.
(`self_test()` re-derives the full matrices symbolically with sympy and checks
this identity, and compares the reduced and the full SDP optimal values.)

NORMALIZATION.  SDP (3.3) of the paper:  max z  s.t.  sum y_k Qhat_k >= z I,
||y||_2 <= 1, y_k >= z.  The authors' Mathematica notebook (public repository
dpapp-github/mono-unstable, homogeneous/) builds the matrices from
QArray[2*quad], i.e. it uses 2*Qhat_k.  `scale` selects the factor (default 2,
i.e. the notebook's normalization; scale=1 is (3.3) as printed).  The sign of
z (certified or not) does not depend on the scale:
z_1 <= z_2 <= 2 z_1 for the optimal values with scale 1 and 2.

CERTIFICATES are verified in exact arithmetic: y is rounded to a dyadic grid
(multiples of 2^-40), and the integer matrix 2^41 * sum_k y_k M_k is checked
to be positive definite by fraction-free (Bareiss) elimination, i.e. all
leading principal minors > 0.  By the decoupling identity this proves
positive definiteness of sum_k y_k Qhat_k, hence (Theorem 8) infeasibility.

Usage:
    python certify.py selftest
    python certify.py repro            # V = 5, 6 + comparison with cases56.csv
    python certify.py run 7 --workers 4 [--level base|b|bc] [--out file]
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time

import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import enumerate_cases as ec  # noqa: E402

LEVELS = ("base", "b", "bc")
DYADIC_BITS = 40


# ---------------------------------------------------------------------------
# Building the systems
# ---------------------------------------------------------------------------

class Ineq:
    """c + sum Gm[p,q] r_p.r_q + sum h[p] r_{p,1} <= 0  (V-vertex form)."""
    __slots__ = ("label", "c", "Gm", "h")

    def __init__(self, label, V):
        self.label = label
        self.c = 0.0
        self.Gm = np.zeros((V, V))
        self.h = np.zeros(V)

    def dot(self, p, q, coef=1.0):
        """add coef * r_p . r_q"""
        self.Gm[p, q] += coef / 2
        self.Gm[q, p] += coef / 2
        return self


def ancestors(parent, i, v):
    """Proper ancestors of i in the tree rooted at v (parent, grandparent, ...,
    v), in order."""
    out = []
    w = parent[i]
    while w != -1:
        out.append(w)
        if w == v:
            break
        w = parent[w]
    return out


def build_inequalities(V, v, parent, face, level="base"):
    """List of Ineq objects for one case.  The shadow inequalities come first
    (non-root vertices in increasing label order), then the CoM inequality,
    then (optionally) the strengthening inequalities."""
    assert level in LEVELS
    others = [i for i in range(V) if i != v]
    ineqs = []
    for i in others:
        j = parent[i]
        ineqs.append(Ineq(f"shadow({i}<-{j})", V).dot(i, i).dot(i, j, -1.0))
    com = Ineq(f"com({v};{face[0]},{face[1]},{face[2]})", V)
    for k in (v, *face):
        com.h[k] += 1.0
    ineqs.append(com)
    if level in ("b", "bc"):
        for i in others:
            ineqs.append(Ineq(f"norm({i}<={v})", V).dot(i, i).dot(v, v, -1.0))
            ineqs.append(Ineq(f"dot({i},{v})<=|{v}|^2", V).dot(i, v).dot(v, v, -1.0))
    if level == "bc":
        for i in others:
            for a in ancestors(parent, i, v):
                if a == v:
                    continue
                ineqs.append(Ineq(f"norm({i}<={a})", V).dot(i, i).dot(a, a, -1.0))
                ineqs.append(Ineq(f"dot({i},{a})<=|{a}|^2", V).dot(i, a).dot(a, a, -1.0))
    return ineqs


def reduced_matrices(ineqs, V, v):
    """Homogenized reduced ("first-coordinate") matrices M_k of size V:
    index 0 = constant 1, indices 1..V-1 = first coordinates of the
    non-root vertices in increasing label order.  Substitutes r_v = e_1."""
    others = [i for i in range(V) if i != v]
    m = len(ineqs)
    M = np.zeros((m, V, V))
    for k, q in enumerate(ineqs):
        a = q.c + q.Gm[v, v] + q.h[v]
        bx = q.h[others] + 2 * q.Gm[v, others]
        A = q.Gm[np.ix_(others, others)]
        M[k, 0, 0] = a
        M[k, 0, 1:] = bx / 2
        M[k, 1:, 0] = bx / 2
        M[k, 1:, 1:] = A
    return M


def full_matrices(ineqs, V, v):
    """The paper's full homogenized matrices Qhat_k of size 1+3(V-1), with
    variables ordered vertex-major: (r_{i,1}, r_{i,2}, r_{i,3}) for the
    non-root vertices i in increasing label order (as in the notebook)."""
    Mred = reduced_matrices(ineqs, V, v)
    m, n = Mred.shape[0], V - 1
    Q = np.zeros((m, 1 + 3 * n, 1 + 3 * n))
    for k in range(m):
        Q[k, 0, 0] = Mred[k, 0, 0]
        Q[k, 0, 1::3] = Mred[k, 0, 1:]       # linear terms: first coords only
        Q[k, 1::3, 0] = Mred[k, 1:, 0]
        Q[k, 1:, 1:] = np.kron(Mred[k, 1:, 1:], np.eye(3))
    return Q


def system_matrices(V, v, parent, face, level="base"):
    ineqs = build_inequalities(V, v, parent, face, level)
    return reduced_matrices(ineqs, V, v), [q.label for q in ineqs]


# ---------------------------------------------------------------------------
# The SDP (3.3), solved directly with Clarabel
# ---------------------------------------------------------------------------

def _svec_index(d):
    """Clarabel's PSDTriangleConeT: upper triangle, column-major, off-diagonal
    entries scaled by sqrt(2)."""
    rows, cols, scal = [], [], []
    for j in range(d):
        for i in range(j + 1):
            rows.append(i)
            cols.append(j)
            scal.append(1.0 if i == j else math.sqrt(2.0))
    return np.array(rows), np.array(cols), np.array(scal)


_SVEC_CACHE = {}


def solve_sdp(M, scale=2.0, tol=1e-9):
    """max z s.t. sum_k y_k (scale*M_k) >= z I, ||y||_2 <= 1, y_k >= z.
    M has shape (m, d, d).  Returns (z, y, status_string)."""
    import clarabel
    m, d, _ = M.shape
    if d not in _SVEC_CACHE:
        _SVEC_CACHE[d] = _svec_index(d)
    ri, ci, sc = _SVEC_CACHE[d]
    nt = len(ri)
    n = m + 1                                   # x = (z, y_1..y_m)
    # --- constraint matrix A and vector b, s = b - A x in cone ---
    blocks = []
    # (1) y_k - z >= 0
    A1 = np.zeros((m, n))
    A1[:, 0] = 1.0
    A1[np.arange(m), 1 + np.arange(m)] = -1.0
    b1 = np.zeros(m)
    # (2) (1, y) in SOC
    A2 = np.zeros((m + 1, n))
    A2[1 + np.arange(m), 1 + np.arange(m)] = -1.0
    b2 = np.zeros(m + 1)
    b2[0] = 1.0
    # (3) svec(sum y_k scale M_k - z I) in PSD cone
    A3 = np.zeros((nt, n))
    A3[:, 0] = (ri == ci).astype(float)          # svec(I)
    A3[:, 1:] = -(scale * M[:, ri, ci] * sc).T
    b3 = np.zeros(nt)
    A = sp.csc_matrix(np.vstack([A1, A2, A3]))
    b = np.concatenate([b1, b2, b3])
    P = sp.csc_matrix((n, n))
    q = np.zeros(n)
    q[0] = -1.0
    cones = [clarabel.NonnegativeConeT(m), clarabel.SecondOrderConeT(m + 1),
             clarabel.PSDTriangleConeT(d)]
    st = clarabel.DefaultSettings()
    st.verbose = False
    st.tol_gap_abs = tol
    st.tol_gap_rel = tol
    st.tol_feas = tol
    st.max_iter = 200
    solver = clarabel.DefaultSolver(P, q, A, b, cones, st)
    sol = solver.solve()
    x = np.array(sol.x)
    return float(x[0]), x[1:], str(sol.status)


# ---------------------------------------------------------------------------
# Exact verification of a certificate
# ---------------------------------------------------------------------------

def _bareiss_pd(N):
    """N: list of lists of Python ints, symmetric.  True iff positive
    definite, via fraction-free Gaussian elimination without pivoting:
    the k-th pivot equals the k-th leading principal minor."""
    N = [row[:] for row in N]
    d = len(N)
    prev = 1
    for k in range(d):
        if N[k][k] <= 0:
            return False
        for i in range(k + 1, d):
            for j in range(k + 1, d):
                N[i][j] = (N[i][j] * N[k][k] - N[i][k] * N[k][j]) // prev
        prev = N[k][k]
    return True


def round_certificate(y):
    """Round y to the dyadic grid 2^-DYADIC_BITS, clipping at 0."""
    Y = np.maximum(np.round(np.asarray(y) * 2.0 ** DYADIC_BITS), 0.0)
    return [int(t) for t in Y]


def verify_certificate(M, y):
    """Exact check that the rounded y (>= 0 by construction) makes
    sum y_k M_k positive definite.  Entries of M_k are multiples of 1/2,
    so 2*M_k is an integer matrix.  Returns (ok, rounded_y_as_floats)."""
    Y = round_certificate(y)
    M2 = np.rint(2 * M).astype(np.int64)
    assert np.array_equal(M2, 2 * M), "matrix entries must be multiples of 1/2"
    d = M.shape[1]
    N = [[sum(Y[k] * int(M2[k, i, j]) for k in range(len(Y)) if M2[k, i, j])
          for j in range(d)] for i in range(d)]
    ok = _bareiss_pd(N)
    return ok, [t / 2.0 ** DYADIC_BITS for t in Y]


# ---------------------------------------------------------------------------
# Running all cases
# ---------------------------------------------------------------------------

def solve_case(V, v, parent, face, level="base", scale=2.0):
    M, labels = system_matrices(V, v, parent, face, level)
    z, y, status = solve_sdp(M, scale=scale)
    ok, yr = verify_certificate(M, y)
    lam = float(np.linalg.eigvalsh(np.tensordot(y, M, axes=1))[0])
    return dict(z=z, y=[float(t) for t in y], status=status,
                certified=bool(ok), lam_min=lam, n_ineq=len(labels))


def _worker(args):
    V, gid, v, trees, faces, level, scale = args
    out = []
    for parent in trees:
        for face in faces:
            r = solve_case(V, v, parent, face, level, scale)
            r.update(graph=gid, root=v, parent=list(parent), face=list(face))
            out.append(r)
    return out


def run_all(V, level="base", scale=2.0, workers=4, gds=None, cases=None,
            chunk=200, verbose=True):
    """Solve every system (or the given list of cases (gid, v, parent, face)).
    Returns the list of result dicts in deterministic order."""
    if gds is None:
        gds = ec.graph_data(V)
    tasks = []
    if cases is None:
        for gd in gds:
            for v in gd.roots:
                T = gd.trees[v]
                for s in range(0, len(T), chunk):
                    tasks.append((V, gd.gid, v, T[s:s + chunk],
                                  gd.decomposition[v], level, scale))
    else:
        for (gid, v, parent, face) in cases:
            tasks.append((V, gid, v, [tuple(parent)], [tuple(face)], level, scale))
    t0 = time.time()
    results = []
    workers = min(workers, 4)          # project policy: at most 4 cores
    if workers > 1:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(workers) as pool:
            for k, res in enumerate(pool.imap(_worker, tasks, chunksize=1 if cases is None else 50)):
                results.extend(res)
                if verbose and cases is None and (k + 1) % 50 == 0:
                    print(f"  {len(results)} systems solved, {time.time() - t0:.0f}s",
                          flush=True)
    else:
        for task in tasks:
            results.extend(_worker(task))
    if verbose:
        print(f"  V={V} level={level}: {len(results)} systems in {time.time() - t0:.1f}s")
    return results


def write_csv(results, path):
    """One row per system.  Columns:
      graph        graph id (index into enumerate_cases.triangulations(V))
      root         candidate unstable vertex v (0-based)
      parent       parent array of the rooted spanning tree (-1 at the root)
      tetrahedron  v followed by the face F (the tetrahedron conv(v, F))
      z            optimal value of SDP (3.3) with matrices scale*Qhat (default
                   scale 2 = the authors' notebook normalization)
      certified    1 iff the (dyadically rounded) y passed the exact check
      lam_min      smallest eigenvalue of sum_k y_k Qhat_k (scale 1, floats)
      status       Clarabel status
      y            multipliers: shadow inequalities of the non-root vertices in
                   increasing label order, then the centroid inequality
    verify_certificates.py re-checks the stored y exactly."""
    cols = ["graph", "root", "parent", "tetrahedron", "z", "certified",
            "lam_min", "status", "y"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in results:
            w.writerow([r["graph"], r["root"], " ".join(map(str, r["parent"])),
                        " ".join(map(str, [r["root"]] + list(r["face"]))),
                        f"{r['z']:.10g}", int(r["certified"]),
                        f"{r['lam_min']:.6g}", r["status"],
                        " ".join(f"{t:.10g}" for t in r["y"])])


def summarize(results, label=""):
    n = len(results)
    cert = sum(r["certified"] for r in results)
    zs = np.array([r["z"] for r in results])
    zc = zs[[r["certified"] for r in results]] if cert else np.array([np.nan])
    print(f"{label}: {n} systems, {cert} certified (exact check), {n - cert} not; "
          f"min z (certified) = {np.min(zc):.6g}, min z overall = {zs.min():.3g}, "
          f"#(z > 1e-6) = {(zs > 1e-6).sum()}")
    return n, cert


# ---------------------------------------------------------------------------
# Reproduction (V = 5, 6) and comparison with the authors' case list
# ---------------------------------------------------------------------------

def parse_papp_cases(path):
    """Cases from cases56.csv of the authors' repository.  Line format (1-based
    labels): u, Vs (V-1 entries), j(i) for i in Vs, T (u followed by face)."""
    out = []
    with open(path) as f:
        for line in f:
            a = [int(t) for t in line.strip().split(",") if t]
            V = (len(a) - 3) // 2
            u = a[0] - 1
            Vs = [t - 1 for t in a[1:V]]
            s = [t - 1 for t in a[V:2 * V - 1]]
            T = [t - 1 for t in a[2 * V - 1:2 * V + 3]]
            assert T[0] == u and len(T) == 4
            parent = [-1] * V
            for i, j in zip(Vs, s):
                parent[i] = j
            out.append((V, u, tuple(parent), tuple(T[1:])))
    return out


def reproduce(papp_csv=None, workers=4):
    report = {}
    for V in (5, 6):
        gds = ec.graph_data(V)
        for scale in (2.0, 1.0):
            res = run_all(V, "base", scale, workers=workers, gds=gds, verbose=False)
            n, c = summarize(res, f"V={V}, scale={scale:g} (ours)")
            zs = sorted(r["z"] for r in res)
            report[(V, scale)] = zs
            if scale == 2.0:
                write_csv(res, os.path.join(HERE, f"results_V{V}.csv"))
    if papp_csv and os.path.exists(papp_csv):
        pc = parse_papp_cases(papp_csv)
        print(f"authors' case list: {len(pc)} cases "
              f"({sum(1 for c in pc if c[0] == 5)} with V=5, "
              f"{sum(1 for c in pc if c[0] == 6)} with V=6)")
        for scale in (2.0, 1.0):
            rs = [solve_case(V, u, p, F, "base", scale) for (V, u, p, F) in pc]
            for V in (5, 6):
                sub = [r for r, c in zip(rs, pc) if c[0] == V]
                zs = sorted(r["z"] for r in sub)
                diff = np.max(np.abs(np.array(zs) - np.array(report[(V, scale)])))
                print(f"  V={V} scale={scale:g}: authors' cases: {len(sub)}, "
                      f"certified {sum(r['certified'] for r in sub)}, "
                      f"min z = {min(zs):.6g}; sorted-z max deviation from our "
                      f"enumeration = {diff:.2e}")
            if scale == 2.0:
                print("  first 3 / last 3 cases (z, y), scale 2 [compare with the "
                      "notebook output Out[16]]:")
                for k in (0, 1, 2, len(pc) - 3, len(pc) - 2, len(pc) - 1):
                    print("   ", k + 1, f"{rs[k]['z']:.8f}",
                          [round(t, 8) for t in rs[k]["y"]])
                print(f"  min z over all 5943 authors' cases (scale 2): "
                      f"{min(r['z'] for r in rs):.10g}  [notebook Out[17]: "
                      f"0.007784953778603651]")
    return report


# ---------------------------------------------------------------------------
# Self tests: decoupling and Clarabel formulation
# ---------------------------------------------------------------------------

def _sympy_full_matrices(V, v, parent, face):
    """Independent construction of the paper's Qhat_k directly from the
    polynomials (sympy), vertex-major variable order."""
    import sympy
    others = [i for i in range(V) if i != v]
    r = {i: sympy.symbols(f"r{i}_1:4") for i in others}
    r[v] = (sympy.Integer(1), sympy.Integer(0), sympy.Integer(0))
    dot = lambda p, q: sum(a * b for a, b in zip(r[p], r[q]))
    polys = [dot(i, i) - dot(i, parent[i]) for i in others]
    polys.append(sum(r[k][0] for k in (v, *face)))
    xs = [x for i in others for x in r[i]]
    n = len(xs)
    out = []
    for p in polys:
        P = sympy.Poly(sympy.expand(p), *xs)
        Q = sympy.zeros(n + 1, n + 1)
        for mon, coef in P.terms():
            idx = [k for k, e in enumerate(mon) for _ in range(e)]
            if len(idx) == 0:
                Q[0, 0] += coef
            elif len(idx) == 1:
                Q[0, idx[0] + 1] += coef / 2
                Q[idx[0] + 1, 0] += coef / 2
            else:
                a, b = idx
                if a == b:
                    Q[a + 1, a + 1] += coef
                else:
                    Q[a + 1, b + 1] += coef / 2
                    Q[b + 1, a + 1] += coef / 2
        out.append(np.array(Q.tolist(), dtype=float))
    return np.array(out)


def solve_sdp_cvxpy(M, scale=2.0):
    import cvxpy as cp
    m, d, _ = M.shape
    y = cp.Variable(m)
    z = cp.Variable()
    S = sum(y[k] * (scale * M[k]) for k in range(m))
    cons = [S - z * np.eye(d) >> 0, cp.norm(y, 2) <= 1, y >= z]
    prob = cp.Problem(cp.Maximize(z), cons)
    prob.solve(solver="CLARABEL")
    return float(z.value), np.array(y.value)


def self_test(n_random=30, seed=1):
    rng = np.random.default_rng(seed)
    print("self test: decoupling identity + Clarabel vs cvxpy (full matrices)")
    worst = 0.0
    for V in (5, 6, 7):
        gds = ec.graph_data(V)
        for t in range(n_random // 3):
            gd = gds[rng.integers(len(gds))]
            v = gd.roots[rng.integers(len(gd.roots))]
            parent = gd.trees[v][rng.integers(len(gd.trees[v]))]
            face = gd.decomposition[v][rng.integers(len(gd.decomposition[v]))]
            ineqs = build_inequalities(V, v, parent, face, "base")
            Qs = _sympy_full_matrices(V, v, parent, face)
            Qf = full_matrices(ineqs, V, v)
            assert np.allclose(Qs, Qf), "full matrices disagree with sympy"
            Mred = reduced_matrices(ineqs, V, v)
            # lambda_min identity for random nonnegative y
            for _ in range(5):
                y = rng.random(len(ineqs))
                l_full = np.linalg.eigvalsh(np.tensordot(y, Qf, 1))[0]
                l_red = np.linalg.eigvalsh(np.tensordot(y, Mred, 1))[0]
                assert abs(l_full - l_red) < 1e-10, (l_full, l_red)
            z_red, y_red, _ = solve_sdp(Mred, 2.0)
            if V <= 6:
                z_full, y_full = solve_sdp_cvxpy(Qf, 2.0)
                worst = max(worst, abs(z_full - z_red))
            z_cvx, _ = solve_sdp_cvxpy(Mred, 2.0)
            worst = max(worst, abs(z_cvx - z_red))
            # strengthened systems also decouple
            for lev in ("b", "bc"):
                iq = build_inequalities(V, v, parent, face, lev)
                Qf2 = full_matrices(iq, V, v)
                M2 = reduced_matrices(iq, V, v)
                y = rng.random(len(iq))
                assert abs(np.linalg.eigvalsh(np.tensordot(y, Qf2, 1))[0]
                           - np.linalg.eigvalsh(np.tensordot(y, M2, 1))[0]) < 1e-10
                # each inequality is valid: check on a random genuine
                # configuration satisfying the shadow constraints exactly
    print(f"  OK; max |z(reduced, Clarabel direct) - z(full or cvxpy)| = {worst:.2e}")
    # validity of the strengthening inequalities on random tree-feasible points
    _check_strengthening_validity(rng)


def _check_strengthening_validity(rng, trials=2000):
    """Sample random points satisfying all shadow inequalities of a random
    tree (sampling r_i inside the Thales ball of its parent) and check that the
    level-bc inequalities hold (they must: they are implied)."""
    gds = ec.graph_data(7)
    worst = -np.inf
    for _ in range(trials):
        gd = gds[rng.integers(len(gds))]
        v = gd.roots[rng.integers(len(gd.roots))]
        parent = gd.trees[v][rng.integers(len(gd.trees[v]))]
        V = gd.V
        r = np.zeros((V, 3))
        r[v] = (1, 0, 0)
        done = {v}
        while len(done) < V:
            for i in range(V):
                if i not in done and parent[i] in done:
                    c = r[parent[i]] / 2
                    rad = np.linalg.norm(c)
                    u = rng.normal(size=3)
                    u *= rng.random() ** (1 / 3) / np.linalg.norm(u)
                    r[i] = c + rad * u
                    done.add(i)
        for q in build_inequalities(V, v, parent, gd.decomposition[v][0], "bc"):
            if q.label.startswith("com"):
                continue
            val = q.c + np.einsum("pq,pk,qk->", q.Gm, r, r) + q.h @ r[:, 0]
            worst = max(worst, val)
    print(f"  strengthening inequalities valid on {trials} random tree-feasible "
          f"configurations: max value = {worst:.2e} (must be <= 0)")
    assert worst <= 1e-12


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["selftest", "repro", "run"])
    ap.add_argument("V", nargs="?", type=int, default=7)
    ap.add_argument("--level", default="base", choices=LEVELS)
    ap.add_argument("--scale", type=float, default=2.0)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default=None)
    ap.add_argument("--papp", default=None, help="path to the authors' cases56.csv")
    a = ap.parse_args()
    if a.cmd == "selftest":
        self_test()
    elif a.cmd == "repro":
        reproduce(a.papp, workers=a.workers)
    else:
        res = run_all(a.V, a.level, a.scale, workers=min(a.workers, 4))
        summarize(res, f"V={a.V} level={a.level}")
        out = a.out or os.path.join(HERE, f"results_V{a.V}.csv")
        write_csv(res, out)
        print("wrote", out)


if __name__ == "__main__":
    main()
