"""
Joint-direction ("mu") certificates for the mono-unstable lower bound.

Setting (see the project notes): P a convex polyhedron with V vertices, uniform
density, centroid at the origin, Gram matrix R = (r_a . r_b) (PSD, nonzero).
If v is the only unstable vertex there is a spanning tree of the (triangulated)
graph G rooted at v with  <S_{i->p(i)}, R> <= 0  for all i != v, where
S_{i->j} = E_ii - (E_ij + E_ji)/2.  The cone decomposition from v into
tetrahedra T = conv(v, F), F a face of G with v not in F, satisfies
sum_T vol(T) * sum_{k in T} r_k = 0.

CERTIFICATE for (G, v, tree): mu in R^V and, for every tetrahedron T of the cone
decomposition from v, y^T in R^V with y^T_i >= 0 (i != v), y^T_v = 0, such that

    M_T = sum_{i != v} y^T_i S_{i->p(i)} + (mu 1_T^T + 1_T mu^T)/2   is positive definite.

(Proof of validity: w = sum_k mu_k r_k; sum_T vol(T) w.(sum_{k in T} r_k) = 0 so
some T* with vol > 0 has <sym(mu 1_T*^T), R> <= 0; adding the shadow inequalities
gives <M_T*, R> <= 0, contradicting M_T* > 0, R >= 0, R != 0.)

We find (mu, y) with ONE SDP per tree:
    max z  s.t.  M_T >= z I  (all T),  y^T >= 0,  sum_i y^T_i <= 1,  |mu|_inf <= 1,
solved by calling Clarabel directly (no modelling layer), then ROUND to integers:
mu_int = round(D mu), y_int = max(0, round(D y)) for a scale D, and check that
every integer matrix 2*M_T(mu_int, y_int) is positive definite EXACTLY
(fraction-free Bareiss elimination: all leading principal minors > 0).
Positive scaling does not affect positive definiteness, so the stored integer
vectors are themselves a certificate.  D is increased until the check passes.

Output: gzip JSON lines,  lb_joint/certs_V{V}.jsonl.gz
  line 1  : {"kind":"header", "V":..., ...}
  graph   : {"kind":"graph", "gid", "edges", "faces", "roots",
             "aut": {u: perm}}   perm[x] = image of x, an automorphism of G with
                                 perm[u] = (orbit representative root of u),
                                 given for every non-root vertex u
  cert    : {"kind":"cert", "gid", "root", "parent", "z", "mu": [ints],
             "Y": [[[a,b,c], [y_0..y_{V-1}]], ...]}   one entry per face not
                                 containing root, y_root = 0
Uncertified trees (if any) are listed in lb_joint/uncertified_V{V}.json.

Usage:  python certify_joint.py run V [--workers 4]
        python certify_joint.py crosscheck      (Clarabel-direct vs cvxpy prototype)
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import sys
import time

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "lowerbound"))
import enumerate_cases as ec  # noqa: E402


# ---------------------------------------------------------------------------
# SDP via Clarabel directly
# ---------------------------------------------------------------------------

_SVEC = {}


def svec_pos(V):
    """Clarabel PSDTriangleConeT layout: upper triangle, column-major,
    off-diagonals scaled by sqrt(2).  Returns dict (a,b)->row (a<=b), size."""
    if V not in _SVEC:
        pos, k = {}, 0
        for j in range(V):
            for i in range(j + 1):
                pos[(i, j)] = k
                k += 1
        _SVEC[V] = (pos, k)
    return _SVEC[V]


def tet_indicator(V, v, face):
    one = np.zeros(V)
    one[list(face)] = 1.0
    one[v] = 1.0
    return one


def solve_joint(V, v, parent, tets, tol=1e-9, max_iter=300, extra=None):
    """max z s.t. M_T >= zI for all T in tets, y^T >= 0, sum y^T <= 1,
    |mu|_inf <= 1.   Returns (z, mu (V,), Y (nT, V) with Y[:, v] = 0, status).
    `extra`: optional list of additional symmetric V x V matrices E_l that may
    be added to every M_T with multipliers >= 0 (valid inequalities <E_l,R> <= 0);
    used only by the strengthening experiments (not for the main certificates)."""
    import clarabel
    pos, nt = svec_pos(V)
    r2 = math.sqrt(2.0)
    others = [i for i in range(V) if i != v]
    m = len(others)
    nT = len(tets)
    ne = 0 if extra is None else len(extra)
    nloc = m + ne                       # multipliers per T
    nx = 1 + V + nT * nloc
    # svec of S_i for the tree, and of the extra matrices
    Ssv = np.zeros((nt, m))
    for c, i in enumerate(others):
        j = parent[i]
        Ssv[pos[(i, i)], c] += 1.0
        a, b = min(i, j), max(i, j)
        Ssv[pos[(a, b)], c] += -0.5 * r2
    if ne:
        Esv = np.zeros((nt, ne))
        for l, E in enumerate(extra):
            for (a, b), k in pos.items():
                Esv[k, l] = E[a, b] * (1.0 if a == b else r2)
    Isv = np.zeros(nt)
    for i in range(V):
        Isv[pos[(i, i)]] = 1.0
    rows_A, rows_b, cones = [], [], []
    import clarabel as cl
    # (1) multipliers >= 0
    A1 = np.zeros((nT * nloc, nx))
    A1[np.arange(nT * nloc), 1 + V + np.arange(nT * nloc)] = -1.0
    rows_A.append(A1); rows_b.append(np.zeros(nT * nloc))
    # (2) sum of shadow multipliers per T <= 1 (extra multipliers also counted)
    A2 = np.zeros((nT, nx))
    for k in range(nT):
        A2[k, 1 + V + k * nloc: 1 + V + (k + 1) * nloc] = 1.0
    rows_A.append(A2); rows_b.append(np.ones(nT))
    # (3) -1 <= mu <= 1
    A3 = np.zeros((2 * V, nx))
    A3[np.arange(V), 1 + np.arange(V)] = 1.0
    A3[V + np.arange(V), 1 + np.arange(V)] = -1.0
    rows_A.append(A3); rows_b.append(np.ones(2 * V))
    cones.append(cl.NonnegativeConeT(nT * nloc + nT + 2 * V))
    # (4) PSD: svec(M_T - zI) = b - A x with b = 0
    for k, f in enumerate(tets):
        one = tet_indicator(V, v, f)
        Ak = np.zeros((nt, nx))
        Ak[:, 0] = Isv
        # sym(mu 1^T): entry (a,b) = (mu_a one_b + mu_b one_a)/2
        for (a, b), q in pos.items():
            if a == b:
                Ak[q, 1 + a] -= one[a]
            else:
                Ak[q, 1 + a] -= r2 * one[b] / 2
                Ak[q, 1 + b] -= r2 * one[a] / 2
        Ak[:, 1 + V + k * nloc: 1 + V + k * nloc + m] = -Ssv
        if ne:
            Ak[:, 1 + V + k * nloc + m: 1 + V + (k + 1) * nloc] = -Esv
        rows_A.append(Ak); rows_b.append(np.zeros(nt))
        cones.append(cl.PSDTriangleConeT(V))
    A = sp.csc_matrix(np.vstack(rows_A))
    b = np.concatenate(rows_b)
    P = sp.csc_matrix((nx, nx))
    q = np.zeros(nx); q[0] = -1.0
    st = cl.DefaultSettings()
    st.verbose = False
    st.tol_gap_abs = tol; st.tol_gap_rel = tol; st.tol_feas = tol
    st.max_iter = max_iter
    sol = cl.DefaultSolver(P, q, A, b, cones, st).solve()
    x = np.array(sol.x)
    z = float(x[0]); mu = x[1:1 + V]
    Y = np.zeros((nT, V))
    Yx = np.zeros((nT, ne))
    for k in range(nT):
        blk = x[1 + V + k * nloc: 1 + V + (k + 1) * nloc]
        Y[k, others] = blk[:m]
        if ne:
            Yx[k] = blk[m:]
    if ne:
        return z, mu, Y, str(sol.status), Yx
    return z, mu, Y, str(sol.status)


# ---------------------------------------------------------------------------
# Exact arithmetic
# ---------------------------------------------------------------------------

def bareiss_pd(N):
    """N symmetric list-of-lists of ints; True iff positive definite (all
    leading principal minors > 0; Bareiss pivots are exactly those minors)."""
    N = [row[:] for row in N]
    d, prev = len(N), 1
    for k in range(d):
        p = N[k][k]
        if p <= 0:
            return False
        Nk = N[k]
        for i in range(k + 1, d):
            Ni = N[i]; nik = Ni[k]
            for j in range(k + 1, d):
                Ni[j] = (Ni[j] * p - nik * Nk[j]) // prev
        prev = p
    return True


def int_matrix_2M(V, v, parent, one_idx, mu, y):
    """2*M_T as an integer matrix for integer mu (len V) and y (len V, y[v]=0);
    one_idx = the 4 vertices of T."""
    N = [[0] * V for _ in range(V)]
    for i in range(V):
        if i == v or y[i] == 0:
            continue
        j = parent[i]
        N[i][i] += 2 * y[i]
        N[i][j] -= y[i]
        N[j][i] -= y[i]
    for a in range(V):
        for b in one_idx:
            N[a][b] += mu[a]
            N[b][a] += mu[a]
    return N


def round_and_check(V, v, parent, tets, mu, Y, z, start_bits=None, max_bits=44):
    """Try integer scalings D = 2^k; returns (mu_int, Y_int, k) or None."""
    if start_bits is None:
        start_bits = max(4, int(math.ceil(math.log2(8.0 / max(z, 1e-12)))))
    for bits in range(start_bits, max_bits + 1, 2):
        D = float(2 ** bits)
        mu_i = [int(t) for t in np.rint(mu * D)]
        Y_i = np.maximum(np.rint(Y * D), 0.0)
        Y_i[:, v] = 0
        Y_l = [[int(t) for t in row] for row in Y_i]
        ok = True
        for f, y in zip(tets, Y_l):
            if not bareiss_pd(int_matrix_2M(V, v, parent, (v,) + tuple(f), mu_i, y)):
                ok = False
                break
        if ok:
            return mu_i, Y_l, bits
    return None


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def solve_tree(V, v, parent, tets):
    z, mu, Y, st = solve_joint(V, v, parent, tets)
    res = None
    if z > 1e-9:
        res = round_and_check(V, v, parent, tets, mu, Y, z)
        if res is None:                      # refine: tighter tolerance, retry
            z2, mu2, Y2, st2 = solve_joint(V, v, parent, tets, tol=1e-11, max_iter=600)
            if z2 > 1e-9:
                res = round_and_check(V, v, parent, tets, mu2, Y2, z2)
                z, st = z2, st2
    return z, st, res


def _worker(task):
    V, gid, v, trees, tets = task
    lines, unc, zs = [], [], []
    for parent in trees:
        z, st, res = solve_tree(V, v, parent, tets)
        zs.append(z)
        if res is None:
            unc.append(dict(graph=gid, root=v, parent=list(parent), z=z, status=st))
            continue
        mu_i, Y_l, bits = res
        lines.append(json.dumps(dict(kind="cert", gid=gid, root=v, parent=list(parent),
                                     z=float(f"{z:.6g}"), mu=mu_i,
                                     Y=[[list(f), y] for f, y in zip(tets, Y_l)]),
                                separators=(",", ":")))
    return gid, v, lines, unc, zs


def automorphism_to_root(gd):
    """For every non-root vertex u an automorphism (as a list perm, perm[x] =
    image of x) with perm[u] = representative root of u's orbit."""
    from networkx.algorithms.isomorphism import GraphMatcher
    G = gd.graph
    rep = {u: o[0] for o in gd.orbits for u in o}
    need = {u for u in range(gd.V) if rep[u] != u}
    out = {}
    for phi in GraphMatcher(G, G).isomorphisms_iter():
        for u in list(need):
            if phi[u] == rep[u]:
                out[u] = [phi[x] for x in range(gd.V)]
                need.discard(u)
        if not need:
            break
    assert not need
    return out


def run(V, workers=4, chunk=100, out=None, only=None):
    t0 = time.time()
    gds = ec.graph_data(V)
    out = out or os.path.join(HERE, f"certs_V{V}.jsonl.gz")
    tasks = []
    for gd in gds:
        for v in gd.roots:
            if only is not None and (gd.gid, v) not in only:
                continue
            T = gd.trees[v]
            for s in range(0, len(T), chunk):
                tasks.append((V, gd.gid, v, T[s:s + chunk], gd.decomposition[v]))
    ntrees = sum(len(t[3]) for t in tasks)
    print(f"V={V}: {len(gds)} graphs, {sum(len(g.roots) for g in gds)} roots, "
          f"{ntrees} trees, {len(tasks)} tasks", flush=True)
    stats = {}
    all_unc = []
    zmin_cert = math.inf
    n_cert = 0
    tmp = out + ".part"
    with gzip.open(tmp, "wt", compresslevel=6) as fo:
        fo.write(json.dumps(dict(kind="header", V=V, n_graphs=len(gds),
                                 description="joint-mu certificates; see certify_joint.py")) + "\n")
        for gd in gds:
            fo.write(json.dumps(dict(kind="graph", gid=gd.gid, edges=[list(e) for e in gd.edges],
                                     faces=[list(f) for f in gd.faces], roots=gd.roots,
                                     aut={str(u): p for u, p in automorphism_to_root(gd).items()}),
                                separators=(",", ":")) + "\n")
        done = 0
        if workers > 1:
            import multiprocessing as mp
            pool = mp.get_context("spawn").Pool(workers)
            it = pool.imap_unordered(_worker, tasks, chunksize=1)
        else:
            pool = None
            it = map(_worker, tasks)
        for k, (gid, v, lines, unc, zs) in enumerate(it):
            for ln in lines:
                fo.write(ln + "\n")
            all_unc.extend(unc)
            s = stats.setdefault((gid, v), dict(trees=0, cert=0, zmin=math.inf))
            s["trees"] += len(zs)
            s["cert"] += len(lines)
            n_cert += len(lines)
            cz = [json.loads(l)["z"] for l in lines]
            if cz:
                s["zmin"] = min(s["zmin"], min(cz))
                zmin_cert = min(zmin_cert, min(cz))
            done += len(zs)
            if (k + 1) % 50 == 0:
                print(f"  {done}/{ntrees} trees, {n_cert} certified, "
                      f"{len(all_unc)} uncertified, {time.time() - t0:.0f}s", flush=True)
        if pool is not None:
            pool.close(); pool.join()
    os.replace(tmp, out)
    print(f"V={V}: {done} trees solved, {n_cert} certified exactly, {len(all_unc)} uncertified; "
          f"min z (certified) = {zmin_cert:.6g}; time {time.time() - t0:.1f}s")
    for (gid, v), s in sorted(stats.items()):
        print(f"   graph {gid} root {v}: trees {s['trees']}, certified {s['cert']}, "
              f"min z {s['zmin']:.4g}")
    with open(os.path.join(HERE, f"uncertified_V{V}.json"), "w") as f:
        json.dump(all_unc, f)
    print("wrote", out)
    return all_unc


def crosscheck(n=40, seed=0):
    """Compare z from the direct Clarabel formulation with the cvxpy prototype
    (explore/joint_mu_test.py) on random V=6,7 trees."""
    sys.path.insert(0, os.path.join(ROOT, "explore"))
    cwd = os.getcwd(); os.chdir(ROOT)
    from joint_mu_test import joint_sdp
    os.chdir(cwd)
    rng = np.random.default_rng(seed)
    worst = 0.0
    for V in (6, 7):
        gds = ec.graph_data(V)
        for _ in range(n // 2):
            gd = gds[rng.integers(len(gds))]
            v = gd.roots[rng.integers(len(gd.roots))]
            p = gd.trees[v][rng.integers(len(gd.trees[v]))]
            t = time.time(); z1, *_ = solve_joint(V, v, p, gd.decomposition[v]); t1 = time.time() - t
            t = time.time(); z2, *_ = joint_sdp(V, v, p, gd.decomposition[v]); t2 = time.time() - t
            worst = max(worst, abs(z1 - z2))
            print(f"V={V} g{gd.gid} r{v}: z direct {z1:.8f} ({t1*1e3:.0f} ms)  cvxpy {z2:.8f} ({t2*1e3:.0f} ms)")
    print("max |dz| =", worst)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "crosscheck"])
    ap.add_argument("V", nargs="?", type=int, default=7)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--chunk", type=int, default=100)
    a = ap.parse_args()
    if a.cmd == "run":
        run(a.V, workers=a.workers, chunk=a.chunk)
    else:
        crosscheck()
