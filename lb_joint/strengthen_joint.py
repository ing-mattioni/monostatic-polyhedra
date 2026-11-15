"""
Strengthenings of the joint-mu certificate, tested on the trees that the plain
joint SDP (certify_joint.py) does not certify, plus analysis tools.

(a) PAIR certificate with a second cone decomposition (from a vertex q != v).
    Validity: the centroid identity holds for EVERY cone decomposition of P
    (triangulate the boundary, cone from any vertex q):
        sum_{T'} vol(T') * sum_{k in T'} r_k = 0,   T' = conv(q, F'), q not in F'.
    Hence, for any mu, mu' in R^V, there are T (cone from v) and T' (cone from q),
    both of positive volume, with
        <sym(mu 1_T^T), R> <= 0   and   <sym(mu' 1_T'^T), R> <= 0
    SIMULTANEOUSLY (the two choices are independent).  Therefore, if for EVERY
    pair (T, T') there are y^{T,T'} >= 0 with
        M_{T,T'} = sum_i y_i S_{i->p(i)} + sym(mu 1_T^T) + sym(mu' 1_T'^T)  > 0,
    then for the actual pair (T*, T'*) we get <M, R> <= 0, contradicting
    M > 0, R >= 0, R != 0.  One convex SDP in (mu, mu', all y), with nT*nT'
    PSD blocks.  (Taking q = v recovers the plain certificate.)
    General form used below: a list of directions d = 1..D with free vectors
    mu_d, and a list of blocks, each block giving one indicator vector per
    direction; M_block = sum y S + sum_d sym(mu_d one_{block,d}^T).

(a') ALL-q certificate restricted to SAME-face pairs is NOT valid (the T and
    T' are independent), so only full products are used.

(c) Feasible Gram matrix for a tree that stays uncertified: R >= 0, tr R = 1,
    <S_i, R> <= 0, R omega = 0 with omega = sum_T lam_T 1_T, lam >= 0,
    sum lam = 1 (T from the cone decomposition from v).  Bilinear; solved by
    alternating convex minimisation of |R omega| (R-step: SDP; lam-step: QP),
    then reported with its eigenvalues (numerical rank).
"""

from __future__ import annotations

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
sys.path.insert(0, HERE)
import certify_joint as cj  # noqa: E402
import enumerate_cases as ec  # noqa: E402


def solve_multi(V, v, parent, blocks, D, tol=1e-8, max_iter=400, want_dual=False):
    """max z s.t. for every block b: sum_i y^b_i S_i + sum_d sym(mu_d one_{b,d}^T) >= zI,
    y^b >= 0, sum y^b <= 1, |mu_d|_inf <= 1.   blocks: list of lists of D
    indicator vectors (numpy arrays of length V, or None)."""
    import clarabel as cl
    pos, nt = cj.svec_pos(V)
    r2 = math.sqrt(2.0)
    others = [i for i in range(V) if i != v]
    m = len(others)
    nB = len(blocks)
    nmu = D * V
    nx = 1 + nmu + nB * m
    Ssv = np.zeros((nt, m))
    for c, i in enumerate(others):
        j = parent[i]
        Ssv[pos[(i, i)], c] += 1.0
        a, b = min(i, j), max(i, j)
        Ssv[pos[(a, b)], c] += -0.5 * r2
    Isv = np.zeros(nt)
    for i in range(V):
        Isv[pos[(i, i)]] = 1.0
    A1 = sp.coo_matrix((-np.ones(nB * m), (np.arange(nB * m), 1 + nmu + np.arange(nB * m))),
                       shape=(nB * m, nx))
    r = np.repeat(np.arange(nB), m)
    A2 = sp.coo_matrix((np.ones(nB * m), (r, 1 + nmu + np.arange(nB * m))), shape=(nB, nx))
    A3 = sp.coo_matrix((np.concatenate([np.ones(nmu), -np.ones(nmu)]),
                        (np.arange(2 * nmu), np.concatenate([1 + np.arange(nmu)] * 2))),
                       shape=(2 * nmu, nx))
    rows = [A1, A2, A3]
    bs = [np.zeros(nB * m), np.ones(nB), np.ones(2 * nmu)]
    cones = [cl.NonnegativeConeT(nB * m + nB + 2 * nmu)]
    # mu-part template for an indicator vector `one` and direction d
    def mu_part(one, d):
        Mk = np.zeros((nt, V))
        for (a, b), q in pos.items():
            if a == b:
                Mk[q, a] -= one[a]
            else:
                Mk[q, a] -= r2 * one[b] / 2
                Mk[q, b] -= r2 * one[a] / 2
        return Mk
    cache = {}
    for k, blk in enumerate(blocks):
        Ak = np.zeros((nt, nx))
        Ak[:, 0] = Isv
        for d, one in enumerate(blk):
            if one is None:
                continue
            key = (d, tuple(one))
            if key not in cache:
                cache[key] = mu_part(one, d)
            Ak[:, 1 + d * V: 1 + (d + 1) * V] += cache[key]
        Ak[:, 1 + nmu + k * m: 1 + nmu + (k + 1) * m] = -Ssv
        rows.append(sp.coo_matrix(Ak))
        bs.append(np.zeros(nt))
        cones.append(cl.PSDTriangleConeT(V))
    A = sp.vstack(rows).tocsc()
    b = np.concatenate(bs)
    P = sp.csc_matrix((nx, nx))
    q = np.zeros(nx); q[0] = -1.0
    st = cl.DefaultSettings()
    st.verbose = False
    st.tol_gap_abs = tol; st.tol_gap_rel = tol; st.tol_feas = tol
    st.max_iter = max_iter
    sol = cl.DefaultSolver(P, q, A, b, cones, st).solve()
    x = np.array(sol.x)
    z = float(x[0])
    mus = x[1:1 + nmu].reshape(D, V)
    Y = np.zeros((nB, V))
    for k in range(nB):
        Y[k, others] = x[1 + nmu + k * m: 1 + nmu + (k + 1) * m]
    out = dict(z=z, mus=mus, Y=Y, status=str(sol.status))
    if want_dual:
        # dual PSD matrices X_b (unscale svec)
        zz = np.array(sol.z)
        off = nB * m + nB + 2 * nmu
        Xs = []
        for k in range(nB):
            s = zz[off + k * nt: off + (k + 1) * nt]
            X = np.zeros((V, V))
            for (a, bb), qq in pos.items():
                if a == bb:
                    X[a, a] = s[qq]
                else:
                    X[a, bb] = X[bb, a] = s[qq] / r2
            Xs.append(X)
        out["X"] = Xs
    return out


def indicators(V, apex, faces):
    return [cj.tet_indicator(V, apex, f) for f in faces]


def plain_blocks(V, v, faces):
    return [[o] for o in indicators(V, v, [f for f in faces if v not in f])]


def pair_blocks(V, v, q, faces):
    A = indicators(V, v, [f for f in faces if v not in f])
    B = indicators(V, q, [f for f in faces if q not in f])
    return [[a, b] for a in A for b in B]


def round_check_multi(V, v, parent, blocks, D, res):
    """Exact rounding as in certify_joint (integers, Bareiss)."""
    z = res["z"]
    if z <= 1e-9:
        return None
    start = max(4, int(math.ceil(math.log2(8.0 * D / z))))
    for bits in range(start, 46, 2):
        Dd = float(2 ** bits)
        mus = [[int(t) for t in np.rint(mu * Dd)] for mu in res["mus"]]
        Yi = np.maximum(np.rint(res["Y"] * Dd), 0.0); Yi[:, v] = 0
        ok = True
        for blk, y in zip(blocks, Yi):
            N = [[0] * V for _ in range(V)]
            for i in range(V):
                if i == v or y[i] == 0:
                    continue
                j = parent[i]; yi = int(y[i])
                N[i][i] += 2 * yi; N[i][j] -= yi; N[j][i] -= yi
            for d, one in enumerate(blk):
                idx = [k for k in range(V) if one[k]]
                for a in range(V):
                    for bb in idx:
                        N[a][bb] += mus[d][a]; N[bb][a] += mus[d][a]
            if not cj.bareiss_pd(N):
                ok = False
                break
        if ok:
            return mus, [[int(t) for t in y] for y in Yi], bits
    return None


# ---------------------------------------------------------------------------
# structural statistics
# ---------------------------------------------------------------------------

def depths(parent, v):
    V = len(parent)
    d = [0] * V
    for i in range(V):
        x, s = i, 0
        while x != v:
            x = parent[x]; s += 1
        d[i] = s
    return d


def active_tets(V, v, parent, faces, thr=1e-3):
    """Faces F (cone from v) whose dual block X_T carries weight in the plain
    joint SDP (relative trace > thr)."""
    blocks = plain_blocks(V, v, faces)
    res = solve_multi(V, v, parent, blocks, 1, want_dual=True)
    tf = [f for f in faces if v not in f]
    tr = np.array([np.trace(X) for X in res["X"]])
    tr = tr / max(tr.sum(), 1e-300)
    return [tuple(f) for f, t in zip(tf, tr) if t > thr], tr, res


# ---------------------------------------------------------------------------
# feasible Gram matrix search (c)
# ---------------------------------------------------------------------------

def gram_search(V, v, parent, faces, lam0=None, iters=60, seed=0):
    import cvxpy as cp
    tf = [f for f in faces if v not in f]
    Ones = np.array(indicators(V, v, tf)).T          # V x nT
    nT = Ones.shape[1]
    rng = np.random.default_rng(seed)
    lam = np.ones(nT) / nT if lam0 is None else np.asarray(lam0, float)
    S = []
    for i in range(V):
        if i == v:
            continue
        M = np.zeros((V, V)); j = parent[i]
        M[i, i] = 1; M[i, j] -= 0.5; M[j, i] -= 0.5
        S.append(M)
    R = cp.Variable((V, V), PSD=True)
    w = cp.Parameter(V)
    consR = [cp.trace(R) == 1] + [cp.trace(Si @ R) <= 0 for Si in S]
    probR = cp.Problem(cp.Minimize(cp.norm(R @ w)), consR)
    L = cp.Variable(nT, nonneg=True)
    Rp = cp.Parameter((V, V))
    probL = cp.Problem(cp.Minimize(cp.norm(Rp @ Ones @ L)), [cp.sum(L) == 1])
    best = (np.inf, None, None)
    for it in range(iters):
        w.value = Ones @ lam
        probR.solve(solver="CLARABEL")
        Rv = R.value
        Rp.value = (Rv + Rv.T) / 2
        probL.solve(solver="CLARABEL")
        lam = np.maximum(L.value, 0); lam /= lam.sum()
        res = np.linalg.norm(Rv @ Ones @ lam)
        if res < best[0]:
            best = (res, Rv.copy(), lam.copy())
        if res < 1e-9:
            break
    return best


# ---------------------------------------------------------------------------
# driver: pair certificates for the uncertified trees
# ---------------------------------------------------------------------------

def q_order(V, v, parent):
    """Heuristic order of the second apex q: increasing tree depth (the deep
    vertices carry the bad tetrahedra, a shallow second apex works best)."""
    d = depths(parent, v)
    return sorted((u for u in range(V) if u != v), key=lambda u: (d[u], u))


def pair_worker(task):
    V, gid, v, parent, faces = task
    tried = []
    for q in q_order(V, v, parent):
        blocks = pair_blocks(V, v, q, faces)
        res = solve_multi(V, v, parent, blocks, 2)
        tried.append((q, res["z"]))
        if res["z"] > 1e-7:
            rc = round_check_multi(V, v, parent, blocks, 2, res)
            if rc is not None:
                mus, Yi, bits = rc
                A = [f for f in faces if v not in f]
                B = [f for f in faces if q not in f]
                prs = [(fa, fb) for fa in A for fb in B]
                line = json.dumps(dict(kind="cert2", gid=gid, root=v, q=q, parent=list(parent),
                                       z=float(f"{res['z']:.6g}"), mu=mus[0], mu2=mus[1],
                                       Y=[[list(fa), list(fb), y] for (fa, fb), y in zip(prs, Yi)]),
                                  separators=(",", ":"))
                return gid, v, parent, line, tried
    return gid, v, parent, None, tried


def run_pairs(V, workers=4, unc_path=None, out=None):
    import gzip
    t0 = time.time()
    unc_path = unc_path or os.path.join(HERE, f"uncertified_V{V}.json")
    unc = json.load(open(unc_path))
    gds = ec.graph_data(V, with_trees=False)
    tasks = [(V, c["graph"], c["root"], tuple(c["parent"]), gds[c["graph"]].faces) for c in unc]
    print(f"V={V}: {len(tasks)} uncertified trees -> pair certificates", flush=True)
    out = out or os.path.join(HERE, f"certs_V{V}_pairs.jsonl.gz")
    remaining = []
    qstat = {}
    n_ok = 0
    import multiprocessing as mp
    with gzip.open(out + ".part", "wt") as fo:
        fo.write(json.dumps(dict(kind="header2", V=V, description="pair certificates")) + "\n")
        with mp.get_context("spawn").Pool(min(workers, 4)) as pool:
            for k, (gid, v, parent, line, tried) in enumerate(
                    pool.imap_unordered(pair_worker, tasks, chunksize=4)):
                if line is None:
                    remaining.append(dict(graph=gid, root=v, parent=list(parent), tried=tried))
                else:
                    fo.write(line + "\n"); n_ok += 1
                    qstat[len(tried)] = qstat.get(len(tried), 0) + 1
                if (k + 1) % 200 == 0:
                    print(f"  {k + 1}/{len(tasks)}: {n_ok} pair-certified, {len(remaining)} not, "
                          f"{time.time() - t0:.0f}s", flush=True)
    os.replace(out + ".part", out)
    print(f"V={V}: {n_ok} of {len(tasks)} trees pair-certified exactly, {len(remaining)} remain; "
          f"#q tried until success: {dict(sorted(qstat.items()))}; {time.time() - t0:.0f}s")
    with open(os.path.join(HERE, f"uncertified_V{V}_after_pairs.json"), "w") as f:
        json.dump(remaining, f)
    return remaining


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["pairs"])
    ap.add_argument("V", type=int)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--unc", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run_pairs(a.V, a.workers, a.unc, a.out)
