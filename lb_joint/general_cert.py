"""
General certificate ("certG") combining two SOUND strengthenings of the joint-mu
certificate, used on the trees the plain joint SDP does not certify.

(S1) Steepest-neighbour tree.  For a configuration we may CHOOSE the parent of
     each i != v as p(i) in argmax_{k ~ i} r_k . r_i (k over G-neighbours).
     Some neighbour weakly shadows i (r_j . r_i >= |r_i|^2), so the argmax does
     too; shadowing strictly increases |r| (|r_i|^2 <= r_i.r_p <= |r_i||r_p|,
     equality would force r_p = r_i), so p is acyclic and every non-root vertex
     has a parent: p is a spanning tree of G rooted at v.  Every configuration
     therefore has a steepest-neighbour tree among the enumerated trees, and for
     THAT tree the extra inequalities
         r_i . r_k <= r_i . r_p(i)   (i != v, k ~ i, k != p(i)),
     i.e. <A_{ik}, R> <= 0 with A_{ik} = sym(E_ik) - sym(E_{i p(i)}), are valid.
     (For non-simplicial P the argmax over the triangulated graph G only enlarges
     the candidate set; still valid.)  Multipliers lam >= 0 may be added in every
     block.

(S2) All cones.  For every vertex q, the cone decomposition from q (T = conv(q,F),
     F a face of G not containing q) gives, for any direction mu^(q), some T_q of
     positive volume with <sym(mu^(q) 1_{T_q}^T), R> <= 0.  Certificate for a set
     Q of apexes: mu^(q), symmetric Z_q with sum_{q in Q} Z_q = 0, and for every
     q in Q and T in cone(q) multipliers y >= 0 (and lam >= 0 if S1) with
         N_q(T) + Z_q > 0,   N_q(T) = sum y_i S_i + sum lam_ik A_ik + sym(mu^(q) 1_T^T).
     Soundness: sum_q N_q(T_q) = sum_q (N_q(T_q) + Z_q) is positive definite but
     <sum_q N_q(T_q), R> <= 0, contradiction with R >= 0, R != 0.
     Q = {v} without S1 is the plain certificate.

Integer form stored (all entries integers; positive scaling is irrelevant):
     2 N_q(T) + W_q  with W_q = 2 Z_q integer symmetric, sum_q W_q = 0.
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


def neighbours(V, edges):
    nb = {u: set() for u in range(V)}
    for a, b in edges:
        nb[a].add(b); nb[b].add(a)
    return nb


def s1_list(V, v, parent, nb):
    return [(i, k) for i in range(V) if i != v for k in sorted(nb[i]) if k != parent[i]]


def solve_general(V, v, parent, nb, faces, Q, use_s1, tol=1e-8, max_iter=400, zbox=10.0,
                  omegas=None):
    import clarabel as cl
    pos, nt = cj.svec_pos(V)
    r2 = math.sqrt(2.0)
    others = [i for i in range(V) if i != v]
    m = len(others)
    s1 = s1_list(V, v, parent, nb) if use_s1 else []
    ns = len(s1)
    nloc = m + ns
    blocks = []                                  # (group index, apex q, face)
    if omegas is not None:                       # S3: arbitrary omega vectors, one cone
        assert len(Q) == 1
        blocks = [(0, Q[0], ("w", j)) for j in range(len(omegas))]
    else:
        for g, q in enumerate(Q):
            for f in faces:
                if q not in f:
                    blocks.append((g, q, tuple(f)))
    nB = len(blocks)
    nQ = len(Q)
    nmu = nQ * V
    nzeta = (nQ - 1) * nt
    off_y = 1 + nmu + nzeta
    nx = off_y + nB * nloc
    # svec columns for S_i and A_ik
    Ssv = np.zeros((nt, nloc))
    for c, i in enumerate(others):
        j = parent[i]
        Ssv[pos[(i, i)], c] += 1.0
        a, b = min(i, j), max(i, j)
        Ssv[pos[(a, b)], c] -= 0.5 * r2
    for c, (i, k) in enumerate(s1):
        j = parent[i]
        a, b = min(i, k), max(i, k); Ssv[pos[(a, b)], m + c] += 0.5 * r2
        a, b = min(i, j), max(i, j); Ssv[pos[(a, b)], m + c] -= 0.5 * r2
    Isv = np.zeros(nt)
    for i in range(V):
        Isv[pos[(i, i)]] = 1.0
    rows, bs = [], []
    nm = nB * nloc
    rows.append(sp.coo_matrix((-np.ones(nm), (np.arange(nm), off_y + np.arange(nm))), shape=(nm, nx)))
    bs.append(np.zeros(nm))
    rows.append(sp.coo_matrix((np.ones(nm), (np.repeat(np.arange(nB), nloc), off_y + np.arange(nm))),
                              shape=(nB, nx)))
    bs.append(np.ones(nB))
    nbox = nmu + nzeta
    idx = 1 + np.arange(nbox)
    rows.append(sp.coo_matrix((np.concatenate([np.ones(nbox), -np.ones(nbox)]),
                               (np.arange(2 * nbox), np.concatenate([idx, idx]))), shape=(2 * nbox, nx)))
    bs.append(np.concatenate([np.ones(nmu), zbox * np.ones(nzeta), np.ones(nmu), zbox * np.ones(nzeta)]))
    cones = [cl.NonnegativeConeT(nm + nB + 2 * nbox)]
    mu_cache = {}

    def mu_part(one):
        key = tuple(one)
        if key not in mu_cache:
            Mk = np.zeros((nt, V))
            for (a, b), qq in pos.items():
                if a == b:
                    Mk[qq, a] -= one[a]
                else:
                    Mk[qq, a] -= r2 * one[b] / 2
                    Mk[qq, b] -= r2 * one[a] / 2
            mu_cache[key] = Mk
        return mu_cache[key]
    for k, (g, q, f) in enumerate(blocks):
        Ak = np.zeros((nt, nx))
        Ak[:, 0] = Isv
        vec = omegas[f[1]] if omegas is not None else cj.tet_indicator(V, q, f)
        Ak[:, 1 + g * V: 1 + (g + 1) * V] = mu_part(vec)
        if nQ > 1:          # + Z_q ; Z_last = -sum of the others (zeta = svec(Z))
            if g < nQ - 1:
                Ak[:, 1 + nmu + g * nt: 1 + nmu + (g + 1) * nt] -= np.eye(nt)
            else:
                for h in range(nQ - 1):
                    Ak[:, 1 + nmu + h * nt: 1 + nmu + (h + 1) * nt] += np.eye(nt)
        Ak[:, off_y + k * nloc: off_y + (k + 1) * nloc] = -Ssv
        rows.append(sp.coo_matrix(Ak)); bs.append(np.zeros(nt))
        cones.append(cl.PSDTriangleConeT(V))
    A = sp.vstack(rows).tocsc()
    b = np.concatenate(bs)
    P = sp.csc_matrix((nx, nx))
    qv = np.zeros(nx); qv[0] = -1.0
    st = cl.DefaultSettings()
    st.verbose = False
    st.tol_gap_abs = tol; st.tol_gap_rel = tol; st.tol_feas = tol
    st.max_iter = max_iter
    sol = cl.DefaultSolver(P, qv, A, b, cones, st).solve()
    x = np.array(sol.x)

    def unsvec(s):
        Z = np.zeros((V, V))
        for (a, bb), qq in pos.items():
            if a == bb:
                Z[a, a] = s[qq]
            else:
                Z[a, bb] = Z[bb, a] = s[qq] / r2
        return Z
    Zs = [unsvec(x[1 + nmu + h * nt: 1 + nmu + (h + 1) * nt]) for h in range(nQ - 1)]
    if nQ > 1:
        Zs.append(-sum(Zs))
    else:
        Zs = [np.zeros((V, V))]
    Ys, Ls = [], []
    for k in range(nB):
        blk = x[off_y + k * nloc: off_y + (k + 1) * nloc]
        y = np.zeros(V); y[others] = blk[:m]
        Ys.append(y); Ls.append(blk[m:])
    return dict(z=float(x[0]), mus=x[1:1 + nmu].reshape(nQ, V), Zs=Zs, Ys=Ys, Ls=Ls,
                blocks=blocks, s1=s1, status=str(sol.status))


def int_block_matrix(V, v, parent, one_idx, mu, W, y, s1, lam):
    N = [[W[a][b] if W is not None else 0 for b in range(V)] for a in range(V)]
    for i in range(V):
        if i != v and y[i]:
            j = parent[i]
            N[i][i] += 2 * y[i]; N[i][j] -= y[i]; N[j][i] -= y[i]
    for (i, k), c in zip(s1, lam):
        if c:
            j = parent[i]
            N[i][k] += c; N[k][i] += c; N[i][j] -= c; N[j][i] -= c
    for a in range(V):
        for b in one_idx:
            N[a][b] += mu[a]; N[b][a] += mu[a]
    return N


def round_general(V, v, parent, Q, res):
    z = res["z"]
    if z <= 1e-9:
        return None
    nQ = len(Q)
    start = max(4, int(math.ceil(math.log2(8.0 * nQ / z))))
    for bits in range(start, 48, 2):
        D = float(2 ** bits)
        mus = [[int(t) for t in np.rint(mu * D)] for mu in res["mus"]]
        Ws = [np.rint(2 * D * Z).astype(np.int64) for Z in res["Zs"][:-1]]
        Ws = [[[int(t) for t in row] for row in W] for W in Ws]
        if nQ > 1:
            last = [[-sum(W[a][b] for W in Ws) for b in range(V)] for a in range(V)]
            Ws.append(last)
        else:
            Ws = [None]
        ok = True
        Yi, Li = [], []
        for (g, q, f), y, lam in zip(res["blocks"], res["Ys"], res["Ls"]):
            yi = [int(t) for t in np.maximum(np.rint(y * D), 0)]; yi[v] = 0
            li = [int(t) for t in np.maximum(np.rint(lam * D), 0)]
            Yi.append(yi); Li.append(li)
            if not cj.bareiss_pd(int_block_matrix(V, v, parent, (q,) + f, mus[g], Ws[g], yi,
                                                  res["s1"], li)):
                ok = False
                break
        if ok:
            groups = []
            for g, q in enumerate(Q):
                bl = []
                for (gg, qq, f), yi, li in zip(res["blocks"], Yi, Li):
                    if gg != g:
                        continue
                    a = [[i, k, c] for (i, k), c in zip(res["s1"], li) if c]
                    bl.append([list(f), yi, a])
                groups.append(dict(q=q, mu=mus[g], W=Ws[g], blocks=bl))
            return groups, bits
    return None


# ---------------------------------------------------------------------------
# (S3) omega-subdivision (branch and bound over the simplex of weights lambda)
# ---------------------------------------------------------------------------
# Soundness: R omega(lam) = 0 with omega(lam) = sum_T lam_T 1_T, lam in the simplex
# Lambda.  If lam lies in a sub-simplex with vertices lam^1..lam^m, then
# omega(lam) = sum_j theta_j omega_j (theta >= 0), so 0 = mu^T R omega(lam) gives
# some j with <sym(mu omega_j^T), R> <= 0; if M_j = sum y S + sym(mu omega_j^T) > 0
# for all vertices j of the sub-simplex, contradiction.  The leaves of a
# bisection tree rooted at Lambda cover Lambda.

def exact_int_check_omega(V, v, parent, s1, mu_i, y_i, lam_i, w_int):
    N = int_block_matrix(V, v, parent, (), mu_i, None, y_i, s1, lam_i)
    for a in range(V):
        for b in range(V):
            N[a][b] += mu_i[a] * w_int[b] + w_int[a] * mu_i[b]
    return cj.bareiss_pd(N)


def s3_bnb(V, v, parent, nb, faces, use_s1, max_nodes=64, tol_z=1e-7):
    """Returns (certified, n_nodes, n_leaves, max_depth)."""
    from fractions import Fraction
    tf = [f for f in faces if v not in f]
    nT = len(tf)
    Ones = np.array([cj.tet_indicator(V, v, f) for f in tf]).T
    root = [tuple(Fraction(int(i == j)) for i in range(nT)) for j in range(nT)]
    stack = [(root, 0)]
    nodes = leaves = maxd = 0
    while stack:
        simp, depth = stack.pop()
        nodes += 1
        maxd = max(maxd, depth)
        if nodes > max_nodes:
            return False, nodes, leaves, maxd
        lam_f = [np.array([float(t) for t in lam]) for lam in simp]
        oms = [Ones @ l for l in lam_f]
        res = solve_general(V, v, parent, nb, faces, [v], use_s1, omegas=oms)
        ok = False
        if res["z"] > tol_z:
            den = [max(t.denominator for t in lam) for lam in simp]
            w_int = [[int(sum(lam[T] * den_j for T in range(nT) if tf[T] and Ones[k, T]))
                      for k in range(V)] for lam, den_j in zip(simp, den)]
            for bits in range(max(4, int(math.ceil(math.log2(8.0 / res["z"])))), 46, 2):
                D = float(2 ** bits)
                mu_i = [int(t) for t in np.rint(res["mus"][0] * D)]
                good = True
                for j in range(len(simp)):
                    y = res["Ys"][j]; lm = res["Ls"][j]
                    yi = [int(t) for t in np.maximum(np.rint(y * D), 0)]; yi[v] = 0
                    li = [int(t) for t in np.maximum(np.rint(lm * D), 0)]
                    # M_j uses omega_j; w_int[j] = den_j * omega_j (positive multiple):
                    # rescale the multipliers by den_j as well to keep the same matrix
                    yi = [t * den[j] for t in yi]; li = [t * den[j] for t in li]
                    mj = [t * 1 for t in mu_i]
                    if not exact_int_check_omega(V, v, parent, res["s1"], mj, yi, li, w_int[j]):
                        good = False
                        break
                if good:
                    ok = True
                    break
        if ok:
            leaves += 1
            continue
        # bisect the longest edge (in lambda coordinates)
        best, pair = -1.0, None
        for a in range(len(simp)):
            for b in range(a + 1, len(simp)):
                d = float(sum((x - y) ** 2 for x, y in zip(simp[a], simp[b])))
                if d > best:
                    best, pair = d, (a, b)
        a, b = pair
        mid = tuple((x + y) / 2 for x, y in zip(simp[a], simp[b]))
        c1 = list(simp); c1[a] = mid
        c2 = list(simp); c2[b] = mid
        stack.append((c1, depth + 1)); stack.append((c2, depth + 1))
    return True, nodes, leaves, maxd


# ---------------------------------------------------------------------------
# experiments on the uncertified trees
# ---------------------------------------------------------------------------

def depths(parent, v):
    d = []
    for i in range(len(parent)):
        x, s = i, 0
        while x != v:
            x = parent[x]; s += 1
        d.append(s)
    return d


def try_tree(task):
    """Returns dict of z-values for the experiments, and an exact certG line
    for the first configuration that certifies (order: S1, S2 pairs, S1+S2 all)."""
    V, gid, v, parent, edges, faces, modes = task
    nb = neighbours(V, edges)
    out = dict(graph=gid, root=v, parent=list(parent))
    line = None

    def attempt(Q, s1, key):
        nonlocal line
        res = solve_general(V, v, parent, nb, faces, Q, s1)
        out[key] = res["z"]
        if res["z"] > 1e-7:
            rc = round_general(V, v, parent, Q, res)
            if rc is not None:
                out[key + "_exact"] = True
                if line is None:
                    groups, bits = rc
                    line = json.dumps(dict(kind="certG", gid=gid, root=v, parent=list(parent),
                                           s1=bool(s1), z=float(f"{res['z']:.6g}"), groups=groups),
                                      separators=(",", ":"))
                return True
            out[key + "_exact"] = False
        return False
    if "S1" in modes:
        attempt([v], True, "S1")
    if "S3" in modes or "S3S1" in modes:
        for key, s1 in (("S3", False), ("S3S1", True)):
            if key in modes:
                ok, nn, nl, md = s3_bnb(V, v, parent, nb, faces, s1)
                out[key + "_exact"] = ok
                out[key + "_nodes"] = nn
    if "S2pair" in modes:
        out["S2pair"] = None
        d = depths(parent, v)
        for q in sorted((u for u in range(V) if u != v), key=lambda u: (d[u], u)):
            if attempt([v, q], False, f"S2_{q}"):
                out["S2pair"] = q
                break
    if "S2all" in modes:
        attempt(list(range(V)), False, "S2all")
    if "S1S2all" in modes:
        attempt(list(range(V)), True, "S1S2all")
    return out, line


def run_experiments(V, unc_path, modes, workers=4, out_prefix=None, limit=None):
    import gzip
    import multiprocessing as mp
    t0 = time.time()
    unc = json.load(open(unc_path))
    if limit:
        unc = unc[:limit]
    gds = ec.graph_data(V, with_trees=False)
    tasks = [(V, c["graph"], c["root"], tuple(c["parent"]), gds[c["graph"]].edges,
              gds[c["graph"]].faces, modes) for c in unc]
    out_prefix = out_prefix or os.path.join(HERE, f"exp_V{V}_{'_'.join(modes)}")
    results = []
    with gzip.open(out_prefix + "_certG.jsonl.gz", "wt") as fo, \
            mp.get_context("spawn").Pool(workers) as pool:
        fo.write(json.dumps(dict(kind="header_certG", V=V, modes=modes)) + "\n")
        for k, (r, line) in enumerate(pool.imap_unordered(try_tree, tasks, chunksize=2)):
            results.append(r)
            if line:
                fo.write(line + "\n")
            if (k + 1) % 250 == 0:
                print(f"  {k + 1}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    json.dump(results, open(out_prefix + "_results.json", "w"))
    n = len(results)
    print(f"V={V}: {n} previously uncertified trees; modes {modes}; {time.time() - t0:.0f}s")
    if "S1" in modes:
        print(f"  S1 alone (Q={{v}}):            {sum(bool(r.get('S1_exact')) for r in results)} certified exactly")
    for key in ("S3", "S3S1"):
        if key in modes:
            print(f"  {key} (omega-subdivision, <=64 nodes): {sum(bool(r.get(key + '_exact')) for r in results)} certified exactly")
    if "S2pair" in modes:
        print(f"  S2 with Q={{v,q}} (some q):    {sum(r['S2pair'] is not None for r in results)} certified exactly")
    if "S2all" in modes:
        print(f"  S2 with Q=all vertices:        {sum(bool(r.get('S2all_exact')) for r in results)} certified exactly")
    if "S1S2all" in modes:
        print(f"  S1+S2 with Q=all vertices:     {sum(bool(r.get('S1S2all_exact')) for r in results)} certified exactly")
    anyc = [r for r in results if any(r.get(k) for k in r if k.endswith("_exact"))]
    print(f"  certified by at least one: {len(anyc)}; remaining: {n - len(anyc)}")
    return results


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("V", type=int)
    ap.add_argument("unc")
    ap.add_argument("--modes", default="S1,S2pair,S1S2all")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run_experiments(a.V, a.unc, a.modes.split(","), a.workers, a.out, a.limit)


# ---------------------------------------------------------------------------
# relaxed Gram matrices for trees that survive every certificate
# ---------------------------------------------------------------------------

def gram_feasible(V, v, parent, nb, faces, Q, use_s1, iters=80, restarts=4, seed=0):
    """Alternating minimisation of sum_q |R omega_q| over R (PSD, tr R = 1,
    shadows, S1 if use_s1) and omega_q in the cone of the 1_T, T in cone(q)
    (lam_q in the simplex).  Returns the best (residual, R, {q: lam_q})."""
    import cvxpy as cp
    rng = np.random.default_rng(seed)
    Ones = {q: np.array([cj.tet_indicator(V, q, f) for f in faces if q not in f]).T for q in Q}
    mats = []
    for i in range(V):
        if i == v:
            continue
        M = np.zeros((V, V)); j = parent[i]
        M[i, i] = 1; M[i, j] -= 0.5; M[j, i] -= 0.5
        mats.append(M)
    if use_s1:
        for (i, k) in s1_list(V, v, parent, nb):
            M = np.zeros((V, V)); j = parent[i]
            M[i, k] += 0.5; M[k, i] += 0.5; M[i, j] -= 0.5; M[j, i] -= 0.5
            mats.append(M)
    R = cp.Variable((V, V), PSD=True)
    W = {q: cp.Parameter(V) for q in Q}
    probR = cp.Problem(cp.Minimize(sum(cp.norm(R @ W[q]) for q in Q)),
                       [cp.trace(R) == 1] + [cp.trace(M @ R) <= 0 for M in mats])
    Rp = cp.Parameter((V, V))
    L = {q: cp.Variable(Ones[q].shape[1], nonneg=True) for q in Q}
    probL = {q: cp.Problem(cp.Minimize(cp.norm(Rp @ Ones[q] @ L[q])), [cp.sum(L[q]) == 1]) for q in Q}
    best = (np.inf, None, None)
    for rs in range(restarts):
        lam = {q: (rng.dirichlet(np.ones(Ones[q].shape[1])) if rs else
                   np.ones(Ones[q].shape[1]) / Ones[q].shape[1]) for q in Q}
        for it in range(iters):
            for q in Q:
                W[q].value = Ones[q] @ lam[q]
            probR.solve(solver="CLARABEL")
            Rv = (R.value + R.value.T) / 2
            Rp.value = Rv
            for q in Q:
                probL[q].solve(solver="CLARABEL")
                lam[q] = np.maximum(L[q].value, 0); lam[q] /= lam[q].sum()
            res = sum(np.linalg.norm(Rv @ Ones[q] @ lam[q]) for q in Q)
            if res < best[0]:
                best = (res, Rv.copy(), {q: lam[q].copy() for q in Q})
            if res < 1e-9:
                break
        if best[0] < 1e-9:
            break
    return best


def hull_check(R, edges, tol=1e-6):
    """If R has numerical rank <= 3: factor R = X X^T and test the point set."""
    from scipy.spatial import ConvexHull
    w, U = np.linalg.eigh(R)
    rank = int((w > tol * w.max()).sum())
    out = dict(eigenvalues=[float(t) for t in w], rank=rank)
    if rank <= 3:
        X = U[:, -3:] * np.sqrt(np.maximum(w[-3:], 0))
        try:
            H = ConvexHull(X)
            hv = set(H.vertices.tolist())
            he = set()
            for s in H.simplices:
                for a in range(3):
                    for b in range(a + 1, 3):
                        he.add(tuple(sorted((int(s[a]), int(s[b])))))
            out["convex_position"] = len(hv) == len(X)
            out["hull_edges_equal_G"] = he == {tuple(sorted(e)) for e in edges}
            # solid centroid of the hull
            c0 = X[H.vertices].mean(axis=0)
            vol, cen = 0.0, np.zeros(3)
            for s in H.simplices:
                P = X[s]
                d = abs(np.linalg.det(np.vstack([P[0] - c0, P[1] - c0, P[2] - c0]))) / 6
                vol += d; cen += d * (c0 + P.sum(axis=0)) / 4
            out["centroid"] = (cen / vol).tolist()
            out["points"] = X.tolist()
        except Exception as e:      # degenerate point sets
            out["hull_error"] = str(e)
    return out
