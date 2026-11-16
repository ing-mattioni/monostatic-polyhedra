"""
fast_pair.py -- faster solver for the two-cone certificate (S2 with Q = {v, q}) of general_cert.py.

Same optimisation problem as general_cert.solve_general(..., Q=[v, q], use_s1=False):
    max z  s.t.  for every block k = (group g, apex Q[g], triangle F):
                 W_g + sum_i y^k_i S_i + sym(mu_g 1_{T_k}^T)  >=  z I,
                 y^k >= 0, sum_i y^k_i <= 1,  |mu_g|_inf <= 1,  |W_0|_entries <= zbox,  W_1 = -W_0.
The only change is the FORMULATION handed to Clarabel: instead of one shared matrix variable W_0 that
appears in every block (45 dense columns, which make the sparse LDL factorisation of the KKT system
fill in completely), every block k gets its own copy W^(k), and consecutive copies are tied by equality
constraints (W^(k) = W^(k+1) inside a group, W^(last of group 0) = -W^(first of group 1)).  The KKT matrix
is then banded and the factorisation is cheap.  The returned dictionary has the format of
solve_general, so general_cert.round_general (exact rounding + exact positive-definiteness check) is
used unchanged; nothing here affects soundness (every certificate is re-checked exactly, and again by
verify_joint.py).
"""
import math

import numpy as np
import scipy.sparse as sp

import certify_joint as cj


def solve_pair_fast(V, v, parent, faces, Q, tol=1e-8, max_iter=400, zbox=10.0):
    import clarabel as cl
    assert len(Q) == 2
    pos, nt = cj.svec_pos(V)
    r2 = math.sqrt(2.0)
    others = [i for i in range(V) if i != v]
    m = len(others)
    blocks = [(g, q, tuple(f)) for g, q in enumerate(Q) for f in faces if q not in f]
    nB = len(blocks)
    # variable layout: [z | mu_0 (V) | mu_1 (V) | W^(0..nB-1) (nt each) | y^(0..nB-1) (m each)]
    off_mu = 1
    off_W = 1 + 2 * V
    off_y = off_W + nB * nt
    nx = off_y + nB * m
    Ssv = np.zeros((nt, m))
    for c, i in enumerate(others):
        j = parent[i]
        Ssv[pos[(i, i)], c] += 1.0
        a, b = min(i, j), max(i, j)
        Ssv[pos[(a, b)], c] -= 0.5 * r2
    Isv = np.zeros(nt)
    for i in range(V):
        Isv[pos[(i, i)]] = 1.0

    def mu_part(one):
        Mk = np.zeros((nt, V))
        for (a, b), qq in pos.items():
            if a == b:
                Mk[qq, a] -= one[a]
            else:
                Mk[qq, a] -= r2 * one[b] / 2
                Mk[qq, b] -= r2 * one[a] / 2
        return Mk

    rows, bs, cones = [], [], []
    # --- equality constraints tying the copies of W (ZeroCone: A x = b)
    eq = []
    I = sp.identity(nt, format="coo")
    for k in range(nB - 1):
        sign = -1.0 if blocks[k][0] == blocks[k + 1][0] else +1.0     # same group: W^k - W^(k+1) = 0
        A = sp.hstack([sp.coo_matrix((nt, off_W + k * nt)), I, sign * I,
                       sp.coo_matrix((nt, nx - off_W - (k + 2) * nt))])
        eq.append(A)
    if eq:
        rows.append(sp.vstack(eq)); bs.append(np.zeros(len(eq) * nt))
        cones.append(cl.ZeroConeT(len(eq) * nt))
    # --- nonnegative constraints
    nm = nB * m
    nn = []
    nn.append(sp.coo_matrix((-np.ones(nm), (np.arange(nm), off_y + np.arange(nm))), shape=(nm, nx)))   # y >= 0
    nb_rows = np.repeat(np.arange(nB), m)
    nn.append(sp.coo_matrix((np.ones(nm), (nb_rows, off_y + np.arange(nm))), shape=(nB, nx)))           # sum y <= 1
    idx = off_mu + np.arange(2 * V)
    nn.append(sp.coo_matrix((np.ones(2 * V), (np.arange(2 * V), idx)), shape=(2 * V, nx)))              # mu <= 1
    nn.append(sp.coo_matrix((-np.ones(2 * V), (np.arange(2 * V), idx)), shape=(2 * V, nx)))             # -mu <= 1
    idw = off_W + np.arange(nt)
    nn.append(sp.coo_matrix((np.ones(nt), (np.arange(nt), idw)), shape=(nt, nx)))                       # W^0 <= zbox
    nn.append(sp.coo_matrix((-np.ones(nt), (np.arange(nt), idw)), shape=(nt, nx)))                      # -W^0 <= zbox
    rows.append(sp.vstack(nn))
    bs.append(np.concatenate([np.zeros(nm), np.ones(nB), np.ones(4 * V), zbox * np.ones(2 * nt)]))
    cones.append(cl.NonnegativeConeT(nm + nB + 4 * V + 2 * nt))
    # --- PSD blocks:  s = b - A x = -z I + sym(mu 1_T) + W^k + sum y S  in PSD
    for k, (g, q, f) in enumerate(blocks):
        Ak = sp.lil_matrix((nt, nx))
        Ak[:, 0] = Isv.reshape(-1, 1)
        Ak[:, off_mu + g * V: off_mu + (g + 1) * V] = mu_part(cj.tet_indicator(V, q, f))
        Ak[:, off_W + k * nt: off_W + (k + 1) * nt] = -np.eye(nt)
        Ak[:, off_y + k * m: off_y + (k + 1) * m] = -Ssv
        rows.append(Ak.tocoo()); bs.append(np.zeros(nt))
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
    W0 = unsvec(x[off_W: off_W + nt])
    Ys = []
    for k in range(nB):
        y = np.zeros(V); y[others] = x[off_y + k * m: off_y + (k + 1) * m]
        Ys.append(y)
    return dict(z=float(x[0]), mus=x[off_mu: off_mu + 2 * V].reshape(2, V), Zs=[W0, -W0], Ys=Ys,
                Ls=[np.zeros(0) for _ in range(nB)], blocks=blocks, s1=[], status=str(sol.status))
