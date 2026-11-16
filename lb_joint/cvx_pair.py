"""
cvx_pair.py -- the two-cone certificate SDP (general_cert.solve_general with Q = [v, q], no S1 terms),
solved with CVXOPT instead of Clarabel.

Same problem:
    max z  s.t.  for every block k = (group g, apex Q[g], triangle F not containing Q[g]):
                 M_k = W_g + sum_i y^k_i S_{i->p(i)} + sym(mu_g 1_{T_k}^T) - z I  >= 0   (PSD),
                 y^k >= 0, sum_i y^k_i <= 1, |mu_g| <= 1 entrywise, |W_0| <= zbox entrywise, W_1 = -W_0.
Clarabel factorises the full KKT system of this problem, and the shared 9x9 matrix W_0 (which enters every
block) makes that factorisation fill in; CVXOPT's interior-point method instead forms the Schur complement
in the ~220 variables, which is cheap here.  The returned dictionary has the format of
general_cert.solve_general, so general_cert.round_general (exact rounding + exact check) is reused
unchanged.  Nothing here affects soundness: every certificate is re-checked in exact arithmetic, and again
by verify_joint.py.
"""
import numpy as np


def solve_pair_cvxopt(V, v, parent, faces, Q, zbox=10.0, tol=1e-8, maxiters=60, kktsolver='chol'):
    from cvxopt import matrix, spmatrix, solvers
    assert len(Q) == 2
    others = [i for i in range(V) if i != v]
    m = len(others)
    blocks = [(g, q, tuple(f)) for g, q in enumerate(Q) for f in faces if q not in f]
    nB = len(blocks)
    tri = [(a, b) for b in range(V) for a in range(b, V)]          # lower-triangular index pairs (a >= b)
    nt = len(tri)
    off_mu, off_W = 1, 1 + 2 * V
    off_y = off_W + nt
    nx = off_y + nB * m
    vidx = lambda a, b: a + b * V                                    # column-major position of entry (a, b)

    # ---- PSD blocks: Gs_k x + s_k = 0,  s_k = M_k(x)  =>  Gs_k column = -vec(coefficient matrix)
    Gs, hs = [], []
    for k, (g, q, f) in enumerate(blocks):
        one = np.zeros(V)
        one[q] = 1.0
        one[list(f)] = 1.0
        I_, J_, X_ = [], [], []

        def add(mat_entries, col, coef):
            for (a, b), val in mat_entries:
                # store both triangles (CVXOPT reads the lower triangle; full storage is harmless)
                I_.append(vidx(a, b)); J_.append(col); X_.append(-coef * val)
                if a != b:
                    I_.append(vidx(b, a)); J_.append(col); X_.append(-coef * val)
        add([((a, a), 1.0) for a in range(V)], 0, -1.0)              # -z I
        for a in range(V):                                           # sym(mu_g one^T)
            ent = {}
            for b in range(V):
                if one[b]:
                    key = (max(a, b), min(a, b))
                    ent[key] = ent.get(key, 0.0) + (1.0 if a == b else 0.5)
            add(list(ent.items()), off_mu + g * V + a, 1.0)
        sign = 1.0 if g == 0 else -1.0                               # W_0 = W, W_1 = -W
        for t, (a, b) in enumerate(tri):
            add([((a, b), 1.0)], off_W + t, sign)
        for c, i in enumerate(others):                               # y_i S_{i->p(i)}
            j = parent[i]
            add([((i, i), 1.0), ((max(i, j), min(i, j)), -0.5)], off_y + k * m + c, 1.0)
        Gs.append(spmatrix(X_, I_, J_, (V * V, nx)))
        hs.append(matrix(0.0, (V, V)))

    # ---- linear inequalities Gl x <= hl
    rows, cols, vals, h = [], [], [], []
    r = 0
    for k in range(nB):
        for c in range(m):                                           # -y <= 0
            rows.append(r); cols.append(off_y + k * m + c); vals.append(-1.0); h.append(0.0); r += 1
        for c in range(m):                                           # sum y <= 1
            rows.append(r); cols.append(off_y + k * m + c); vals.append(1.0)
        h.append(1.0); r += 1
    for t in range(2 * V):                                           # |mu| <= 1
        rows += [r, r + 1]; cols += [off_mu + t, off_mu + t]; vals += [1.0, -1.0]; h += [1.0, 1.0]; r += 2
    for t in range(nt):                                              # |W entries| <= zbox
        rows += [r, r + 1]; cols += [off_W + t, off_W + t]; vals += [1.0, -1.0]; h += [zbox, zbox]; r += 2
    Gl = spmatrix(vals, rows, cols, (r, nx))
    hl = matrix(h)
    c = matrix(0.0, (nx, 1))
    c[0] = -1.0
    solvers.options.update(dict(show_progress=False, abstol=tol, reltol=tol, feastol=tol, maxiters=maxiters))
    sol = solvers.sdp(c, Gl=Gl, hl=hl, Gs=Gs, hs=hs, kktsolver=kktsolver)
    x = np.array(sol["x"]).ravel() if sol["x"] is not None else np.zeros(nx)
    W = np.zeros((V, V))
    for t, (a, b) in enumerate(tri):
        W[a, b] = W[b, a] = x[off_W + t]
    Ys = []
    for k in range(nB):
        y = np.zeros(V)
        y[others] = x[off_y + k * m: off_y + (k + 1) * m]
        Ys.append(y)
    return dict(z=float(x[0]), mus=x[off_mu: off_mu + 2 * V].reshape(2, V), Zs=[W, -W], Ys=Ys,
                Ls=[np.zeros(0) for _ in range(nB)], blocks=blocks, s1=[], status=sol["status"])
