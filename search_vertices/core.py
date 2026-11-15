"""Fast evaluator for homogeneous mono-unstable polyhedra (numba).

X : (V,3) float64 points, P = conv(X), c = centroid of the SOLID P.
Angular margin of point i:  m_i = max_{j!=i} (x_j-x_i).(x_i-c) / (|x_j-x_i| |x_i-c|).
F(X) = min_{i != farthest} m_i     (F > 0  <=>  mono-unstable).
Points that are not hull vertices automatically have m_i > 0 (if inside) so F is
continuous and needs no vertex identification.

Hull: brute force over triples on a deterministically joggled copy of X (joggle
1e-10 * scale), so coplanar groups are triangulated consistently; volume/centroid
are then summed over the triangles using the ORIGINAL coordinates.  A closure
check (sum of area vectors == 0, Euler count) triggers a retry with a new joggle.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from numba import njit

_JOG = np.sin(np.arange(1, 3 * 64 * 8 + 1) * 12.9898) * 43758.5453
_JOG = (_JOG - np.floor(_JOG)) - 0.5           # deterministic pseudo-random in [-0.5,0.5)
_JOG = _JOG.reshape(8, 64, 3).copy()
R0_REL = 0.02     # parents closer than R0_REL*diam are penalised in the optimisation objective
KAPPA = 0.1
RPEN = 0.03       # global min pairwise distance (relative to diameter) below which a slope-2 penalty applies
R2MAX = 0.98      # second-largest radius / largest radius above which a slope-2 penalty applies (twin-T trap)


@njit(cache=True)
def _hull_tris(X, J, tris):
    """Brute-force facet triangles of conv(X) (X assumed in general position).
    Writes outward-oriented triangles into tris, returns count (or -1 if flat)."""
    V = X.shape[0]
    nt = 0
    for i in range(V):
        for j in range(i + 1, V):
            ax = X[j, 0] - X[i, 0]; ay = X[j, 1] - X[i, 1]; az = X[j, 2] - X[i, 2]
            for k in range(j + 1, V):
                bx = X[k, 0] - X[i, 0]; by = X[k, 1] - X[i, 1]; bz = X[k, 2] - X[i, 2]
                nx = ay * bz - az * by; ny = az * bx - ax * bz; nz = ax * by - ay * bx
                pos = False; neg = False
                for l in range(V):
                    if l == i or l == j or l == k:
                        continue
                    s = nx * (X[l, 0] - X[i, 0]) + ny * (X[l, 1] - X[i, 1]) + nz * (X[l, 2] - X[i, 2])
                    if s > 0.0:
                        pos = True
                        if neg:
                            break
                    elif s < 0.0:
                        neg = True
                        if pos:
                            break
                if pos and neg:
                    continue
                if not pos and not neg:
                    return -1
                if nt >= tris.shape[0]:
                    return -2
                if neg:      # all others on negative side: n outward, (i,j,k) ccw from outside
                    tris[nt, 0] = i; tris[nt, 1] = j; tris[nt, 2] = k
                else:
                    tris[nt, 0] = i; tris[nt, 1] = k; tris[nt, 2] = j
                nt += 1
    return nt


@njit(cache=True)
def hull_centroid(X, JOG):
    """Returns (volume, c[3], ntri, isvert[V], ok).  ok=False on failure (flat/degenerate)."""
    V = X.shape[0]
    c = np.zeros(3)
    isv = np.zeros(V, dtype=np.bool_)
    # scale
    mx = np.zeros(3)
    for d in range(3):
        s = 0.0
        for i in range(V):
            s += X[i, d]
        mx[d] = s / V
    scale = 0.0
    for i in range(V):
        r = 0.0
        for d in range(3):
            r += (X[i, d] - mx[d]) ** 2
        if r > scale:
            scale = r
    scale = np.sqrt(scale)
    if scale <= 0.0:
        return 0.0, c, 0, isv, False
    tris = np.empty((4 * V, 3), dtype=np.int64)
    Xj = np.empty_like(X)
    for attempt in range(JOG.shape[0]):
        for i in range(V):
            for d in range(3):
                Xj[i, d] = X[i, d] + 1e-10 * scale * JOG[attempt, i % JOG.shape[1], d]
        nt = _hull_tris(Xj, JOG, tris)
        if nt < 4:
            if nt == -1:
                return 0.0, c, 0, isv, False
            continue
        # closure + Euler check
        for i in range(V):
            isv[i] = False
        sx = 0.0; sy = 0.0; sz = 0.0; atot = 0.0
        for t in range(nt):
            a = tris[t, 0]; b = tris[t, 1]; e = tris[t, 2]
            isv[a] = True; isv[b] = True; isv[e] = True
            ux = Xj[b, 0] - Xj[a, 0]; uy = Xj[b, 1] - Xj[a, 1]; uz = Xj[b, 2] - Xj[a, 2]
            wx = Xj[e, 0] - Xj[a, 0]; wy = Xj[e, 1] - Xj[a, 1]; wz = Xj[e, 2] - Xj[a, 2]
            nx = uy * wz - uz * wy; ny = uz * wx - ux * wz; nz = ux * wy - uy * wx
            sx += nx; sy += ny; sz += nz
            atot += np.sqrt(nx * nx + ny * ny + nz * nz)
        nv = 0
        for i in range(V):
            if isv[i]:
                nv += 1
        if nt != 2 * nv - 4 or np.sqrt(sx * sx + sy * sy + sz * sz) > 1e-9 * atot:
            continue
        # volume & centroid with original coordinates
        ox = 0.0; oy = 0.0; oz = 0.0
        for i in range(V):
            ox += X[i, 0]; oy += X[i, 1]; oz += X[i, 2]
        ox /= V; oy /= V; oz /= V
        vol = 0.0; cx = 0.0; cy = 0.0; cz = 0.0
        for t in range(nt):
            a = tris[t, 0]; b = tris[t, 1]; e = tris[t, 2]
            ax_ = X[a, 0] - ox; ay_ = X[a, 1] - oy; az_ = X[a, 2] - oz
            bx_ = X[b, 0] - ox; by_ = X[b, 1] - oy; bz_ = X[b, 2] - oz
            ex_ = X[e, 0] - ox; ey_ = X[e, 1] - oy; ez_ = X[e, 2] - oz
            v6 = ax_ * (by_ * ez_ - bz_ * ey_) + ay_ * (bz_ * ex_ - bx_ * ez_) + az_ * (bx_ * ey_ - by_ * ex_)
            vol += v6
            cx += v6 * (ax_ + bx_ + ex_); cy += v6 * (ay_ + by_ + ey_); cz += v6 * (az_ + bz_ + ez_)
        if vol <= 0.0:
            return 0.0, c, nt, isv, False
        c[0] = ox + cx / (4.0 * vol); c[1] = oy + cy / (4.0 * vol); c[2] = oz + cz / (4.0 * vol)
        return vol / 6.0, c, nt, isv, True
    return 0.0, c, 0, isv, False


@njit(cache=True)
def margins_c(X, c):
    """Angular margins w.r.t. reference point c.  Returns (m[V], parent[V], far)."""
    V = X.shape[0]
    m = np.empty(V)
    par = np.empty(V, dtype=np.int64)
    d = np.empty(V)
    far = 0
    for i in range(V):
        d[i] = np.sqrt((X[i, 0] - c[0]) ** 2 + (X[i, 1] - c[1]) ** 2 + (X[i, 2] - c[2]) ** 2)
        if d[i] > d[far]:
            far = i
    for i in range(V):
        ux = X[i, 0] - c[0]; uy = X[i, 1] - c[1]; uz = X[i, 2] - c[2]
        best = -2.0; bj = -1
        di = d[i] if d[i] > 1e-300 else 1e-300
        for j in range(V):
            if j == i:
                continue
            wx = X[j, 0] - X[i, 0]; wy = X[j, 1] - X[i, 1]; wz = X[j, 2] - X[i, 2]
            nw = np.sqrt(wx * wx + wy * wy + wz * wz)
            if nw <= 1e-300:
                continue
            v = (wx * ux + wy * uy + wz * uz) / (nw * di)
            if v > best:
                best = v
                bj = j
        m[i] = best
        par[i] = bj
    return m, par, far


@njit(cache=True)
def margins_pen(X, c, r0, kappa):
    """Like margins_c but a parent j closer than r0 to i loses kappa*(1-r/r0) of its cosine
    (kills the coincident-pair attractor).  Returns (Fpen, far)."""
    V = X.shape[0]
    far = 0; dfar = -1.0
    for i in range(V):
        di = np.sqrt((X[i, 0] - c[0]) ** 2 + (X[i, 1] - c[1]) ** 2 + (X[i, 2] - c[2]) ** 2)
        if di > dfar:
            dfar = di; far = i
    F = 10.0
    for i in range(V):
        if i == far:
            continue
        ux = X[i, 0] - c[0]; uy = X[i, 1] - c[1]; uz = X[i, 2] - c[2]
        di = np.sqrt(ux * ux + uy * uy + uz * uz)
        if di < 1e-300:
            di = 1e-300
        best = -2.0
        for j in range(V):
            if j == i:
                continue
            wx = X[j, 0] - X[i, 0]; wy = X[j, 1] - X[i, 1]; wz = X[j, 2] - X[i, 2]
            nw = np.sqrt(wx * wx + wy * wy + wz * wz)
            if nw <= 1e-300:
                continue
            v = (wx * ux + wy * uy + wz * uz) / (nw * di)
            if nw < r0:
                v -= kappa * (1.0 - nw / r0)
            if v > best:
                best = v
        if best < F:
            F = best
    return F, far


@njit(cache=True)
def evalF(X, JOG):
    """Raw F (no penalties). Returns (F, vol, c, far, worst_index)."""
    vol, c, nt, isv, ok = hull_centroid(X, JOG)
    if not ok:
        return -10.0, 0.0, c, -1, -1
    m, par, far = margins_c(X, c)
    F = 10.0; wi = -1
    for i in range(X.shape[0]):
        if i != far and m[i] < F:
            F = m[i]; wi = i
    return F, vol, c, far, wi


@njit(cache=True)
def shape_stats(X):
    """(min pairwise distance / diameter, diameter)."""
    V = X.shape[0]
    dmin = 1e300; dmax = 0.0
    for i in range(V):
        for j in range(i + 1, V):
            r = np.sqrt((X[i, 0] - X[j, 0]) ** 2 + (X[i, 1] - X[j, 1]) ** 2 + (X[i, 2] - X[j, 2]) ** 2)
            if r < dmin:
                dmin = r
            if r > dmax:
                dmax = r
    return dmin / dmax, dmax


@njit(cache=True)
def trap_pen(Q, c, dim):
    """Penalty against degenerate traps: near-coincident points and twin farthest points.
    Q: (V,dim) points, c: reference (dim,)."""
    V = Q.shape[0]
    dmin = 1e300; dmax = 0.0
    for i in range(V):
        for j in range(i + 1, V):
            r = 0.0
            for k in range(dim):
                r += (Q[i, k] - Q[j, k]) ** 2
            r = np.sqrt(r)
            if r < dmin: dmin = r
            if r > dmax: dmax = r
    r1 = 0.0; r2 = 0.0
    for i in range(V):
        r = 0.0
        for k in range(dim):
            r += (Q[i, k] - c[k]) ** 2
        r = np.sqrt(r)
        if r > r1:
            r2 = r1; r1 = r
        elif r > r2:
            r2 = r
    pen = 0.0
    if dmin < RPEN * dmax:
        pen += 2.0 * (RPEN - dmin / dmax)
    if r2 > R2MAX * r1:
        pen += 2.0 * (r2 / r1 - R2MAX)
    return pen


@njit(cache=True)
def objective(X, JOG, dmin_rel, vol_rel):
    """Penalised objective (to MAXIMISE)."""
    vol, c, nt, isv, ok = hull_centroid(X, JOG)
    if not ok:
        return -10.0
    r, diam = shape_stats(X)
    F, far = margins_pen(X, c, R0_REL * diam, KAPPA)
    pen = trap_pen(X, c, 3)
    if r < dmin_rel:
        pen += (dmin_rel - r) / dmin_rel * 0.05
    vr = vol / diam ** 3
    if vr < vol_rel:
        pen += (vol_rel - vr) / vol_rel * 0.05
    return F - pen


@njit(cache=True)
def objective_batch(Z, V, JOG, dmin_rel, vol_rel):
    n = Z.shape[0]
    out = np.empty(n)
    for t in range(n):
        out[t] = objective(Z[t].reshape(V, 3), JOG, dmin_rel, vol_rel)
    return out


# ---------------------------------------------------------------- flat-limit surrogate
@njit(cache=True)
def evalF_flat(P, Y, JOG, r0rel=0.0, kappa=0.0):
    """Flat limit: points (p_i, eps*y_i), eps->0.  P: (V,2) planar (x,z), Y: (V,) heights.
    c' = xz-projection of the centroid of conv{(p_i, y_i)} (independent of eps).
    Returns (F_planar, cprime[2], far)."""
    V = P.shape[0]
    X = np.empty((V, 3))
    for i in range(V):
        X[i, 0] = P[i, 0]; X[i, 1] = Y[i]; X[i, 2] = P[i, 1]
    vol, c, nt, isv, ok = hull_centroid(X, JOG)
    cp = np.zeros(2)
    if not ok:
        return -10.0, cp, -1
    cp[0] = c[0]; cp[1] = c[2]
    d = np.empty(V)
    far = 0
    for i in range(V):
        d[i] = np.sqrt((P[i, 0] - cp[0]) ** 2 + (P[i, 1] - cp[1]) ** 2)
        if d[i] > d[far]:
            far = i
    F = 10.0
    for i in range(V):
        if i == far:
            continue
        ux = P[i, 0] - cp[0]; uz = P[i, 1] - cp[1]
        best = -2.0
        for j in range(V):
            if j == i:
                continue
            wx = P[j, 0] - P[i, 0]; wz = P[j, 1] - P[i, 1]
            nw = np.sqrt(wx * wx + wz * wz)
            if nw <= 1e-12 * d[far]:
                continue
            v = (wx * ux + wz * uz) / (nw * d[i])
            if nw < r0rel * d[far]:
                v -= kappa * (1.0 - nw / (r0rel * d[far]))
            if v > best:
                best = v
        if best < F:
            F = best
    return F, cp, far


@njit(cache=True)
def interior_depth(P, n_out):
    """Sum over i < n_out of the depth of P[i] inside the planar convex hull of all P (0 if on hull)."""
    V = P.shape[0]
    idx = np.argsort(P[:, 0] + 1e-12 * P[:, 1])
    H = np.empty(2 * V + 1, dtype=np.int64)
    k = 0
    for t in range(V):
        i = idx[t]
        while k >= 2 and ((P[H[k-1], 0] - P[H[k-2], 0]) * (P[i, 1] - P[H[k-2], 1]) - (P[H[k-1], 1] - P[H[k-2], 1]) * (P[i, 0] - P[H[k-2], 0])) <= 0:
            k -= 1
        H[k] = i; k += 1
    lo = k + 1
    for t in range(V - 2, -1, -1):
        i = idx[t]
        while k >= lo and ((P[H[k-1], 0] - P[H[k-2], 0]) * (P[i, 1] - P[H[k-2], 1]) - (P[H[k-1], 1] - P[H[k-2], 1]) * (P[i, 0] - P[H[k-2], 0])) <= 0:
            k -= 1
        H[k] = i; k += 1
    nh = k - 1          # H[0..nh-1] ccw hull, H[nh] == H[0]
    onh = np.zeros(V, dtype=np.bool_)
    for t in range(nh):
        onh[H[t]] = True
    tot = 0.0
    for i in range(n_out):
        if onh[i]:
            continue
        dm = 1e300
        for t in range(nh):
            a = H[t]; b = H[t + 1]
            ex = P[b, 0] - P[a, 0]; ez = P[b, 1] - P[a, 1]
            ne = np.sqrt(ex * ex + ez * ez) + 1e-300
            dd = (ex * (P[i, 1] - P[a, 1]) - ez * (P[i, 0] - P[a, 0])) / ne
            if dd < dm:
                dm = dd
        if dm > 0:
            tot += dm
    return tot


@njit(cache=True)
def objective_flat_batch(Z, V, JOG, dmin_rel, n_out=0):
    n = Z.shape[0]
    out = np.empty(n)
    for t in range(n):
        z = Z[t]
        P = z[:2 * V].reshape(V, 2)
        Y = z[2 * V:3 * V]
        F, cp, far = evalF_flat(P, Y, JOG, 2.0 * R0_REL, KAPPA)
        X = np.empty((V, 3))
        for i in range(V):
            X[i, 0] = P[i, 0]; X[i, 1] = Y[i]; X[i, 2] = P[i, 1]
        r, diam = shape_stats(X)
        pen = trap_pen(P, cp, 2)
        if n_out > 0:
            pen += 2.0 * interior_depth(P, n_out) / diam
        if r < dmin_rel:
            pen += (dmin_rel - r) / dmin_rel * 0.05
        out[t] = F - pen
    return out


# ---------------------------------------------------------------- python helpers
JOG = _JOG


def F(X):
    return evalF(np.ascontiguousarray(X, dtype=np.float64), JOG)[0]


def info(X):
    X = np.ascontiguousarray(X, dtype=np.float64)
    vol, c, nt, isv, ok = hull_centroid(X, JOG)
    m, par, far = margins_c(X, c)
    Fv = min(m[i] for i in range(len(X)) if i != far)
    r, diam = shape_stats(X)
    nvert = true_vertex_count(X)
    return dict(F=float(Fv), vol=float(vol), c=c.copy(), far=int(far), parent=par.copy(), m=m.copy(),
                dmin_rel=float(r), diam=float(diam), nvert=nvert, unstable=[int(i) for i in range(len(X)) if m[i] <= 0])


def true_vertex_count(X, rel=1e-9):
    """Count points whose removal shrinks the hull volume by more than rel*vol."""
    X = np.ascontiguousarray(X, dtype=np.float64)
    vol = hull_centroid(X, JOG)[0]
    cnt = 0
    for i in range(len(X)):
        Y = np.delete(X, i, axis=0)
        v2 = hull_centroid(Y, JOG)[0]
        if vol - v2 > rel * vol:
            cnt += 1
    return cnt


def normalize(X):
    """Translate centroid to 0 and scale so that max |x-c| = 1."""
    X = np.asarray(X, dtype=np.float64)
    vol, c, nt, isv, ok = hull_centroid(np.ascontiguousarray(X), JOG)
    Y = X - c
    return Y / np.linalg.norm(Y, axis=1).max()
