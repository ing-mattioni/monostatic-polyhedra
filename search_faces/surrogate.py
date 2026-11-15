"""Planar long-spindle surrogate (limit tau -> 0) of the 3D stability problem.

Given planes in the p-parametrisation (P_i . x <= |P_i|^2, x = (X, y, z), long axis X):
  c   = solid centroid (exact for this polyhedron; invariant under X-rescaling of the tilts),
  B   = cross-section at X = c_X,  c' = (c_y, c_z).
Face i is stable in the limit iff edge i is present in B and the foot of c' on its line lies inside it.
G2 = -(second largest planar signed margin)/sqrt(area B).  G2 > 0 => for small enough tilt scale the
3D spindle is mono-stable (limit argument), so G2 is an elongation-free surrogate.
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from numba import njit
from core import G_planes, HDR


@njit(cache=True)
def planar_G(P):
    F = P.shape[0]
    b = np.empty(F)
    for i in range(F):
        b[i] = P[i, 0] ** 2 + P[i, 1] ** 2 + P[i, 2] ** 2
    out = G_planes(P, b)
    if out[0] != 0:
        return -10.0 - out[0], out[5], -1
    cX, cy, cz = out[2], out[3], out[4]
    nu = np.empty((F, 2))
    d = np.empty(F)
    for i in range(F):
        nn = np.sqrt(P[i, 1] ** 2 + P[i, 2] ** 2)
        if nn < 1e-12:
            return -9.0, out[5], -1
        nu[i, 0] = P[i, 1] / nn
        nu[i, 1] = P[i, 2] / nn
        d[i] = (b[i] - P[i, 0] * cX) / nn - (nu[i, 0] * cy + nu[i, 1] * cz)   # relative to c'
    # vertices by pairs
    L = 0.0
    for i in range(F):
        if abs(d[i]) > L:
            L = abs(d[i])
    tol = 1e-10 * L
    maxv = F * (F - 1) // 2
    V = np.empty((maxv, 2))
    nv = 0
    for i in range(F):
        for j in range(i + 1, F):
            det = nu[i, 0] * nu[j, 1] - nu[i, 1] * nu[j, 0]
            if abs(det) < 1e-13:
                continue
            x = (d[i] * nu[j, 1] - d[j] * nu[i, 1]) / det
            y = (nu[i, 0] * d[j] - nu[j, 0] * d[i]) / det
            ok = True
            for l in range(F):
                if nu[l, 0] * x + nu[l, 1] * y > d[l] + tol:
                    ok = False
                    break
            if not ok:
                continue
            dup = False
            for m in range(nv):
                if abs(V[m, 0] - x) + abs(V[m, 1] - y) < 10 * tol:
                    dup = True
                    break
            if not dup:
                V[nv, 0] = x
                V[nv, 1] = y
                nv += 1
    if nv < 3:
        return -8.0, out[5], -1
    # area (sort by angle around mean)
    mx = 0.0
    my = 0.0
    for m in range(nv):
        mx += V[m, 0]
        my += V[m, 1]
    mx /= nv
    my /= nv
    ang = np.empty(nv)
    for m in range(nv):
        ang[m] = np.arctan2(V[m, 1] - my, V[m, 0] - mx)
    o = np.argsort(ang)
    A = 0.0
    for t in range(nv):
        p = o[t]
        q = o[(t + 1) % nv]
        A += 0.5 * (V[p, 0] * V[q, 1] - V[q, 0] * V[p, 1])
    if A <= 0:
        return -7.0, out[5], -1
    sc = np.sqrt(A)
    s = np.full(F, -1e9)
    npres = 0
    for i in range(F):
        # endpoints on line i: params along tangent t = (-nu1, nu0)
        tmin = 1e300
        tmax = -1e300
        cnt = 0
        for m in range(nv):
            if abs(nu[i, 0] * V[m, 0] + nu[i, 1] * V[m, 1] - d[i]) <= 1e-8 * L:
                tt = -nu[i, 1] * V[m, 0] + nu[i, 0] * V[m, 1]
                cnt += 1
                if tt < tmin:
                    tmin = tt
                if tt > tmax:
                    tmax = tt
        if cnt < 2 or tmax - tmin < 1e-9 * L:
            continue
        npres += 1
        # foot of c' (origin) on line i has tangent param 0
        if tmin < 0 and tmax > 0:
            s[i] = min(-tmin, tmax) / sc
        else:
            s[i] = -min(abs(tmin), abs(tmax)) / sc
    i1 = -1
    for i in range(F):
        if s[i] > -1e8 and (i1 < 0 or s[i] > s[i1]):
            i1 = i
    i2 = -1
    for i in range(F):
        if s[i] > -1e8 and i != i1 and (i2 < 0 or s[i] > s[i2]):
            i2 = i
    if i2 < 0:
        return -6.0, out[5], npres
    return -s[i2], out[5], npres


@njit(cache=True)
def batch_planar(X, F):
    n = X.shape[0]
    J = np.empty(n)
    NF = np.empty(n)
    NP = np.empty(n)
    for k in range(n):
        g, nf, npres = planar_G(X[k].reshape(F, 3).copy())
        J[k] = g
        NF[k] = nf
        NP[k] = npres
    return J, NF, NP
