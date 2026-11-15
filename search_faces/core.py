"""Fast evaluator for mono-stability of a convex polyhedron given by planes a_i . x <= b_i.

Main entry points
-----------------
evaluate(A, b)            -> dict with vertices, faces, volume, centroid, s_i, G, margins ...  (python wrapper)
eval_core(A, b, out)      -> numba kernel, fills `out` (status, vol, c, nfaces, G, margins, s_i, area_i)
G_from_p(P)               -> numba: objective for the "Reshetov" parametrisation (plane i: p_i . x = |p_i|^2)
batch_G(X, F)             -> numba: objective for a population (rows = flattened P)

Definitions
-----------
c      : centroid of the solid.
q_i    : orthogonal projection of c on the plane of face i.
s_i    : signed in-plane distance of q_i to the boundary of face i:
           q_i inside  ->  +min over edges of distance to the edge LINES (= distance to the boundary)
           q_i outside ->  -Euclidean distance from q_i to the polygon (continuous, also under
                           combinatorial changes of the face polygon).
         (all divided by vol^(1/3) in the scaled outputs)
viol_i : Reshetov's measure for an outside point = max over edge lines of the violation
         (his Mathematica "Distance": min over edges of signed squared distance to edge lines).
G      : -(second largest s_i)/vol^(1/3)  (continuous everywhere).  G > 0 <=> exactly one face
         stable (the one with largest s_i, which is then the face closest to c).
M_res  : Reshetov-style safety margin = min( s_{i*}, min_{i != i*} viol_i ) / vol^(1/3)
M_euc  : min( s_{i*}, G*vol^(1/3) ) / vol^(1/3)
A face is "present" when its polygon has >= 3 distinct vertices and area > AREA_TOL * L^2
(L = max plane distance from the vertex centroid).  `nfaces` counts present faces.
status: 0 ok; 1 too few vertices; 2 unbounded / not closed; 3 bad volume; 4 bad input.
"""
import os
import numpy as np
from numba import njit

AREA_TOL = 1e-12

# output layout of eval_core: out[0]=status, out[1]=vol, out[2:5]=c, out[5]=nfaces, out[6]=G,
# out[7]=M_res, out[8]=M_euc, out[9]=istar, out[10]=nstable, out[11]=min rel. face area,
# out[12]=min rel. edge length, out[13]=nvert, out[14]=w=sqrt(vol/diam) (transverse scale),
# out[15]=Gw = G*vol^(1/3)/w (elongation-invariant objective) ; then s[F], viol[F], area[F], present[F]
HDR = 16


@njit(cache=True)
def _cross(a0, a1, a2, b0, b1, b2):
    return a1 * b2 - a2 * b1, a2 * b0 - a0 * b2, a0 * b1 - a1 * b0


@njit(cache=True)
def eval_core(A, b, out, verts_out, face_vidx, face_len):
    """A: (F,3), b: (F,). Fills out (len >= HDR+4F), verts_out (maxV,3), face_vidx (F,maxV), face_len (F,).
    Returns number of vertices."""
    F = A.shape[0]
    for k in range(out.shape[0]):
        out[k] = 0.0
    out[6] = -1e9
    out[7] = -1e9
    out[8] = -1e9
    # unit normals / distances
    U = np.empty((F, 3))
    h = np.empty(F)
    for i in range(F):
        nn = np.sqrt(A[i, 0] ** 2 + A[i, 1] ** 2 + A[i, 2] ** 2)
        if nn < 1e-300:
            out[0] = 4
            return 0
        U[i, 0] = A[i, 0] / nn
        U[i, 1] = A[i, 1] / nn
        U[i, 2] = A[i, 2] / nn
        h[i] = b[i] / nn
    L = 0.0
    for i in range(F):
        if abs(h[i]) > L:
            L = abs(h[i])
    if L == 0.0:
        L = 1.0
    tol = 1e-10 * L
    maxV = verts_out.shape[0]
    nv = 0
    for i in range(F):
        for j in range(i + 1, F):
            for k in range(j + 1, F):
                # det = u_i . (u_j x u_k)
                jk0, jk1, jk2 = _cross(U[j, 0], U[j, 1], U[j, 2], U[k, 0], U[k, 1], U[k, 2])
                det = U[i, 0] * jk0 + U[i, 1] * jk1 + U[i, 2] * jk2
                if abs(det) < 1e-13:
                    continue
                ki0, ki1, ki2 = _cross(U[k, 0], U[k, 1], U[k, 2], U[i, 0], U[i, 1], U[i, 2])
                ij0, ij1, ij2 = _cross(U[i, 0], U[i, 1], U[i, 2], U[j, 0], U[j, 1], U[j, 2])
                x0 = (h[i] * jk0 + h[j] * ki0 + h[k] * ij0) / det
                x1 = (h[i] * jk1 + h[j] * ki1 + h[k] * ij1) / det
                x2 = (h[i] * jk2 + h[j] * ki2 + h[k] * ij2) / det
                ok = True
                for l in range(F):
                    if U[l, 0] * x0 + U[l, 1] * x1 + U[l, 2] * x2 > h[l] + tol:
                        ok = False
                        break
                if not ok:
                    continue
                dup = False
                for m in range(nv):
                    d = abs(verts_out[m, 0] - x0) + abs(verts_out[m, 1] - x1) + abs(verts_out[m, 2] - x2)
                    if d < 10 * tol:
                        dup = True
                        break
                if dup:
                    continue
                if nv >= maxV:
                    out[0] = 4
                    return nv
                verts_out[nv, 0] = x0
                verts_out[nv, 1] = x1
                verts_out[nv, 2] = x2
                nv += 1
    out[13] = nv
    if nv < 4:
        out[0] = 1
        return nv
    # interior point
    o0 = 0.0
    o1 = 0.0
    o2 = 0.0
    for m in range(nv):
        o0 += verts_out[m, 0]
        o1 += verts_out[m, 1]
        o2 += verts_out[m, 2]
    o0 /= nv
    o1 /= nv
    o2 /= nv
    # scale: max distance of planes from interior point
    L2 = 0.0
    for i in range(F):
        dd = h[i] - (U[i, 0] * o0 + U[i, 1] * o1 + U[i, 2] * o2)
        if dd > L2:
            L2 = dd
    for m in range(nv):
        dd = np.sqrt((verts_out[m, 0] - o0) ** 2 + (verts_out[m, 1] - o1) ** 2 + (verts_out[m, 2] - o2) ** 2)
        if dd > L2:
            L2 = dd
    itol = 1e-9 * L2
    area = np.zeros(F)
    present = np.zeros(F, dtype=np.bool_)
    ang = np.empty(nv)
    clos0 = 0.0
    clos1 = 0.0
    clos2 = 0.0
    sumA = 0.0
    nf = 0
    minedge = 1e300
    for l in range(F):
        cnt = 0
        for m in range(nv):
            dd = U[l, 0] * verts_out[m, 0] + U[l, 1] * verts_out[m, 1] + U[l, 2] * verts_out[m, 2] - h[l]
            if abs(dd) <= itol:
                face_vidx[l, cnt] = m
                cnt += 1
        face_len[l] = cnt
        if cnt < 3:
            face_len[l] = 0
            continue
        # order by angle around mean
        m0 = 0.0
        m1 = 0.0
        m2 = 0.0
        for t in range(cnt):
            m = face_vidx[l, t]
            m0 += verts_out[m, 0]
            m1 += verts_out[m, 1]
            m2 += verts_out[m, 2]
        m0 /= cnt
        m1 /= cnt
        m2 /= cnt
        # basis e1 perpendicular to U[l]
        ux, uy, uz = U[l, 0], U[l, 1], U[l, 2]
        if abs(ux) < 0.9:
            e10, e11, e12 = _cross(ux, uy, uz, 1.0, 0.0, 0.0)
        else:
            e10, e11, e12 = _cross(ux, uy, uz, 0.0, 1.0, 0.0)
        ne = np.sqrt(e10 * e10 + e11 * e11 + e12 * e12)
        e10 /= ne
        e11 /= ne
        e12 /= ne
        e20, e21, e22 = _cross(ux, uy, uz, e10, e11, e12)
        for t in range(cnt):
            m = face_vidx[l, t]
            dx = verts_out[m, 0] - m0
            dy = verts_out[m, 1] - m1
            dz = verts_out[m, 2] - m2
            ang[t] = np.arctan2(e20 * dx + e21 * dy + e22 * dz, e10 * dx + e11 * dy + e12 * dz)
        # insertion sort
        for t in range(1, cnt):
            ka = ang[t]
            kv = face_vidx[l, t]
            r = t - 1
            while r >= 0 and ang[r] > ka:
                ang[r + 1] = ang[r]
                face_vidx[l, r + 1] = face_vidx[l, r]
                r -= 1
            ang[r + 1] = ka
            face_vidx[l, r + 1] = kv
        # area (vector area along U)
        ar = 0.0
        for t in range(cnt):
            p = face_vidx[l, t]
            q = face_vidx[l, (t + 1) % cnt]
            c0, c1, c2 = _cross(verts_out[p, 0] - m0, verts_out[p, 1] - m1, verts_out[p, 2] - m2,
                                verts_out[q, 0] - m0, verts_out[q, 1] - m1, verts_out[q, 2] - m2)
            ar += 0.5 * (c0 * ux + c1 * uy + c2 * uz)
        if ar <= AREA_TOL * L2 * L2:
            face_len[l] = 0
            continue
        area[l] = ar
        present[l] = True
        nf += 1
        sumA += ar
        clos0 += ar * ux
        clos1 += ar * uy
        clos2 += ar * uz
        for t in range(cnt):
            p = face_vidx[l, t]
            q = face_vidx[l, (t + 1) % cnt]
            el = np.sqrt((verts_out[p, 0] - verts_out[q, 0]) ** 2 + (verts_out[p, 1] - verts_out[q, 1]) ** 2
                         + (verts_out[p, 2] - verts_out[q, 2]) ** 2)
            if el < minedge:
                minedge = el
    out[5] = nf
    if nf < 4:
        out[0] = 1
        return nv
    if np.sqrt(clos0 ** 2 + clos1 ** 2 + clos2 ** 2) > 1e-7 * sumA:
        out[0] = 2
        return nv
    # volume & centroid (fan from o)
    vol = 0.0
    g0 = 0.0
    g1 = 0.0
    g2 = 0.0
    for l in range(F):
        if not present[l]:
            continue
        cnt = face_len[l]
        w = face_vidx[l, 0]
        for t in range(1, cnt - 1):
            p = face_vidx[l, t]
            q = face_vidx[l, t + 1]
            a0 = verts_out[w, 0] - o0
            a1 = verts_out[w, 1] - o1
            a2 = verts_out[w, 2] - o2
            c0, c1, c2 = _cross(verts_out[p, 0] - o0, verts_out[p, 1] - o1, verts_out[p, 2] - o2,
                                verts_out[q, 0] - o0, verts_out[q, 1] - o1, verts_out[q, 2] - o2)
            v6 = a0 * c0 + a1 * c1 + a2 * c2
            vol += v6
            g0 += v6 * (verts_out[w, 0] + verts_out[p, 0] + verts_out[q, 0] + o0)
            g1 += v6 * (verts_out[w, 1] + verts_out[p, 1] + verts_out[q, 1] + o1)
            g2 += v6 * (verts_out[w, 2] + verts_out[p, 2] + verts_out[q, 2] + o2)
    if vol <= 1e-14 * L2 * L2 * L2:
        out[0] = 3
        return nv
    g0 /= 4 * vol
    g1 /= 4 * vol
    g2 /= 4 * vol
    vol /= 6.0
    out[1] = vol
    out[2] = g0
    out[3] = g1
    out[4] = g2
    sc = vol ** (1.0 / 3.0)
    s = out[HDR:HDR + F]
    viol = out[HDR + F:HDR + 2 * F]
    minarea = 1e300
    nst = 0
    for l in range(F):
        out[HDR + 2 * F + l] = area[l] / (sc * sc)
        out[HDR + 3 * F + l] = 1.0 if present[l] else 0.0
        if not present[l]:
            s[l] = -1e9
            viol[l] = 1e9
            continue
        if area[l] / (sc * sc) < minarea:
            minarea = area[l] / (sc * sc)
        ux, uy, uz = U[l, 0], U[l, 1], U[l, 2]
        dpl = ux * g0 + uy * g1 + uz * g2 - h[l]
        q0 = g0 - dpl * ux
        q1 = g1 - dpl * uy
        q2 = g2 - dpl * uz
        cnt = face_len[l]
        minline = 1e300
        minseg = 1e300
        inside = True
        for t in range(cnt):
            p = face_vidx[l, t]
            r = face_vidx[l, (t + 1) % cnt]
            t0 = verts_out[r, 0] - verts_out[p, 0]
            t1 = verts_out[r, 1] - verts_out[p, 1]
            t2 = verts_out[r, 2] - verts_out[p, 2]
            tl = np.sqrt(t0 * t0 + t1 * t1 + t2 * t2)
            if tl < 1e-300:
                continue
            n0, n1, n2 = _cross(ux, uy, uz, t0 / tl, t1 / tl, t2 / tl)
            w0 = q0 - verts_out[p, 0]
            w1 = q1 - verts_out[p, 1]
            w2 = q2 - verts_out[p, 2]
            dline = n0 * w0 + n1 * w1 + n2 * w2
            if dline < minline:
                minline = dline
            if dline < 0:
                inside = False
            # distance to segment
            tt = (w0 * t0 + w1 * t1 + w2 * t2) / (tl * tl)
            if tt < 0:
                tt = 0.0
            elif tt > 1:
                tt = 1.0
            e0 = w0 - tt * t0
            e1 = w1 - tt * t1
            e2 = w2 - tt * t2
            ds = np.sqrt(e0 * e0 + e1 * e1 + e2 * e2)
            if ds < minseg:
                minseg = ds
        if inside:
            s[l] = minline / sc
            viol[l] = -minline / sc
            nst += 1
        else:
            s[l] = -minseg / sc
            viol[l] = -minline / sc  # max violation (positive)
    # G = -(second largest s)
    i1 = -1
    for l in range(F):
        if present[l] and (i1 < 0 or s[l] > s[i1]):
            i1 = l
    i2 = -1
    for l in range(F):
        if present[l] and l != i1 and (i2 < 0 or s[l] > s[i2]):
            i2 = l
    G = -s[i2]
    out[6] = G
    diam2 = 0.0
    for m in range(nv):
        for m2 in range(m + 1, nv):
            dd = ((verts_out[m, 0] - verts_out[m2, 0]) ** 2 + (verts_out[m, 1] - verts_out[m2, 1]) ** 2
                  + (verts_out[m, 2] - verts_out[m2, 2]) ** 2)
            if dd > diam2:
                diam2 = dd
    wsc = np.sqrt(vol / np.sqrt(diam2))
    out[14] = wsc
    out[15] = G * sc / wsc
    out[9] = i1
    out[10] = nst
    out[11] = minarea
    out[12] = minedge / sc
    # Reshetov margin
    mv = 1e300
    for l in range(F):
        if present[l] and l != i1:
            if viol[l] < mv:
                mv = viol[l]
    if nst == 1:
        out[7] = min(s[i1], mv)
        out[8] = min(s[i1], G)
    else:
        out[7] = -mv if mv < 0 else -abs(G)
        out[8] = G
    return nv


@njit(cache=True)
def G_planes(A, b):
    F = A.shape[0]
    maxV = F * (F - 1) * (F - 2) // 6 + 1
    out = np.zeros(HDR + 4 * F)
    verts = np.zeros((maxV, 3))
    fv = np.zeros((F, maxV), dtype=np.int64)
    fl = np.zeros(F, dtype=np.int64)
    eval_core(A, b, out, verts, fv, fl)
    return out


@njit(cache=True)
def planes_from_p(P):
    F = P.shape[0]
    b = np.empty(F)
    for i in range(F):
        b[i] = P[i, 0] ** 2 + P[i, 1] ** 2 + P[i, 2] ** 2
    return P, b


@njit(cache=True)
def obj_from_x(x, F, lam, norm=0):
    """Objective (to MAXIMISE) for flattened p-parametrisation.
    J = G + lam * sum_{i != i*} min(0, -s_i)  ; penalised for degenerate bodies.
    norm=0: lengths / vol^(1/3);  norm=1: lengths / sqrt(vol/diam) (elongation invariant)."""
    P = x.reshape(F, 3).copy()
    A, b = planes_from_p(P)
    out = G_planes(A, b)
    st = out[0]
    if st != 0:
        return -10.0 - st, out
    G = out[6]
    i1 = int(out[9])
    pen = 0.0
    for l in range(F):
        if out[HDR + 3 * F + l] > 0 and l != i1:
            v = -out[HDR + l]
            if v < 0:
                pen += v
    if norm == 1:
        fac = (out[1] ** (1.0 / 3.0)) / out[14]
        return (G + lam * pen) * fac, out
    return G + lam * pen, out


@njit(cache=True)
def batch_obj(X, F, lam, minarea, wF, norm=0, amax=0.0):
    """Population evaluation. Returns (J, G, nfaces). Faces with relative area < minarea
    are penalised (J -= (minarea - area)*100) to keep them present (set minarea=0 to allow vanishing).
    wF: reward per missing face (0 = none)."""
    n = X.shape[0]
    J = np.empty(n)
    Gs = np.empty(n)
    NF = np.empty(n)
    for k in range(n):
        j, out = obj_from_x(X[k], F, lam, norm)
        Gs[k] = out[6] if out[0] == 0 else -10.0
        NF[k] = out[5]
        if out[0] == 0:
            if minarea > 0:
                for l in range(F):
                    a = out[HDR + 2 * F + l]
                    if a < minarea:
                        j -= (minarea - a) * 10.0
            j += wF * (F - out[5])
            if amax > 0:
                asp = out[1] ** (2.0 / 3.0) / (out[14] * out[14])   # diam / vol^(1/3)
                if asp > amax:
                    j -= 0.01 * (asp - amax)
        J[k] = j
        Gs[k] = Gs[k]
    return J, Gs, NF


def evaluate(A, b):
    """Python wrapper returning a dict (vertices, faces as ordered vertex index lists, etc.)."""
    A = np.ascontiguousarray(A, dtype=np.float64)
    b = np.ascontiguousarray(b, dtype=np.float64)
    F = A.shape[0]
    maxV = F * (F - 1) * (F - 2) // 6 + 1
    out = np.zeros(HDR + 4 * F)
    verts = np.zeros((maxV, 3))
    fv = np.zeros((F, maxV), dtype=np.int64)
    fl = np.zeros(F, dtype=np.int64)
    nv = eval_core(A, b, out, verts, fv, fl)
    res = dict(status=int(out[0]), vol=out[1], c=out[2:5].copy(), nfaces=int(out[5]), G=out[6],
               M_res=out[7], M_euc=out[8], istar=int(out[9]), nstable=int(out[10]),
               min_area=out[11], min_edge=out[12], nvert=int(nv), w=out[14], Gw=out[15],
               s=out[HDR:HDR + F].copy(), viol=out[HDR + F:HDR + 2 * F].copy(),
               area=out[HDR + 2 * F:HDR + 3 * F].copy(), present=out[HDR + 3 * F:HDR + 4 * F] > 0,
               vertices=verts[:nv].copy(),
               faces=[list(fv[l, :fl[l]]) if fl[l] > 0 else [] for l in range(F)])
    return res


def evaluate_p(P):
    P = np.asarray(P, dtype=np.float64).reshape(-1, 3)
    return evaluate(P, (P ** 2).sum(1))


# ----------------------------------------------------------------------------------------------
# Slow reference implementation (scipy)
# ----------------------------------------------------------------------------------------------
def reference(A, b):
    from scipy.spatial import HalfspaceIntersection, ConvexHull
    from scipy.optimize import linprog
    A = np.asarray(A, float)
    b = np.asarray(b, float)
    F = len(b)
    # interior point via Chebyshev center
    nrm = np.linalg.norm(A, axis=1)
    res = linprog(np.r_[0, 0, 0, -1], A_ub=np.c_[A, nrm], b_ub=b, bounds=[(None, None)] * 3 + [(0, None)])
    if res.status != 0:
        return None
    ip = res.x[:3]
    hs = HalfspaceIntersection(np.c_[A, -b], ip)
    V = hs.intersections
    hull = ConvexHull(V)
    # volume / centroid by tetrahedra from ip
    vol = 0.0
    g = np.zeros(3)
    for simp in hull.simplices:
        a, bb, cc = V[simp]
        v6 = abs(np.dot(a - ip, np.cross(bb - ip, cc - ip)))
        vol += v6
        g += v6 * (a + bb + cc + ip)
    g /= 4 * vol
    vol /= 6
    sc = vol ** (1 / 3)
    s = np.full(F, -1e9)
    for i in range(F):
        u = A[i] / nrm[i]
        hi = b[i] / nrm[i]
        on = np.abs(V @ u - hi) < 1e-9 * np.abs(V).max()
        W = V[on]
        if len(W) < 3:
            continue
        e1 = np.cross(u, [1, 0, 0] if abs(u[0]) < 0.9 else [0, 1, 0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(u, e1)
        W2 = np.c_[W @ e1, W @ e2]
        try:
            h2 = ConvexHull(W2)
        except Exception:
            continue
        if h2.volume < 1e-12 * sc * sc:
            continue
        q = g - (u @ g - hi) * u
        q2 = np.array([q @ e1, q @ e2])
        # h2.equations: n.x + d <= 0 inside
        dl = -(h2.equations[:, :2] @ q2 + h2.equations[:, 2])
        if dl.min() >= 0:
            s[i] = dl.min() / sc
        else:
            P2 = W2[h2.vertices]
            dmin = 1e300
            for k in range(len(P2)):
                p0, p1 = P2[k], P2[(k + 1) % len(P2)]
                t = np.clip(np.dot(q2 - p0, p1 - p0) / np.dot(p1 - p0, p1 - p0), 0, 1)
                dmin = min(dmin, np.linalg.norm(q2 - p0 - t * (p1 - p0)))
            s[i] = -dmin / sc
    return dict(vol=vol, c=g, s=s)


# ----------------------------------------------------------------------------------------------
# Reshetov's examples (plane i: n_i . x = |n_i|^2)
# ----------------------------------------------------------------------------------------------
RESHETOV = {
    14: [[19, 70, -410], [0, 289, -186], [-26, -389, -2], [-6, 207, -331], [13, -329, 101], [1, -167, 233],
         [-5, -263, 181], [0, -73, 253], [-27, -85, -443], [0, 107, 201], [0, 181, 103], [0, 287, -27],
         [32, -394, -134], [0, 18, 244]],
    15: [[-64, 0, -37], [80, 0, 0], [-46, 0, -52], [29, 0, -47], [-74, -2, -18], [-84, 0, 6], [-90, 7, 40],
         [-23, 0, -61], [2, 0, -61], [63, 0, -30], [80, -1, 35], [-75, -7, 50], [40, -6, 87], [66, -2, 67],
         [24, 8, 107]],
    16: [[6, -7, -60], [8, -29, -46], [8, -42, -20], [-23, 34, 88], [-5, 3, 80], [-3, 79, 0], [-18, 82, -25],
         [-1, -23, 72], [6, -56, 29], [8, -57, 5], [3, 14, -63], [1, -44, 57], [0, 47, -50], [-3, 33, 74],
         [0, 32, -60], [-1, 60, -37]],
    17: [[0, -11, 92], [49, 0, -39], [-79, 4, -41], [-89, 1, -21], [70, 0, 3], [-52, 2, -63], [69, 0, 20],
         [63, -1, 41], [64, 0, -18], [20, 0, -52], [-95, 1, 4], [-91, -11, 36], [29, -8, 79], [28, 8, 85],
         [50, -3, 62], [-94, 12, 35], [-21, 1, -71]],
}


def reshetov_p(F):
    return np.array(RESHETOV[F], dtype=np.float64)
