"""Exploratory flat-limit frustum model. Outline K = tight double spiral (arms a, b; reference O = origin, T=(0,1)).
Candidate thickness for pushing the planar centroid in direction -u: frustum t = cone_E, E = K ∩ {u.p <= zeta}
(plateau 1 on E, minimal concave extension elsewhere).  Body = conv(K x {0}  ∪  E x {1}).
We compute how far 'below' O the planar centroid can be pushed:  val = max over (zeta) of  -(u . c')  with u=(0,1)
for symmetric outlines, and report the best found over the spiral angles. The numerical searches do not
certify global optimality, optimality over arbitrary thickness functions, or a vertex-count bound."""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.spatial import ConvexHull
from scipy.optimize import minimize, minimize_scalar
sys.path.insert(0, 'explore')
from flat_affine import build

def clip_below(P, zeta):
    """vertices of K ∩ {z <= zeta} (P ccw polygon)"""
    out = []
    n = len(P)
    for i in range(n):
        a, b = P[i], P[(i+1) % n]
        ia, ib = a[1] <= zeta, b[1] <= zeta
        if ia: out.append(a)
        if ia != ib:
            s = (zeta - a[1])/(b[1] - a[1]); out.append(a + s*(b - a))
    return np.array(out)

def centroid3(X):
    hull = ConvexHull(X)
    p0 = X[hull.vertices].mean(0); T = X[hull.simplices] - p0
    vol = np.abs(np.linalg.det(T)); return p0 + ((T.sum(1)/4)*vol[:, None]).sum(0)/vol.sum()

def zbar(P, zeta):
    E = clip_below(P, zeta)
    if len(E) < 3: return 1e9
    X = np.vstack([np.c_[P[:,0], np.zeros(len(P)), P[:,1]], np.c_[E[:,0], np.ones(len(E)), E[:,1]]])
    try:
        return centroid3(X)[2]
    except Exception:
        return 1e9

def best_zbar(dR, dL):
    P = build(dR, dL)
    zlo, zhi = P[:,1].min(), P[:,1].max()
    # coarse scan then refine
    zs = np.linspace(zlo + 0.02*(zhi-zlo), zhi - 0.02*(zhi-zlo), 25)
    vals = [zbar(P, z) for z in zs]
    k = int(np.argmin(vals))
    lo, hi = zs[max(k-1,0)], zs[min(k+1,len(zs)-1)]
    r = minimize_scalar(lambda z: zbar(P, z), bounds=(lo, hi), method="bounded", options={"xatol":1e-6})
    return min(r.fun, vals[k]), r.x

if __name__ == "__main__":
    rng = np.random.default_rng(0)
    for m in range(3, 8):
        best = (9, None, None)
        for trial in range(12):
            v = rng.uniform(0.7, 1.3, size=m); x0 = v*np.radians(rng.uniform(140, 175))/v.sum()
            def obj(x):
                x = np.abs(x)
                if x.sum() > np.pi or x.max() > 1.5: return 5.0
                return best_zbar(x, x)[0]
            res = minimize(obj, x0, method="Nelder-Mead", options={"xatol":1e-6, "fatol":1e-8, "maxfev":1500})
            if res.fun < best[0]: best = (res.fun, np.abs(res.x), best_zbar(np.abs(res.x), np.abs(res.x))[1])
        zb, ang, zeta = best
        print(f"symmetric m={m} (outline k={2*m+1}): best sampled frustum centroid z = {zb:+.5f} (need < 0); plateau level zeta={zeta:+.4f}; angles {np.degrees(ang).round(1)}", flush=True)
