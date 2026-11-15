"""Core evaluator for homogeneous mono-unstable polyhedra.
X: (V,3) points. P = conv(X), c = centroid of the solid P.
Vertex i is NOT an unstable equilibrium iff exists j: (x_j - x_i).(x_i - c) > 0.
margin_i = max_j (u_j.u_i)/(u_i.u_i) - 1,  u = x - c.
f(X) = min over hull vertices i != farthest of margin_i   (mono-unstable iff f > 0)."""
import numpy as np
from scipy.spatial import ConvexHull, QhullError

def centroid(X):
    hull = ConvexHull(X)
    p0 = X[hull.vertices].mean(axis=0)
    T = X[hull.simplices] - p0            # (nf,3,3)
    vol = np.abs(np.linalg.det(T))        # 6*vol of each tetra (p0 interior => all same sign)
    cen = T.sum(axis=1) / 4.0
    V = vol.sum()
    c = p0 + (cen * vol[:, None]).sum(axis=0) / V
    return c, V / 6.0, hull

def margins(X):
    c, vol, hull = centroid(X)
    U = X - c
    G = U @ U.T
    d = np.diag(G).copy()
    M = G / d[:, None] - 1.0
    np.fill_diagonal(M, -np.inf)
    m = M.max(axis=1)
    return m, d, hull, c

def f(X):
    try:
        m, d, hull, c = margins(X)
    except QhullError:
        return -1e3
    hv = hull.vertices
    apex = hv[np.argmax(d[hv])]
    mm = [m[i] for i in hv if i != apex]
    return min(mm)

def info(X):
    m, d, hull, c = margins(X)
    hv = hull.vertices
    apex = hv[np.argmax(d[hv])]
    unstable = [int(i) for i in hv if m[i] <= 0]
    return dict(nvert=len(hv), apex=int(apex), unstable=unstable,
                fmin=min(m[i] for i in hv if i != apex), c=c)
