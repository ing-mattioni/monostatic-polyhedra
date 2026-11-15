"""Flat-pyramid family: symmetric double Conway spiral polygon B with reference point at origin.
Top vertex (0,1). Spiral vertices at angles phi_1<...<phi_m (from the top), r_i = r_{i-1} cos(dphi_i).
Polygon vertices: top, P_{+-1..m}  (k = 2m+1), bottom edge horizontal iff ... (general).
Need: -3 g in B  (then apex A' = -3g gives pyramid centroid (3g+A')/4 = 0).
We maximize slack s = d_bot - 3 z_g  (by symmetry g on the axis)."""
import numpy as np
from scipy.optimize import minimize

def polygon(dphi):
    phis = np.cumsum(dphi)
    r = np.cumprod(np.cos(dphi))
    xs = r*np.sin(phis); zs = r*np.cos(phis)
    # right side from top going down, then left side going up
    P = [(0.0, 1.0)] + list(zip(xs, zs)) + list(zip(-xs[::-1], zs[::-1]))
    return np.array(P)

def area_centroid(P):
    x, z = P[:,0], P[:,1]
    x1, z1 = np.roll(x,-1), np.roll(z,-1)
    cr = x*z1 - x1*z
    A = cr.sum()/2
    cx = ((x+x1)*cr).sum()/(6*A); cz = ((z+z1)*cr).sum()/(6*A)
    return A, cx, cz

def slack(dphi):
    P = polygon(dphi)
    A, cx, cz = area_centroid(P)
    # bottom depth along axis: the polygon's lowest point on the axis x=0
    # bottom edge connects P_m and P_{-m}: at height z_m (horizontal) -> depth = -z_m
    zs = P[:,1]
    m = len(dphi)
    d_bot = -P[m,1]
    return d_bot - 3*cz, P, cz, d_bot

def main():
    for m in range(2, 9):
        best = None
        rng = np.random.default_rng(0)
        for trial in range(200):
            x0 = rng.uniform(0.1, 1.0, size=m); x0 *= rng.uniform(1.5, 3.0)/x0.sum()
            obj = lambda d: -slack(np.abs(d))[0] + 1e3*max(0, np.abs(d).sum()-np.pi)**2 + 1e3*sum(max(0, a-1.55)**2 for a in np.abs(d))
            res = minimize(obj, x0, method="Nelder-Mead", options={"xatol":1e-10,"fatol":1e-12,"maxiter":20000})
            d = np.abs(res.x)
            if d.sum() <= np.pi and (best is None or -res.fun > best[0]):
                best = (-res.fun, d)
        s, P, cz, d_bot = slack(best[1])
        print(f"m={m}: k={2*m+1} polygon vertices, V={2*m+2}: best slack={s:+.5f}  z_g={cz:.5f} d_bot={d_bot:.5f}  angles(deg)={np.degrees(best[1]).round(2)}")

if __name__ == '__main__':
    main()
