"""Flat pyramid family with tight (Conway) spirals, asymmetric arms.
Reference point O = origin, top T = (0,1).  Right arm: polar angles (from top, clockwise) cumulative
sums of dR, radii r_i = r_{i-1} cos(dR_i).  Left arm similarly (counter-clockwise).
Pyramid apex A' in B gives body centroid (3g + A')/4; need it to be O  <=>  -3g in B.
slack = signed distance of -3g to the boundary of B (positive inside)."""
import numpy as np, sys, os
os.environ.setdefault("OMP_NUM_THREADS", "1")
from scipy.optimize import minimize

def build(dR, dL):
    def arm(d, sgn):
        ph = np.cumsum(d); r = np.cumprod(np.cos(d))
        return np.stack([sgn*r*np.sin(ph), r*np.cos(ph)], axis=1)
    R = arm(dR, 1.0); L = arm(dL, -1.0)
    # ccw order: T, left arm going down, right arm going up
    P = np.vstack([[0.0, 1.0], L, R[::-1]])
    return P

def slack(dR, dL):
    P = build(dR, dL)
    x, z = P[:,0], P[:,1]; x1, z1 = np.roll(x,-1), np.roll(z,-1)
    cr = x*z1 - x1*z
    A = cr.sum()/2
    if A <= 0: return -10.0, P, None
    g = np.array([((x+x1)*cr).sum(), ((z+z1)*cr).sum()])/(6*A)
    # convexity: all turns ccw
    ex, ez = x1-x, z1-z
    turn = ex*np.roll(ez,-1) - ez*np.roll(ex,-1)
    q = -3*g
    el = np.hypot(ex, ez)
    s = (ex*(q[1]-z) - ez*(q[0]-x))/el          # signed distance to each edge line (inside > 0)
    conv = min(0.0, turn.min())
    return min(s.min(), 0 if conv >= 0 else 1e9) + 50*conv, P, g

def opt(a, b, ntrial=300, seed=0):
    rng = np.random.default_rng(seed)
    best = (-1e9, None)
    for t in range(ntrial):
        def mk(n):
            v = rng.uniform(0.6, 1.4, size=n); return v*np.radians(rng.uniform(140, 178))/v.sum()
        x0 = np.concatenate([mk(a), mk(b)])
        def obj(x):
            x = np.abs(x)
            if x.max() > 1.5 or x[:a].sum() > np.pi or x[a:].sum() > np.pi: return 10.0
            return -slack(x[:a], x[a:])[0]
        res = minimize(obj, x0, method="Nelder-Mead", options={"xatol":1e-11,"fatol":1e-13,"maxiter":40000,"maxfev":40000})
        res = minimize(obj, res.x, method="Nelder-Mead", options={"xatol":1e-12,"fatol":1e-14,"maxiter":40000,"maxfev":40000})
        if -res.fun > best[0]: best = (-res.fun, np.abs(res.x))
    return best

def job(ab):
    a, b = ab
    s, x = opt(a, b)
    return a, b, s, x

if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor
    pairs = [(7,7),(7,6),(8,6),(6,6),(8,5),(9,5),(7,5),(9,4),(10,4),(8,7),(10,3),(11,3),(12,2)]
    with ProcessPoolExecutor(12) as ex:
        for a, b, s, x in ex.map(job, pairs):
            print(f"arms {a}+{b}: k={a+b+1}, V={a+b+2}: best slack = {s:+.6f}   R angles {np.degrees(x[:a]).round(2)}  L angles {np.degrees(x[a:]).round(2)}", flush=True)
            np.save(f"explore/flat_asym_{a}_{b}.npy", x)
