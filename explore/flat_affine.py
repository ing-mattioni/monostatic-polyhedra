"""Tight double-spiral polygon B (reference O = origin, top T=(0,1)).  Two weightings:
 cone   : body centroid (3g + A')/4 for a pyramid apex A' in B         -> feasible iff -3g in B
 affine : density l(p) = 1 + beta.(p-g) >= 0 on B (oblique prism)      -> centroid g + Sigma beta = O
          feasible iff l(v) = 1 - g^T Sigma^{-1} (v - g) >= 0 at all vertices v."""
import numpy as np, sys, os
os.environ.setdefault("OMP_NUM_THREADS", "1")
from scipy.optimize import minimize

def build(dR, dL):
    def arm(d, sgn):
        ph = np.cumsum(d); r = np.cumprod(np.cos(d))
        return np.stack([sgn*r*np.sin(ph), r*np.cos(ph)], axis=1)
    R = arm(dR, 1.0); L = arm(dL, -1.0)
    return np.vstack([[0.0, 1.0], L, R[::-1]])   # ccw

def moments(P):
    x, z = P[:,0], P[:,1]; x1, z1 = np.roll(x,-1), np.roll(z,-1)
    cr = x*z1 - x1*z
    A = cr.sum()/2
    gx = ((x+x1)*cr).sum()/(6*A); gz = ((z+z1)*cr).sum()/(6*A)
    Ixx = ((x*x + x*x1 + x1*x1)*cr).sum()/12
    Izz = ((z*z + z*z1 + z1*z1)*cr).sum()/12
    Ixz = ((x*z1 + 2*x*z + 2*x1*z1 + x1*z)*cr).sum()/24
    S = np.array([[Ixx/A - gx*gx, Ixz/A - gx*gz],[Ixz/A - gx*gz, Izz/A - gz*gz]])
    return A, np.array([gx, gz]), S

def convex_pen(P):
    x, z = P[:,0], P[:,1]; x1, z1 = np.roll(x,-1), np.roll(z,-1)
    ex, ez = x1-x, z1-z
    turn = ex*np.roll(ez,-1) - ez*np.roll(ex,-1)
    return min(0.0, turn.min())

def slack_affine(dR, dL):
    P = build(dR, dL)
    A, g, S = moments(P)
    if A <= 0: return -10.0
    pen = convex_pen(P)
    l = 1 - (P - g) @ np.linalg.solve(S, g)
    return l.min() + 50*pen

def opt(a, b, fn, ntrial=200, seed=0):
    rng = np.random.default_rng(seed)
    best = (-1e9, None)
    for t in range(ntrial):
        def mk(n):
            v = rng.uniform(0.6, 1.4, size=n); return v*np.radians(rng.uniform(130, 178))/v.sum()
        x0 = np.concatenate([mk(a), mk(b)])
        def obj(x):
            x = np.abs(x)
            if x.max() > 1.5 or x[:a].sum() > np.pi or x[a:].sum() > np.pi: return 10.0
            return -fn(x[:a], x[a:])
        res = minimize(obj, x0, method="Nelder-Mead", options={"xatol":1e-11,"fatol":1e-13,"maxiter":40000,"maxfev":40000})
        res = minimize(obj, res.x, method="Nelder-Mead", options={"xatol":1e-12,"fatol":1e-14,"maxiter":40000,"maxfev":40000})
        if -res.fun > best[0]: best = (-res.fun, np.abs(res.x))
    return best

def job(ab):
    a, b = ab
    s, x = opt(a, b, slack_affine)
    return a, b, s, x

if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor
    pairs = [(3,3),(4,3),(4,4),(5,4),(5,5),(6,5),(6,6),(7,6),(7,7),(6,4),(7,4),(7,5),(8,4),(8,5)]
    with ProcessPoolExecutor(10) as ex:
        for a, b, s, x in ex.map(job, pairs):
            print(f"affine arms {a}+{b}: k={a+b+1} (prism F={a+b+3}): best slack = {s:+.6f}   R {np.degrees(x[:a]).round(2)}  L {np.degrees(x[a:]).round(2)}", flush=True)
            np.save(f"explore/flat_affine_{a}_{b}.npy", x)
