"""Validate core.py against the scipy reference explore/mu3d.py and time it."""
import os, sys, time
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "explore"))
import core
from mu3d import centroid as ref_centroid
from scipy.spatial import ConvexHull


def ref_F(X):
    c, vol, hull = ref_centroid(X)
    U = X - c
    d = np.linalg.norm(U, axis=1)
    far = np.argmax(d)
    m = np.full(len(X), np.inf)
    for i in range(len(X)):
        best = -2
        for j in range(len(X)):
            if j == i: continue
            w = X[j] - X[i]
            best = max(best, w @ U[i] / (np.linalg.norm(w) * d[i]))
        m[i] = best
    return min(m[i] for i in range(len(X)) if i != far), c, vol, set(hull.vertices.tolist())


rng = np.random.default_rng(0)
maxdc = 0; maxdF = 0; maxdv = 0; nbad_v = 0
tests = []
for t in range(300):
    V = rng.integers(5, 25)
    kind = t % 3
    if kind == 0:
        X = rng.normal(size=(V, 3))
    elif kind == 1:
        X = rng.normal(size=(V, 3)); X /= np.linalg.norm(X, axis=1)[:, None]
    else:   # flat-ish + coplanar subsets
        X = rng.normal(size=(V, 3)); X[: V // 2, 1] = 0.0; X[:, 1] *= 0.05
    tests.append(X)
for name in ("seed_V16.npy", "seed_V18.npy", "reshetov_dual14.npy"):
    tests.append(np.load(os.path.join(HERE, "..", "explore", name)))
for X in tests:
    X = np.ascontiguousarray(X)
    Fr, cr, vr, hv = ref_F(X)
    Fc, vol, c, far, wi = core.evalF(X, core.JOG)
    vol2, c2, nt, isv, ok = core.hull_centroid(X, core.JOG)
    maxdc = max(maxdc, np.abs(c - cr).max()); maxdv = max(maxdv, abs(vol - vr) / vr)
    maxdF = max(maxdF, abs(Fc - Fr))
    if set(np.nonzero(isv)[0].tolist()) != hv: nbad_v += 1
print(f"{len(tests)} tests: max|dc|={maxdc:.2e}  max rel dvol={maxdv:.2e}  max|dF|={maxdF:.2e}  vertex-set mismatches={nbad_v}")
for name in ("seed_V16.npy", "seed_V18.npy", "reshetov_dual14.npy"):
    X = np.load(os.path.join(HERE, "..", "explore", name))
    i = core.info(X)
    print(name, "F=%.4e" % i["F"], "nvert", i["nvert"], "far", i["far"], "unstable", i["unstable"])
# timing
for name in ("seed_V16.npy", "reshetov_dual14.npy"):
    X = np.ascontiguousarray(np.load(os.path.join(HERE, "..", "explore", name)))
    Z = np.repeat(X.ravel()[None, :], 2000, axis=0) + 1e-6 * rng.normal(size=(2000, X.size))
    core.objective_batch(Z[:2], len(X), core.JOG, 1e-3, 1e-7)
    t0 = time.perf_counter(); core.objective_batch(Z, len(X), core.JOG, 1e-3, 1e-7); dt = time.perf_counter() - t0
    print(f"{name}: {dt/len(Z)*1e6:.1f} us / eval (batch)")
    X = rng.normal(size=(16, 3)); Z = np.repeat(X.ravel()[None, :], 2000, axis=0) + 1e-3 * rng.normal(size=(2000, 48))
    t0 = time.perf_counter(); core.objective_batch(Z, 16, core.JOG, 1e-3, 1e-7); dt = time.perf_counter() - t0
print(f"random V=16: {dt/len(Z)*1e6:.1f} us / eval (batch)")
