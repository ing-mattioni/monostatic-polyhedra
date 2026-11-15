"""Maximise the angular margin of the symmetric flat pyramid (7+7 outline + apex), true 3D evaluation."""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np, cma, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, 'explore')
from mu3d import centroid
from scipy.spatial import QhullError

def build(x, m=7):
    dphi = np.abs(x[:m]); rr = np.abs(x[m:2*m]); zA, eps = x[2*m], abs(x[2*m+1])
    ph = np.cumsum(dphi)
    xs = rr*np.sin(ph); zs = rr*np.cos(ph)
    pts = [(0.0, 0.0, 1.0)] + [(xs[i], 0.0, zs[i]) for i in range(m)] + [(-xs[i], 0.0, zs[i]) for i in range(m)] + [(0.0, eps, zA)]
    return np.array(pts)

def angmargin(X):
    try:
        c, vol, hull = centroid(X)
    except QhullError:
        return -1.0, None
    if len(hull.vertices) < len(X): return -1.0, None
    U = X - c
    d = np.linalg.norm(U, axis=1)
    D = X[None, :, :] - X[:, None, :]            # D[i,j] = x_j - x_i
    num = np.einsum('ijk,ik->ij', D, U)
    den = np.linalg.norm(D, axis=2) * d[:, None]
    np.fill_diagonal(den, 1.0)
    M = num/den; np.fill_diagonal(M, -np.inf)
    m = M.max(axis=1)
    apex = np.argmax(d); m[apex] = np.inf
    return m.min(), m

ang7 = np.radians([34.85,20.71,18.86,18.47,18.88,20.05,22.21])
r0 = np.cumprod(np.cos(ang7))
X0 = np.load('explore/seed_V16.npy')
x0 = np.concatenate([ang7, r0*0.9995, [X0[-1,2], 0.02]])
print("start margin:", angmargin(build(x0))[0])
best = (-1, None)
for seed in range(6):
    es = cma.CMAEvolutionStrategy(x0, 2e-3, {"seed": seed+1, "verbose": -9, "maxfevals": 60000, "tolfun": 1e-15, "tolx": 1e-13})
    while not es.stop():
        xs = es.ask(); es.tell(xs, [-angmargin(build(np.asarray(x)))[0] for x in xs])
    val = -es.result.fbest
    print(f"seed {seed}: margin {val:.6e}", flush=True)
    if val > best[0]: best = (val, es.result.xbest)
    x0 = best[1]
X = build(best[1])
np.save('explore/V16_polished.npy', X)
np.set_printoptions(precision=6, suppress=True)
print("best angular margin (cos of angle slack):", best[0]); print(X)
print("params:", best[1])
