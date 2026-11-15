"""Mono-stable pentahedra?  Affine shape x metric formulation.
Simple pentahedron = triangular frustum: apex O=0, A_i = t_i e_i, B_i = e_i (0<t_i<1), faces Q1,Q2,Q3 (x_i>=0), B (sum x<=1),
A (sum x_i/t_i >= 1).  Centroid g_i = (1 - tau t_i)/(4(1-tau)), tau = t1 t2 t3.
With origin at the centroid the planes are l_i . x <= 1.  For an inner product with inverse Gram W (3x3 PD) the foot of
the centroid on plane i is W l_i/(l_i W l_i); face i is NOT stable iff  l_i W (l_j - l_i) > 0 for some adjacent j.
For a rooted spanning tree pi of the face-adjacency graph:  SDP  max m : l_i W (l_pi(i) - l_i) >= m * |l_i| |l_pi(i)-l_i|, tr W = 1, W >= 0."""
import os, sys, itertools, time
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np, cvxpy as cp

FACES = ['Q1', 'Q2', 'Q3', 'A', 'B']
ADJ = {0: [1, 2, 3, 4], 1: [0, 2, 3, 4], 2: [0, 1, 3, 4], 3: [0, 1, 2], 4: [0, 1, 2]}

def covectors(t):
    t = np.asarray(t, float); tau = t.prod()
    g = (1 - tau*t)/(4*(1 - tau))
    L = np.zeros((5, 3))
    for i in range(3): L[i, i] = -1/g[i]
    L[3] = -(1/t)/((g/t).sum() - 1)
    L[4] = np.ones(3)/(1 - g.sum())
    return L, g

def rooted_trees():
    out = []
    for s in range(5):
        others = [i for i in range(5) if i != s]
        for choice in itertools.product(*[ADJ[i] for i in others]):
            par = {i: j for i, j in zip(others, choice)}
            ok = True
            for i in others:
                x, k = i, 0
                while x != s and k < 6: x = par[x]; k += 1
                if x != s: ok = False; break
            if ok: out.append((s, par))
    return out

TREES = rooted_trees()

W = cp.Variable((3, 3), PSD=True); m = cp.Variable()
Cs = [cp.Parameter((3, 3), symmetric=True) for _ in range(4)]
prob = cp.Problem(cp.Maximize(m), [cp.trace(W) == 1] + [cp.trace(C @ W) >= m for C in Cs])

def best_margin(L, tree):
    s, par = tree
    k = 0
    for i, j in par.items():
        d = L[j] - L[i]
        C = np.outer(L[i], d); C = (C + C.T)/2
        C = C/(np.linalg.norm(L[i])*np.linalg.norm(d))
        Cs[k].value = C; k += 1
    prob.solve(solver="CLARABEL")
    return m.value

if __name__ == "__main__":
    print(len(TREES), "rooted spanning trees")
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    worst = (-9, None, None)
    t0 = time.time()
    for it in range(n):
        if it % 3 == 0: t = rng.uniform(0.02, 0.98, size=3)
        elif it % 3 == 1: t = 1 - rng.uniform(0.001, 0.2, size=3)**1.5
        else: t = rng.uniform(0.001, 0.2, size=3)**1.5 * rng.choice([1, 1, 5], size=3); t = np.clip(t, 1e-4, 0.999)
        L, g = covectors(t)
        vals = [best_margin(L, tr) for tr in TREES]
        k = int(np.argmax(vals))
        if vals[k] > worst[0]: worst = (vals[k], t.copy(), TREES[k])
        if it < 12 or it % 10 == 0:
            print(f"t={np.round(t,4)}: max over trees of SDP margin = {vals[k]:+.3e}  (tree root {FACES[TREES[k][0]]}, {dict((FACES[i],FACES[j]) for i,j in TREES[k][1].items())})", flush=True)
    print("WORST (largest) margin overall:", worst[0], "at t =", worst[1], "tree", worst[2], f" time {time.time()-t0:.0f}s")
