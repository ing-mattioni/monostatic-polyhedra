"""Quick test: for uncertified V=7 trees, is the Gram relaxation with EXACT centroid weights infeasible?
 R >= 0 (7x7), tr R = 1, <S_i,R> <= t (shadow), R w = 0 with w = sum_T lam_T 1_T (cone triangulation from v).
 psi(lam) = min t.  Relaxation infeasible at lam iff psi(lam) > 0.  We minimise psi over the simplex of lam."""
import sys, json, os, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, 'lowerbound')
import numpy as np, cvxpy as cp
from enumerate_cases import graph_data
from scipy.optimize import minimize

V = 7
gds = graph_data(V, with_trees=False)
unc = json.load(open('lowerbound/uncertified_V7.json'))

def make_problem(v, parent):
    R = cp.Variable((V, V), PSD=True); t = cp.Variable()
    w = cp.Parameter(V)
    cons = [cp.trace(R) == 1, R @ w == 0]
    for i in range(V):
        if i == v: continue
        j = parent[i]
        cons.append(R[i, i] - R[i, j] <= t)
    prob = cp.Problem(cp.Minimize(t), cons)
    return prob, w, R

def psi_min(case, nsamp=60, seed=0):
    gd = gds[case['graph']]; v = case['root']; parent = case['parent']
    tets = gd.decomposition[v]
    ind = np.zeros((len(tets), V))
    for k, f in enumerate(tets):
        ind[k, list(f)] = 1; ind[k, v] = 1
    prob, w, R = make_problem(v, parent)
    def psi(lam):
        w.value = lam @ ind
        prob.solve(solver="CLARABEL")
        return prob.value
    rng = np.random.default_rng(seed)
    best = (1e9, None)
    cands = [np.eye(len(tets))[k] for k in range(len(tets))] + [np.ones(len(tets))/len(tets)]
    cands += [rng.dirichlet(np.ones(len(tets))*a) for a in (0.3, 1.0) for _ in range(nsamp//2)]
    vals = []
    for lam in cands:
        p = psi(lam); vals.append(p)
        if p < best[0]: best = (p, lam)
    # local refinement (softmax parametrisation)
    def f(z):
        e = np.exp(z - z.max()); return psi(e/e.sum())
    z0 = np.log(np.maximum(best[1], 1e-6))
    res = minimize(f, z0, method="Nelder-Mead", options={"maxfev": 400, "xatol": 1e-4, "fatol": 1e-9})
    if res.fun < best[0]:
        e = np.exp(res.x - res.x.max()); best = (res.fun, e/e.sum())
    return best, len(tets), min(vals)

if __name__ == "__main__":
    idx = [int(a) for a in sys.argv[1:]] or [0, 100, 300, 500, 700, 900, 1100, 1300, 1425]
    t0 = time.time()
    for k in idx:
        c = unc[k]
        (p, lam), nt, _ = psi_min(c)
        print(f"case {k}: graph {c['graph']} root {c['root']} parent {c['parent']} badface {c['face']}  ntets={nt}  min_lam psi = {p:+.3e}  at lam={np.round(lam,3)}", flush=True)
    print("time", time.time()-t0)
