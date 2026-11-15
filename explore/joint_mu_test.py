"""Joint certificate test: for a (graph, root v, tree) find mu in R^V and, for every tetrahedron T of the cone
triangulation from v, multipliers y^T >= 0 such that  M_T = sum_i y^T_i S_i + sym(mu 1_T^T)  is positive definite.
(S_i = E_ii - sym(E_{i,parent(i)}).)  Then no PSD Gram matrix R != 0 satisfies the shadow inequalities <S_i,R> <= 0
together with R w = 0 for any w in the cone spanned by the 1_T  ==> no mono-unstable polyhedron with this data."""
import sys, json, os, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, 'lowerbound')
import numpy as np, cvxpy as cp
from enumerate_cases import graph_data

def joint_sdp(V, v, parent, tets):
    nT = len(tets)
    mu = cp.Variable(V); z = cp.Variable()
    ys = [cp.Variable(V, nonneg=True) for _ in range(nT)]
    cons = []
    I = np.eye(V)
    for k, f in enumerate(tets):
        one = np.zeros(V); one[list(f)] = 1; one[v] = 1
        M = 0
        for i in range(V):
            if i == v:
                cons.append(ys[k][i] == 0); continue
            j = parent[i]
            S = np.zeros((V, V)); S[i, i] = 1; S[i, j] -= 0.5; S[j, i] -= 0.5
            M = M + ys[k][i]*S
        O = np.outer(np.ones(V), one)      # (mu 1_T^T)_{ab} = mu_a one_b
        A = cp.reshape(mu, (V, 1), order='F') @ one.reshape(1, V)
        M = M + (A + A.T)/2
        cons.append(M >> z*I)
        cons.append(cp.sum(ys[k]) <= 1)
    cons.append(cp.norm(mu, 'inf') <= 1)
    prob = cp.Problem(cp.Maximize(z), cons)
    prob.solve(solver="CLARABEL")
    return z.value, mu.value, [y.value for y in ys]

if __name__ == "__main__":
    V = 7
    gds = graph_data(V, with_trees=False)
    unc = json.load(open('lowerbound/uncertified_V7.json'))
    seen = {}
    t0 = time.time()
    zs = []
    for k, c in enumerate(unc):
        key = (c['graph'], c['root'], tuple(c['parent']))
        if key in seen: continue
        gd = gds[c['graph']]
        zval, mu, ys = joint_sdp(V, c['root'], c['parent'], gd.decomposition[c['root']])
        seen[key] = zval; zs.append(zval)
        if len(zs) <= 5 or len(zs) % 200 == 0:
            print(f"tree {len(zs)}: graph {key[0]} root {key[1]} parent {key[2]}: z = {zval:+.4e}  mu = {np.round(mu,3)}", flush=True)
    zs = np.array(zs)
    print(f"{len(zs)} distinct trees; min z = {zs.min():+.4e}; # with z <= 1e-7: {(zs <= 1e-7).sum()}; time {time.time()-t0:.0f}s")
