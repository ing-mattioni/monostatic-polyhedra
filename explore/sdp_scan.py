"""Scan all rooted trees on n nodes: is the 'falling pattern' realizable by a simplex Gram matrix?
G PSD (n x n), G 1 = 0, tr G = 1, G[i,parent(i)] - G[i,i] >= t for non-root i.  maximize t."""
import sys, time
import numpy as np
import cvxpy as cp

def rooted_trees(n):
    """Beyer-Hedetniemi: generate level sequences of all rooted trees with n nodes."""
    L = list(range(n))  # path
    while True:
        yield tuple(L)
        # find last position p with L[p] > 1 (levels start at 0 for root; here use 0-based: need L[p]>=2)
        p = n - 1
        while p > 0 and L[p] == 1:
            p -= 1
        if p == 0:
            return
        # find q < p with L[q] == L[p]-1 (parent of p)
        q = p - 1
        while L[q] != L[p] - 1:
            q -= 1
        for i in range(p, n):
            L[i] = L[i - (p - q)]

def parents(L):
    n = len(L); par = [-1]*n; last = {}
    for i, l in enumerate(L):
        if l > 0: par[i] = last[l-1]
        last[l] = i
    return par

def solve(par, solver="CLARABEL"):
    n = len(par)
    G = cp.Variable((n, n), PSD=True)
    t = cp.Variable()
    cons = [cp.sum(G, axis=1) == 0, cp.trace(G) == 1]
    for i in range(n):
        if par[i] >= 0:
            cons.append(G[i, par[i]] - G[i, i] >= t)
    prob = cp.Problem(cp.Maximize(t), cons)
    prob.solve(solver=solver)
    return t.value, G.value

if __name__ == "__main__":
    n = int(sys.argv[1])
    t0 = time.time(); best = []
    cnt = 0
    for L in rooted_trees(n):
        cnt += 1
        par = parents(L)
        tv, G = solve(par)
        best.append((tv, L))
    best.sort(reverse=True)
    print(f"n={n} (dimension d={n-1}): {cnt} rooted trees, time {time.time()-t0:.1f}s")
    for tv, L in best[:15]:
        print(f"  t*={tv:+.3e}  levels={L}")
    print("  number with t*>1e-7:", sum(1 for tv,_ in best if tv>1e-7))
