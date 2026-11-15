"""Flat-pyramid family, no symmetry assumed.  Outline polygon B = conv(p_1..p_k) in the plane,
apex A' inside B (lifted by eps in 3D).  Reference point c' = (3 g + A')/4, g = area centroid of B.
f = min over outline vertices (except farthest) of the planar domination margin (dominators: outline
vertices and A').  f > 0  ==> flat pyramid with k+1 vertices is mono-unstable."""
import numpy as np, sys, os, time
os.environ.setdefault("OMP_NUM_THREADS", "1")

def hull2d(P):
    idx = np.lexsort((P[:,1], P[:,0]))
    def cross(o, a, b): return (P[a,0]-P[o,0])*(P[b,1]-P[o,1]) - (P[a,1]-P[o,1])*(P[b,0]-P[o,0])
    lower = []
    for i in idx:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], i) <= 0: lower.pop()
        lower.append(i)
    upper = []
    for i in idx[::-1]:
        while len(upper) >= 2 and cross(upper[-2], upper[-1], i) <= 0: upper.pop()
        upper.append(i)
    return lower[:-1] + upper[:-1]     # ccw

def fobj(x, k):
    P = x[:2*k].reshape(k, 2); a = x[2*k:2*k+2]
    h = hull2d(P)
    H = P[h]
    x0, z0 = H[:,0], H[:,1]; x1, z1 = np.roll(x0,-1), np.roll(z0,-1)
    cr = x0*z1 - x1*z0; A = cr.sum()/2
    if A <= 1e-12: return -1e3, None
    g = np.array([((x0+x1)*cr).sum(), ((z0+z1)*cr).sum()])/(6*A)
    # apex must be inside B: compute violation
    ex, ez = x1-x0, z1-z0
    s = ex*(a[1]-z0) - ez*(a[0]-x0)      # >=0 inside (ccw)
    pen = min(0.0, s.min()/np.sqrt(A))
    c = (3*g + a)/4
    U = np.vstack([H - c, (a - c)[None,:]])
    G = U @ U.T
    d = np.diag(G).copy()
    nh = len(h)
    # angular margin: cos of angle at vertex i between (v_j - v_i) and (v_i - c)
    num = G[:nh,:] - d[:nh,None]                       # (u_j - u_i).u_i
    dist = np.sqrt(np.maximum(d[:nh,None] + d[None,:] - 2*G[:nh,:], 1e-300))
    M = num/(dist*np.sqrt(d[:nh,None]))
    M[np.arange(nh), np.arange(nh)] = -np.inf
    m = M.max(axis=1)
    apex = np.argmax(d[:nh])
    m[apex] = np.inf
    # the pyramid apex A' is itself a vertex and must be dominated by an outline vertex
    numa = G[nh,:nh] - d[nh]
    dista = np.sqrt(np.maximum(d[nh] + d[:nh] - 2*G[nh,:nh], 1e-300))
    ma = (numa/(dista*np.sqrt(max(d[nh], 1e-300)))).max()
    val = min(m.min(), ma)
    if pen < 0: val = min(val, 0.0) + 10*pen      # apex must project inside B
    return val - 1e-3*(k-nh), (h, c, g, m)

def run(args):
    k, seed, x0, sigma, maxevals = args
    import cma, warnings; warnings.filterwarnings("ignore")
    es = cma.CMAEvolutionStrategy(x0, sigma, {"seed": seed+1, "verbose": -9, "maxfevals": maxevals,
                                              "tolfun":1e-14, "tolx":1e-13, "tolfunhist":1e-14})
    while not es.stop():
        xs = es.ask(); es.tell(xs, [-fobj(np.asarray(x), k)[0] for x in xs])
    return -es.result.fbest, es.result.xbest, seed

def sym_seed(dphi_deg_right, dphi_deg_left):
    def arm(dd, sgn):
        d = np.radians(dd); ph = np.cumsum(d); r = np.cumprod(np.cos(d))
        return [(sgn*r[i]*np.sin(ph[i]), r[i]*np.cos(ph[i])) for i in range(len(d))]
    P = [(0.0,1.0)] + arm(dphi_deg_right, 1) + arm(dphi_deg_left, -1)
    P = np.array(P)
    return P

if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor
    kR, kL = int(sys.argv[1]), int(sys.argv[2])
    nrest = int(sys.argv[3]) if len(sys.argv) > 3 else 32
    ang7 = [34.85,20.71,18.86,18.47,18.88,20.05,22.21]
    rng = np.random.default_rng(1)
    tasks = []
    for s in range(nrest):
        def mk(n):
            a = rng.uniform(0.6, 1.4, size=n); a *= rng.uniform(140, 175)/a.sum(); return list(a)
        P = sym_seed(mk(kR), mk(kL))
        a0 = np.array([0.0, P[:,1].min()*0.8])
        x0 = np.concatenate([P.ravel(), a0])
        tasks.append((len(P), s, x0, 0.05, 300000))
    k = 1 + kR + kL
    t0 = time.time(); best = None
    with ProcessPoolExecutor(16) as ex:
        for val, xb, seed in ex.map(run, tasks):
            if best is None or val > best[0]: best = (val, xb)
    val, xb = best
    v, (h, c, g, m) = fobj(xb, k)
    print(f"arms {kR}+{kL}: k={k}, V={k+1}: best f={val:+.5e}  hull size={len(h)}  time={time.time()-t0:.0f}s")
    np.save(f"explore/flat2d_k{k}_{kR}_{kL}.npy", xb)
