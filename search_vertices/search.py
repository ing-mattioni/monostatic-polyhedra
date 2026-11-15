"""Optimisation drivers: lean CMA-ES (batch numba evaluation), tree-fixed SLSQP polish,
structured seeds, vertex-removal ladder.  Results are logged to results.jsonl."""
import os, sys, json, time, argparse
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import core
from numba import njit

ROOT = os.path.dirname(HERE)
SOLDIR = os.path.join(HERE, "sols")
LOG = os.path.join(HERE, "results.jsonl")
DMIN_REL = 1e-3
VOL_REL = 1e-7


# ------------------------------------------------------------------ CMA-ES
def cma_run(fbatch, x0, sigma0, maxevals, seed=0, popsize=None, tolx=1e-13, stall=400, target=None):
    """Maximise fbatch (maps (lam,n) -> (lam,)).  Returns (fbest, xbest, nevals)."""
    rng = np.random.default_rng(seed)
    x0 = np.asarray(x0, dtype=np.float64)
    n = x0.size
    lam = popsize or (4 + int(3 * np.log(n)))
    mu = lam // 2
    w = np.log(mu + 0.5) - np.log(np.arange(1, mu + 1)); w /= w.sum(); mueff = 1 / np.sum(w ** 2)
    cc = (4 + mueff / n) / (n + 4 + 2 * mueff / n); cs = (mueff + 2) / (n + mueff + 5)
    c1 = 2 / ((n + 1.3) ** 2 + mueff); cmu = min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((n + 2) ** 2 + mueff))
    damps = 1 + 2 * max(0, np.sqrt((mueff - 1) / (n + 1)) - 1) + cs
    chiN = np.sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n * n))
    m = x0.copy(); sigma = sigma0
    pc = np.zeros(n); ps = np.zeros(n); C = np.eye(n); B = np.eye(n); D = np.ones(n); invsqrtC = np.eye(n)
    eigeneval = 0; counteval = 0; g = 0; last_imp = 0
    fbest = fbatch(x0[None, :])[0]; xbest = x0.copy()
    while counteval < maxevals:
        Zs = rng.standard_normal((lam, n))
        Ys = (Zs * D) @ B.T
        Xs = m + sigma * Ys
        f = fbatch(Xs)
        counteval += lam; g += 1
        idx = np.argsort(-f)
        if f[idx[0]] > fbest:
            if f[idx[0]] > fbest + 1e-12 * max(1e-3, abs(fbest)):
                last_imp = g
            fbest = f[idx[0]]; xbest = Xs[idx[0]].copy()
        if target is not None and fbest >= target:
            break
        Ysel = Ys[idx[:mu]]
        ymean = w @ Ysel
        m = m + sigma * ymean
        ps = (1 - cs) * ps + np.sqrt(cs * (2 - cs) * mueff) * (invsqrtC @ ymean)
        hsig = np.linalg.norm(ps) / np.sqrt(1 - (1 - cs) ** (2 * counteval / lam)) / chiN < 1.4 + 2 / (n + 1)
        pc = (1 - cc) * pc + hsig * np.sqrt(cc * (2 - cc) * mueff) * ymean
        C = (1 - c1 - cmu) * C + c1 * (np.outer(pc, pc) + (1 - hsig) * cc * (2 - cc) * C) + cmu * (Ysel.T * w) @ Ysel
        sigma *= np.exp((cs / damps) * (np.linalg.norm(ps) / chiN - 1))
        if counteval - eigeneval > lam / (c1 + cmu) / n / 10:
            eigeneval = counteval
            C = np.triu(C) + np.triu(C, 1).T
            D2, B = np.linalg.eigh(C)
            D = np.sqrt(np.maximum(D2, 1e-30))
            invsqrtC = (B / D) @ B.T
        if sigma * D.max() < tolx or g - last_imp > stall or D.max() > 1e7 * D.min() or not np.isfinite(sigma):
            break
    return fbest, xbest, counteval


def fb3d(V):
    def f(Z):
        return core.objective_batch(np.ascontiguousarray(Z), V, core.JOG, DMIN_REL, VOL_REL)
    return f


def fbflat(V, n_out=0):
    def f(Z):
        return core.objective_flat_batch(np.ascontiguousarray(Z), V, core.JOG, DMIN_REL, n_out)
    return f


# ------------------------------------------------------------------ polish
@njit(cache=True)
def cons_fixed(X, par, far, JOG):
    vol, c, nt, isv, ok = core.hull_centroid(X, JOG)
    V = X.shape[0]
    out = np.empty(V - 1)
    k = 0
    for i in range(V):
        if i == far:
            continue
        if not ok:
            out[k] = -1.0; k += 1; continue
        j = par[i]
        ux = X[i, 0] - c[0]; uy = X[i, 1] - c[1]; uz = X[i, 2] - c[2]
        wx = X[j, 0] - X[i, 0]; wy = X[j, 1] - X[i, 1]; wz = X[j, 2] - X[i, 2]
        out[k] = (wx * ux + wy * uy + wz * uz) / (np.sqrt(wx * wx + wy * wy + wz * wz) * np.sqrt(ux * ux + uy * uy + uz * uz) + 1e-300)
        k += 1
    return out


def slp_polish(X, iters=300, delta=0.01, h=1e-7, min_delta=1e-9, multi=2):
    """Sequential LP with trust region on the max-min of fixed-parent cosine margins.
    For each i != far up to `multi` best parents are linearised; the LP maximises t with
    t <= g_ij + grad.d for the CURRENT best parent of each i (others ignored)."""
    from scipy.optimize import linprog
    X = core.normalize(X); V = len(X); n = 3 * V
    Fcur = core.F(X)
    for it in range(iters):
        inf = core.info(X); far = inf["far"]; par = inf["parent"].astype(np.int64)
        g0 = cons_fixed(np.ascontiguousarray(X), par, far, core.JOG)
        G = np.empty((len(g0), n))
        x0 = X.ravel()
        for k in range(n):
            xp = x0.copy(); xp[k] += h
            G[:, k] = (cons_fixed(np.ascontiguousarray(xp.reshape(V, 3)), par, far, core.JOG) - g0) / h
        # variables: d (n), t ; maximise t
        c = np.zeros(n + 1); c[-1] = -1.0
        A = np.hstack([-G, np.ones((len(g0), 1))]); b = g0
        bounds = [(-delta, delta)] * n + [(None, None)]
        res = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method="highs")
        if res.status != 0:
            delta *= 0.5
            if delta < min_delta: break
            continue
        Xn = core.normalize((x0 + res.x[:n]).reshape(V, 3))
        r, diam = core.shape_stats(Xn)
        Fn = core.F(Xn) if r > DMIN_REL else -1.0
        pred = res.x[-1] - Fcur
        if Fn > Fcur + 1e-15:
            ratio = (Fn - Fcur) / max(pred, 1e-18)
            X, Fcur = Xn, Fn
            if ratio > 0.5: delta = min(delta * 1.5, 0.05)
        else:
            delta *= 0.4
            if delta < min_delta: break
    return X, Fcur


def polish(X, rounds=4, box=0.02, maxiter=300):
    Xs, Fs = slp_polish(X, iters=150, delta=min(box, 0.01))
    return Xs, Fs


def polish_slsqp(X, rounds=4, box=0.02, maxiter=300):
    """Tree-fixed SLSQP: maximise t s.t. cos-margin(i, parent(i)) >= t.  Returns improved X."""
    from scipy.optimize import minimize
    X = core.normalize(X)
    V = len(X)
    bestF = core.F(X); bestX = X.copy()
    for r in range(rounds):
        inf = core.info(bestX)
        far = inf["far"]; par = inf["parent"].astype(np.int64)
        X0 = bestX.copy()
        z0 = np.append(X0.ravel(), inf["F"])

        def cons(z):
            Xz = np.ascontiguousarray(z[:-1].reshape(V, 3))
            g = cons_fixed(Xz, par, far, core.JOG) - z[-1]
            rr, diam = core.shape_stats(Xz)
            return np.append(g, rr - DMIN_REL)
        bounds = [(v - box, v + box) for v in X0.ravel()] + [(-1, 1)]
        try:
            res = minimize(lambda z: -z[-1], z0, jac=lambda z: np.r_[np.zeros(3 * V), -1.0],
                           constraints=[{"type": "ineq", "fun": cons}], bounds=bounds, method="SLSQP",
                           options=dict(maxiter=maxiter, ftol=1e-15))
        except Exception:
            break
        Xn = res.x[:-1].reshape(V, 3)
        Fn = core.F(Xn)
        if Fn > bestF + 1e-14:
            bestF = Fn; bestX = core.normalize(Xn)
        else:
            break
    return bestX, bestF


# ------------------------------------------------------------------ seeds
def spiral_polygon(dR, dL):
    """Tight double spiral, reference O=0, top T=(0,1); returns planar (k,2) as (x,z): T, left arm, right arm reversed (ccw)."""
    def arm(d, sgn):
        ph = np.cumsum(d); r = np.cumprod(np.cos(d))
        return np.stack([sgn * r * np.sin(ph), r * np.cos(ph)], axis=1)
    R = arm(np.asarray(dR), 1.0); L = arm(np.asarray(dL), -1.0)
    return np.vstack([[0.0, 1.0], L, R[::-1]])


def random_arm_angles(n, rng, total=None):
    v = rng.uniform(0.6, 1.4, size=n)
    tot = np.radians(total if total is not None else rng.uniform(140, 176))
    return v * tot / v.sum()


def arm_angles(a, b, rng):
    """Use the optimised flat-pyramid spiral angles when available, else random."""
    f = os.path.join(ROOT, "explore", f"flat_asym_{a}_{b}.npy")
    if os.path.exists(f) and rng.random() < 0.6:
        x = np.load(f)
        return x[:a] * rng.uniform(0.97, 1.03, a), x[a:] * rng.uniform(0.97, 1.03, b)
    f = os.path.join(ROOT, "explore", f"flat_asym_{b}_{a}.npy")
    if os.path.exists(f) and rng.random() < 0.6:
        x = np.load(f)
        return x[b:] * rng.uniform(0.97, 1.03, a), x[:b] * rng.uniform(0.97, 1.03, b)
    return random_arm_angles(a, rng), random_arm_angles(b, rng)


HEIGHT_PATTERNS = ("alt", "armpm", "random", "grow", "growalt", "zero_top_alt", "bottom_pm")


def heights(pattern, a, b, rng):
    """Heights for polygon order [T, L_1..L_b, R_a..R_1]."""
    k = 1 + a + b
    iL = np.arange(1, b + 1)                 # left arm, index along arm = 1..b
    iR = np.arange(a, 0, -1)                 # right arm reversed
    along = np.r_[0, iL / b, iR / a]         # 0 at T, 1 at bottom
    side = np.r_[0, -np.ones(b), np.ones(a)]
    idx = np.r_[0, iL, iR]
    if pattern == "alt":
        y = (-1.0) ** idx
    elif pattern == "armpm":
        y = side.copy()
    elif pattern == "random":
        y = rng.normal(size=k)
    elif pattern == "grow":
        y = side * along ** 1.5
    elif pattern == "growalt":
        y = (-1.0) ** idx * along
    elif pattern == "zero_top_alt":
        y = (-1.0) ** idx * np.minimum(1, 2 * along)
    elif pattern == "bottom_pm":
        y = np.zeros(k); y[b] = 1.0; y[b + 1] = -1.0
    y = y + 0.15 * rng.normal(size=k)
    y[0] = 0.0
    return y


def seed_spiral(V, rng, a=None, pattern=None, eps=None, apex=False):
    k = V - 1 if apex else V
    if a is None:
        a = int(rng.integers((k - 1) // 2, k - 1 - 1))   # right arm >= left arm
    b = k - 1 - a
    dR, dL = arm_angles(a, b, rng)
    P = spiral_polygon(dR, dL)
    pattern = pattern or HEIGHT_PATTERNS[rng.integers(len(HEIGHT_PATTERNS))]
    eps = eps if eps is not None else [0.02, 0.1, 0.3, 1.0][rng.integers(4)]
    y = heights(pattern, a, b, rng)
    X = np.stack([P[:, 0], eps * y, P[:, 1]], axis=1)
    if apex:
        zb = P[:, 1].min()
        X = np.vstack([X, [rng.normal() * 0.05, eps * rng.uniform(0.5, 1.5), zb * rng.uniform(0.6, 0.95)]])
    return X, dict(a=a, b=b, pattern=pattern, eps=eps, apex=apex)


def seed_flatpyr(V, rng):
    return seed_spiral(V, rng, apex=True)


def seed_v16_minus(V, rng):
    X = np.load(os.path.join(ROOT, "explore", "seed_V16.npy"))
    drops = rng.choice(16, size=16 - V, replace=False)
    X = np.delete(X, drops, axis=0)
    X = X + 1e-3 * rng.normal(size=X.shape)
    return X, dict(drop=[int(d) for d in drops])


def seed_random(V, rng):
    X = rng.normal(size=(V, 3)); X /= np.linalg.norm(X, axis=1)[:, None]
    X *= rng.uniform(0.6, 1.2, size=(V, 1))
    return X, {}


def seed_reshetov(V, rng):
    """Reshetov-dual template: 14 points (x thin). For V != 14 delete / insert points."""
    X = np.load(os.path.join(ROOT, "explore", "reshetov_dual14.npy"))
    info = {}
    if V < 14:
        drops = rng.choice(14, size=14 - V, replace=False)
        drops = [d for d in drops if d != 10]
        while len(drops) < 14 - V:
            d = int(rng.integers(14))
            if d != 10 and d not in drops: drops.append(d)
        X = np.delete(X, drops, axis=0); info["drop"] = [int(d) for d in drops]
    elif V > 14:
        hv = np.arange(14)
        for t in range(V - 14):   # insert midpoint of a random pair of nearby points, pushed outward
            i = int(rng.integers(len(X))); d = np.linalg.norm(X - X[i], axis=1); d[i] = 1e9; j = int(np.argmin(d))
            p = 0.5 * (X[i] + X[j]); p[1:] *= 1.02
            X = np.vstack([X, p])
        info["insert"] = V - 14
    X = X + rng.uniform(0, 3e-3) * rng.normal(size=X.shape)
    return X, info


SEEDS = dict(spiral=seed_spiral, flatpyr=seed_flatpyr, v16minus=seed_v16_minus, random=seed_random, reshetov=seed_reshetov)


# ------------------------------------------------------------------ tasks
def log(rec):
    with open(LOG, "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def save_sol(X, tag):
    os.makedirs(SOLDIR, exist_ok=True)
    p = os.path.join(SOLDIR, tag + ".npy")
    np.save(p, X)
    return p


def run_3d(X0, sigmas, maxevals, seed, restarts=1, popsize=None, do_polish=True):
    """CMA from X0 (several sigmas / restarts), keep best, polish."""
    V = len(X0)
    f = fb3d(V)
    X0 = core.normalize(X0)
    best = (-1e9, X0.ravel()); tot = 0
    for r in range(restarts):
        for s in sigmas:
            start = best[1] if (r > 0 and best[0] > -1e8) else X0.ravel()
            fv, xv, ne = cma_run(f, start, s, maxevals, seed=seed * 1000 + r * 10 + int(s * 1e4) % 97, popsize=popsize)
            tot += ne
            if fv > best[0]:
                best = (fv, xv)
    X = core.normalize(best[1].reshape(V, 3))
    if do_polish:
        Xp, Fp = polish(X)
        if Fp > core.F(X):
            X = Xp
    return X, tot


def task(args):
    """Generic worker. args: dict(kind, V, seed, ...)."""
    t0 = time.time()
    kind = args["kind"]; V = args["V"]; seed = args["seed"]
    rng = np.random.default_rng(seed)
    meta = {}
    if kind == "flat":       # flat-limit surrogate first, then 3D at several eps
        Xs, meta = SEEDS[args.get("family", "spiral")](V, rng)
        P = Xs[:, [0, 2]].copy(); Y = Xs[:, 1] / max(1e-9, meta.get("eps", 1.0))
        z0 = np.r_[P.ravel(), Y]
        fv, zb, ne = cma_run(fbflat(V), z0, args.get("sigma", 0.05), args.get("maxevals", 150000), seed=seed)
        P = zb[:2 * V].reshape(V, 2); Y = zb[2 * V:]
        Fflat = float(core.evalF_flat(np.ascontiguousarray(P), np.ascontiguousarray(Y), core.JOG)[0])
        best = None
        for eps in (0.003, 0.01, 0.03, 0.1):
            X = np.stack([P[:, 0], eps * Y * np.ptp(P) / max(1e-9, np.ptp(Y)), P[:, 1]], axis=1)
            Fx = core.F(X)
            if best is None or Fx > best[0]:
                best = (Fx, X, eps)
        X, ne2 = run_3d(best[1], [args.get("sigma3", 0.003)], args.get("maxevals3", 100000), seed)
        meta.update(Fflat=Fflat, eps_best=best[2], F_at_eps=float(best[0]))
        ne += ne2
    elif kind == "seed3d":
        Xs, meta = SEEDS[args["family"]](V, rng)
        X, ne = run_3d(Xs, args.get("sigmas", [0.05]), args.get("maxevals", 200000), seed, restarts=args.get("restarts", 1))
    elif kind == "refine":   # from a given npy file
        Xs = np.load(args["file"])
        if args.get("drop") is not None:
            Xs = np.delete(Xs, args["drop"], axis=0)
        Xs = core.normalize(Xs) + args.get("jitter", 0.0) * rng.normal(size=Xs.shape)
        X, ne = run_3d(Xs, args.get("sigmas", [0.003]), args.get("maxevals", 200000), seed, restarts=args.get("restarts", 1))
        meta = dict(file=os.path.basename(args["file"]), drop=args.get("drop"))
    else:
        raise ValueError(kind)
    inf = core.info(X)
    tag = f"V{V}_{kind}_{args.get('family', '')}_{seed}"
    path = save_sol(X, tag)
    rec = dict(kind=kind, family=args.get("family"), V=V, seed=seed, F=inf["F"], nvert=inf["nvert"],
               dmin_rel=inf["dmin_rel"], evals=int(ne), time=round(time.time() - t0, 1), file=path,
               meta={k: (v if not isinstance(v, np.generic) else v.item()) for k, v in meta.items()})
    log(rec)
    return rec


def run_pool(tasks, nproc=6, verbose=True):
    from concurrent.futures import ProcessPoolExecutor, as_completed
    out = []
    with ProcessPoolExecutor(nproc) as ex:
        futs = [ex.submit(task, t) for t in tasks]
        for fu in as_completed(futs):
            try:
                r = fu.result()
            except Exception as e:
                print("task failed:", e, flush=True); continue
            out.append(r)
            if verbose:
                print(f"V={r['V']} {r['kind']}/{r['family']} seed={r['seed']}: F={r['F']:+.4e} nvert={r['nvert']} "
                      f"evals={r['evals']} t={r['time']}s meta={r['meta']}", flush=True)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode")
    ap.add_argument("--V", type=int, default=15)
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--family", default="spiral")
    ap.add_argument("--file", default=None)
    ap.add_argument("--sigma", type=float, default=0.003)
    ap.add_argument("--maxevals", type=int, default=200000)
    ap.add_argument("--restarts", type=int, default=1)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--nproc", type=int, default=6)
    ap.add_argument("--jitter", type=float, default=0.0)
    a = ap.parse_args()
    t0 = time.time()
    if a.mode == "refine":
        tasks = [dict(kind="refine", V=len(np.load(a.file)), seed=a.seed0 + s, file=a.file, jitter=a.jitter,
                      sigmas=[a.sigma * (1 + (s % 3))], maxevals=a.maxevals, restarts=a.restarts) for s in range(a.n)]
    elif a.mode == "ladder":
        X = np.load(a.file); Vn = len(X) - 1
        tasks = [dict(kind="refine", V=Vn, seed=a.seed0 + d, file=a.file, drop=d, jitter=a.jitter,
                      sigmas=[a.sigma], maxevals=a.maxevals, restarts=a.restarts) for d in range(len(X))]
    elif a.mode == "seed3d":
        tasks = [dict(kind="seed3d", V=a.V, seed=a.seed0 + s, family=a.family, sigmas=[a.sigma],
                      maxevals=a.maxevals, restarts=a.restarts) for s in range(a.n)]
    elif a.mode == "flat":
        tasks = [dict(kind="flat", V=a.V, seed=a.seed0 + s, family=a.family, sigma=0.05, maxevals=a.maxevals,
                      sigma3=a.sigma, maxevals3=a.maxevals) for s in range(a.n)]
    else:
        raise SystemExit("unknown mode")
    res = run_pool(tasks, a.nproc)
    res.sort(key=lambda r: -r["F"])
    print(f"done {len(res)} tasks in {time.time()-t0:.0f}s; best F={res[0]['F']:+.4e} file={res[0]['file']}")
