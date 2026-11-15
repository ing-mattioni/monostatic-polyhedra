"""Flat-limit (generalised flat family) search.
Variables: planar positions p_i (x,z) and heights y_i; c' = xz-projection of centroid of conv{(p_i,y_i)}.
Objective: planar angular margin w.r.t. c' (F_flat > 0 => mono-unstable for all small eps).
Then realise in 3D at small eps and polish with the true 3D objective."""
import os, sys, json, time, argparse
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import core, search
from search import cma_run, fbflat, spiral_polygon, arm_angles, heights, HEIGHT_PATTERNS, log, save_sol, polish

EPS_LIST = (3e-4, 1e-3, 3e-3, 1e-2, 3e-2)


def flat_eval(z, V):
    P = np.ascontiguousarray(z[:2 * V].reshape(V, 2)); Y = np.ascontiguousarray(z[2 * V:3 * V])
    return core.evalF_flat(P, Y, core.JOG)[0]


def to3d(z, V, eps):
    P = z[:2 * V].reshape(V, 2); Y = z[2 * V:3 * V]
    Y = Y / max(1e-12, np.abs(Y).max())
    s = np.ptp(P[:, 0]) + np.ptp(P[:, 1])
    return np.stack([P[:, 0], eps * s * Y, P[:, 1]], axis=1)


def from3d(X):
    """principal axes: thin direction -> heights"""
    X = core.normalize(X)
    U, S, Vt = np.linalg.svd(X - X.mean(0))
    Y = (X - X.mean(0)) @ Vt.T
    h = Y[:, 2] / max(1e-12, np.abs(Y[:, 2]).max())
    return np.r_[Y[:, :2].ravel(), h]


def flat_opt(z0, V, seed, sigmas=(0.05, 0.01, 0.002), maxevals=150000, n_out=0):
    f = fbflat(V, n_out)
    best = (f(z0[None, :])[0], z0); ne = 0
    for k, s in enumerate(sigmas):
        start = best[1] if k > 0 else z0
        fv, zb, n = cma_run(f, start, s, maxevals, seed=seed * 7 + k, stall=300)
        ne += n
        if fv > best[0]:
            best = (fv, zb)
    return best[0], best[1], ne


def realise(z, V, do_polish=True):
    best = None
    for eps in EPS_LIST:
        X = to3d(z, V, eps)
        Fx = core.F(X)
        if best is None or Fx > best[0]:
            best = (Fx, X, eps)
    Fx, X, eps = best
    if do_polish and Fx > -0.05:
        Xp, Fp = polish(X, rounds=3, box=0.01)
        if Fp > Fx:
            X, Fx = Xp, Fp
    return X, Fx, eps


def seed_flat(V, rng, a=None, pattern=None, apex=False):
    k = V - 1 if apex else V
    if a is None:
        a = int(rng.integers((k - 1 + 1) // 2, k - 1 - 1 + 1))   # a >= b >= 1
    b = k - 1 - a
    dR, dL = arm_angles(a, b, rng)
    P = spiral_polygon(dR, dL)
    if pattern is None:
        pattern = HEIGHT_PATTERNS[rng.integers(len(HEIGHT_PATTERNS))]
    y = heights(pattern, a, b, rng)
    if rng.random() < 0.5:
        P = P * np.array([[-1.0, 1.0]])    # mirror (arms swap sides; harmless)
    if apex:
        P = np.vstack([P, [0.0, P[:, 1].min() * rng.uniform(0.6, 0.95)]]); y = np.r_[y, rng.uniform(0.5, 2)]
    return np.r_[P.ravel(), y], dict(a=a, b=b, pattern=pattern, apex=apex)


def task(args):
    t0 = time.time()
    V = args["V"]; seed = args["seed"]
    rng = np.random.default_rng(seed)
    if args["kind"] == "flatseed":
        z0, meta = seed_flat(V, rng, a=args.get("a"), pattern=args.get("pattern"), apex=args.get("apex", False))
    elif args["kind"] == "flatfile":
        X = np.load(args["file"])
        if args.get("drop") is not None:
            X = np.delete(X, args["drop"], axis=0)
        z0 = from3d(X)
        z0 = z0 + args.get("jitter", 0.0) * rng.normal(size=z0.shape)
        meta = dict(file=os.path.basename(args["file"]), drop=args.get("drop"))
    n_out = (V - 1 if args.get("apex") else V) if args.get("outline", True) else 0
    Ff, zb, ne = flat_opt(z0, V, seed, sigmas=args.get("sigmas", (0.05, 0.01, 0.002)), maxevals=args.get("maxevals", 150000), n_out=n_out)
    X, F3, eps = realise(zb, V)
    inf = core.info(X)
    tag = f"V{V}_{args['kind']}_{seed}"
    path = save_sol(X, tag)
    np.save(os.path.join(search.SOLDIR, tag + "_flat.npy"), zb)
    rec = dict(kind=args["kind"], family="flat", V=V, seed=seed, F_flat=float(Ff), F=inf["F"], nvert=inf["nvert"],
               dmin_rel=inf["dmin_rel"], eps=eps, evals=int(ne), time=round(time.time() - t0, 1), file=path, meta=meta)
    log(rec)
    return rec


def run_pool(tasks, nproc=6):
    from concurrent.futures import ProcessPoolExecutor, as_completed
    out = []
    with ProcessPoolExecutor(nproc) as ex:
        futs = [ex.submit(task, t) for t in tasks]
        for fu in as_completed(futs):
            try:
                r = fu.result()
            except Exception as e:
                import traceback; traceback.print_exc(); continue
            out.append(r)
            print(f"V={r['V']} {r['kind']} seed={r['seed']}: F_flat={r['F_flat']:+.4e} F3d={r['F']:+.4e} nvert={r['nvert']} "
                  f"eps={r['eps']} evals={r['evals']} t={r['time']}s meta={r['meta']}", flush=True)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode")      # seeds | file | ladder
    ap.add_argument("--V", type=int, default=13)
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--file")
    ap.add_argument("--a", type=int, default=None)
    ap.add_argument("--pattern", default=None)
    ap.add_argument("--apex", action="store_true")
    ap.add_argument("--no-outline", action="store_true")
    ap.add_argument("--maxevals", type=int, default=150000)
    ap.add_argument("--sigmas", default="0.05,0.01,0.002")
    ap.add_argument("--jitter", type=float, default=0.0)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--nproc", type=int, default=6)
    a = ap.parse_args()
    sig = tuple(float(s) for s in a.sigmas.split(","))
    t0 = time.time()
    if a.mode == "seeds":
        tasks = [dict(kind="flatseed", V=a.V, seed=a.seed0 + s, a=a.a, pattern=a.pattern, apex=a.apex, outline=not a.no_outline,
                      sigmas=sig, maxevals=a.maxevals) for s in range(a.n)]
    elif a.mode == "file":
        V = len(np.load(a.file))
        tasks = [dict(kind="flatfile", V=V, seed=a.seed0 + s, file=a.file, jitter=a.jitter, sigmas=sig,
                      maxevals=a.maxevals) for s in range(a.n)]
    elif a.mode == "ladder":
        X = np.load(a.file)
        tasks = [dict(kind="flatfile", V=len(X) - 1, seed=a.seed0 + d, file=a.file, drop=d, jitter=a.jitter,
                      sigmas=sig, maxevals=a.maxevals) for d in range(len(X)) for r in range(max(1, a.n // len(X)))]
        for i, t in enumerate(tasks):
            t["seed"] = a.seed0 + i
    res = run_pool(tasks, a.nproc)
    res.sort(key=lambda r: -r["F_flat"])
    print(f"done {len(res)} in {time.time()-t0:.0f}s; best F_flat={res[0]['F_flat']:+.4e} F3d={res[0]['F']:+.4e} {res[0]['file']}")
