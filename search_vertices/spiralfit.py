"""Reduced search on the tight-double-spiral manifold: params = arm angles (a+b) + heights (k).
The polygon is a tight Conway double spiral about O; c' = flat-limit centroid from the heights.
Objective = penalised flat F (planar angular margins w.r.t. c').  Then full flat CMA + 3D realisation."""
import os, sys, time, argparse, json
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import core, search, flat
from search import cma_run, spiral_polygon, arm_angles, heights, HEIGHT_PATTERNS, log, save_sol


def build_z(p, a, b):
    dR = np.abs(p[:a]); dL = np.abs(p[a:a + b]); y = p[a + b:]
    P = spiral_polygon(dR, dL)
    return np.r_[P.ravel(), y]


def make_f(a, b):
    k = a + b + 1
    fb = search.fbflat(k, k)

    def f(Ps):
        Z = np.empty((len(Ps), 3 * k))
        bad = np.zeros(len(Ps), dtype=bool)
        for t, p in enumerate(Ps):
            dR = np.abs(p[:a]); dL = np.abs(p[a:a + b])
            if dR.sum() > np.pi or dL.sum() > np.pi or dR.max() > 1.5 or dL.max() > 1.5:
                bad[t] = True
            Z[t] = build_z(p, a, b)
        out = fb(Z)
        out[bad] = -5.0
        return out
    return f


def task(args):
    t0 = time.time()
    a, b, seed = args["a"], args["b"], args["seed"]
    rng = np.random.default_rng(seed)
    k = a + b + 1
    dR, dL = arm_angles(a, b, rng)
    pattern = args.get("pattern") or ["growalt", "growalt", "alt", "zero_top_alt", "random", "grow"][seed % 6]
    y = heights(pattern, a, b, rng)
    p0 = np.r_[dR, dL, y]
    f = make_f(a, b)
    best = (f(p0[None, :])[0], p0); ne = 0
    for s in (0.1, 0.03, 0.01):
        fv, pv, n = cma_run(f, best[1], s, 60000, seed=seed * 3 + int(s * 100), stall=200)
        ne += n
        if fv > best[0]: best = (fv, pv)
    Fred = best[0]
    z = build_z(best[1], a, b)
    Ff, zb, n2 = flat.flat_opt(z, k, seed, sigmas=(0.01, 0.003, 0.001), maxevals=100000, n_out=k)
    ne += n2
    X, F3, eps = flat.realise(zb, k)
    inf = core.info(X)
    tag = f"V{k}_spiral_{a}_{b}_{seed}"
    path = save_sol(X, tag); np.save(os.path.join(search.SOLDIR, tag + "_flat.npy"), zb)
    rec = dict(kind="spiralfit", family="flat", V=k, seed=seed, F_red=float(Fred), F_flat=float(Ff), F=inf["F"],
               nvert=inf["nvert"], eps=eps, evals=int(ne), time=round(time.time() - t0, 1), file=path,
               meta=dict(a=a, b=b, pattern=pattern))
    log(rec)
    return rec


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--V", type=int, default=13)
    ap.add_argument("--n", type=int, default=6, help="runs per split")
    ap.add_argument("--splits", default=None)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--nproc", type=int, default=6)
    a = ap.parse_args()
    k = a.V
    if a.splits:
        splits = [tuple(int(v) for v in s.split("+")) for s in a.splits.split(",")]
    else:
        splits = [(aa, k - 1 - aa) for aa in range((k - 1 + 1) // 2, k - 1) if k - 1 - aa >= 2]
    tasks = [dict(a=sa, b=sb, seed=a.seed0 + 100 * i + r) for i, (sa, sb) in enumerate(splits) for r in range(a.n)]
    from concurrent.futures import ProcessPoolExecutor, as_completed
    t0 = time.time(); res = []
    with ProcessPoolExecutor(a.nproc) as ex:
        for fu in as_completed([ex.submit(task, t) for t in tasks]):
            try:
                r = fu.result()
            except Exception:
                import traceback; traceback.print_exc(); continue
            res.append(r)
            print(f"V={r['V']} {r['meta']['a']}+{r['meta']['b']} {r['meta']['pattern']} seed={r['seed']}: F_red={r['F_red']:+.4e} "
                  f"F_flat={r['F_flat']:+.4e} F3d={r['F']:+.4e} nvert={r['nvert']} t={r['time']}s", flush=True)
    res.sort(key=lambda r: -r["F"])
    print(f"done {len(res)} in {time.time()-t0:.0f}s; best F3d={res[0]['F']:+.4e} {res[0]['file']}")
