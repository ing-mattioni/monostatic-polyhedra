"""CMA-ES search for homogeneous mono-unstable polyhedra with V vertices."""
import sys, os, time, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mu3d import f, info

def run_one(args):
    V, seed, x0, sigma0, maxevals = args
    import cma, warnings
    warnings.filterwarnings("ignore")
    rng = np.random.default_rng(seed)
    if x0 is None:
        # random start: points on a sphere-ish blob
        X = rng.normal(size=(V, 3))
        X /= np.linalg.norm(X, axis=1)[:, None]
        X *= rng.uniform(0.5, 1.5, size=(V, 1))
        x0 = X.ravel()
    def obj(x):
        X = x.reshape(V, 3)
        val = f(X)
        if not np.isfinite(val): return 1e3
        return -val
    es = cma.CMAEvolutionStrategy(x0, sigma0, {"seed": int(seed) + 1, "verbose": -9, "maxfevals": maxevals,
                                               "tolfun": 1e-13, "tolx": 1e-12, "tolfunhist": 1e-13})
    while not es.stop():
        xs = es.ask()
        es.tell(xs, [obj(x) for x in xs])
    xb = es.result.xbest
    return -es.result.fbest, xb, seed, es.result.evaluations

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("V", type=int)
    ap.add_argument("--seedfile", default=None)
    ap.add_argument("--drop", type=int, default=None, help="drop this vertex index from seed")
    ap.add_argument("--restarts", type=int, default=16)
    ap.add_argument("--sigma", type=float, default=0.3)
    ap.add_argument("--maxevals", type=int, default=200000)
    ap.add_argument("--out", default=None)
    ap.add_argument("--jitter", type=float, default=0.0)
    a = ap.parse_args()
    from concurrent.futures import ProcessPoolExecutor
    tasks = []
    base = None
    if a.seedfile:
        base = np.load(a.seedfile)
        if a.drop is not None:
            base = np.delete(base, a.drop, axis=0)
        assert base.shape[0] == a.V, base.shape
    rng = np.random.default_rng(12345)
    for s in range(a.restarts):
        x0 = None
        if base is not None:
            x0 = (base + a.jitter * rng.normal(size=base.shape)).ravel()
        tasks.append((a.V, s, x0, a.sigma, a.maxevals))
    t0 = time.time()
    best = None
    with ProcessPoolExecutor(max_workers=min(16, a.restarts)) as ex:
        for val, xb, seed, nev in ex.map(run_one, tasks):
            print(f"  restart {seed}: f={val:+.4e} evals={nev}", flush=True)
            if best is None or val > best[0]:
                best = (val, xb)
    X = best[1].reshape(a.V, 3)
    print(f"V={a.V} best f={best[0]:+.6e}  info={info(X)}  time={time.time()-t0:.0f}s")
    if a.out:
        np.save(a.out, X)
