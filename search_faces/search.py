"""Search drivers for mono-stable polyhedra with few faces.

Parametrisation: plane i = {x : p_i . x = |p_i|^2}, origin interior; x = P.ravel() (3F reals).
Objective (maximise): J = G + lam * sum_{i != i*} min(0, -s_i)   (G = -(second largest s_i)/vol^(1/3)).

modes
  improve  F                      : CMA-ES polish of Reshetov's F-face example
  ladder   --seedfile f.npy        : delete each face in turn, re-optimise (many restarts / sigmas)
  spindle  F                      : spindle seeds from the long-spindle limit theory
  file     --seedfile f.npy        : re-optimise a given solution
All runs append one JSON line per CMA run to runs.jsonl.
"""
import os, sys, json, time, argparse
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np

LOG = os.path.join(HERE, "runs.jsonl")
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)


def normalize(P):
    P = np.asarray(P, float).reshape(-1, 3)
    return P / np.median(np.linalg.norm(P, axis=1))


def cma_run(task):
    """task: dict(x0, F, sigma, seed, maxevals, lam, popsize, minarea, wF, tag, ...)"""
    import cma, warnings
    warnings.filterwarnings("ignore")
    from core import batch_obj, evaluate_p
    F = task["F"]
    x0 = np.asarray(task["x0"], float)
    lam = task.get("lam", 0.1)
    minarea = task.get("minarea", 0.0)
    wF = task.get("wF", 0.0)
    norm = task.get("norm", 1)
    opts = {"seed": int(task["seed"]) + 1, "verbose": -9, "maxfevals": task["maxevals"],
            "tolfun": 1e-14, "tolx": 1e-13, "tolfunhist": 1e-14, "tolstagnation": 2000}
    if task.get("popsize"):
        opts["popsize"] = task["popsize"]
    t0 = time.time()
    es = cma.CMAEvolutionStrategy(x0, task["sigma"], opts)
    bestJ, bestx = -1e18, x0
    nev = 0
    planar = task.get("obj") == "planar"
    if planar:
        from surrogate import batch_planar, planar_G
    while not es.stop():
        X = np.array(es.ask())
        if planar:
            J, NF, NP = batch_planar(X, F)
        else:
            J, G, NF = batch_obj(X, F, lam, minarea, wF, norm, task.get("amax", 0.0))
        es.tell(list(X), list(-J))
        nev += len(X)
        k = int(np.argmax(J))
        if J[k] > bestJ:
            bestJ, bestx = J[k], X[k].copy()
    r = evaluate_p(bestx.reshape(F, 3))
    out = dict(tag=task.get("tag", ""), F=F, seed=int(task["seed"]), sigma=task["sigma"], lam=lam,
               popsize=task.get("popsize"), J=float(bestJ), G=float(r["G"]), Gw=float(r["Gw"]), M_res=float(r["M_res"]), norm=norm,
               nfaces=int(r["nfaces"]), nstable=int(r["nstable"]), status=int(r["status"]),
               min_area=float(r["min_area"]), min_edge=float(r["min_edge"]), evals=nev,
               secs=round(time.time() - t0, 1), aspect=float(r["vol"] ** (2 / 3) / r["w"] ** 2))
    out.update({k: v for k, v in task.items() if k in ("drop", "split", "elong", "pattern", "parent", "obj")})
    if planar:
        g2, nf, npres = planar_G(bestx.reshape(F, 3).copy())
        out["G2"] = float(g2)
        out["npresent_cs"] = int(npres)
    return out, bestx


def log(rec):
    with open(LOG, "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def run_pool(tasks, nworkers=6, save_prefix=None, keep_best=True):
    from concurrent.futures import ProcessPoolExecutor, as_completed
    best = None
    results = []
    with ProcessPoolExecutor(max_workers=nworkers) as ex:
        futs = {ex.submit(cma_run, t): t for t in tasks}
        for fu in as_completed(futs):
            rec, x = fu.result()
            log(rec)
            results.append((rec, x))
            print(json.dumps(rec), flush=True)
            if save_prefix and rec["status"] == 0:
                key = (rec.get("G2", rec["G"]) > 0, rec.get("G2", rec["Gw"]))
                if best is None or key > best[0]:
                    best = (key, rec, x)
                    np.save(f"{save_prefix}_best.npy", x.reshape(rec["F"], 3))
                if rec.get("G2", rec["G"]) > 0:
                    np.save(f"{save_prefix}_F{rec['nfaces']}_G{rec['G']:.6f}_s{rec['seed']}.npy",
                            x.reshape(rec["F"], 3))
    return results


# ------------------------------------------------------------------------------------------
def mode_improve(a):
    from core import reshetov_p
    P0 = normalize(reshetov_p(a.F)).ravel()
    tasks = []
    for s in range(a.restarts):
        sig = a.sigma * (1 if s % 2 == 0 else 3)
        tasks.append(dict(x0=P0, F=a.F, sigma=sig, seed=s, maxevals=a.maxevals, lam=a.lam, norm=a.norm,
                          popsize=a.popsize, tag=f"improve{a.F}"))
    run_pool(tasks, a.workers, save_prefix=os.path.join(RES, f"improve{a.F}"))


def mode_ladder(a):
    P = normalize(np.load(a.seedfile))
    F = P.shape[0]
    drops = range(F) if a.drop is None else [int(t) for t in a.drop.split(",")]
    rng = np.random.default_rng(a.rseed)
    tasks = []
    sigmas = [float(t) for t in a.sigmas.split(",")]
    for d in drops:
        Q = np.delete(P, d, axis=0)
        for s in range(a.restarts):
            sig = sigmas[s % len(sigmas)]
            x0 = (Q + a.jitter * rng.normal(size=Q.shape)).ravel()
            tasks.append(dict(x0=x0, F=F - 1, sigma=sig, seed=1000 * d + s, maxevals=a.maxevals, lam=a.lam, norm=a.norm, obj=a.obj, amax=a.amax,
                              popsize=a.popsize, tag=a.tag or f"ladder{F}to{F-1}", drop=d,
                              parent=os.path.basename(a.seedfile)))
    # interleave drops so that all faces get early attention
    tasks.sort(key=lambda t: (t["seed"] % 1000, t["drop"]))
    run_pool(tasks, a.workers, save_prefix=os.path.join(RES, (a.tag or f"ladder{F}to{F-1}")))


def mode_cap(a):
    """(F-2) faces of a good F-face body + one small cap face truncating an acute end -> F-1 faces."""
    from core import evaluate_p
    P = normalize(np.load(a.seedfile))
    F = P.shape[0]
    rng = np.random.default_rng(a.rseed)
    V = evaluate_p(P)["vertices"]
    pairs = [(i, j) for i in range(F) for j in range(i + 1, F)]
    rng.shuffle(pairs)
    tasks = []
    for (i, j) in pairs[:a.restarts]:
        Q = np.delete(P, [i, j], axis=0)
        end = 1 if rng.random() < 0.5 else -1
        f = rng.uniform(0.6, 0.92)
        xe = V[:, 0].max() if end > 0 else V[:, 0].min()
        cap = np.array([[f * xe, 0.0, 0.0]])
        x0 = np.vstack([Q, cap]).ravel()
        tasks.append(dict(x0=x0, F=F - 1, sigma=[0.003, 0.03][len(tasks) % 2], seed=len(tasks), maxevals=a.maxevals,
                          lam=a.lam, norm=0, amax=40.0, tag=a.tag or f"cap{F - 1}", drop=f"{i},{j}",
                          pattern=f"cap{'+' if end > 0 else '-'}{f:.2f}"))
    run_pool(tasks, a.workers, save_prefix=os.path.join(RES, a.tag or f"cap{F - 1}"))


def mode_file(a):
    P = normalize(np.load(a.seedfile))
    F = P.shape[0]
    rng = np.random.default_rng(a.rseed)
    tasks = []
    sigmas = [float(t) for t in a.sigmas.split(",")]
    for s in range(a.restarts):
        x0 = (P + a.jitter * rng.normal(size=P.shape)).ravel()
        tasks.append(dict(x0=x0, F=F, sigma=sigmas[s % len(sigmas)], seed=s, maxevals=a.maxevals, lam=a.lam, norm=a.norm,
                          popsize=a.popsize, tag=a.tag or "file", parent=os.path.basename(a.seedfile)))
    run_pool(tasks, a.workers, save_prefix=os.path.join(RES, a.tag or "file"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode")
    ap.add_argument("F", type=int, nargs="?", default=14)
    ap.add_argument("--seedfile")
    ap.add_argument("--drop")
    ap.add_argument("--restarts", type=int, default=12)
    ap.add_argument("--sigma", type=float, default=0.01)
    ap.add_argument("--sigmas", default="0.003,0.01,0.03,0.1")
    ap.add_argument("--maxevals", type=int, default=300000)
    ap.add_argument("--lam", type=float, default=0.1)
    ap.add_argument("--popsize", type=int, default=0)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--jitter", type=float, default=0.0)
    ap.add_argument("--rseed", type=int, default=0)
    ap.add_argument("--tag", default="")
    ap.add_argument("--norm", type=int, default=1)
    ap.add_argument("--obj", default="3d")
    ap.add_argument("--amax", type=float, default=40.0)
    a = ap.parse_args()
    a.popsize = a.popsize or None
    {"improve": mode_improve, "ladder": mode_ladder, "file": mode_file, "cap": mode_cap}[a.mode](a)
