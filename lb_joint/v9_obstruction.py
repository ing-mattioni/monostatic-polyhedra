"""
v9_obstruction.py -- where the Gram-matrix certificates stop: a random sample of the 751,513 nine-vertex
(sphere, root, tree) cases that the single-cone certificate does not settle (stage 1, run_V9_stage1.log).

Each sampled tree is attacked with (1) hybrid two-cone certificates and (2) coupled two-cone certificates
(stage2_hybrid.work, every certificate checked exactly); trees that resist both are tested with the
strongest relaxation available here, all nine cone decompositions at once plus the steepest-neighbour
inequalities (general_cert.solve_general with Q = all vertices, S1). A failed search leaves the tree
uncertified. It does not prove feasibility of the relaxation or non-existence of these certificates;
no exact feasible Gram matrix is constructed or checked for such failures here.

    python v9_obstruction.py N_SAMPLE WORKERS SEED    ->  v9_obstruction.json, printed summary
"""
import collections
import json
import os
import random
import sys
import time

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "lowerbound"))
V = 9


def classify(task):
    import general_cert as gc
    import stage2_hybrid as s2
    gid, v, parent = task
    t0 = time.time()
    _, _, _, line, how, attempts = s2.work(task)
    out = dict(graph=gid, root=v, parent=list(parent), how=how, attempts=attempts,
               depths=gc.depths(list(parent), v))
    if how == "none":
        gd = s2._graphs()[gid]
        nb = gc.neighbours(V, gd.edges)
        r = gc.solve_general(V, v, list(parent), nb, gd.faces, list(range(V)), True)
        out["z_allcones_S1"] = r["z"]
        if r["z"] > 1e-7 and gc.round_general(V, v, list(parent), list(range(V)), r) is not None:
            out["how"] = "allcones+S1"
    out["time"] = time.time() - t0
    return out


if __name__ == "__main__":
    import multiprocessing as mp
    n, workers, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    t0 = time.time()
    unc = json.load(open(os.path.join(HERE, "uncertified_V9.json")))
    rng = random.Random(seed)
    sample = rng.sample(unc, n)
    total_unc = len(unc)
    del unc
    tasks = [(c["graph"], c["root"], tuple(c["parent"])) for c in sample]
    res = []
    with mp.get_context("spawn").Pool(workers) as pool:
        for k, r in enumerate(pool.imap_unordered(classify, tasks, chunksize=4)):
            res.append(r)
            if (k + 1) % 250 == 0:
                print(f"  {k + 1}/{n}  {dict(collections.Counter(x['how'] for x in res))}  {time.time() - t0:.0f}s",
                      flush=True)
    cnt = collections.Counter(r["how"] for r in res)
    surv = [r for r in res if r["how"] == "none"]
    frac = len(surv) / n
    print(f"sample of {n} of the {total_unc} stage-1-uncertified trees (seed {seed}): {dict(cnt)}")
    print(f"estimated number of trees left uncertified by the tested searches: {frac * total_unc:.0f} "
          f"({100 * frac:.2f}% of stage 2, {100 * frac * total_unc / 9915284:.3f}% of all 9,915,284)")
    md = collections.Counter(max(r["depths"]) for r in surv)
    md_all = collections.Counter(max(r["depths"]) for r in res)
    print("max tree depth, survivors:", dict(sorted(md.items())), "  all sampled:", dict(sorted(md_all.items())))
    json.dump(dict(n=n, seed=seed, total_stage2=total_unc, counts=dict(cnt), results=res),
              open(os.path.join(HERE, "v9_obstruction.json"), "w"))
    print(f"done in {time.time() - t0:.0f}s")
