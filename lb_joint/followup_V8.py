"""Follow-up pipeline for V = 8 after certify_joint.py has finished:
  1. S1 alone and S2 with Q = {v, q} (q over all vertices) on all uncertified trees
  2. S2 (all cones) and S1+S2 (all cones) on the survivors of step 1
  3. S3 / S3+S1 (omega-subdivision, <= 64 nodes) on a random sample of 100
     uncertified trees (diagnostic: is the single-cone relaxation the problem?)
  4. relaxed Gram matrices for trees surviving step 2 (all cones + S1)
Writes exp*_results.json, exp*_certG.jsonl.gz and survivors_V8.json in lb_joint/."""
import json
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import general_cert as gc  # noqa: E402
import enumerate_cases as ec  # noqa: E402


def certified(r):
    return any(r.get(k) for k in r if k.endswith("_exact"))


def main():
    V = 8
    unc_path = os.path.join(HERE, "uncertified_V8.json")
    t0 = time.time()
    print("=== step 1: S1 alone, S2 with Q={v,q} ===", flush=True)
    r1 = gc.run_experiments(V, unc_path, ["S1", "S2pair"], 4, os.path.join(HERE, "exp1_V8"))
    surv = [dict(graph=r["graph"], root=r["root"], parent=r["parent"]) for r in r1 if not certified(r)]
    sp = os.path.join(HERE, "survivors1_V8.json")
    json.dump(surv, open(sp, "w"))
    print(f"step 1 survivors: {len(surv)}  ({time.time() - t0:.0f}s)", flush=True)
    surv2 = []
    if surv:
        print("=== step 2: S2 all cones; S1+S2 all cones ===", flush=True)
        r2 = gc.run_experiments(V, sp, ["S2all", "S1S2all"], 4, os.path.join(HERE, "exp2_V8"))
        surv2 = [dict(graph=r["graph"], root=r["root"], parent=r["parent"]) for r in r2 if not certified(r)]
    json.dump(surv2, open(os.path.join(HERE, "survivors_V8.json"), "w"))
    print(f"step 2 survivors: {len(surv2)}  ({time.time() - t0:.0f}s)", flush=True)
    print("=== step 3: S3 / S3+S1 diagnostic on 100 random uncertified trees ===", flush=True)
    unc = json.load(open(unc_path))
    random.seed(1)
    sample = random.sample(unc, min(100, len(unc)))
    spath = os.path.join(HERE, "sample100_V8.json")
    json.dump(sample, open(spath, "w"))
    gc.run_experiments(V, spath, ["S3", "S3S1"], 4, os.path.join(HERE, "exp3_V8"))
    if surv2:
        print("=== step 4: relaxed Gram matrices (all cones + S1) for survivors ===", flush=True)
        gds = ec.graph_data(V, with_trees=False)
        out = []
        for c in surv2[:20]:
            gd = gds[c["graph"]]
            nb = gc.neighbours(V, gd.edges)
            res, R, lam = gc.gram_feasible(V, c["root"], c["parent"], nb, gd.faces,
                                           list(range(V)), True, iters=60, restarts=2)
            hc = gc.hull_check(R, gd.edges) if R is not None else {}
            out.append(dict(c, residual=float(res), R=R.tolist() if R is not None else None,
                            lam={str(q): l.tolist() for q, l in lam.items()} if lam else None, **hc))
            print(f"  graph {c['graph']} root {c['root']} parent {c['parent']}: residual {res:.2e}, "
                  f"rank {hc.get('rank')}, eig {[round(t, 5) for t in hc.get('eigenvalues', [])]}", flush=True)
        json.dump(out, open(os.path.join(HERE, "survivors_gram_V8.json"), "w"))
    print(f"done ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
