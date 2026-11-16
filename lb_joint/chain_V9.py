"""Driver for the V = 9 lower-bound computation (all 9,915,284 (sphere, root, tree) cases):
  1. single-cone certificates (certify_joint.run)            -> certs_V9.jsonl.gz, uncertified_V9.json
  2. two-cone certificates for the rest (S2 with Q = {v,q})  -> exp1_V9_certG.jsonl.gz
  3. all cones / all cones + S1 for survivors of step 2       -> exp2_V9_certG.jsonl.gz
  4. exact verification                                       -> verify_V9.log
Run:  python chain_V9.py [workers]      (several hours; logs: run_V9.log is this script's stdout)"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

if __name__ == "__main__":
    import certify_joint as cj
    import general_cert as gc
    V = 9
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    t0 = time.time()
    unc = cj.run(V, workers=workers, chunk=200)
    print(f"=== step 1 done: {len(unc)} uncertified ({time.time() - t0:.0f}s)", flush=True)
    files = ["certs_V9.jsonl.gz"]
    surv = []
    if unc:
        r1 = gc.run_experiments(V, os.path.join(HERE, "uncertified_V9.json"), ["S2pair"], workers,
                                os.path.join(HERE, "exp1_V9"))
        files.append("exp1_V9_certG.jsonl.gz")
        surv = [dict(graph=r["graph"], root=r["root"], parent=r["parent"]) for r in r1
                if not any(r.get(k) for k in r if k.endswith("_exact"))]
        json.dump(surv, open(os.path.join(HERE, "survivors1_V9.json"), "w"))
        print(f"=== step 2 done: {len(surv)} survivors ({time.time() - t0:.0f}s)", flush=True)
    surv2 = []
    if surv:
        r2 = gc.run_experiments(V, os.path.join(HERE, "survivors1_V9.json"), ["S2all", "S1S2all"], workers,
                                os.path.join(HERE, "exp2_V9"))
        files.append("exp2_V9_certG.jsonl.gz")
        surv2 = [dict(graph=r["graph"], root=r["root"], parent=r["parent"]) for r in r2
                 if not any(r.get(k) for k in r if k.endswith("_exact"))]
        print(f"=== step 3 done: {len(surv2)} survivors ({time.time() - t0:.0f}s)", flush=True)
    json.dump(surv2, open(os.path.join(HERE, "survivors_V9.json"), "w"))
    if surv2:
        print("SOME TREES REMAIN UNCERTIFIED: no theorem for V = 9.", flush=True)
    else:
        with open(os.path.join(HERE, "verify_V9.log"), "w") as fh:
            rc = subprocess.call([sys.executable, os.path.join(HERE, "verify_joint.py"), str(V)] +
                                 [os.path.join(HERE, f) for f in files], stdout=fh, stderr=subprocess.STDOUT)
            fh.write(f"exit={rc}\n")
        print(f"=== verification exit code {rc} ({time.time() - t0:.0f}s)", flush=True)
    print("CHAIN_DONE", flush=True)
