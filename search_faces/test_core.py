"""Validation of core.py against the scipy reference and on Reshetov's examples."""
import os, sys, time
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from core import evaluate, evaluate_p, reference, reshetov_p, batch_obj, RESHETOV


def random_poly(rng, F):
    U = rng.normal(size=(F, 3))
    U /= np.linalg.norm(U, axis=1)[:, None]
    d = rng.uniform(0.5, 1.5, F)
    M = np.diag(rng.uniform(0.3, 3, 3))
    A = U @ np.linalg.inv(M)   # anisotropic
    return A, d


def main():
    rng = np.random.default_rng(1)
    nbad = 0
    ntest = 0
    maxerr = 0.0
    for t in range(300):
        F = int(rng.integers(5, 19))
        A, b = random_poly(rng, F)
        r = evaluate(A, b)
        ref = reference(A, b)
        if ref is None or r["status"] != 0:
            if r["status"] == 0:
                nbad += 1
            continue
        ntest += 1
        err = max(abs(r["vol"] - ref["vol"]) / ref["vol"], np.abs(r["c"] - ref["c"]).max() / ref["vol"] ** (1 / 3))
        pres = r["present"]
        refp = ref["s"] > -1e8
        if not np.array_equal(pres, refp):
            print("face presence mismatch", t, pres.astype(int), refp.astype(int))
            nbad += 1
            continue
        err = max(err, np.abs(r["s"][pres] - ref["s"][pres]).max())
        maxerr = max(maxerr, err)
        if err > 1e-8:
            nbad += 1
            print("mismatch", t, err)
    # unbounded test
    A = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, -1, 0.0]])
    print("unbounded status (expect 2 or 1):", evaluate(A, np.ones(4))["status"])
    print(f"random tests: {ntest} compared, {nbad} bad, max err {maxerr:.2e}")

    for F in (14, 15, 16, 17):
        P = reshetov_p(F)
        r = evaluate_p(P)
        ref = reference(P, (P ** 2).sum(1))
        print(f"Reshetov F={F}: status={r['status']} nfaces={r['nfaces']} nvert={r['nvert']} nstable={r['nstable']} "
              f"G(euclid)={r['G']:.6f} M_res={r['M_res']:.6f} M_euc={r['M_euc']:.6f} stable face={r['istar']} "
              f"s*={r['s'][r['istar']]:.6f}  ref nstable={(ref['s'] > 0).sum()} max|ds|={np.abs(ref['s'] - r['s']).max():.1e}")
    # throughput
    for F in (13, 14, 17):
        P = reshetov_p(14)[:F] if F <= 14 else reshetov_p(F)
        X = np.tile(P.ravel(), (2000, 1)) * (1 + 1e-3 * rng.normal(size=(2000, 3 * F)))
        batch_obj(X[:2], F, 0.1, 0.0, 0.0)
        t0 = time.perf_counter()
        batch_obj(X, F, 0.1, 0.0, 0.0)
        dt = time.perf_counter() - t0
        print(f"throughput F={F}: {2000 / dt:.0f} evals/s per core")


if __name__ == "__main__":
    main()
