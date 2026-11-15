"""Spindle seeds from the long-spindle limit theory + planar surrogate.

Cross-section: tight double Conway spiral (explore/flat_affine.build convention) with F edges,
arms a + b = F - 1.  Edge lines nu_i . p = d_i (O = origin).  3D planes (tau_i, nu_i).(X, p) <= d_i.
Tilt patterns: 'alt'  - the k faces farthest from O get alternating-sign tilts (angular order), rest 0
               'rand' - random signs on the k top faces, small random tilts elsewhere
               'altall' - alternating signs on all faces, magnitude growing with distance
usage: python spindle.py F [--n N] [--maxevals M] [--tag T]
"""
import os, sys, json, argparse
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np


def build(dR, dL):
    def arm(d, sgn):
        ph = np.cumsum(d)
        r = np.cumprod(np.cos(d))
        return np.stack([sgn * r * np.sin(ph), r * np.cos(ph)], axis=1)
    R = arm(dR, 1.0)
    L = arm(dL, -1.0)
    return np.vstack([[0.0, 1.0], L, R[::-1]])   # ccw


def random_spiral(rng, F, a):
    b = F - 1 - a
    tot = np.radians(rng.uniform(250, 315))
    fa = (a + 0.6) / (a + b + 1.2) * rng.uniform(0.85, 1.15)
    tR = min(np.radians(179), tot * fa)
    tL = min(np.radians(179), tot - tR)

    def steps(n, T):
        w = rng.uniform(0.8, 1.2, n)
        w[0] *= rng.uniform(1.3, 2.3)
        return w * T / w.sum()
    dR = steps(a, tR) if a > 0 else np.zeros(0)
    dL = steps(b, tL) if b > 0 else np.zeros(0)
    return build(dR, dL)


def edge_lines(Pg):
    """ccw polygon -> outward unit normals nu (k,2), offsets d (k,) with nu.p <= d."""
    E = np.roll(Pg, -1, axis=0) - Pg
    nu = np.stack([E[:, 1], -E[:, 0]], axis=1)
    nu /= np.linalg.norm(nu, axis=1)[:, None]
    d = (nu * Pg).sum(1)
    return nu, d


def make_planes(nu, d, tau, drop_bottom=0.0):
    """p-parametrisation: normal n = (tau, nu) ; plane n.x = d  ->  p = n d / |n|^2 ."""
    n = np.c_[tau, nu]
    return n * (d / (n ** 2).sum(1))[:, None]


def tilt_pattern(rng, d, nu, pattern, k, eps):
    F = len(d)
    ang = np.arctan2(nu[:, 1], nu[:, 0])
    order = np.argsort(-d)            # farthest first
    tau = np.zeros(F)
    top = order[:k]
    if pattern == "alt":
        # alternate signs in angular order among the top faces
        topang = top[np.argsort(ang[top])]
        sg = 1.0 if rng.random() < 0.5 else -1.0
        for t, i in enumerate(topang):
            tau[i] = sg * (-1) ** t * eps * rng.uniform(0.5, 1.3) * (d[i] / d.max()) ** 2
    elif pattern == "rand":
        sg = rng.choice([-1.0, 1.0], size=k)
        if np.all(sg > 0) or np.all(sg < 0):
            sg[0] = -sg[0]
        tau[top] = sg * eps * rng.uniform(0.5, 1.3, k)
        rest = order[k:]
        tau[rest] = eps * 0.1 * rng.normal(size=len(rest))
    elif pattern == "altall":
        oa = np.argsort(ang)
        for t, i in enumerate(oa):
            tau[i] = (-1) ** t * eps * ((d[i] - d.min()) / (d.max() - d.min())) ** 1.5 * rng.uniform(0.6, 1.2)
    return tau


def make_seed(rng, F, a, pattern, k, eps, lower=0.03):
    Pg = random_spiral(rng, F, a)
    Pg = Pg - np.array([0.0, -lower])     # move reference point O down: O' = (0,-lower) -> shift polygon up
    nu, d = edge_lines(Pg)
    tau = tilt_pattern(rng, d, nu, pattern, k, eps)
    return make_planes(nu, d, tau)


def gen_tasks(F, n, maxevals, tag, rseed=0, sigmas=(0.01, 0.03), lam=0.1):
    from core import evaluate_p
    rng = np.random.default_rng(rseed)
    tasks = []
    tries = 0
    while len(tasks) < n and tries < 50 * n:
        tries += 1
        a = int(rng.integers(max(1, (F - 1) // 2 - 1), F - 2))   # a >= b-ish (asymmetric allowed)
        pattern = rng.choice(["alt", "alt", "rand", "altall"])
        k = int(rng.integers(3, min(8, F - 3) + 1))
        eps = float(rng.choice([0.03, 0.06, 0.1, 0.15]))
        P = make_seed(rng, F, a, pattern, k, eps, lower=float(rng.uniform(0.0, 0.08)))
        r = evaluate_p(P)
        if r["status"] != 0 or r["nfaces"] != F:
            continue
        P = P / np.median(np.linalg.norm(P, axis=1))
        tasks.append(dict(x0=P.ravel(), F=F, sigma=float(sigmas[len(tasks) % len(sigmas)]), seed=rseed * 100000 + len(tasks), norm=0, amax=40.0,
                          maxevals=maxevals, lam=lam, tag=tag, split=f"{a}+{F - 1 - a}", pattern=str(pattern),
                          elong=eps))
    return tasks


if __name__ == "__main__":
    from search import run_pool, RES
    ap = argparse.ArgumentParser()
    ap.add_argument("F", type=int)
    ap.add_argument("--n", type=int, default=24)
    ap.add_argument("--maxevals", type=int, default=400000)
    ap.add_argument("--tag", default="")
    ap.add_argument("--rseed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    tag = a.tag or f"spindle{a.F}"
    tasks = gen_tasks(a.F, a.n, a.maxevals, tag, a.rseed)
    print(f"{len(tasks)} seeds generated", flush=True)
    run_pool(tasks, a.workers, save_prefix=os.path.join(RES, tag))
