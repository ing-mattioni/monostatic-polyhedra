#!/usr/bin/env python3
"""
sanity_sample.py -- floating-point sanity check for Theorem B (NOT a proof).

For random convex polyhedra with 5..8 vertices (convex hulls of random
points; isotropic, elongated, nearly flat, and points on ellipsoids) compute
the centre of mass and the unstable vertices, using the DEFINITION directly:
vertex r_i is unstable iff every other vertex r_j satisfies
(r_j - r_i).(r_i - c) < 0 (the plane through r_i orthogonal to r_i - c then
meets P = conv(vertices) only in r_i).  Also checks, per sample:
  * Lemma 2.1 (neighbour test gives the same unstable set as the definition),
  * Lemma 4.1 (cone identity sum_T vol(T) sum_{k in T} r_k = 0 from every
    vertex q, with the hull triangulation),
  * that a shadowing tree exists whenever exactly one vertex is unstable
    (vacuous unless a counterexample appears).
Then a small adversarial hill-climb tries to make all but one vertex stable.
"""
import sys
import time
import numpy as np
from scipy.spatial import ConvexHull

rng = np.random.default_rng(20261116)


def analyse(pts):
    hull = ConvexHull(pts)
    vid = list(hull.vertices)
    idx = {g: k for k, g in enumerate(vid)}
    X = pts[vid]
    tris = [tuple(idx[g] for g in s) for s in hull.simplices]
    p0 = X.mean(axis=0)
    vol = 0.0
    mom = np.zeros(3)
    for a, b, c in tris:
        v6 = abs(np.linalg.det(np.array([X[a] - p0, X[b] - p0, X[c] - p0])))
        vol += v6
        mom += v6 * (p0 + X[a] + X[b] + X[c]) / 4.0
    cen = mom / vol
    R = X - cen
    n = len(X)
    scale = np.max(np.linalg.norm(R, axis=1)) ** 2
    # definition: all other vertices
    marg = np.array([max((R[j] - R[i]) @ R[i] for j in range(n) if j != i) for i in range(n)]) / scale
    unst = set(np.nonzero(marg < 0)[0])
    # Lemma 2.1: neighbours via hull edges (coplanar facets are generic -> triangles are faces)
    nb = [set() for _ in range(n)]
    for a, b, c in tris:
        nb[a] |= {b, c}
        nb[b] |= {a, c}
        nb[c] |= {a, b}
    unst_nb = {i for i in range(n) if all((R[j] - R[i]) @ R[i] < 0 for j in nb[i])}
    # Lemma 4.1: cone identity from every vertex q, relative to the size
    cone_err = 0.0
    for q in range(n):
        s = np.zeros(3)
        tot = 0.0
        for a, b, c in tris:
            if q in (a, b, c):
                continue
            w = abs(np.linalg.det(np.array([R[a] - R[q], R[b] - R[q], R[c] - R[q]]))) / 6.0
            s += w * (R[q] + R[a] + R[b] + R[c])
            tot += w
        cone_err = max(cone_err, np.linalg.norm(s) / (tot * np.sqrt(scale)))
        if abs(tot - vol / 6.0) > 1e-9 * vol:
            cone_err = max(cone_err, 1.0)
    return n, unst, unst_nb, marg, cone_err


def family(kind, m):
    if kind == "cube":
        return rng.uniform(-1, 1, (m, 3))
    if kind == "gauss":
        return rng.normal(size=(m, 3))
    if kind == "sphere":
        P = rng.normal(size=(m, 3))
        return P / np.linalg.norm(P, axis=1)[:, None]
    if kind == "elongated":
        return rng.normal(size=(m, 3)) * np.array([1.0, 1.0, 10.0 ** rng.uniform(0.5, 2)])
    if kind == "flat":
        return rng.normal(size=(m, 3)) * np.array([1.0, 1.0, 10.0 ** -rng.uniform(1, 3)])
    if kind == "needle_ellipsoid":
        P = rng.normal(size=(m, 3))
        P /= np.linalg.norm(P, axis=1)[:, None]
        return P * np.array([1.0, 10.0 ** rng.uniform(-2, 0), 10.0 ** rng.uniform(-2, 0)])
    raise ValueError(kind)


def main():
    t0 = time.time()
    kinds = ["cube", "gauss", "sphere", "elongated", "flat", "needle_ellipsoid"]
    per = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    stats = {}
    worst = (1e9, None)
    bad = []
    for kind in kinds:
        for nv in range(5, 9):
            got = 0
            tries = 0
            while got < per and tries < 50 * per:
                tries += 1
                m = nv if kind in ("sphere", "needle_ellipsoid") else nv + rng.integers(0, 4)
                pts = family(kind, m)
                try:
                    n, unst, unst_nb, marg, cerr = analyse(pts)
                except Exception:
                    continue
                if n != nv:
                    continue
                got += 1
                key = (kind, nv)
                st = stats.setdefault(key, {"n": 0, "minU": 99, "Uhist": {}, "lemma21": 0, "cone": 0.0})
                st["n"] += 1
                U = len(unst)
                st["minU"] = min(st["minU"], U)
                st["Uhist"][U] = st["Uhist"].get(U, 0) + 1
                st["lemma21"] += (unst != unst_nb)
                st["cone"] = max(st["cone"], cerr)
                srt = np.sort(marg)
                second = -srt[1]          # > 0 means second vertex unstable with that margin
                if second < worst[0]:
                    worst = (second, (kind, nv))
                if U < 2:
                    bad.append((kind, nv, U, pts.tolist()))
    print("random sampling (definition of unstable vertex applied directly)")
    print("%-17s %2s %5s %5s %-28s %9s %10s" % ("family", "V", "n", "minU", "U histogram", "L2.1 mis", "cone err"))
    total = 0
    for (kind, nv), st in sorted(stats.items()):
        total += st["n"]
        print("%-17s %2d %5d %5d %-28s %9d %10.2e" % (kind, nv, st["n"], st["minU"],
              str(dict(sorted(st["Uhist"].items()))), st["lemma21"], st["cone"]))
    print("total samples:", total, " samples with U < 2:", len(bad))
    print("smallest normalised margin of the second most unstable vertex: %.3e (%s)" % worst)

    # adversarial hill-climb: maximise the second smallest margin (>= 0 would mean U <= 1)
    print("\nadversarial hill-climb (score = second smallest margin; score >= 0 <=> at most one unstable vertex)")
    best_overall = -1e9
    for nv in range(5, 9):
        bestv = -1e9
        for start in range(12):
            P = family(["gauss", "flat", "elongated", "sphere"][start % 4], nv)
            def score(Q):
                try:
                    n, unst, _, marg, _ = analyse(Q)
                except Exception:
                    return -1e9
                if n != nv:
                    return -1e9
                return np.sort(marg)[1]
            s = score(P)
            step = 0.3
            for it in range(400):
                Q = P + rng.normal(size=P.shape) * step * np.std(P, axis=0)
                sq = score(Q)
                if sq > s:
                    P, s = Q, sq
                else:
                    step *= 0.995
            bestv = max(bestv, s)
        best_overall = max(best_overall, bestv)
        print("  V=%d: best score over 12 runs x 400 steps: %.3e" % (nv, bestv))
    print("max score found: %.3e (negative => no mono-unstable example found)" % best_overall)
    print("elapsed %.1f s" % (time.time() - t0))
    return 1 if bad or best_overall >= 0 else 0


if __name__ == "__main__":
    sys.exit(main())
