"""
sanity.py -- independent Monte-Carlo sanity check of the pentahedron theorem (no shared code with the proof).

Random convex pentahedra are generated directly in Euclidean space:
  * frusta:   apex O, three random rays, cut by two planes (0 < alpha_i < beta_i on ray i);
  * prisms:   random triangle, random common direction, random positive lengths;
  * pyramids: frusta/prisms with one lateral edge of length zero;
each additionally distorted by a random linear map with widely spread singular values.
For every sample the centroid of the solid is computed (three tetrahedra), then for each of the five
faces the foot of the perpendicular from the centroid and its signed distance to the boundary of the face
(positive = strictly inside = stable face).  Reported: histogram of the number of stable faces
(the theorem says: never fewer than 2) and the smallest second-largest margin encountered.
"""
import sys
import time
import numpy as np


def stable_margins(A, B):
    """A, B: arrays (n, 3, 3) of the vertices A_i, B_i (lateral edges A_i B_i).  Returns (n, 5) margins,
    normalised by the cube root of the volume."""
    n = A.shape[0]
    A1, A2, A3 = A[:, 0], A[:, 1], A[:, 2]
    B1, B2, B3 = B[:, 0], B[:, 1], B[:, 2]
    tets = [(B1, B2, B3, A1), (A1, B2, B3, A3), (A1, B2, A3, A2)]
    vol = np.zeros(n)
    mom = np.zeros((n, 3))
    for (a, b, c, d) in tets:
        det = np.abs(np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a))
        vol += det
        mom += det[:, None] * (a + b + c + d) / 4
    cen = mom / vol[:, None]
    scale = (vol / 6) ** (1 / 3)
    faces = [(A2, A3, B3, B2), (A3, A1, B1, B3), (A1, A2, B2, B1), (A1, A2, A3), (B1, B2, B3)]
    out = np.zeros((n, 5))
    for f, poly in enumerate(faces):
        m = len(poly)
        # robust normal: sum of cross products (Newell)
        nrm = np.zeros((n, 3))
        for k in range(m):
            nrm += np.cross(poly[k], poly[(k + 1) % m])
        nn = np.linalg.norm(nrm, axis=1)
        nn[nn == 0] = 1e-300
        nrm = nrm / nn[:, None]
        foot = cen - np.einsum("ij,ij->i", cen - poly[0], nrm)[:, None] * nrm
        marg = np.full(n, np.inf)
        for k in range(m):
            e = poly[(k + 1) % m] - poly[k]
            le = np.linalg.norm(e, axis=1)
            good = le > 1e-14 * scale          # zero-length edge (pyramid): no constraint
            s = np.einsum("ij,ij->i", np.cross(e, foot - poly[k]), nrm) / np.where(good, le, 1.0)
            marg = np.where(good, np.minimum(marg, s), marg)
        out[:, f] = marg / scale
    return out


def random_maps(rng, n):
    Q1, _ = np.linalg.qr(rng.normal(size=(n, 3, 3)))
    Q2, _ = np.linalg.qr(rng.normal(size=(n, 3, 3)))
    sv = np.exp(rng.uniform(-3.5, 3.5, size=(n, 3)))
    return np.einsum("nij,nj,nkj->nik", Q1, sv, Q2)


def sample(rng, n, kind):
    if kind in ("frustum", "frustum_pyramid"):
        D = rng.normal(size=(n, 3, 3))
        beta = np.exp(rng.uniform(-1, 1, size=(n, 3)))
        mode = rng.integers(0, 3, size=(n, 1))
        t = np.where(mode == 0, rng.uniform(0.01, 0.99, size=(n, 3)),
                     np.where(mode == 1, 1 - 10 ** rng.uniform(-4, -0.3, size=(n, 3)), 10 ** rng.uniform(-4, -0.3, size=(n, 3))))
        if kind == "frustum_pyramid":
            t[:, 2] = 1.0
        B = beta[:, :, None] * D
        A = (t * beta)[:, :, None] * D
    else:
        Bt = rng.normal(size=(n, 3, 3))
        d = rng.normal(size=(n, 3))
        s = 10 ** rng.uniform(-3, 0.5, size=(n, 3))
        if kind == "prism_pyramid":
            s[:, 2] = 0.0
        B = Bt
        A = Bt + s[:, :, None] * d[:, None, :]
    T = random_maps(rng, n)
    A = np.einsum("nij,nkj->nki", T, A)
    B = np.einsum("nij,nkj->nki", T, B)
    return A, B


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
    rng = np.random.default_rng(2026)
    t0 = time.time()
    total_hist = np.zeros(6, dtype=np.int64)
    worst = (np.inf, None)
    for kind in ("frustum", "prism", "frustum_pyramid", "prism_pyramid"):
        hist = np.zeros(6, dtype=np.int64)
        second_min = np.inf
        done = 0
        while done < N // 4:
            n = min(200_000, N // 4 - done)
            A, B = sample(rng, n, kind)
            M = stable_margins(A, B)
            if kind.endswith("pyramid"):
                pass
            cnt = (M > 0).sum(axis=1)
            hist += np.bincount(cnt, minlength=6)
            srt = np.sort(M, axis=1)
            sec = srt[:, -2]
            k = int(np.argmin(sec))
            if sec[k] < second_min:
                second_min = sec[k]
            done += n
        total_hist += hist
        print(f"{kind:16s}: {done} samples, #stable faces histogram (0..5) = {hist.tolist()}, "
              f"smallest 2nd-best margin = {second_min:.3e}", flush=True)
        worst = min(worst, (second_min, kind))
    print(f"TOTAL histogram of number of stable faces: {total_hist.tolist()}  "
          f"(samples with fewer than 2 stable faces: {int(total_hist[:2].sum())});  time {time.time() - t0:.0f}s")
