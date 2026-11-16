"""
falsify.py -- INDEPENDENT attempt to find a homogeneous convex pentahedron with < 2 stable faces (audit).
Own generators (directly in Euclidean space), own centroid (cone decomposition from an interior point over
the triangulated faces), own stability margin.  No code shared with penta/.

Families: frustum (apex + 3 rays, cut at a_i < b_i), prism (triangle + common direction + lengths),
pyramid (convex quadrilateral base + apex; base drawn directly, incl. near-parallelograms / trapezoids /
near-triangles), plus extreme versions (random linear maps with singular values up to e^{+-4}).

margin(face) = signed distance (in the face plane) from the foot of the perpendicular from the centroid to
the face boundary, positive inside; normalised by vol^(1/3).  Stable face <=> margin > 0.

  python falsify.py random <N> <seed>
  python falsify.py opt <restarts> <seed> <workers>
"""
import sys
import time
import numpy as np


# ------------------------------------------------------------------ geometry (vectorised over n shapes)
def margins(V, faces):
    """V: (n, m, 3) vertices; faces: list of vertex-index tuples (any orientation, convex polygons).
    Returns (n, 5) normalised margins and (n,) volume."""
    n = V.shape[0]
    ip = V.mean(axis=1)                                  # interior point (vertex average)
    vol = np.zeros(n)
    mom = np.zeros((n, 3))
    for f in faces:
        for k in range(1, len(f) - 1):
            a, b, c = V[:, f[0]], V[:, f[k]], V[:, f[k + 1]]
            dt = np.abs(np.einsum("ij,ij->i", np.cross(b - a, c - a), ip - a))   # 6 * vol of tet (ip,a,b,c)
            vol += dt
            mom += dt[:, None] * (a + b + c + ip) / 4
    cen = mom / vol[:, None]
    scale = (vol / 6) ** (1 / 3)
    out = np.zeros((n, len(faces)))
    for fi, f in enumerate(faces):
        P = [V[:, k] for k in f]
        # normal from the polygon (Newell), oriented outward w.r.t. ip
        nr = np.zeros((n, 3))
        for k in range(len(P)):
            nr += np.cross(P[k], P[(k + 1) % len(P)])
        nn = np.linalg.norm(nr, axis=1)
        nr /= np.where(nn > 0, nn, 1)[:, None]
        sgn = np.sign(np.einsum("ij,ij->i", nr, P[0] - ip))
        nr *= sgn[:, None]
        foot = cen - np.einsum("ij,ij->i", cen - P[0], nr)[:, None] * nr
        # polygon centroid (inside) to orient edge normals
        pc = sum(P) / len(P)
        m = np.full(n, np.inf)
        for k in range(len(P)):
            a, b = P[k], P[(k + 1) % len(P)]
            e = b - a
            le = np.linalg.norm(e, axis=1)
            ok = le > 1e-12 * scale
            en = np.cross(nr, e)                              # in-plane normal to the edge
            en /= np.where(np.linalg.norm(en, axis=1) > 0, np.linalg.norm(en, axis=1), 1)[:, None]
            s = np.sign(np.einsum("ij,ij->i", en, pc - a))   # make it point inward
            en *= s[:, None]
            dist = np.einsum("ij,ij->i", en, foot - a)
            m = np.where(ok, np.minimum(m, dist), m)
        out[:, fi] = m / scale
    return out, vol


PRISM_FACES = [(1, 2, 5, 4), (0, 2, 5, 3), (0, 1, 4, 3), (0, 1, 2), (3, 4, 5)]   # verts: B1 B2 B3 A1 A2 A3
PYR_FACES = [(0, 1, 2, 3), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)]          # verts: base q0..q3, apex


def rand_linear(rng, n, spread):
    Q1, _ = np.linalg.qr(rng.normal(size=(n, 3, 3)))
    Q2, _ = np.linalg.qr(rng.normal(size=(n, 3, 3)))
    sv = np.exp(rng.uniform(-spread, spread, size=(n, 3)))
    return np.einsum("nij,nj,nkj->nik", Q1, sv, Q2)


def gen(rng, n, kind):
    if kind == "frustum":
        O = rng.normal(size=(n, 1, 3))
        d = rng.normal(size=(n, 3, 3))
        b = np.exp(rng.uniform(-1.5, 1.5, size=(n, 3)))
        mode = rng.integers(0, 4, size=(n, 1))
        a = np.where(mode == 0, rng.uniform(0, 1, (n, 3)),
                     np.where(mode == 1, 1 - 10 ** rng.uniform(-5, 0, (n, 3)),
                              np.where(mode == 2, 10 ** rng.uniform(-5, 0, (n, 3)), rng.beta(0.3, 0.3, (n, 3)))))
        a = np.clip(a, 1e-9, 1 - 1e-9)
        B = O + b[:, :, None] * d
        A = O + (a * b)[:, :, None] * d
        V = np.concatenate([B, A], axis=1)
        faces = PRISM_FACES
    elif kind == "prism":
        B = rng.normal(size=(n, 3, 3))
        d = rng.normal(size=(n, 1, 3))
        s = 10 ** rng.uniform(-4, 1, size=(n, 3))
        A = B + s[:, :, None] * d
        V = np.concatenate([B, A], axis=1)
        faces = PRISM_FACES
    else:   # pyramid: convex quadrilateral base drawn directly
        ang = np.sort(rng.uniform(0, 2 * np.pi, size=(n, 4)), axis=1)
        mode = rng.integers(0, 3, size=(n, 1))
        rad = np.where(mode == 0, np.exp(rng.uniform(-3, 1, (n, 4))), 1.0)
        base = np.stack([rad * np.cos(ang), rad * np.sin(ang), np.zeros((n, 4))], axis=2)
        # mode 2: near-parallelogram / trapezoid
        par = rng.normal(size=(n, 2, 2))
        lam = np.where(rng.random((n, 1)) < 0.5, 1.0, rng.uniform(0, 1, (n, 1)))
        q0 = np.zeros((n, 2)); q1 = par[:, 0]; q3 = par[:, 1]; q2 = q3 + lam * par[:, 0]
        trap = np.stack([q0, q1, q2, q3], axis=1)
        trap = np.concatenate([trap, np.zeros((n, 4, 1))], axis=2) + 1e-3 * rng.normal(size=(n, 4, 3)) * np.array([1, 1, 0])
        base = np.where((mode == 2)[:, :, None], trap, base)
        apex = np.concatenate([rng.normal(size=(n, 2)) * 2, np.exp(rng.uniform(-3, 2, (n, 1)))], axis=1)[:, None, :]
        V = np.concatenate([base, apex], axis=1)
        faces = PYR_FACES
    T = rand_linear(rng, n, 4.0)
    V = np.einsum("nij,nkj->nki", T, V)
    return V, faces


def convex_ok(V, faces):
    """discard non-convex / degenerate samples: every vertex must be weakly inside every face plane,
    vertices of a face coplanar."""
    n = V.shape[0]
    ip = V.mean(axis=1)
    good = np.ones(n, bool)
    sc = np.abs(V - ip[:, None, :]).max(axis=(1, 2))
    for f in faces:
        P = [V[:, k] for k in f]
        nr = np.zeros((n, 3))
        for k in range(len(P)):
            nr += np.cross(P[k], P[(k + 1) % len(P)])
        nn = np.linalg.norm(nr, axis=1)
        good &= nn > 1e-14 * sc ** 2
        nr /= np.where(nn > 0, nn, 1)[:, None]
        sg = np.sign(np.einsum("ij,ij->i", nr, P[0] - ip))
        nr *= sg[:, None]
        h = np.einsum("nkj,nj->nk", V - P[0][:, None, :], nr)
        good &= (h <= 1e-9 * sc[:, None]).all(axis=1)
        inplane = np.abs(np.stack([np.einsum("ij,ij->i", Pk - P[0], nr) for Pk in P], 1))
        good &= (inplane <= 1e-9 * sc[:, None]).all(axis=1)
        # convex polygon in the face: all turns in the same direction (strictly)
        turns = []
        for k in range(len(P)):
            a, b, c = P[k], P[(k + 1) % len(P)], P[(k + 2) % len(P)]
            turns.append(np.einsum("ij,ij->i", np.cross(b - a, c - b), nr))
        turns = np.stack(turns, 1)
        good &= (turns > 0).all(1) | (turns < 0).all(1)
    return good


def run_random(N, seed):
    rng = np.random.default_rng(seed)
    t0 = time.time()
    overall = []
    for kind in ("frustum", "prism", "pyramid"):
        hist = np.zeros(6, np.int64)
        best = (np.inf, None)
        done = 0
        while done < N // 3:
            n = min(100_000, N // 3 - done)
            V, faces = gen(rng, n, kind)
            good = convex_ok(V, faces)
            M, vol = margins(V[good], faces)
            cnt = (M > 0).sum(1)
            hist += np.bincount(cnt, minlength=6)
            sec = np.sort(M, 1)[:, -2]
            k = int(np.argmin(sec))
            if sec[k] < best[0]:
                best = (sec[k], V[good][k])
            done += n
        overall.append((kind, best[0]))
        print(f"{kind:8s}: {hist.sum()} valid samples, #stable-faces histogram 0..5 = {hist.tolist()}, "
              f"min second-best margin = {best[0]:.3e}", flush=True)
    print(f"random sampling done in {time.time() - t0:.0f}s")


# ------------------------------------------------------------------ local optimisation
def shape_from_params(z, kind):
    z = np.asarray(z, float)
    if kind == "frustum":
        d = z[:9].reshape(3, 3)
        a = 1 / (1 + np.exp(-z[9:12]))
        V = np.concatenate([d, a[:, None] * d], 0)
        return V, PRISM_FACES
    if kind == "prism":
        B = z[:9].reshape(3, 3)
        dv = z[9:12]
        s = np.exp(z[12:15])
        V = np.concatenate([B, B + s[:, None] * dv[None, :]], 0)
        return V, PRISM_FACES
    if kind == "pyramid":   # frustum with a_3 = 1 exactly (apex = B3), general base via rays
        d = z[:9].reshape(3, 3)
        a = 1 / (1 + np.exp(-z[9:11]))
        A = np.stack([a[0] * d[0], a[1] * d[1]])
        # vertices: B1 B2 B3 A1 A2 ; base quad = (B1, B2, A2, A1) ; apex B3
        V = np.stack([d[0], d[1], A[1], A[0], d[2]])
        return V, PYR_FACES
    if kind == "pyr_par":   # pyramid over a parallelogram / trapezoid (parallel sides ratio lam)
        P = z[:9].reshape(3, 3)          # q0 = 0, q1 = P0, q3 = P1, apex = P2
        lam = 1 / (1 + np.exp(-z[9]))
        q0 = np.zeros(3); q1 = P[0]; q3 = P[1]; q2 = q3 + lam * P[0]
        V = np.stack([q0, q1, q2, q3, P[2]])
        return V, PYR_FACES
    raise ValueError(kind)


NPAR = {"frustum": 12, "prism": 15, "pyramid": 11, "pyr_par": 10}


def objective(z, kind):
    V, faces = shape_from_params(z, kind)
    V = V[None]
    if not convex_ok(V, faces)[0]:
        return 10.0
    M, vol = margins(V, faces)
    if not np.isfinite(M).all() or vol[0] <= 1e-12 * np.abs(V).max() ** 3:
        return 10.0
    srt = np.sort(M[0])
    return srt[-2]


def opt_worker(args):
    from scipy.optimize import minimize
    seed, kind = args
    rng = np.random.default_rng(seed)
    best = (np.inf, None)
    for _ in range(3):
        z0 = rng.normal(size=NPAR[kind]) * rng.choice([0.5, 1.0, 2.0])
        if objective(z0, kind) >= 10:
            continue
        r = minimize(objective, z0, args=(kind,), method="Nelder-Mead",
                     options={"maxiter": 4000, "maxfev": 6000, "xatol": 1e-10, "fatol": 1e-12, "adaptive": True})
        r = minimize(objective, r.x, args=(kind,), method="Nelder-Mead",
                     options={"maxiter": 4000, "maxfev": 6000, "xatol": 1e-12, "fatol": 1e-14, "adaptive": True})
        if r.fun < best[0]:
            best = (r.fun, r.x)
    return kind, seed, best[0], (None if best[1] is None else best[1].tolist())


def describe(z, kind):
    V, faces = shape_from_params(np.array(z), kind)
    M, vol = margins(V[None], faces)
    sv = np.linalg.svd(V - V.mean(0), compute_uv=False)
    diam = max(np.linalg.norm(V[i] - V[j]) for i in range(len(V)) for j in range(len(V)))
    return (f"margins={np.round(np.sort(M[0]), 5).tolist()}, singular values of vertex cloud "
            f"(aspect) {np.round(sv / sv[0], 4).tolist()}, vol/diam^3={vol[0] / 6 / diam ** 3:.3e}")


def run_opt(restarts, seed, workers):
    from concurrent.futures import ProcessPoolExecutor
    t0 = time.time()
    kinds = ["frustum", "prism", "pyramid", "pyr_par"]
    tasks = [(seed * 100000 + i, kinds[i % 4]) for i in range(restarts)]
    res = {k: [] for k in kinds}
    with ProcessPoolExecutor(workers) as ex:
        for kind, sd, val, z in ex.map(opt_worker, tasks, chunksize=2):
            res[kind].append((val, sd, z))
    for kind in kinds:
        vals = sorted(res[kind], key=lambda t: t[0])
        nneg = sum(1 for v in vals if v[0] <= 0)
        print(f"{kind:8s}: {len(vals)} restarts x 3 NM runs; best second-best margin {vals[0][0]:.4e}; "
              f"#runs ending <= 0: {nneg}; quartiles {[round(vals[len(vals) * q // 4][0], 5) for q in range(4)]}")
        print("     best shape:", describe(vals[0][2], kind))
    print(f"optimisation done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    if sys.argv[1] == "random":
        run_random(int(sys.argv[2]), int(sys.argv[3]))
    elif sys.argv[1] == "opt":
        run_opt(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))


# ------------------------------------------------------------------ vectorised (1+1)-ES with many parallel chains
def shapes_batch(Z, kind):
    n = Z.shape[0]
    sig = lambda t: 1 / (1 + np.exp(-t))
    if kind == "frustum":
        d = Z[:, :9].reshape(n, 3, 3); a = sig(Z[:, 9:12])
        return np.concatenate([d, a[:, :, None] * d], 1), PRISM_FACES
    if kind == "prism":
        B = Z[:, :9].reshape(n, 3, 3); dv = Z[:, 9:12]; s = np.exp(np.clip(Z[:, 12:15], -30, 30))
        return np.concatenate([B, B + s[:, :, None] * dv[:, None, :]], 1), PRISM_FACES
    if kind == "pyramid":
        d = Z[:, :9].reshape(n, 3, 3); a = sig(Z[:, 9:11])
        return np.stack([d[:, 0], d[:, 1], a[:, 1:2] * d[:, 1], a[:, 0:1] * d[:, 0], d[:, 2]], 1), PYR_FACES
    if kind == "pyr_par":
        P = Z[:, :9].reshape(n, 3, 3); lam = sig(Z[:, 9:10])
        q0 = np.zeros((n, 3))
        return np.stack([q0, P[:, 0], P[:, 1] + lam * P[:, 0], P[:, 1], P[:, 2]], 1), PYR_FACES


def obj_batch(Z, kind):
    V, faces = shapes_batch(Z, kind)
    with np.errstate(all="ignore"):
        good = convex_ok(V, faces)
        M, vol = margins(V, faces)
        sec = np.sort(M, 1)[:, -2]
    diam = np.max(np.linalg.norm(V[:, :, None, :] - V[:, None, :, :], axis=3), axis=(1, 2))
    fat = vol / 6 / np.where(diam > 0, diam, 1) ** 3
    bad = ~good | ~np.isfinite(sec) | (vol <= 1e-10 * np.abs(V).max(axis=(1, 2)) ** 3) | (fat < FATMIN)
    return np.where(bad, 10.0, sec)


def run_es(chains, steps, seed):
    rng = np.random.default_rng(seed)
    t0 = time.time()
    summary = {}
    for kind in ("frustum", "prism", "pyramid", "pyr_par"):
        k = NPAR[kind]
        Z = rng.normal(size=(chains, k)) * rng.choice([0.5, 1.0, 2.0], size=(chains, 1))
        f = obj_batch(Z, kind)
        for _ in range(20):          # re-draw invalid starts
            bad = f >= 10
            if not bad.any():
                break
            Z[bad] = rng.normal(size=(bad.sum(), k)); f[bad] = obj_batch(Z[bad], kind)
        step = np.full(chains, 0.3)
        for it in range(steps):
            Zn = Z + step[:, None] * rng.normal(size=Z.shape)
            fn = obj_batch(Zn, kind)
            acc = fn < f
            Z[acc] = Zn[acc]; f[acc] = fn[acc]
            step = np.where(acc, step * 1.5, step * 0.9)
            step = np.clip(step, 1e-9, 2.0)
            if it % 500 == 499:      # restart stagnated chains far from the current best value
                stag = step < 1e-7
                if stag.any():
                    Z[stag] = rng.normal(size=(stag.sum(), k)); f[stag] = obj_batch(Z[stag], kind); step[stag] = 0.3
        order = np.argsort(f)
        # Nelder-Mead polish of the 10 best chains
        from scipy.optimize import minimize
        pol = []
        for j in order[:6]:
            r = minimize(lambda z: float(obj_batch(z[None], kind)[0]), Z[j], method="Nelder-Mead",
                         options={"maxfev": 2000, "xatol": 1e-12, "fatol": 1e-14, "adaptive": True})
            pol.append((r.fun, r.x))
        pol.sort(key=lambda t: t[0])
        nneg = int((f <= 0).sum())
        summary[kind] = pol[0][0]
        print(f"{kind:8s}: {chains} chains x {steps} ES steps; best second-best margin after ES {f[order[0]]:.4e}, "
              f"after NM polish {pol[0][0]:.4e}; chains ending <= 0: {nneg}; median final {np.median(f):.4e}", flush=True)
        print("     best shape:", describe(pol[0][1], kind), flush=True)
    print(f"ES done in {time.time() - t0:.0f}s")


FATMIN = 0.0
if __name__ == "__main__" and sys.argv[1] == "es":
    if len(sys.argv) > 5:
        FATMIN = float(sys.argv[5])
        print("non-degeneracy constraint: vol/diam^3 >=", FATMIN)
    run_es(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
