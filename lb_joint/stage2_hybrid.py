"""
stage2_hybrid.py -- second stage of the V = 9 computation: certificates for the trees that the single-cone
certificate (certify_joint.py, stage 1) does not certify.

For every tree in uncertified_V9.json, in this order:
  1. hybrid certificate (hybrid_cert.py, record "certH"): the "bad" triangles are chosen from the dual of
     the (failing) single-cone SDP (largest dual blocks), then by tree depth; second apex q = the shallow
     vertices of the tree;
  2. fallback: the two-cone certificate with a coupling matrix (general_cert.py, record "certG");
  3. otherwise the tree is reported as a survivor (no theorem for this V then).
Writes  exp_hybrid_V9.jsonl.gz  (certificates) and  survivors_hybrid_V9.json.  All certificates are
re-checked exactly here and, independently, by verify_joint.py.

    python stage2_hybrid.py [workers]
"""
import gzip
import json
import os
import sys
import time

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "lowerbound"))

V = 9
_GD = None


def _graphs():
    global _GD
    if _GD is None:
        import enumerate_cases as ec
        _GD = ec.graph_data(V, with_trees=False)
    return _GD


def work(task):
    import numpy as np
    import certify_joint as cj
    import general_cert as gc
    import hybrid_cert as hc
    import strengthen_joint as sj
    gid, v, parent = task
    gd = _graphs()[gid]
    faces = gd.faces
    d = hc.depths(parent, v)
    qs = sorted((u for u in range(V) if u != v), key=lambda u: (d[u], u))
    Cv = [tuple(f) for f in faces if v not in f]
    tried = set()
    # 1a. dual-guided bad sets
    res = sj.solve_multi(V, v, parent, [[cj.tet_indicator(V, v, f)] for f in Cv], 1, want_dual=True)
    tr = np.array([np.trace(X) for X in res["X"]])
    candidates = []
    for frac in (0.2, 0.05, 0.01):
        candidates.append(frozenset(Cv[i] for i in range(len(Cv)) if tr[i] >= frac * tr.max()))
    # 1b. depth-ranked bad sets
    ranked = sorted(Cv, key=lambda f: (-min(d[x] for x in f), -sum(d[x] for x in f), f))
    for k in (1, 2, 3, 4):
        candidates.append(frozenset(ranked[:k]))
    attempts = 0
    for bad in candidates:
        for q in qs[:3]:
            if (bad, q) in tried or not bad:
                continue
            tried.add((bad, q))
            attempts += 1
            z, rc, labels = hc.try_hybrid(V, v, parent, faces, q, set(bad))
            if rc is not None:
                mus, Y, bits = rc
                single = [[list(f), y] for (f, g), y in zip(labels, Y) if g is None]
                pairs = [[list(f), list(g), y] for (f, g), y in zip(labels, Y) if g is not None]
                line = json.dumps(dict(kind="certH", gid=gid, root=v, parent=list(parent), q=q,
                                       mu=mus[0], mu2=mus[1], z=float(f"{z:.6g}"),
                                       single=single, pairs=pairs), separators=(",", ":"))
                return gid, v, parent, line, "H", attempts
    # 2. fallback: coupled two-cone certificate (certG)
    nb = gc.neighbours(V, gd.edges)
    for q in qs:
        r = gc.solve_general(V, v, parent, nb, faces, [v, q], False)
        attempts += 1
        if r["z"] > 1e-7:
            rc = gc.round_general(V, v, parent, [v, q], r)
            if rc is not None:
                groups, bits = rc
                line = json.dumps(dict(kind="certG", gid=gid, root=v, parent=list(parent), s1=False,
                                       z=float(f"{r['z']:.6g}"), groups=groups), separators=(",", ":"))
                return gid, v, parent, line, "G", attempts
    return gid, v, parent, None, "none", attempts


if __name__ == "__main__":
    import multiprocessing as mp
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    t0 = time.time()
    unc = json.load(open(os.path.join(HERE, "uncertified_V9.json")))
    tasks = [(c["graph"], c["root"], tuple(c["parent"])) for c in unc]
    del unc
    print(f"stage 2 (hybrid): {len(tasks)} trees, {workers} workers", flush=True)
    out = os.path.join(HERE, "exp_hybrid_V9.jsonl.gz")
    tmp = out + ".part"
    counts = {"H": 0, "G": 0, "none": 0}
    survivors = []
    att = 0
    with gzip.open(tmp, "wt", compresslevel=6) as fo, mp.get_context("spawn").Pool(workers) as pool:
        fo.write(json.dumps(dict(kind="header_certH", V=V)) + "\n")
        for k, (gid, v, parent, line, how, a) in enumerate(pool.imap_unordered(work, tasks, chunksize=16)):
            counts[how] += 1
            att += a
            if line:
                fo.write(line + "\n")
            else:
                survivors.append(dict(graph=gid, root=v, parent=list(parent)))
            if (k + 1) % 5000 == 0:
                el = time.time() - t0
                print(f"  {k + 1}/{len(tasks)}  {counts}  mean attempts {att / (k + 1):.2f}  {el:.0f}s  "
                      f"(eta {el / (k + 1) * (len(tasks) - k - 1) / 3600:.1f} h)", flush=True)
    os.replace(tmp, out)
    json.dump(survivors, open(os.path.join(HERE, "survivors_hybrid_V9.json"), "w"))
    print(f"stage 2 done: {counts}; survivors {len(survivors)}; {time.time() - t0:.0f}s; wrote {out}", flush=True)
