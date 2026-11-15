"""Probe: are V=9 trees certifiable by the Gram certificates (single cone / two cones / all cones)?
Random (graph, root, tree) samples, plus deep 'path-like' trees which are the hard ones."""
import os, sys, random, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, 'lowerbound'); sys.path.insert(0, 'lb_joint')
import numpy as np, networkx as nx
import enumerate_cases as ec
import certify_joint as cj
import general_cert as gc

def random_tree(G, v, rng, deep=False):
    """random spanning tree rooted at v (parent array); deep=True: DFS-like long paths"""
    V = G.number_of_nodes(); parent = [-1]*V
    if not deep:   # random BFS/Prim-like growth
        intree = {v}; frontier = [(v, w) for w in G[v]]
        while len(intree) < V:
            a, b = frontier.pop(rng.randrange(len(frontier)))
            if b in intree: continue
            parent[b] = a; intree.add(b); frontier += [(b, w) for w in G[b] if w not in intree]
    else:          # random DFS (long paths)
        seen = {v}; stack = [v]
        while stack:
            a = stack[-1]
            nxt = [w for w in G[a] if w not in seen]
            if not nxt: stack.pop(); continue
            b = rng.choice(nxt); parent[b] = a; seen.add(b); stack.append(b)
    return parent

def probe(task):
    V, gid, edges, faces, v, parent = task
    tets = [f for f in faces if v not in f]
    z, mu, Y, st = cj.solve_joint(V, v, parent, tets)
    out = dict(gid=gid, root=v, parent=parent, z1=z, depth=max(gc.depths(parent, v)))
    if z > 1e-7: out["level"] = "single"; return out
    nb = gc.neighbours(V, edges)
    d = gc.depths(parent, v)
    best = -1
    for q in sorted((u for u in range(V) if u != v), key=lambda u: (d[u], u)):
        res = gc.solve_general(V, v, parent, nb, faces, [v, q], False)
        best = max(best, res["z"])
        if res["z"] > 1e-7: out["level"] = "pair"; out["z2"] = res["z"]; out["q"] = q; return out
    out["z2"] = best
    res = gc.solve_general(V, v, parent, nb, faces, list(range(V)), False)
    out["zall"] = res["z"]
    if res["z"] > 1e-7: out["level"] = "all"; return out
    res = gc.solve_general(V, v, parent, nb, faces, list(range(V)), True)
    out["zall_s1"] = res["z"]
    out["level"] = "all+S1" if res["z"] > 1e-7 else "NONE"
    return out

if __name__ == "__main__":
    V = int(sys.argv[1]); n_rand = int(sys.argv[2]); n_deep = int(sys.argv[3]); workers = int(sys.argv[4])
    rng = random.Random(1)
    graphs = ec.triangulations(V)
    print(V, "vertices:", len(graphs), "triangulations", flush=True)
    tasks = []
    for k in range(n_rand + n_deep):
        gid = rng.randrange(len(graphs)); G = graphs[gid]
        v = rng.randrange(V)
        parent = random_tree(G, v, rng, deep=(k >= n_rand))
        tasks.append((V, gid, sorted(G.edges()), ec.faces_of(G), v, parent))
    from concurrent.futures import ProcessPoolExecutor
    import collections
    t0 = time.time(); cnt = collections.Counter(); bad = []
    with ProcessPoolExecutor(workers) as ex:
        for i, out in enumerate(ex.map(probe, tasks, chunksize=4)):
            kind = "rand" if i < n_rand else "deep"
            cnt[(kind, out["level"])] += 1
            if out["level"] in ("NONE", "all", "all+S1"): bad.append(out)
            if (i+1) % 200 == 0: print(i+1, dict(cnt), f"{time.time()-t0:.0f}s", flush=True)
    print("FINAL", dict(cnt), f"{time.time()-t0:.0f}s")
    json.dump(bad, open(f'explore/probe_V{V}_hard.json', 'w'))
    for b in bad[:10]: print(b)
