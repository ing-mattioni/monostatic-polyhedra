#!/usr/bin/env python3
"""
enumerate_spheres.py -- self-contained check that the list of simplicial 2-spheres used by the
certificates is COMPLETE, without relying on the literature count (OEIS A000109).

    python enumerate_spheres.py V certs_V{V}.jsonl.gz        (V <= 8; standard library only)

Method (exhaustive, no topology beyond the definition):
  1. Enumerate ALL labelled simplicial 2-spheres on the vertex set {0..V-1}: every set E of 3V-6 edges
     (minimum degree >= 3, every edge in at least two 3-cliques) is generated, and for each E all sets F
     of 2V-4 triangles (3-cliques of E) such that every edge lies in exactly two triangles of F are found
     by exhaustive backtracking; (E, F) is counted iff it is a simplicial 2-sphere: the link of every
     vertex is a single cycle through all its neighbours, the graph is connected and V - E + F = 2.
     (Any simplicial 2-sphere has exactly 3V-6 edges and 2V-4 triangles, every edge in exactly two
     triangles, and minimum degree >= 3, so nothing is missed.)  Let N be the number found.
  2. For each sphere listed in the certificate file compute, by brute force over all V! permutations, the
     number a_k of permutations preserving its edges and triangles; it has V!/a_k labelled copies.
  3. The listed spheres are pairwise non-isomorphic (checked by the main verifier and again here), so their
     labelled copies are distinct; if  sum_k V!/a_k == N  then every labelled sphere is a copy of a listed
     one, i.e. the list is complete.
"""
import gzip
import itertools
import json
import math
import sys
import time
from multiprocessing import Pool


def is_sphere(V, E, F):
    """E: set of frozenset pairs, F: list of frozenset triples (each edge in exactly two of them)."""
    if V - len(E) + len(F) != 2:
        return False
    nbrs = {u: set() for u in range(V)}
    for e in E:
        a, b = tuple(e)
        nbrs[a].add(b)
        nbrs[b].add(a)
    for u in range(V):
        link = [tuple(f - {u}) for f in F if u in f]
        deg = {x: 0 for x in nbrs[u]}
        for a, b in link:
            if a not in deg or b not in deg:
                return False
            deg[a] += 1
            deg[b] += 1
        if len(nbrs[u]) < 3 or any(d != 2 for d in deg.values()):
            return False
        start = next(iter(nbrs[u]))
        seen, stack = {start}, [start]
        while stack:
            x = stack.pop()
            for a, b in link:
                for p, q in ((a, b), (b, a)):
                    if p == x and q not in seen:
                        seen.add(q)
                        stack.append(q)
        if seen != nbrs[u]:
            return False
    seen, stack = {0}, [0]
    while stack:
        x = stack.pop()
        for y in nbrs[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen) == V


def count_structures(V, adj):
    """adj: list of neighbour bitmasks of a graph with 3V-6 edges.  Number of simplicial 2-sphere
    structures (sets of triangles) on this graph, by exhaustive exact-2-cover backtracking."""
    edges = [(a, b) for a in range(V) for b in range(a + 1, V) if adj[a] >> b & 1]
    eidx = {e: k for k, e in enumerate(edges)}
    tris = [(a, b, c) for a in range(V) for b in range(a + 1, V) if adj[a] >> b & 1
            for c in range(b + 1, V) if adj[a] >> c & 1 and adj[b] >> c & 1]
    tri_edges = [(eidx[(a, b)], eidx[(a, c)], eidx[(b, c)]) for a, b, c in tris]
    need = 2 * V - 4
    by_edge = [[] for _ in edges]
    for t, es in enumerate(tri_edges):
        for e in es:
            by_edge[e].append(t)
    if any(len(l) < 2 for l in by_edge):
        return 0
    E = {frozenset(e) for e in edges}
    cover = [0] * len(edges)
    chosen = []
    found = set()

    def rec2(e):
        while e < len(edges) and cover[e] == 2:
            e += 1
        if e == len(edges):
            if len(chosen) == need:
                found.add(frozenset(chosen))
            return
        if len(chosen) >= need:
            return
        for t in by_edge[e]:
            if t in chosen or any(cover[x] >= 2 for x in tri_edges[t]):
                continue
            for x in tri_edges[t]:
                cover[x] += 1
            chosen.append(t)
            rec2(e)
            chosen.pop()
            for x in tri_edges[t]:
                cover[x] -= 1

    rec2(0)
    return sum(1 for S in found if is_sphere(V, E, [frozenset(tris[t]) for t in S]))


def worker(args):
    """Enumerate all graphs on V vertices with exactly 3V-6 edges whose first non-edges (in the fixed
    order of vertex pairs) are the given prefix; count sphere structures."""
    V, prefix = args
    pairs = [(a, b) for a in range(V) for b in range(a + 1, V)]
    n_non = len(pairs) - (3 * V - 6)
    max_removed = (V - 1) - 3                      # each vertex keeps degree >= 3
    full = [((1 << V) - 1) & ~(1 << a) for a in range(V)]
    adj = list(full)
    removed = [0] * V
    count_graphs = 0
    count_spheres = 0

    def remove(k):
        a, b = pairs[k]
        adj[a] &= ~(1 << b)
        adj[b] &= ~(1 << a)
        removed[a] += 1
        removed[b] += 1

    def restore(k):
        a, b = pairs[k]
        adj[a] |= 1 << b
        adj[b] |= 1 << a
        removed[a] -= 1
        removed[b] -= 1

    def finish():
        nonlocal count_graphs, count_spheres
        # every edge must lie in at least two 3-cliques
        for a in range(V):
            m = adj[a] >> (a + 1)
            b = a + 1
            while m:
                if m & 1 and bin(adj[a] & adj[b]).count("1") < 2:
                    return
                m >>= 1
                b += 1
        count_graphs += 1
        count_spheres += count_structures(V, adj)

    def rec(k_start, left):
        if left == 0:
            finish()
            return
        for k in range(k_start, len(pairs) - left + 1):
            a, b = pairs[k]
            if removed[a] >= max_removed or removed[b] >= max_removed:
                continue
            remove(k)
            rec(k + 1, left - 1)
            restore(k)

    ok = True
    for k in prefix:
        a, b = pairs[k]
        if removed[a] >= max_removed or removed[b] >= max_removed:
            ok = False
            break
        remove(k)
    if ok:
        rec((prefix[-1] + 1) if prefix else 0, n_non - len(prefix))
    return count_graphs, count_spheres


def automorphisms(V, E, F):
    n = 0
    for perm in itertools.permutations(range(V)):
        if all(frozenset(perm[x] for x in e) in E for e in E) and all(frozenset(perm[x] for x in f) in F for f in F):
            n += 1
    return n


def main(V, path, workers=12):
    t0 = time.time()
    pairs = [(a, b) for a in range(V) for b in range(a + 1, V)]
    n_non = len(pairs) - (3 * V - 6)
    depth = min(3, n_non)
    tasks = [(V, list(pref)) for pref in itertools.combinations(range(len(pairs)), depth)]
    with Pool(workers) as pool:
        res = pool.map(worker, tasks, chunksize=8)
    n_graphs = sum(r[0] for r in res)
    N = sum(r[1] for r in res)
    print(f"V = {V}: exhaustive enumeration: {n_graphs} labelled graphs with {3 * V - 6} edges pass the filters; "
          f"{N} labelled simplicial 2-spheres  [{time.time() - t0:.0f}s]")
    listed = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            rec = json.loads(line)
            if rec.get("kind") == "graph":
                E = {frozenset(e) for e in rec["edges"]}
                F = {frozenset(f) for f in rec["faces"]}
                # explicit checks (not `assert`: they must survive `python -O`)
                if not (is_sphere(V, E, list(F)) and len(E) == 3 * V - 6 and len(F) == 2 * V - 4):
                    raise RuntimeError("a listed complex is not a simplicial 2-sphere")
                listed.append((E, F))
    # pairwise non-isomorphic (as graphs)
    for (E1, _), (E2, _) in itertools.combinations(listed, 2):
        d1 = sorted(sum(1 for e in E1 if u in e) for u in range(V))
        d2 = sorted(sum(1 for e in E2 if u in e) for u in range(V))
        if d1 == d2:
            if any(all(frozenset(p[x] for x in e) in E2 for e in E1) for p in itertools.permutations(range(V))):
                raise RuntimeError("two listed spheres are isomorphic")
    total = 0
    for k, (E, F) in enumerate(listed):
        a = automorphisms(V, E, F)
        if math.factorial(V) % a != 0:
            raise RuntimeError("automorphism count does not divide V!")
        total += math.factorial(V) // a
    print(f"listed spheres: {len(listed)}, pairwise non-isomorphic; labelled copies sum_k V!/|Aut_k| = {total}  "
          f"[{time.time() - t0:.0f}s]")
    if total == N:
        print(f"COMPLETE: every simplicial 2-sphere on {V} vertices is isomorphic to a listed one "
              f"(no external count used).")
        return 0
    print("MISMATCH: the list is NOT shown to be complete.")
    return 1


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]), sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 12))
