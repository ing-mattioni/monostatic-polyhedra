#!/usr/bin/env python3
"""
flip_closure.py -- completeness of the list of simplicial 2-spheres in a certificate file, proved from
Wagner's flip theorem instead of a literature count.  Standard library only.

    python flip_closure.py V certs_V{V}.jsonl.gz

WAGNER'S THEOREM (K. Wagner, "Bemerkungen zum Vierfarbenproblem", Jber. DMV 46 (1936) 26-32): for n >= 4,
any two triangulations of the 2-sphere with n vertices (maximal planar simple graphs) are connected by a
finite sequence of edge flips.  (A flip replaces the edge ab of the two triangles abc, abd by the edge cd,
allowed when cd is not already an edge.)

CHECKS: the listed complexes are simplicial 2-spheres with pairwise non-isomorphic graphs; for every listed
sphere and every flippable edge, the flipped graph is isomorphic to a listed one.  The list is non-empty
and closed under flips, so by Wagner's theorem it contains (a copy of) every triangulation with V
vertices: the list is complete.  (For a triangulation the graph determines the triangles -- the faces are
exactly the triangles of the graph that are not separating -- so graph isomorphism is the right notion;
the flipped triangle set is computed explicitly and re-checked to be a 2-sphere.)
"""
import gzip
import json
import sys
import time


def need(cond, msg):
    if not cond:
        raise RuntimeError(msg)


def check_sphere(V, E, F):
    """E: set of frozenset edges, F: set of frozenset triangles"""
    need(len(E) == 3 * V - 6 and len(F) == 2 * V - 4, "wrong numbers of edges/triangles")
    nbrs = {u: set() for u in range(V)}
    for e in E:
        a, b = tuple(e)
        nbrs[a].add(b)
        nbrs[b].add(a)
    for e in E:
        need(sum(1 for f in F if e <= f) == 2, "edge not in exactly two triangles")
    for f in F:
        a, b, c = tuple(f)
        need(frozenset((a, b)) in E and frozenset((a, c)) in E and frozenset((b, c)) in E, "triangle side not an edge")
    for u in range(V):
        link = [tuple(f - {u}) for f in F if u in f]
        deg = {x: 0 for x in nbrs[u]}
        for a, b in link:
            need(a in deg and b in deg, "link vertex not a neighbour")
            deg[a] += 1
            deg[b] += 1
        need(len(nbrs[u]) >= 3 and all(d == 2 for d in deg.values()), "vertex link is not a cycle")
        start = next(iter(nbrs[u]))
        seen, stack = {start}, [start]
        while stack:
            x = stack.pop()
            for a, b in link:
                for p, q in ((a, b), (b, a)):
                    if p == x and q not in seen:
                        seen.add(q)
                        stack.append(q)
        need(seen == nbrs[u], "vertex link is not a single cycle")
    need(V - len(E) + len(F) == 2, "Euler characteristic is not 2")
    return nbrs


def isomorphic(V, n1, n2):
    """backtracking graph isomorphism test (adjacency dicts of sets), with degree pruning"""
    d1 = sorted(len(n1[u]) for u in range(V))
    d2 = sorted(len(n2[u]) for u in range(V))
    if d1 != d2:
        return False
    order = sorted(range(V), key=lambda u: -len(n1[u]))
    img = {}
    used = set()

    def rec(k):
        if k == V:
            return True
        u = order[k]
        for w in range(V):
            if w in used or len(n2[w]) != len(n1[u]):
                continue
            if all(((x in n1[u]) == (img[x] in n2[w])) for x in img):
                img[u] = w
                used.add(w)
                if rec(k + 1):
                    return True
                del img[u]
                used.discard(w)
        return False
    return rec(0)


def flips(V, E, F):
    """all triangle sets obtained by one flip: (E', F')"""
    out = []
    for e in E:
        a, b = tuple(e)
        tri = [f for f in F if e <= f]
        c = next(iter(tri[0] - e))
        d = next(iter(tri[1] - e))
        if frozenset((c, d)) in E:
            continue
        E2 = (E - {e}) | {frozenset((c, d))}
        F2 = (F - set(tri)) | {frozenset((a, c, d)), frozenset((b, c, d))}
        out.append((E2, F2))
    return out


def main(V, path):
    t0 = time.time()
    listed = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            rec = json.loads(line)
            if rec.get("kind") == "graph":
                E = {frozenset(e) for e in rec["edges"]}
                F = {frozenset(f) for f in rec["faces"]}
                listed.append((E, F, check_sphere(V, E, F)))
            elif rec.get("kind", "").startswith("cert"):
                break                                    # graph records come first
    need(len(listed) >= 1, "no spheres listed")
    for i in range(len(listed)):
        for j in range(i):
            need(not isomorphic(V, listed[i][2], listed[j][2]), f"listed spheres {j} and {i} are isomorphic")
    n_flips = 0
    for k, (E, F, nb) in enumerate(listed):
        for E2, F2 in flips(V, E, F):
            n_flips += 1
            nb2 = check_sphere(V, E2, F2)
            need(any(isomorphic(V, nb2, other[2]) for other in listed),
                 f"a flip of sphere {k} is not isomorphic to any listed sphere: the list is NOT closed")
    print(f"V = {V}: {len(listed)} listed simplicial 2-spheres, pairwise non-isomorphic; all {n_flips} flips "
          f"stay in the list (up to isomorphism)  [{time.time() - t0:.1f}s]")
    print(f"COMPLETE (by Wagner's flip theorem): every triangulation of the sphere with {V} vertices is "
          f"isomorphic to a listed one.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main(int(sys.argv[1]), sys.argv[2]))
    except RuntimeError as e:
        print("CHECK FAILED:", e)
        sys.exit(1)
