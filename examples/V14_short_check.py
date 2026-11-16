#!/usr/bin/env python3
"""
V14_short_check.py -- a short, self-contained exact check that the 14-point polyhedron of Theorem C is
mono-unstable.  Standard library only, integers only; about one second.  It does not compute a convex
hull: the 24 boundary triangles are given below and the script only CHECKS them.

    python examples/V14_short_check.py

THE ARGUMENT (each step is checked below; numbers refer to the checks).
 (1) For every listed triangle abc (counter-clockwise seen from outside), n = (b-a) x (c-a) satisfies
     n.(r_k - a) < 0 for each of the other 11 points r_k.  So the plane of abc supports the point set and
     meets it in a, b, c only: abc is a facet of the convex hull P, and no point outside abc lies on it.
 (2) Every edge of a listed triangle lies in exactly two listed triangles, with opposite orientations,
     and every point is a corner of some triangle.  Each edge of P lies in exactly two facets of P and
     the facets of P are connected through edges, so the listed triangles are ALL the facets of P: P has
     14 vertices, 36 edges and 24 triangular faces, and the neighbours of a vertex are the vertices it
     shares an edge with.
 (3) With the outward orientation, D = sum det(a,b,c) = 6 vol(P) > 0 and C = sum det(a,b,c)(a+b+c) =
     24 vol(P) c, where c is the centre of mass (divergence theorem; the cones from the origin over
     the oriented facets add up to P with signs).  So 4D (r_i - c) = 4D r_i - C =: g_i is an integer vector.
 (4) Lemma 2.1 of the paper (origin at c): vertex i is NOT an unstable equilibrium iff some neighbour j
     has (r_j - r_i).(r_i - c) >= 0, i.e. s(i,j) := (r_j - r_i).g_i >= 0.  The table lists one neighbour
     with s > 0 for each of the vertices 1..13, and for vertex 0 every neighbour has s < 0.
 Hence vertex 0 is the only unstable vertex: P is mono-unstable.
"""
import sys

POINTS = [(-4, 1000, -1), (-477, 638, -3), (-633, 344, -2), (-657, 62, -4), (-577, -182, 5),
          (-432, -355, -11), (-242, -454, 13), (181, -490, -12), (359, -426, 12), (518, -295, -6),
          (630, -103, 5), (670, 128, 3), (626, 374, 2), (474, 636, 2)]
TRIANGLES = [(1, 0, 5), (1, 6, 0), (5, 0, 7), (6, 13, 0), (0, 9, 7), (0, 13, 9), (3, 2, 1), (2, 4, 1),
             (3, 1, 5), (4, 6, 1), (3, 4, 2), (3, 5, 4), (4, 5, 6), (5, 7, 6), (6, 7, 8), (6, 8, 13),
             (7, 9, 8), (8, 9, 10), (8, 10, 11), (8, 11, 12), (8, 12, 13), (9, 11, 10), (9, 12, 11),
             (13, 12, 9)]
UNSTABLE = 0
WITNESS = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 8, 8: 9, 9: 10, 10: 11, 11: 12, 12: 13, 13: 0}


def need(cond, msg):
    if not cond:
        print("CHECK FAILED:", msg)
        sys.exit(1)


def sub(p, q):
    return tuple(a - b for a, b in zip(p, q))


def dot(p, q):
    return sum(a * b for a, b in zip(p, q))


def cross(p, q):
    return (p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0])


def main():
    n = len(POINTS)
    # (1) every triangle is a facet, strictly supporting
    for t in TRIANGLES:
        a, b, c = (POINTS[k] for k in t)
        nrm = cross(sub(b, a), sub(c, a))
        for k in range(n):
            if k not in t:
                need(dot(nrm, sub(POINTS[k], a)) < 0, f"point {k} is not strictly inside the plane of {t}")
    # (2) closed, consistently oriented surface through all points
    directed = [(t[i], t[(i + 1) % 3]) for t in TRIANGLES for i in range(3)]
    need(len(set(directed)) == len(directed), "a directed edge occurs twice")
    need(all((b, a) in set(directed) for a, b in directed), "an edge lies in only one triangle")
    edges = {frozenset(e) for e in directed}
    need({k for t in TRIANGLES for k in t} == set(range(n)), "some point is not a vertex")
    need(n - len(edges) + len(TRIANGLES) == 2, "Euler characteristic is not 2")
    nbrs = {i: {j for e in edges if i in e for j in e if j != i} for i in range(n)}
    # (3) volume and centre of mass
    D, C = 0, (0, 0, 0)
    for t in TRIANGLES:
        a, b, c = (POINTS[k] for k in t)
        d = dot(a, cross(b, c))
        D += d
        C = tuple(C[m] + d * (a[m] + b[m] + c[m]) for m in range(3))
    need(D > 0, "orientation is not outward")
    print(f"V = {n}, E = {len(edges)}, F = {len(TRIANGLES)};  6 vol = D = {D};  24 vol c = C = {C}")
    print(f"centre of mass c = C/(4D) = ({C[0]}/{4 * D}, {C[1]}/{4 * D}, {C[2]}/{4 * D})")
    # (4) equilibria at the vertices
    print(" vertex  witness j   s(i,j) = (r_j - r_i).(4D r_i - C)")
    for i in range(n):
        g = tuple(4 * D * POINTS[i][m] - C[m] for m in range(3))
        s = {j: dot(sub(POINTS[j], POINTS[i]), g) for j in nbrs[i]}
        if i == UNSTABLE:
            need(all(val < 0 for val in s.values()), f"vertex {i} is not unstable")
            print(f"  {i:4d}   (all {len(s)} neighbours below)  max s = {max(s.values())}")
        else:
            j = WITNESS.get(i)
            need(j in nbrs[i], f"witness {j} is not a neighbour of {i}")
            need(s[j] > 0, f"vertex {i} is not shadowed by {j}")
            print(f"  {i:4d}   {j:4d}        {s[j]:>16d}")
    print("CHECKED: exactly one unstable vertex (vertex 0) -- the polyhedron is mono-unstable.")


if __name__ == "__main__":
    main()
