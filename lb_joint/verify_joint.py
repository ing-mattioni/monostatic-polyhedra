#!/usr/bin/env python3
"""
verify_joint.py -- standalone exact verifier for the "joint mu" certificates.

    python verify_joint.py V certs_V{V}.jsonl.gz [more_certs.jsonl.gz ...]

Uses only the Python standard library (exact integer arithmetic throughout;
no floating point is used in any check).  It does not import the code that
produced the certificates.

WHAT IS PROVED (for the given V) if the verifier exits with status 0:

    No homogeneous convex polyhedron with V vertices has exactly one
    unstable equilibrium vertex.

LOGICAL STRUCTURE.  Let P be such a polyhedron, centroid at the origin, vertex
vectors r_1..r_V, Gram matrix R = (r_a . r_b): R is positive semidefinite and
R != 0.  For a matrix A write <A,R> = sum_ab A_ab R_ab.  S_{i->j} denotes
E_ii - (E_ij + E_ji)/2, so <S_{i->j}, R> = (r_i - r_j) . r_i.

 (1) [geometry, Bozoki-Domokos-Papp-Regos, arXiv:2401.17906]  If v is the only
     unstable vertex, every vertex i != v is "weakly shadowed" by an adjacent
     vertex p(i):  <S_{i->p(i)}, R> <= 0, and the map p has no cycles, so it is
     a spanning tree of the graph rooted at v.  Triangulating the faces of P
     (no new vertices) gives a simplicial 2-sphere on V vertices whose graph
     contains the tree.  The tetrahedra T = conv(v, F), F a triangle of that
     sphere not containing v, have disjoint interiors and fill P, so
     sum_T vol(T) * (sum_{k in T} r_k) = 4 vol(P) * centroid = 0.

 (2) [certificate]  Suppose mu in Z^V and, for every such T, y^T in Z^V with
     y^T >= 0 and y^T_v = 0 satisfy:
         M_T = sum_{i != v} y^T_i S_{i->p(i)} + (mu 1_T^T + 1_T mu^T)/2  is positive definite,
     where 1_T is the 0/1 indicator vector of the four vertices of T.
     Put w = sum_k mu_k r_k.  Dotting (1) with w gives
     sum_T vol(T) * <(mu 1_T^T + 1_T mu^T)/2, R> = 0, hence some T* with
     vol(T*) > 0 has <(mu 1_T*^T + 1_T* mu^T)/2, R> <= 0, and with the shadow
     inequalities <M_T*, R> <= 0.  But M_T* is positive definite and R is PSD
     and nonzero, so <M_T*, R> > 0.  Contradiction: no polyhedron has this
     (sphere triangulation, root v, tree p).

 (3) [completeness, checked here]  Every simplicial 2-sphere on V vertices is
     isomorphic (by a vertex bijection carrying triangles to triangles) to one
     of the listed ones; every vertex is equivalent under a listed
     automorphism to a listed root; for every listed (graph, root) EVERY
     rooted spanning tree has a certificate, and every certificate is valid.

 (2G) [general certificate, record kind "certG", used only where (2) fails]
     Two further sound ingredients may be used:
     (S1) steepest-neighbour tree: a configuration may be assigned the tree
          with p(i) in argmax_{k ~ i} r_k . r_i.  That neighbour shadows i too
          (some neighbour does, and the argmax has a larger inner product), and
          shadowing strictly increases |r| (|r_i|^2 <= r_i.r_p <= |r_i||r_p|,
          with equality only if r_p = r_i), so p is acyclic, i.e. one of the
          enumerated trees.  For that tree, r_i . r_k <= r_i . r_p(i) for every
          neighbour k != p(i) of i != v:  <A_ik, R> <= 0 with
          A_ik = sym(E_ik) - sym(E_i,p(i)).  Since EVERY configuration has such
          a tree, it suffices to refute, for each enumerated tree p, the
          configurations whose steepest-neighbour tree is p; so these
          inequalities may be added (with multipliers c >= 0).
     (S2) several cones: the identity in (1) holds for the cone decomposition
          from ANY vertex q (tetrahedra conv(q,F), F a triangle not containing
          q).  For a set Q of apexes with directions mu^(q) and integer
          symmetric W_q with sum_q W_q = 0, suppose for every q in Q and every
          T in the cone decomposition from q the matrix
             W_q + sum y_i (2 S_i) + sum c_ik (2 A_ik) + (mu^(q) 1_T^T + 1_T mu^(q)T)
          is positive definite.  Picking T_q as in (2) for each q and summing,
          the W_q cancel and we get a positive definite matrix whose pairing
          with R is <= 0: contradiction.  (Q = {v}, W = 0, no c: this is (2).)

 External fact used (not verified here): the number of combinatorial types of
 triangulations of the 2-sphere with V vertices (equivalently, of maximal
 planar simple graphs) is 1, 1, 2, 5, 14, 50 for V = 4, 5, 6, 7, 8, 9
 (OEIS A000109; e.g. Bowen-Fisk 1967, Brinkmann-McKay "plantri").
 Since the listed spheres are verified to be pairwise non-isomorphic AS GRAPHS
 (hence as complexes) and their number equals that count, the list is complete.

CHECKS performed:
 (i)   each listed (edges, faces) is a simplicial 2-sphere on {0..V-1}: simple
       graph with 3V-6 edges; 2V-4 distinct triangles, each of whose sides is an
       edge; each edge lies in exactly two triangles; the link of every vertex
       is a single cycle through all its neighbours (so the complex is a closed
       surface); the graph is connected and V - E + F = 2 (so it is the sphere);
       the graphs are pairwise non-isomorphic (degree sequence or exhaustive
       search over all V! bijections) and their number is A000109(V);
 (ii)  every vertex u is a listed root, or the listed permutation for u is an
       automorphism (edges onto edges, triangles onto triangles) mapping u to
       a listed root;
 (iii) for every (graph, root v) all rooted spanning trees are enumerated
       independently (all maps p: V\\{v} -> neighbours without directed
       cycles), their number is compared with Kirchhoff's matrix-tree
       determinant (exact), and each tree must have a certificate;
 (iv)  each certificate: integer mu, integer y^T >= 0 with y^T_v = 0, one y^T
       for each triangle not containing v (all of them, each once), and every
       integer matrix 2*M_T positive definite (fraction-free Gaussian
       elimination = all leading principal minors positive, Sylvester).
       For a general certificate (2G): distinct apexes q, integer symmetric
       W_q with sum_q W_q = 0 exactly (W = null for a single apex), for each
       apex all triangles avoiding it (each once), y >= 0 with y_v = 0, S1
       terms (i, k, c) only if the record's s1 flag is set, with i != v,
       {i,k} an edge of G, k != p(i), c >= 0, and every block matrix exactly
       positive definite.
"""

import gzip
import itertools
import json
import sys
import time

A000109 = {4: 1, 5: 1, 6: 2, 7: 5, 8: 14, 9: 50}


class Fail(Exception):
    pass


def need(cond, msg):
    if not cond:
        raise Fail(msg)


# ---------------------------------------------------------------- exact algebra

def is_positive_definite(N):
    """Symmetric integer matrix -> True iff positive definite.  Bareiss
    fraction-free elimination without pivoting; the k-th pivot equals the k-th
    leading principal minor, and all divisions are exact."""
    N = [list(row) for row in N]
    d, prev = len(N), 1
    for k in range(d):
        piv = N[k][k]
        if piv <= 0:
            return False
        for i in range(k + 1, d):
            for j in range(k + 1, d):
                num = N[i][j] * piv - N[i][k] * N[k][j]
                q, r = divmod(num, prev)
                assert r == 0
                N[i][j] = q
        prev = piv
    return True


def determinant(N):
    """Exact determinant of an integer matrix (Bareiss with row pivoting)."""
    N = [list(row) for row in N]
    d, prev, sign = len(N), 1, 1
    for k in range(d):
        if N[k][k] == 0:
            for i in range(k + 1, d):
                if N[i][k] != 0:
                    N[k], N[i] = N[i], N[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, d):
            for j in range(k + 1, d):
                N[i][j] = (N[i][j] * N[k][k] - N[i][k] * N[k][j]) // prev
        prev = N[k][k]
    return sign * N[d - 1][d - 1]


# ---------------------------------------------------------------- (i) spheres

def check_sphere(V, edges, faces):
    E = set()
    for e in edges:
        need(len(e) == 2, "edge must have 2 ends")
        a, b = e
        need(isinstance(a, int) and isinstance(b, int) and 0 <= a < V and 0 <= b < V and a != b,
             f"bad edge {e}")
        need(frozenset(e) not in E, f"repeated edge {e}")
        E.add(frozenset(e))
    need(len(E) == 3 * V - 6, f"{len(E)} edges, expected 3V-6 = {3 * V - 6}")
    F = set()
    for f in faces:
        need(len(f) == 3 and len(set(f)) == 3 and all(isinstance(x, int) and 0 <= x < V for x in f),
             f"bad triangle {f}")
        need(frozenset(f) not in F, f"repeated triangle {f}")
        F.add(frozenset(f))
        for a, b in itertools.combinations(f, 2):
            need(frozenset((a, b)) in E, f"side {a}{b} of triangle {f} is not an edge")
    need(len(F) == 2 * V - 4, f"{len(F)} triangles, expected 2V-4")
    for e in E:
        k = sum(1 for f in F if e <= f)
        need(k == 2, f"edge {sorted(e)} lies in {k} triangles")
    nbrs = {u: {x for e in E if u in e for x in e if x != u} for u in range(V)}
    for u in range(V):
        link = [tuple(f - {u}) for f in F if u in f]
        # link = graph on nbrs[u] with these edges; must be one cycle through all nbrs
        ldeg = {x: 0 for x in nbrs[u]}
        for a, b in link:
            need(a in ldeg and b in ldeg, "link vertex not a neighbour")
            ldeg[a] += 1
            ldeg[b] += 1
        need(len(nbrs[u]) >= 3 and all(d == 2 for d in ldeg.values()),
             f"link of vertex {u} is not a cycle")
        start = next(iter(nbrs[u]))
        seen, stack = {start}, [start]
        while stack:
            x = stack.pop()
            for a, b in link:
                for p, q in ((a, b), (b, a)):
                    if p == x and q not in seen:
                        seen.add(q)
                        stack.append(q)
        need(seen == nbrs[u], f"link of vertex {u} is not a single cycle")
    seen, stack = {0}, [0]
    while stack:
        x = stack.pop()
        for y in nbrs[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    need(len(seen) == V, "graph not connected")
    need(V - len(E) + len(F) == 2, "Euler characteristic is not 2")
    return E, F, nbrs


def isomorphic(V, E1, E2):
    d1 = sorted(sum(1 for e in E1 if u in e) for u in range(V))
    d2 = sorted(sum(1 for e in E2 if u in e) for u in range(V))
    if d1 != d2:
        return False
    for perm in itertools.permutations(range(V)):
        if all(frozenset(perm[x] for x in e) in E2 for e in E1):
            return True
    return False


# ---------------------------------------------------------------- (iii) trees

def rooted_spanning_trees(V, nbrs, v):
    """All maps p on V\\{v} with p(i) a neighbour of i and no directed cycle
    (i.e. iterating p from any vertex reaches v): exactly the spanning trees
    oriented towards v.  Returned as parent tuples with p[v] = -1."""
    others = [i for i in range(V) if i != v]
    out = []
    for choice in itertools.product(*[sorted(nbrs[i]) for i in others]):
        p = [-1] * V
        for i, j in zip(others, choice):
            p[i] = j
        ok = True
        for i in others:
            x, steps = i, 0
            while x != v and steps <= V:
                x = p[x]
                steps += 1
            if x != v:
                ok = False
                break
        if ok:
            out.append(tuple(p))
    return out


def kirchhoff(V, nbrs):
    L = [[(len(nbrs[a]) if a == b else (-1 if b in nbrs[a] else 0)) for b in range(V)]
         for a in range(V)]
    return determinant([row[1:] for row in L[1:]])


# ---------------------------------------------------------------- (iv) certificates

def check_certificate(V, F, nbrs, v, parent, mu, Y):
    need(isinstance(mu, list) and len(mu) == V and all(isinstance(t, int) for t in mu),
         "mu must be a list of V integers")
    need(parent[v] == -1, "parent of root must be -1")
    want = {f for f in F if v not in f}
    got = set()
    for face, y in Y:
        T = frozenset(face)
        need(T in want, f"{face} is not a triangle avoiding the root")
        need(T not in got, f"triangle {face} listed twice")
        got.add(T)
        need(isinstance(y, list) and len(y) == V and all(isinstance(t, int) for t in y),
             "y must be a list of V integers")
        need(all(t >= 0 for t in y), "negative multiplier")
        need(y[v] == 0, "multiplier of the root must be 0")
        one = [1 if (k == v or k in T) else 0 for k in range(V)]
        # N = 2 M_T = sum_i y_i (2E_ii - E_ij - E_ji) + mu one^T + one mu^T
        N = [[mu[a] * one[b] + one[a] * mu[b] for b in range(V)] for a in range(V)]
        for i in range(V):
            if i != v and y[i]:
                j = parent[i]
                need(j in nbrs[i], "parent not adjacent")
                N[i][i] += 2 * y[i]
                N[i][j] -= y[i]
                N[j][i] -= y[i]
        need(is_positive_definite(N), f"M_T not positive definite for T = {v}+{face}")
    need(got == want, "some triangle avoiding the root has no multipliers")


def int_list(x, n, what):
    need(isinstance(x, list) and len(x) == n and all(isinstance(t, int) for t in x),
         f"{what} must be a list of {n} integers")


def check_general_certificate(V, F, nbrs, v, parent, rec):
    """Record kind "certG" (see (2G) in the module docstring)."""
    need(parent[v] == -1, "parent of root must be -1")
    use_s1 = rec["s1"]
    need(use_s1 in (True, False), "s1 flag must be boolean")
    groups = rec["groups"]
    Q = [g["q"] for g in groups]
    need(len(Q) >= 1 and len(set(Q)) == len(Q) and all(isinstance(q, int) and 0 <= q < V for q in Q),
         "bad apex set Q")
    if len(groups) == 1:
        need(groups[0]["W"] is None, "single cone: W must be null (= 0)")
        Ws = [[[0] * V for _ in range(V)]]
    else:
        Ws = []
        for g in groups:
            W = g["W"]
            need(isinstance(W, list) and len(W) == V, "W must be V x V")
            for row in W:
                int_list(row, V, "row of W")
            need(all(W[a][b] == W[b][a] for a in range(V) for b in range(V)), "W not symmetric")
            Ws.append(W)
        need(all(sum(W[a][b] for W in Ws) == 0 for a in range(V) for b in range(V)),
             "the W_q do not sum to zero")
    for g, W in zip(groups, Ws):
        q, mu = g["q"], g["mu"]
        int_list(mu, V, "mu")
        want = {f for f in F if q not in f}
        got = set()
        for face, y, extra in g["blocks"]:
            T = frozenset(face)
            need(T in want, f"{face} is not a triangle avoiding the apex {q}")
            need(T not in got, f"triangle {face} listed twice for apex {q}")
            got.add(T)
            int_list(y, V, "y")
            need(all(t >= 0 for t in y), "negative multiplier")
            need(y[v] == 0, "multiplier of the root must be 0")
            need(use_s1 or not extra, "steepest-neighbour terms used without s1 flag")
            one = [1 if (k == q or k in T) else 0 for k in range(V)]
            N = [[W[a][b] + mu[a] * one[b] + one[a] * mu[b] for b in range(V)] for a in range(V)]
            for i in range(V):
                if i != v and y[i]:
                    j = parent[i]
                    N[i][i] += 2 * y[i]
                    N[i][j] -= y[i]
                    N[j][i] -= y[i]
            for i, k, c in extra:
                need(isinstance(c, int) and c >= 0, "negative steepest-neighbour multiplier")
                need(isinstance(i, int) and 0 <= i < V and i != v, "bad vertex in S1 term")
                j = parent[i]
                need(k in nbrs[i] and k != j,
                     f"S1 term ({i},{k}): not a neighbour of {i} other than its parent")
                N[i][k] += c
                N[k][i] += c
                N[i][j] -= c
                N[j][i] -= c
            need(is_positive_definite(N), f"matrix not positive definite (apex {q}, triangle {face})")
        need(got == want, f"some triangle avoiding apex {q} has no multipliers")


# ---------------------------------------------------------------- main

def main(V, paths):
    t0 = time.time()
    print(f"verify_joint.py: V = {V}, files {paths}")
    print("External fact used: number of triangulations of the 2-sphere with V vertices "
          f"(OEIS A000109) = {A000109[V]} for V = {V}.")
    graphs = {}
    verified = {}          # (gid, root) -> set of parent tuples with a valid certificate
    n_cert = n_dup = n_gen = 0
    for nfile, path in enumerate(paths):
      with gzip.open(path, "rt") as fh:
        header = json.loads(fh.readline())
        need(header.get("kind", "").startswith("header") and header.get("V") == V, "bad header")
        for line in fh:
            rec = json.loads(line)
            if rec["kind"] == "graph":
                need(nfile == 0, "graph records must be in the first file")
                gid = rec["gid"]
                need(gid not in graphs, "repeated graph id")
                E, F, nbrs = check_sphere(V, rec["edges"], rec["faces"])
                roots = rec["roots"]
                need(len(set(roots)) == len(roots) and all(0 <= r < V for r in roots), "bad roots")
                for u in range(V):
                    if u in roots:
                        continue
                    perm = rec["aut"].get(str(u))
                    need(perm is not None, f"graph {gid}: vertex {u} neither root nor mapped")
                    need(sorted(perm) == list(range(V)), "not a permutation")
                    need({frozenset(perm[x] for x in e) for e in E} == E,
                         f"graph {gid}: map for vertex {u} is not an automorphism")
                    need({frozenset(perm[x] for x in f) for f in F} == F,
                         f"graph {gid}: map for vertex {u} does not preserve triangles")
                    need(perm[u] in roots, f"graph {gid}: vertex {u} not mapped to a root")
                graphs[gid] = (E, F, nbrs, roots)
                for r in roots:
                    verified[(gid, r)] = set()
            elif rec["kind"] == "cert":
                gid, v = rec["gid"], rec["root"]
                need((gid, v) in verified, "certificate for unknown graph/root (graphs must come first)")
                E, F, nbrs, roots = graphs[gid]
                parent = tuple(rec["parent"])
                need(len(parent) == V, "bad parent array")
                check_certificate(V, F, nbrs, v, parent, rec["mu"], rec["Y"])
                if parent in verified[(gid, v)]:
                    n_dup += 1
                verified[(gid, v)].add(parent)
                n_cert += 1
                if n_cert % 50000 == 0:
                    print(f"   {n_cert} certificates verified ({time.time() - t0:.0f}s)", flush=True)
            elif rec["kind"] == "certG":
                gid, v = rec["gid"], rec["root"]
                need((gid, v) in verified, "certificate for unknown graph/root (graphs must come first)")
                E, F, nbrs, roots = graphs[gid]
                parent = tuple(rec["parent"])
                need(len(parent) == V, "bad parent array")
                check_general_certificate(V, F, nbrs, v, parent, rec)
                if parent in verified[(gid, v)]:
                    n_dup += 1
                verified[(gid, v)].add(parent)
                n_cert += 1
                n_gen += 1
            else:
                raise Fail("unknown record kind")
    t1 = time.time()
    print(f"(iv) {n_cert} certificates read ({n_gen} of the general kind (2G)), each verified "
          f"exactly (duplicates: {n_dup}); {t1 - t0:.1f}s")

    # (i) completeness of the list of spheres
    gids = sorted(graphs)
    for a, b in itertools.combinations(gids, 2):
        need(not isomorphic(V, graphs[a][0], graphs[b][0]), f"graphs {a} and {b} are isomorphic")
    need(len(gids) == A000109[V], f"{len(gids)} spheres listed, A000109 gives {A000109[V]}")
    print(f"(i)  {len(gids)} simplicial 2-spheres, pairwise non-isomorphic = A000109({V}): complete")
    print(f"(ii) every vertex is a root or mapped to one by a verified automorphism")

    # (iii) every rooted spanning tree has a verified certificate
    total = missing = 0
    for gid in gids:
        E, F, nbrs, roots = graphs[gid]
        tau = kirchhoff(V, nbrs)
        for v in roots:
            trees = rooted_spanning_trees(V, nbrs, v)
            need(len(trees) == tau, f"graph {gid} root {v}: {len(trees)} trees vs Kirchhoff {tau}")
            miss = [p for p in trees if p not in verified[(gid, v)]]
            total += len(trees)
            missing += len(miss)
            print(f"     graph {gid} root {v}: {len(trees)} rooted spanning trees "
                  f"(= Kirchhoff {tau}), {len(trees) - len(miss)} certified, {len(miss)} missing")
    print(f"(iii) {total} (graph, root, tree) cases, {missing} without certificate; "
          f"total time {time.time() - t0:.1f}s")
    need(missing == 0, f"{missing} trees have no certificate")
    print(f"ALL CHECKS PASSED: no homogeneous convex polyhedron with {V} vertices "
          f"is mono-unstable (given the geometric lemma (1) and A000109).")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    try:
        main(int(sys.argv[1]), sys.argv[2:])
    except Fail as e:
        print("VERIFICATION FAILED:", e)
        sys.exit(1)
    sys.exit(0)
