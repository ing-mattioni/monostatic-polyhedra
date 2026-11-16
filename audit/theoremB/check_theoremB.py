#!/usr/bin/env python3
"""
check_theoremB.py -- independent referee re-check of the finite part of
Theorem B of paper/paper.md:

    every homogeneous convex polyhedron with at most 8 vertices has at least
    two unstable vertices.

Written from the mathematical statements of paper/paper.md (Sections 2 and 4)
and from the file-format description in the docstring of
lb_joint/verify_joint.py.  It does NOT import, call or copy any repository
code.  Standard library only; every check is exact (Python integers and
fractions.Fraction); no floating point is used anywhere in a decision.

Usage:
    python check_theoremB.py V CERTFILE [CERTFILE ...] [--procs N] [--quiet]
    (exit status 0 iff every check passes)

WHAT IS CHECKED (for one V):

 (A) Spheres.  Every listed (edges, faces) is a simplicial 2-sphere on
     {0..V-1}: distinct triangles, edge list = set of triangle sides, every
     edge in exactly two triangles, the link of every vertex is one cycle
     through all its neighbours (=> closed combinatorial surface), connected,
     V - E + F = 2 (=> the 2-sphere, by the classification of surfaces).

     Completeness is proved by an own GENERATION of all simplicial 2-spheres,
     starting from the boundary of the tetrahedron, by the operation
     "cone a fan disk":  pick a disk D in the sphere that is
       (3) one triangle abc, or
       (4) two triangles abc, acd sharing the edge ac, or
       (5) three consecutive triangles a b0 b1, a b1 b2, a b2 b3 around a
           vertex a of degree >= 4 (so a, b0, b1, b2, b3 are distinct),
     delete the triangles of D (and its interior edges) and join a new vertex
     x to the boundary cycle of D.
     Lemma.  Every simplicial 2-sphere K with n >= 5 vertices arises in this
     way from a simplicial 2-sphere with n-1 vertices.
     Proof.  2E = 3F and V - E + F = 2 give sum(deg) = 6n - 12 < 6n, so some
     vertex x has degree d in {3,4,5}.  Removing the open star of x leaves a
     triangulated disk D' whose boundary is the link cycle C = (a_1..a_d) of
     x.  An edge a_i a_j between non-consecutive vertices of C is an interior
     chord of D'; two such chords cannot have interleaved endpoints (they
     would cross in the disk D').  d = 3: the triangle a1a2a3 is not a face
     of K (else D' would be that single triangle and n = 4).  d = 4: at most
     one of a1a3, a2a4 is an edge; use the other.  d = 5: two disjoint
     diagonals of a pentagon always cross, and no vertex has more than two
     diagonals, so the chords present lie in a fan {a1a3, a1a4} after
     cyclic relabelling; the diagonals a2a4, a2a5 are not in that set, hence
     are not edges, and we use the fan from a2.  Filling C by the fan of new
     diagonals glues a disk to D' along C: a simplicial 2-sphere K0 with n-1
     vertices (new edges are non-edges, new triangles contain a new edge, so
     no simplex is duplicated).  K is obtained from K0 by coning that fan
     disk (of type (3), (4) or (5)).  QED.
     So the generated set is complete; every generated complex is also
     re-validated as a sphere.  Isomorphism classes are compared by an exact
     canonical form (lexicographically least relabelled triangle list over
     all relabellings that sort the vertices by an isomorphism invariant).
     The listed spheres must be pairwise non-isomorphic and their canonical
     forms must be exactly the generated set.  No literature count is used.

 (B) Roots.  The full automorphism group of every listed sphere is computed
     by brute force (all invariant-preserving bijections, test triangles);
     every orbit must contain a listed root.  The listed "aut" entries are
     also checked: perm is a bijection, maps triangles onto triangles, and
     perm[u] is a listed root.

 (C) Trees.  For every listed sphere all spanning trees are enumerated as
     (V-1)-subsets of the edge set that are acyclic (union-find); the count
     is compared with the matrix-tree theorem (exact rational determinant of
     the reduced Laplacian).  Each tree is oriented towards each listed root.
     Every such (sphere, root, tree) must carry a VALID certificate record.

 (D) Certificates (every record must be valid, not only one per tree):
     kind "cert" (Prop. 4.2): integer mu (length V); for each triangle F of G
       not containing the root v exactly one block [F, y] with integer y >= 0
       (length V) and y_v = 0; parent array = a spanning tree of G oriented
       to v (parent[v] = -1, every {i, parent[i]} an edge of G, acyclic);
       for T = F + {v} the integer matrix
          2 M_T = sum_i y_i (2 E_ii - E_ip(i) - E_p(i)i) + mu 1_T^T + 1_T mu^T
       must be positive definite.
     kind "certG" (Prop. 4.3): as above, but for every apex q of a set Q of
       distinct vertices an integer symmetric W_q, an integer mu^(q), and for
       each triangle F not containing q exactly one block [F, y, c] (y >= 0,
       y_v = 0 where v is the ROOT); sum_q W_q = 0 exactly; c must be empty
       and the s1 flag false (the paper's Prop. 4.3 has no S1 terms, so this
       checker accepts none); for T = F + {q} the matrix
          W_q + sum_i y_i (2 E_ii - E_ip(i) - E_p(i)i) + mu 1_T^T + 1_T mu^T
       must be positive definite.
     Positive definiteness: exact LDL^T (Schur-complement recursion) over
     fractions.Fraction: A > 0 iff A_00 > 0 and the Schur complement
     A_11 - a a^T / A_00 > 0.  This is NOT fraction-free Bareiss.
"""
import sys
import os
import gzip
import json
import time
import itertools
from fractions import Fraction
from multiprocessing import Pool


def isint(x):
    return isinstance(x, int) and not isinstance(x, bool)


def ntri(a, b, c):
    return tuple(sorted((a, b, c)))


def tri_edges(tris):
    E = set()
    for a, b, c in tris:
        for p, q in ((a, b), (b, c), (a, c)):
            E.add((min(p, q), max(p, q)))
    return E


def neighbours(V, E):
    nb = [set() for _ in range(V)]
    for a, b in E:
        nb[a].add(b)
        nb[b].add(a)
    return nb


# ---------------------------------------------------------------------------
# (A) simplicial 2-spheres
# ---------------------------------------------------------------------------

def sphere_errors(V, faces, edges=None):
    """Reasons why (faces[, edges]) is not a simplicial 2-sphere on 0..V-1."""
    err = []
    if not isint(V) or V < 4:
        return ["V < 4"]
    tris = []
    for t in faces:
        if (not isinstance(t, (list, tuple)) or len(t) != 3
                or not all(isint(x) and 0 <= x < V for x in t)
                or len(set(t)) != 3):
            return ["bad triangle %r" % (t,)]
        tris.append(ntri(*t))
    if len(set(tris)) != len(tris):
        err.append("duplicate triangle")
    tris = sorted(set(tris))
    E = tri_edges(tris)
    if edges is not None:
        given = []
        for e in edges:
            if (not isinstance(e, (list, tuple)) or len(e) != 2
                    or not all(isint(x) and 0 <= x < V for x in e)
                    or e[0] == e[1]):
                return ["bad edge %r" % (e,)]
            given.append((min(e), max(e)))
        if len(set(given)) != len(given):
            err.append("duplicate edge")
        if set(given) != E:
            err.append("edge list differs from triangle sides")
    cnt = {e: 0 for e in E}
    for a, b, c in tris:
        for p, q in ((a, b), (b, c), (a, c)):
            cnt[(min(p, q), max(p, q))] += 1
    if any(k != 2 for k in cnt.values()):
        err.append("some edge not in exactly two triangles")
    nb = neighbours(V, E)
    for v in range(V):
        if len(nb[v]) < 3:
            err.append("vertex %d has degree < 3" % v)
            continue
        ladj = {u: [] for u in nb[v]}
        nl = 0
        for t in tris:
            if v in t:
                a, b = [x for x in t if x != v]
                ladj[a].append(b)
                ladj[b].append(a)
                nl += 1
        if nl != len(nb[v]) or any(len(ladj[u]) != 2 for u in nb[v]):
            err.append("link of %d is not a cycle" % v)
            continue
        start = next(iter(nb[v]))
        seen = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for w in ladj[u]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        if seen != nb[v]:
            err.append("link of %d is not a single cycle" % v)
    seen = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for w in nb[u]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    if len(seen) != V:
        err.append("not connected")
    if V - len(E) + len(tris) != 2:
        err.append("Euler characteristic != 2")
    return err


def _groups(V, nb):
    """Vertices grouped by an isomorphism invariant, groups in sorted key order."""
    deg = [len(nb[v]) for v in range(V)]
    key = [(deg[v], tuple(sorted(deg[u] for u in nb[v]))) for v in range(V)]
    keys = sorted(set(key))
    return [[v for v in range(V) if key[v] == k] for k in keys]


def canon(V, tris):
    """Exact canonical form of a simplicial complex: least sorted triangle
    tuple over all relabellings listing the vertices in invariant order."""
    nb = neighbours(V, tri_edges(tris))
    groups = _groups(V, nb)
    best = None
    lab = [0] * V
    for combo in itertools.product(*(itertools.permutations(g) for g in groups)):
        i = 0
        for perm in combo:
            for x in perm:
                lab[x] = i
                i += 1
        code = tuple(sorted(ntri(lab[a], lab[b], lab[c]) for a, b, c in tris))
        if best is None or code < best:
            best = code
    return best


def automorphisms(V, tris):
    nb = neighbours(V, tri_edges(tris))
    groups = _groups(V, nb)
    triset = set(tris)
    auts = []
    for combo in itertools.product(*(itertools.permutations(g) for g in groups)):
        phi = [0] * V
        for g, perm in zip(groups, combo):
            for a, b in zip(g, perm):
                phi[a] = b
        if all(ntri(phi[a], phi[b], phi[c]) in triset for a, b, c in tris):
            auts.append(tuple(phi))
    return auts


def link_cycle(v, tris):
    ladj = {}
    for t in tris:
        if v in t:
            a, b = [x for x in t if x != v]
            ladj.setdefault(a, []).append(b)
            ladj.setdefault(b, []).append(a)
    start = min(ladj)
    cyc = [start]
    prev, cur = None, start
    while True:
        nxt = [w for w in ladj[cur] if w != prev]
        w = nxt[0] if prev is not None else ladj[cur][0]
        if w == start:
            break
        cyc.append(w)
        prev, cur = cur, w
    return cyc


def fan_extensions(V, tris):
    """All complexes on V+1 vertices obtained by coning a fan disk."""
    x = V
    T = set(tris)
    out = []
    for t in tris:                                   # type (3)
        a, b, c = t
        out.append((T - {t}) | {ntri(x, a, b), ntri(x, b, c), ntri(x, a, c)})
    for (a, c) in tri_edges(tris):                   # type (4)
        two = [t for t in tris if a in t and c in t]
        b = [y for y in two[0] if y not in (a, c)][0]
        d = [y for y in two[1] if y not in (a, c)][0]
        out.append((T - set(two)) | {ntri(x, a, b), ntri(x, b, c),
                                     ntri(x, c, d), ntri(x, d, a)})
    for a in range(V):                               # type (5)
        cyc = link_cycle(a, tris)
        k = len(cyc)
        if k < 4:
            continue
        for i in range(k):
            b0, b1, b2, b3 = (cyc[(i + j) % k] for j in range(4))
            rem = {ntri(a, b0, b1), ntri(a, b1, b2), ntri(a, b2, b3)}
            add = {ntri(x, a, b0), ntri(x, b0, b1), ntri(x, b1, b2),
                   ntri(x, b2, b3), ntri(x, b3, a)}
            out.append((T - rem) | add)
    return [sorted(s) for s in out]


def generate_spheres(Vmax):
    """dict n -> set of canonical forms of all simplicial 2-spheres, 4<=n<=Vmax."""
    res = {4: {canon(4, [ntri(0, 1, 2), ntri(0, 1, 3), ntri(0, 2, 3), ntri(1, 2, 3)])}}
    nbad = 0
    for n in range(5, Vmax + 1):
        s = set()
        for rep in res[n - 1]:
            for K in fan_extensions(n - 1, list(rep)):
                if sphere_errors(n, K):
                    nbad += 1
                    continue
                s.add(canon(n, K))
        res[n] = s
    return res, nbad


# ---------------------------------------------------------------------------
# (C) spanning trees and the matrix-tree theorem
# ---------------------------------------------------------------------------

def spanning_trees(V, E):
    E = sorted(E)
    trees = []
    for sub in itertools.combinations(E, V - 1):
        par = list(range(V))

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a
        ok = True
        for a, b in sub:
            ra, rb = find(a), find(b)
            if ra == rb:
                ok = False
                break
            par[ra] = rb
        if ok:
            trees.append(sub)
    return trees


def orient(V, tree_edges, root):
    nb = neighbours(V, tree_edges)
    parent = [None] * V
    parent[root] = -1
    stack = [root]
    while stack:
        u = stack.pop()
        for w in nb[u]:
            if parent[w] is None:
                parent[w] = u
                stack.append(w)
    return tuple(parent)


def det_fraction(A):
    A = [[Fraction(x) for x in row] for row in A]
    n = len(A)
    det = Fraction(1)
    for k in range(n):
        piv = next((i for i in range(k, n) if A[i][k] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != k:
            A[k], A[piv] = A[piv], A[k]
            det = -det
        det *= A[k][k]
        for i in range(k + 1, n):
            f = A[i][k] / A[k][k]
            if f:
                for j in range(k, n):
                    A[i][j] -= f * A[k][j]
    return det


def kirchhoff(V, E):
    L = [[0] * V for _ in range(V)]
    for a, b in E:
        L[a][a] += 1
        L[b][b] += 1
        L[a][b] -= 1
        L[b][a] -= 1
    return det_fraction([row[1:] for row in L[1:]])


# ---------------------------------------------------------------------------
# (D) certificates
# ---------------------------------------------------------------------------

def is_pd_ldl(M):
    """Exact LDL^T over the rationals (no pivoting, no Bareiss): True iff the
    symmetric matrix M is positive definite."""
    n = len(M)
    A = [[Fraction(M[i][j]) for j in range(i + 1)] for i in range(n)]  # lower
    for k in range(n):
        d = A[k][k]
        if d <= 0:
            return False
        col = [A[i][k] for i in range(k + 1, n)]
        for ii, i in enumerate(range(k + 1, n)):
            aik = col[ii]
            if aik:
                f = aik / d
                row = A[i]
                for jj, j in enumerate(range(k + 1, i + 1)):
                    ajk = col[jj]
                    if ajk:
                        row[j] -= f * ajk
    return True


_G = None      # gid -> dict
_V = None


def make_graph(V, r):
    """Data of one listed sphere (graph record r, already validated)."""
    tris = sorted(set(ntri(*t) for t in r["faces"]))
    return {"T": set(tris), "E": tri_edges(tris), "tris": tris,
            "roots": list(r.get("roots", [])),
            "aut": r.get("aut", {}),
            "avoid": [set(t for t in tris if q not in t) for q in range(V)]}


def _init(V, graphs):
    global _G, _V
    _V = V
    _G = graphs


def int_vector(x, V):
    return isinstance(x, list) and len(x) == V and all(isint(t) for t in x)


def tree_error(g, v, parent):
    V = _V
    if not int_vector(parent, V):
        return "parent not an integer vector of length V"
    if parent[v] != -1:
        return "parent[root] != -1"
    for i in range(V):
        if i == v:
            continue
        p = parent[i]
        if not (0 <= p < V) or p == i:
            return "bad parent value"
        if (min(i, p), max(i, p)) not in g["E"]:
            return "parent edge %d-%d not an edge of G" % (i, p)
    for i in range(V):
        x, steps = i, 0
        while x != v:
            x = parent[x]
            steps += 1
            if steps > V:
                return "parent map has a cycle"
    return None


def block_matrix(V, v, parent, y, mu, Tset, W):
    M = [list(r) for r in W] if W is not None else [[0] * V for _ in range(V)]
    for i in range(V):
        yi = y[i]
        if yi and i != v:
            j = parent[i]
            M[i][i] += 2 * yi
            M[i][j] -= yi
            M[j][i] -= yi
    for a in range(V):
        ma = mu[a]
        if ma:
            for b in Tset:
                M[a][b] += ma
                M[b][a] += ma
    return M


def check_blocks(g, v, q, parent, mu, W, blocks, general):
    """Blocks of one apex q.  Returns (error or None, number of PD tests)."""
    V = _V
    if not int_vector(mu, V):
        return "mu not an integer vector of length V", 0
    if not isinstance(blocks, list):
        return "blocks not a list", 0
    seen = set()
    todo = []
    for blk in blocks:
        if general:
            if not (isinstance(blk, list) and len(blk) == 3):
                return "bad block shape", 0
            F, y, c = blk
            if c != []:
                return "S1 terms present (not part of Prop. 4.3)", 0
        else:
            if not (isinstance(blk, list) and len(blk) == 2):
                return "bad block shape", 0
            F, y = blk
        if not (isinstance(F, list) and len(F) == 3 and all(isint(t) for t in F)):
            return "bad triangle in block", 0
        Ft = ntri(*F)
        if Ft not in g["T"]:
            return "block triangle %r is not a triangle of G" % (F,), 0
        if q in Ft:
            return "block triangle %r contains the apex" % (F,), 0
        if Ft in seen:
            return "triangle %r has two blocks" % (F,), 0
        seen.add(Ft)
        if not int_vector(y, V):
            return "y not an integer vector of length V", 0
        if any(t < 0 for t in y):
            return "negative multiplier y", 0
        if y[v] != 0:
            return "y_root != 0", 0
        todo.append((Ft, y))
    if seen != g["avoid"][q]:
        return "blocks do not cover all triangles avoiding apex %d" % q, 0
    npd = 0
    for Ft, y in todo:
        M = block_matrix(V, v, parent, y, mu, set(Ft) | {q}, W)
        npd += 1
        if not is_pd_ldl(M):
            return "block %r (apex %d) not positive definite" % (Ft, q), npd
    return None, npd


def check_record(r):
    """Returns (kind, key or None, error or None, n_pd_tests)."""
    V = _V
    kind = r.get("kind") if isinstance(r, dict) else None
    if kind not in ("cert", "certG"):
        return (kind, None, "unexpected record kind %r in certificate section" % kind, 0)
    gid = r.get("gid")
    if not isint(gid) or gid not in _G:
        return (kind, None, "unknown gid", 0)
    g = _G[gid]
    v = r.get("root")
    if not isint(v) or v not in g["roots"]:
        return (kind, None, "root not a listed root", 0)
    parent = r.get("parent")
    e = tree_error(g, v, parent)
    if e:
        return (kind, None, e, 0)
    key = (gid, v, tuple(parent))
    if kind == "cert":
        err, n = check_blocks(g, v, v, parent, r.get("mu"), None, r.get("Y"), False)
        return (kind, key, err, n)
    # certG
    if r.get("s1", False) is not False:
        return (kind, key, "s1 flag set (S1 is not part of Prop. 4.3)", 0)
    groups = r.get("groups")
    if not isinstance(groups, list) or not groups:
        return (kind, key, "no groups", 0)
    qs = [grp.get("q") if isinstance(grp, dict) else None for grp in groups]
    if not all(isint(q) and 0 <= q < V for q in qs) or len(set(qs)) != len(qs):
        return (kind, key, "apexes not distinct vertices", 0)
    Ws = []
    for grp in groups:
        W = grp.get("W")
        if W is None and len(groups) == 1:
            W = [[0] * V for _ in range(V)]
        if not (isinstance(W, list) and len(W) == V and all(int_vector(row, V) for row in W)):
            return (kind, key, "W not an integer VxV matrix", 0)
        if any(W[a][b] != W[b][a] for a in range(V) for b in range(a)):
            return (kind, key, "W not symmetric", 0)
        Ws.append(W)
    if any(sum(W[a][b] for W in Ws) != 0 for a in range(V) for b in range(V)):
        return (kind, key, "sum of W_q is not 0", 0)
    total = 0
    for grp, W in zip(groups, Ws):
        err, n = check_blocks(g, v, grp["q"], parent, grp.get("mu"), W, grp.get("blocks"), True)
        total += n
        if err:
            return (kind, key, err, total)
    return (kind, key, None, total)


def _work(lines):
    out = []
    for ln in lines:
        try:
            r = json.loads(ln)
        except ValueError:
            out.append(("?", None, "unparsable line", 0))
            continue
        if isinstance(r, dict) and r.get("kind") in ("header", "header_certG"):
            out.append(("header", None, None if r.get("V") == _V else "header V mismatch", 0))
            continue
        try:
            out.append(check_record(r))
        except Exception as ex:   # malformed record -> rejected, never accepted
            out.append((r.get("kind") if isinstance(r, dict) else "?", None,
                        "exception %s: %s" % (type(ex).__name__, ex), 0))
    return out


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    procs = os.cpu_count() or 1
    for a in argv[1:]:
        if a.startswith("--procs="):
            procs = int(a.split("=", 1)[1])
    V = int(args[0])
    files = args[1:]
    t0 = time.time()
    fails = []

    def fail(msg):
        fails.append(msg)
        print("FAIL:", msg, flush=True)

    import hashlib
    with open(os.path.abspath(__file__), "rb") as fh:
        me = hashlib.sha256(fh.read()).hexdigest()
    print("=== independent Theorem B check, V = %d, files = %s" % (V, files), flush=True)
    for fn in files:
        with open(fn, "rb") as fh:
            print("    sha256(%s) = %s" % (os.path.basename(fn), hashlib.sha256(fh.read()).hexdigest()))
    print("    sha256(check_theoremB.py) = %s, workers = %d" % (me, procs), flush=True)

    # ---- stream the files: header/graph records first, then certificates
    def lines():
        for fn in files:
            with gzip.open(fn, "rt") as f:
                for ln in f:
                    if ln.strip():
                        yield ln
    it = lines()
    graphs_raw = {}
    first_cert = None
    for ln in it:
        r = json.loads(ln)
        k = r.get("kind")
        if k in ("header", "header_certG"):
            if r.get("V") != V:
                fail("header V mismatch")
            continue
        if k == "graph":
            gid = r.get("gid")
            if gid in graphs_raw:
                fail("duplicate gid %r" % gid)
            graphs_raw[gid] = r
            continue
        first_cert = ln
        break

    # ---- (A) spheres
    gen, nbad = generate_spheres(V)
    print("(A) generated sphere classes by fan-disk coning: %s ; invalid intermediate complexes: %d"
          % (", ".join("V=%d:%d" % (n, len(gen[n])) for n in sorted(gen)), nbad), flush=True)
    if nbad:
        fail("generation produced a non-sphere (bug)")
    graphs = {}
    canon_of = {}
    for gid, r in sorted(graphs_raw.items()):
        errs = sphere_errors(V, r.get("faces", []), r.get("edges"))
        if errs:
            fail("gid %r is not a simplicial 2-sphere: %s" % (gid, errs))
            continue
        tris = sorted(set(ntri(*t) for t in r["faces"]))
        c = canon(V, tris)
        if c in canon_of.values():
            fail("gid %r isomorphic to an earlier listed sphere" % gid)
        canon_of[gid] = c
        graphs[gid] = make_graph(V, r)
    listed = set(canon_of.values())
    print("(A) listed spheres: %d valid, pairwise non-isomorphic: %s; generated classes: %d; "
          "listed == generated: %s"
          % (len(canon_of), len(listed) == len(canon_of), len(gen[V]), listed == gen[V]), flush=True)
    if listed != gen[V]:
        fail("listed spheres are not exactly the generated classes (missing %d, extra %d)"
             % (len(gen[V] - listed), len(listed - gen[V])))

    # ---- (B) roots
    for gid, g in sorted(graphs.items()):
        auts = automorphisms(V, g["tris"])
        orbit = list(range(V))

        def find(a):
            while orbit[a] != a:
                a = orbit[a]
            return a
        for phi in auts:
            for x in range(V):
                ra, rb = find(x), find(phi[x])
                if ra != rb:
                    orbit[ra] = rb
        orbits = {}
        for x in range(V):
            orbits.setdefault(find(x), []).append(x)
        orbits = sorted(orbits.values())
        roots = g["roots"]
        if not all(isint(x) and 0 <= x < V for x in roots) or len(set(roots)) != len(roots):
            fail("gid %d: bad root list" % gid)
        uncovered = [o for o in orbits if not set(o) & set(roots)]
        if uncovered:
            fail("gid %d: orbits without listed root: %s" % (gid, uncovered))
        redundant = [o for o in orbits if len(set(o) & set(roots)) > 1]
        # stated automorphisms
        autset = set(auts)
        for u in range(V):
            if u in roots:
                continue
            perm = g["aut"].get(str(u))
            if (not isinstance(perm, list) or sorted(perm) != list(range(V))
                    or tuple(perm) not in autset or perm[u] not in roots):
                fail("gid %d: stated automorphism for vertex %d invalid" % (gid, u))
        print("(B) gid %2d: |Aut| = %3d, orbits = %d, listed roots = %s, orbits w/o root = %d, "
              "orbits with >1 root = %d, stated auts ok"
              % (gid, len(auts), len(orbits), roots, len(uncovered), len(redundant)), flush=True)

    # ---- (C) trees
    required = set()
    per_pair = {}
    for gid, g in sorted(graphs.items()):
        trees = spanning_trees(V, g["E"])
        kh = kirchhoff(V, g["E"])
        if kh != len(trees):
            fail("gid %d: %d enumerated trees but Kirchhoff gives %s" % (gid, len(trees), kh))
        for v in g["roots"]:
            ks = set((gid, v, orient(V, t, v)) for t in trees)
            if len(ks) != len(trees):
                fail("gid %d root %d: orientation not injective (bug)" % (gid, v))
            per_pair[(gid, v)] = len(ks)
            required |= ks
        print("(C) gid %2d: spanning trees = %d (Kirchhoff %s)" % (gid, len(trees), kh), flush=True)
    print("(C) (sphere, root) pairs = %d, rooted trees required = %d" % (len(per_pair), len(required)),
          flush=True)

    # ---- (D) certificates
    covered = set()
    stats = {"cert": [0, 0], "certG": [0, 0]}
    npd = 0
    nbadrec = 0
    shown = 0

    def chunks():
        buf = [first_cert] if first_cert is not None else []
        for ln in it:
            buf.append(ln)
            if len(buf) >= 400:
                yield buf
                buf = []
        if buf:
            yield buf
    with Pool(procs, initializer=_init, initargs=(V, graphs)) as pool:
        for res in pool.imap_unordered(_work, chunks()):
            for kind, key, err, n in res:
                npd += n
                if kind == "header":
                    if err:
                        fail(err)
                    continue
                if kind not in stats:
                    stats[kind] = [0, 0]
                if err is None:
                    stats[kind][0] += 1
                    if key not in required:
                        fail("certificate for a case that is not a required tree: %r" % (key,))
                    covered.add(key)
                else:
                    stats[kind][1] += 1
                    nbadrec += 1
                    if shown < 25:
                        print("FAIL: invalid %s record %r: %s" % (kind, key, err), flush=True)
                        shown += 1
    if nbadrec:
        fail("%d invalid certificate records" % nbadrec)
    missing = required - covered
    for k in sorted(missing)[:25]:
        print("FAIL: no valid certificate for (gid, root, parent) = %r" % (k,), flush=True)
    if missing:
        fail("%d required trees without valid certificate" % len(missing))
    for kind, (ok, bad) in sorted(stats.items()):
        print("(D) %-6s records: valid %d, invalid %d" % (kind, ok, bad), flush=True)
    print("(D) exact positive-definiteness tests performed: %d" % npd, flush=True)
    print("(D) required trees: %d, covered by a valid certificate: %d, missing: %d"
          % (len(required), len(required & covered), len(missing)), flush=True)
    print("elapsed %.1f s" % (time.time() - t0), flush=True)
    if fails:
        print("RESULT V=%d: FAIL (%d problems)" % (V, len(fails)), flush=True)
        return 1
    print("RESULT V=%d: PASS" % V, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
