"""
Combinatorial case generation for "Algorithm 2" of

    S. Bozoki, G. Domokos, D. Papp, K. Regos,
    "The smallest mono-unstable, homogeneous convex polyhedron has at least
     7 vertices", arXiv:2401.17906.

For a given number of vertices V this module produces

  * all maximal planar graphs (triangulations of the sphere) on V vertices,
    up to isomorphism (1, 1, 2, 5, 14 graphs for V = 4, 5, 6, 7, 8),
  * the triangular faces of each graph (from its unique planar embedding),
  * the vertex orbits of the automorphism group (one representative root v
    per orbit is enough, exactly as in the paper),
  * all spanning trees rooted at v, encoded as parent arrays
    (parent[v] = -1, parent[i] = j(i) is a graph neighbour of i and following
    parents always ends at v), i.e. the acyclic "shadowing maps" of the paper,
  * the tetrahedral decomposition of the polyhedron w.r.t. v: the faces not
    containing v (each gives the tetrahedron conv(v, face)).

A "case" (= one quadratic system of the paper) is a tuple
(graph id, root v, parent array, face not containing v).

All vertex labels are 0-based.  Graph ids are 0-based indices into the list
returned by `triangulations(V)`, which is sorted deterministically by the
sorted degree sequence (descending), so ids are reproducible.

Enumeration of triangulations: breadth-first search over the edge-flip graph
(Wagner's theorem: any two triangulations with the same number of vertices are
connected by flips), deduplicated up to isomorphism.  For V <= 7 this is
cross-checked against a brute-force filter of all graphs with 3V-6 edges
(`brute_force_triangulations`).

Run as a script to print the summary table (counts of graphs, trees, systems).
"""

from __future__ import annotations

import itertools
import sys
from dataclasses import dataclass, field

import networkx as nx
import numpy as np
from networkx.algorithms.isomorphism import GraphMatcher


# ---------------------------------------------------------------------------
# Triangulations
# ---------------------------------------------------------------------------

def bipyramid(V: int) -> nx.Graph:
    """Bipyramid over a (V-2)-gon: a triangulation on V >= 5 vertices
    (for V = 4 we return K4, the tetrahedron)."""
    if V == 4:
        return nx.complete_graph(4)
    G = nx.cycle_graph(V - 2)
    for apex in (V - 2, V - 1):
        G.add_edges_from((apex, k) for k in range(V - 2))
    return G


def faces_of(G: nx.Graph) -> list[tuple[int, int, int]]:
    """Faces of a maximal planar graph from its planar embedding (unique up
    to reflection since the graph is 3-connected).  Returned as sorted
    triples, sorted lexicographically."""
    ok, emb = nx.check_planarity(G)
    if not ok:
        raise ValueError("graph is not planar")
    seen: set[tuple[int, int]] = set()
    faces = []
    for u, v in emb.edges():
        if (u, v) in seen:
            continue
        f = emb.traverse_face(u, v, mark_half_edges=seen)
        if len(f) != 3:
            raise ValueError("graph is not a triangulation (face of length %d)" % len(f))
        faces.append(tuple(sorted(f)))
    faces = sorted(set(faces))
    assert len(faces) == 2 * G.number_of_nodes() - 4
    return faces


def is_triangulation(G: nx.Graph) -> bool:
    V = G.number_of_nodes()
    return G.number_of_edges() == 3 * V - 6 and nx.check_planarity(G)[0]


def _flips(G: nx.Graph):
    """All graphs obtained from triangulation G by one edge flip."""
    faces = faces_of(G)
    by_edge: dict[tuple[int, int], list[int]] = {}
    for f in faces:
        for a, b in itertools.combinations(f, 2):
            third = [x for x in f if x != a and x != b][0]
            by_edge.setdefault((a, b), []).append(third)
    for (a, b), (c, d) in by_edge.items():   # faces abc and abd -> flip ab to cd
        if G.has_edge(c, d):
            continue                     # flip would create a multi-edge
        H = G.copy()
        H.remove_edge(a, b)
        H.add_edge(c, d)
        yield H


def _dedupe_add(H: nx.Graph, reps: list[nx.Graph]) -> bool:
    """Append H to reps unless isomorphic to a member; return True if new."""
    inv = sorted(d for _, d in H.degree())
    for R in reps:
        if sorted(d for _, d in R.degree()) == inv and nx.is_isomorphic(H, R):
            return False
    reps.append(H)
    return True


def _sort_key(G: nx.Graph):
    deg = sorted((d for _, d in G.degree()), reverse=True)
    # secondary key: sorted list of (deg u, deg v) over edges, to break ties
    ed = sorted(tuple(sorted((G.degree(u), G.degree(v)), reverse=True))
                for u, v in G.edges())
    return (deg, ed)


def triangulations(V: int) -> list[nx.Graph]:
    """All triangulations of the sphere with V >= 4 vertices up to
    isomorphism, via BFS on the flip graph.  Deterministic order."""
    start = nx.convert_node_labels_to_integers(bipyramid(V))
    reps = [start]
    queue = [start]
    while queue:
        G = queue.pop()
        for H in _flips(G):
            if _dedupe_add(H, reps):
                queue.append(H)
    reps.sort(key=_sort_key, reverse=True)
    # Relabel each representative canonically-ish: vertices sorted by
    # (degree descending, original label) so that labels are reproducible.
    out = []
    for G in reps:
        order = sorted(G.nodes(), key=lambda u: (-G.degree(u), u))
        mapping = {u: k for k, u in enumerate(order)}
        H = nx.relabel_nodes(G, mapping)
        H2 = nx.Graph()
        H2.add_nodes_from(range(V))
        H2.add_edges_from(sorted(tuple(sorted(e)) for e in H.edges()))
        out.append(H2)
    return out


def brute_force_triangulations(V: int) -> list[nx.Graph]:
    """Independent check (practical for V <= 7): filter all edge subsets of
    size 3V-6 of K_V (with min degree >= 3) for planarity, dedupe."""
    allE = list(itertools.combinations(range(V), 2))
    reps: list[nx.Graph] = []
    for E in itertools.combinations(allE, 3 * V - 6):
        deg = [0] * V
        for a, b in E:
            deg[a] += 1
            deg[b] += 1
        if min(deg) < 3:
            continue
        G = nx.Graph(E)
        if nx.check_planarity(G)[0]:
            _dedupe_add(G, reps)
    return reps


# ---------------------------------------------------------------------------
# Automorphisms, spanning trees, decompositions
# ---------------------------------------------------------------------------

def vertex_orbits(G: nx.Graph) -> list[list[int]]:
    """Orbits of Aut(G) on the vertices, each sorted, sorted by min element."""
    orbit_of = {u: {u} for u in G.nodes()}
    for phi in GraphMatcher(G, G).isomorphisms_iter():
        for u, w in phi.items():
            orbit_of[u].add(w)
    orbits = sorted({tuple(sorted(o)) for o in orbit_of.values()})
    return [list(o) for o in orbits]


def num_automorphisms(G: nx.Graph) -> int:
    return sum(1 for _ in GraphMatcher(G, G).isomorphisms_iter())


def rooted_spanning_trees(G: nx.Graph, root: int) -> list[tuple[int, ...]]:
    """All spanning trees of G as parent arrays rooted at `root`
    (parent[root] = -1).  Equivalently: all maps j from V\\{root} to
    neighbours without directed cycles.  Enumerated by depth-first choice of
    parents with early cycle rejection."""
    V = G.number_of_nodes()
    others = [u for u in range(V) if u != root]
    nbrs = {u: sorted(G.neighbors(u)) for u in others}
    parent = [-2] * V
    parent[root] = -1
    out: list[tuple[int, ...]] = []

    def creates_cycle(u: int) -> bool:
        # The parents assigned before u form a forest, so a new cycle must
        # pass through u: follow parents from u until we reach the root, an
        # unassigned vertex, or u itself.
        w = parent[u]
        while True:
            if w == u:
                return True
            if w == root or parent[w] == -2:
                return False
            w = parent[w]

    def rec(k: int):
        if k == len(others):
            out.append(tuple(parent))
            return
        u = others[k]
        for w in nbrs[u]:
            parent[u] = w
            if not creates_cycle(u):
                rec(k + 1)
        parent[u] = -2

    rec(0)
    # sanity: every vertex reaches the root
    for p in out:
        for u in others:
            w, s = u, 0
            while w != root:
                w = p[w]
                s += 1
                assert s <= V
    return out


def kirchhoff_count(G: nx.Graph) -> int:
    """Number of spanning trees by the matrix-tree theorem (exact)."""
    import sympy
    L = sympy.Matrix(nx.laplacian_matrix(G, nodelist=sorted(G.nodes())).toarray())
    return int(L[1:, 1:].det())


def tetrahedral_decomposition(faces, v: int) -> list[tuple[int, int, int]]:
    """Faces not containing v; each spans the tetrahedron conv(v, face)."""
    return [f for f in faces if v not in f]


# ---------------------------------------------------------------------------
# Case lists
# ---------------------------------------------------------------------------

@dataclass
class GraphData:
    gid: int
    V: int
    edges: list[tuple[int, int]]
    faces: list[tuple[int, int, int]]
    degrees: list[int]
    orbits: list[list[int]]
    n_aut: int
    roots: list[int] = field(default_factory=list)          # orbit representatives
    trees: dict[int, list[tuple[int, ...]]] = field(default_factory=dict)
    decomposition: dict[int, list[tuple[int, int, int]]] = field(default_factory=dict)

    @property
    def graph(self) -> nx.Graph:
        G = nx.Graph()
        G.add_nodes_from(range(self.V))
        G.add_edges_from(self.edges)
        return G


def graph_data(V: int, with_trees: bool = True) -> list[GraphData]:
    """All per-graph data needed to generate the systems for V vertices."""
    out = []
    for gid, G in enumerate(triangulations(V)):
        faces = faces_of(G)
        orbits = vertex_orbits(G)
        gd = GraphData(gid=gid, V=V, edges=sorted(G.edges()), faces=faces,
                       degrees=[G.degree(u) for u in range(V)], orbits=orbits,
                       n_aut=num_automorphisms(G))
        gd.roots = [o[0] for o in orbits]
        for v in gd.roots:
            gd.decomposition[v] = tetrahedral_decomposition(faces, v)
            if with_trees:
                gd.trees[v] = rooted_spanning_trees(G, v)
        out.append(gd)
    return out


def iter_cases(gds: list[GraphData]):
    """Yield (gid, v, parent, face) for every system of Algorithm 2."""
    for gd in gds:
        for v in gd.roots:
            for parent in gd.trees[v]:
                for face in gd.decomposition[v]:
                    yield gd.gid, v, parent, face


def count_systems(gds: list[GraphData]) -> int:
    return sum(len(gd.trees[v]) * len(gd.decomposition[v])
               for gd in gds for v in gd.roots)


# ---------------------------------------------------------------------------
# Self-checks and summary
# ---------------------------------------------------------------------------

KNOWN_NUM_TRIANGULATIONS = {4: 1, 5: 1, 6: 2, 7: 5, 8: 14, 9: 50}


def summary(V: int, check_brute_force: bool | None = None) -> list[GraphData]:
    gds = graph_data(V)
    assert len(gds) == KNOWN_NUM_TRIANGULATIONS[V], (V, len(gds))
    if check_brute_force is None:
        check_brute_force = V <= 7
    if check_brute_force:
        bf = brute_force_triangulations(V)
        assert len(bf) == len(gds)
        for G in bf:
            assert any(nx.is_isomorphic(G, gd.graph) for gd in gds)
    print(f"V = {V}: {len(gds)} triangulations"
          + ("  (brute-force check OK)" if check_brute_force else ""))
    for gd in gds:
        G = gd.graph
        tau = kirchhoff_count(G)
        for v in gd.roots:
            assert len(gd.trees[v]) == tau, (len(gd.trees[v]), tau)
        print(f"  graph {gd.gid}: degrees {gd.degrees}, |Aut| = {gd.n_aut}, "
              f"spanning trees = {tau} (matrix-tree theorem agrees)")
        print(f"     edges: {gd.edges}")
        print(f"     faces: {gd.faces}")
        for v in gd.roots:
            orb = [o for o in gd.orbits if v in o][0]
            print(f"     root {v} (deg {gd.degrees[v]}, orbit {orb}): "
                  f"{len(gd.decomposition[v])} tetrahedra -> "
                  f"{len(gd.trees[v]) * len(gd.decomposition[v])} systems")
    print(f"  total systems for V = {V}: {count_systems(gds)}")
    return gds


if __name__ == "__main__":
    Vs = [int(a) for a in sys.argv[1:]] or [5, 6, 7]
    for V in Vs:
        summary(V)
