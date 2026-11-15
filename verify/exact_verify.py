#!/usr/bin/env python3
r"""
exact_verify.py -- an exact-arithmetic certificate for the static equilibria
of a homogeneous convex polyhedron with rational vertex coordinates.

Usage
-----
    python exact_verify.py polyhedron.json [--report report.json] [--verbose]
    python exact_verify.py --from-npy points.npy --denom 1000000 --out poly.json

Input format (JSON)::

    {"name": "...", "vertices": [["x", "y", "z"], ...]}

Each coordinate is a string accepted by ``fractions.Fraction`` ("3/7", "-12",
"0.125", "1e-3") or a JSON integer.  JSON decimal literals such as ``0.1`` are
also accepted and are read *exactly as written* (0.1 -> 1/10), never through a
binary float.

Mathematical setting
--------------------
Let P = conv{p_1, ..., p_n} be a full-dimensional convex polyhedron in R^3 and
let c be the centre of mass of the *solid* P with uniform density.  A boundary
point q of P is an equilibrium iff the plane through q perpendicular to q - c
supports P.  In the nondegenerate situation the equilibria are

  * stable   (one per face F):    the orthogonal projection of c onto the plane
                                  of F lies in the relative interior of F;
  * saddle   (one per edge [a,b]): the foot q of the perpendicular from c to the
                                  line ab lies in the open segment (a,b), and the
                                  plane through q perpendicular to q - c meets P
                                  exactly in the edge [a,b];
  * unstable (one per vertex v):  the plane through v perpendicular to v - c
                                  meets P only in v, i.e. (w - v).(v - c) < 0
                                  for every other vertex w.

For a nondegenerate polyhedron Poincare-Hopf gives S - H + U = 2.

How decisions are made
----------------------
Every decision taken by this program is the *sign of an exact rational number*
computed with ``fractions.Fraction`` and Python integers.  Floating-point
numbers appear only in the human-readable output (approximate values and
"margins" printed for information); they never influence a decision.

Each candidate (face, edge, vertex) is described by a finite list of rational
"margins" that must ALL be strictly positive for it to be a (nondegenerate)
equilibrium:

    all margins > 0                     ->  equilibrium
    at least one margin < 0             ->  not an equilibrium
    otherwise (all >= 0, some == 0)     ->  DEGENERATE

A polyhedron is *certified* when no candidate is degenerate, all internal
self-consistency checks pass and S - H + U = 2.

Algorithm outline
-----------------
1. Convex hull, by brute force over triples of input points.  The plane
   spanned by three non-collinear points is a facet plane iff all points lie
   (weakly) on one side of it.  All input points on such a plane are collected
   and their planar convex hull is computed with an exact Jarvis march, so
   faces with any number of coplanar vertices (e.g. a 17-gon base) are found as
   single polygons.  The hull computation is done on integer coordinates
   obtained by multiplying all coordinates by their common denominator L; this
   is a homothety with ratio L > 0 and does not change any orientation sign.
2. Volume and centroid: every face is fan-triangulated, and the signed volumes
   and centroids of the tetrahedra (O, a, b, c) are summed, where O is the
   average of the vertices (an interior point).
3. Classification of every face, edge and vertex as described above, together
   with self-consistency checks based on elementary lemmas of convex geometry
   (see the comments in the classification functions).

The program uses only the Python standard library; ``numpy`` is imported
lazily, and only by the ``--from-npy`` conversion path.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import numbers
import os
import sys
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

__all__ = [
    "InputError", "SelfCheckError",
    "parse_coordinate", "parse_points", "load_polyhedron_json",
    "convex_hull", "volume_and_centroid", "verify",
    "format_report", "rationalize", "write_polyhedron_json", "npy_to_json",
    "main",
]

STABLE = "stable"
SADDLE = "saddle"
UNSTABLE = "unstable"
NONEQ = "non-equilibrium"
DEGENERATE = "degenerate"


# =============================================================================
# Errors and the self-check primitive
# =============================================================================

class InputError(ValueError):
    """The input does not describe a valid full-dimensional set of points."""


class SelfCheckError(RuntimeError):
    """An internal consistency check failed.

    This can only happen through a bug in this program (or an input situation
    that it was not designed for); a result is never reported in that case.
    """


def _check(condition: bool, message: str) -> None:
    """Raise SelfCheckError unless `condition` holds.

    Used instead of ``assert`` so that the checks cannot be disabled by
    running Python with -O.
    """
    if not condition:
        raise SelfCheckError(message)


# =============================================================================
# Exact 3-vector arithmetic.  Entries are Python ints or Fractions; no floats.
# =============================================================================

def vsub(p, q):
    return (p[0] - q[0], p[1] - q[1], p[2] - q[2])


def vadd(p, q):
    return (p[0] + q[0], p[1] + q[1], p[2] + q[2])


def vscale(s, p):
    return (s * p[0], s * p[1], s * p[2])


def dot(p, q):
    return p[0] * q[0] + p[1] * q[1] + p[2] * q[2]


def cross(p, q):
    return (p[1] * q[2] - p[2] * q[1],
            p[2] * q[0] - p[0] * q[2],
            p[0] * q[1] - p[1] * q[0])


def det3(a, b, c):
    """Determinant of the 3x3 matrix with rows a, b, c (= a . (b x c))."""
    return dot(a, cross(b, c))


def sign(x) -> int:
    return (x > 0) - (x < 0)


def _norm_float(p) -> float:
    """Euclidean norm as a float -- for DISPLAY purposes only."""
    return math.sqrt(float(dot(p, p)))


# =============================================================================
# Input parsing
# =============================================================================

Point = Tuple[Fraction, Fraction, Fraction]


def parse_coordinate(value: Any) -> Fraction:
    """Convert one coordinate to an exact Fraction.

    Accepted: int, Fraction, or a string accepted by ``Fraction`` (e.g. "3/7",
    "-12", "0.125", "1e-3").  Binary floats are rejected on purpose: their
    exact value is rarely what was meant (use :func:`rationalize`).
    """
    if isinstance(value, bool):
        raise InputError(f"invalid coordinate {value!r}")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, numbers.Integral):          # int, numpy integers, ...
        return Fraction(int(value))
    if isinstance(value, str):
        try:
            return Fraction(value.strip())
        except (ValueError, ZeroDivisionError) as exc:
            raise InputError(f"cannot parse coordinate {value!r} as a rational number") from exc
    if isinstance(value, float):
        raise InputError(f"floating-point coordinate {value!r}: give coordinates as strings "
                         f"or integers (see rationalize())")
    raise InputError(f"invalid coordinate {value!r} of type {type(value).__name__}")


def parse_points(raw: Iterable[Sequence[Any]]) -> List[Point]:
    """Convert a sequence of 3-sequences of coordinates into exact points."""
    points: List[Point] = []
    for i, row in enumerate(raw):
        if isinstance(row, (str, bytes)) or not hasattr(row, "__len__") or len(row) != 3:
            raise InputError(f"point {i} is not a list of three coordinates: {row!r}")
        points.append(tuple(parse_coordinate(x) for x in row))  # type: ignore[misc]
    return points


def load_polyhedron_json(path: str) -> Tuple[str, List[Point]]:
    """Read a polyhedron file; return (name, list of exact points)."""
    try:
        with open(path, encoding="utf-8") as fh:
            # parse_float=Fraction: a decimal literal such as 0.1 is read exactly
            # as the rational number it denotes (1/10), not as a binary float.
            data = json.load(fh, parse_float=Fraction)
    except (OSError, json.JSONDecodeError) as exc:
        raise InputError(f"cannot read {path}: {exc}") from exc
    default_name = os.path.splitext(os.path.basename(path))[0]
    if isinstance(data, list):
        return default_name, parse_points(data)
    if not isinstance(data, dict) or "vertices" not in data:
        raise InputError(f"{path}: expected an object with a 'vertices' list")
    return str(data.get("name", default_name)), parse_points(data["vertices"])


def input_fingerprint(points: Sequence[Point]) -> str:
    """SHA-256 of a canonical text form of the input coordinates.

    Lets a printed certificate be tied unambiguously to its input.  The
    canonical form is the compact JSON list of lists of reduced fractions,
    e.g. [["0","0","1"],["-3/7","1","0"]].
    """
    canon = json.dumps([[str(x) for x in p] for p in points], separators=(",", ":"))
    return hashlib.sha256(canon.encode("ascii")).hexdigest()


# =============================================================================
# Exact convex hull
# =============================================================================

@dataclass
class Face:
    vertices: Tuple[int, ...]         # input indices, counter-clockwise seen from outside
    normal: Tuple[int, int, int]      # outward normal, primitive integer vector
    points_on_plane: Tuple[int, ...]  # all input points on the face plane (incl. redundant ones)


@dataclass
class Hull:
    num_points: int
    common_denominator: int                       # L: integer frame is L * (input coordinates)
    vertices: List[int]                           # input indices of the extreme points
    faces: List[Face]
    edges: List[Tuple[int, int]]                  # (a, b) with a < b
    edge_faces: Dict[Tuple[int, int], Tuple[int, int]]  # edge -> its two incident faces
    neighbours: Dict[int, List[int]]              # vertex -> adjacent vertices
    redundant: Dict[int, str]                     # input index -> "interior" | "face" | "edge"


def _common_denominator(points: Sequence[Point]) -> int:
    L = 1
    for p in points:
        for x in p:
            L = math.lcm(L, x.denominator)
    return L


def _primitive(v: Tuple[int, int, int]) -> Tuple[int, int, int]:
    """Divide an integer vector by the gcd of its entries (direction unchanged)."""
    g = math.gcd(math.gcd(v[0], v[1]), v[2])
    return (v[0] // g, v[1] // g, v[2] // g)


def _find_spanning_tetrahedron(Z: Sequence[Tuple[int, int, int]]) -> Tuple[int, int, int, int]:
    """Return indices of four affinely independent points, or raise InputError."""
    i0, i1 = 0, 1                       # distinct: duplicates were rejected before
    e1 = vsub(Z[i1], Z[i0])
    i2 = next((k for k in range(len(Z)) if cross(e1, vsub(Z[k], Z[i0])) != (0, 0, 0)), None)
    if i2 is None:
        raise InputError("degenerate input: all points are collinear")
    nrm = cross(e1, vsub(Z[i2], Z[i0]))
    i3 = next((k for k in range(len(Z)) if dot(nrm, vsub(Z[k], Z[i0])) != 0), None)
    if i3 is None:
        raise InputError("degenerate input: all points are coplanar")
    return i0, i1, i2, i3


def _supporting_planes(Z: Sequence[Tuple[int, int, int]]) -> List[Tuple[List[int], Tuple[int, int, int]]]:
    """All facet planes of conv(Z), found by brute force over triples.

    A plane spanned by three non-collinear points of Z is the plane of a facet
    iff Z does not have points strictly on both sides of it.  Conversely every
    facet plane contains three non-collinear points of Z, so every facet is
    found.  Returns a list of (indices of the points on the
    plane, outward normal).  Z must not be coplanar.
    """
    n = len(Z)
    covered = set()   # triples inside an already-found facet plane: nothing new to learn
    planes = []
    for i in range(n - 2):
        for j in range(i + 1, n - 1):
            e1 = vsub(Z[j], Z[i])
            for k in range(j + 1, n):
                if (i, j, k) in covered:
                    continue
                nrm = cross(e1, vsub(Z[k], Z[i]))
                if nrm == (0, 0, 0):
                    continue                       # collinear triple spans no plane
                offset = dot(nrm, Z[i])
                above = below = False
                on_plane = []
                for m in range(n):
                    s = dot(nrm, Z[m]) - offset
                    if s > 0:
                        above = True
                    elif s < 0:
                        below = True
                    else:
                        on_plane.append(m)
                    if above and below:
                        break                      # points on both sides: not a facet
                if above and below:
                    continue
                # All points satisfy s <= 0 (then nrm points outward) or s >= 0.
                outward = nrm if below else vscale(-1, nrm)
                planes.append((on_plane, _primitive(outward)))
                covered.update(itertools.combinations(on_plane, 3))
    return planes


def _orient(Z, a: int, b: int, x: int, normal) -> int:
    """((Z[b]-Z[a]) x (Z[x]-Z[a])) . normal  for three points of one face plane.

    Positive iff x lies to the LEFT of the directed line a -> b when the plane
    is viewed from outside (from the side the outward normal points to).
    """
    return det3(normal, vsub(Z[b], Z[a]), vsub(Z[x], Z[a]))


def _face_polygon(Z, on_plane: Sequence[int], normal) -> Tuple[List[int], List[int]]:
    """Exact Jarvis march on coplanar points.

    Returns (polygon, on_edges): the extreme points of conv(on_plane) in
    counter-clockwise order seen from outside (no three consecutive collinear),
    and the points lying in the relative interior of a polygon edge.
    """
    def dist2(a, b):
        d = vsub(Z[a], Z[b])
        return dot(d, d)

    # The lexicographically smallest point of a finite set is an extreme point
    # of its convex hull, hence a valid starting vertex.
    start = min(on_plane, key=lambda i: Z[i])
    polygon = [start]
    on_edges: List[int] = []
    current = start
    while True:
        # Choose the next vertex: no point may lie strictly to its right; among
        # collinear candidates take the farthest one.
        cand = None
        for x in on_plane:
            if x == current:
                continue
            if cand is None:
                cand = x
                continue
            o = _orient(Z, current, cand, x, normal)
            if o < 0 or (o == 0 and dist2(current, x) > dist2(current, cand)):
                cand = x
        # Certify the choice independently of how it was found: every other
        # point is strictly left of current->cand, or strictly inside the
        # segment [current, cand].  Hence current->cand is a hull edge.
        d = vsub(Z[cand], Z[current])
        dd = dot(d, d)
        for x in on_plane:
            if x in (current, cand):
                continue
            o = _orient(Z, current, cand, x, normal)
            _check(o >= 0, "Jarvis march: a point lies outside a face polygon edge")
            if o == 0:
                s = dot(vsub(Z[x], Z[current]), d)
                _check(0 < s < dd, "Jarvis march: collinear point outside the edge segment")
                on_edges.append(x)
        if cand == start:
            break
        _check(len(polygon) < len(on_plane), "Jarvis march did not close")
        polygon.append(cand)
        current = cand
    _check(len(polygon) >= 3, "face polygon with fewer than 3 vertices")
    return polygon, on_edges


def convex_hull(points: Sequence[Point]) -> Hull:
    """Exact convex hull of a finite set of rational points in R^3.

    Raises InputError for fewer than 4 points, duplicate points, or a set that
    is contained in a plane.
    """
    n = len(points)
    if n < 4:
        raise InputError(f"need at least 4 points, got {n}")

    first_index: Dict[Point, int] = {}
    duplicates = []
    for i, p in enumerate(points):
        if p in first_index:
            duplicates.append((first_index[p], i))
        else:
            first_index[p] = i
    if duplicates:
        raise InputError("duplicate points (index pairs): "
                         + ", ".join(f"{a}={b}" for a, b in duplicates))

    # Integer frame: Z = L * points.  A homothety of ratio L > 0 preserves all
    # orientation signs used below, so the combinatorics of the hull are those
    # of the original points.
    L = _common_denominator(points)
    Z = [tuple(int(x * L) for x in p) for p in points]
    _find_spanning_tetrahedron(Z)

    planes = _supporting_planes(Z)
    _check(len({frozenset(on) for on, _ in planes}) == len(planes), "a facet plane was found twice")
    faces: List[Face] = []
    on_some_edge = set()
    for on_plane, normal in planes:
        polygon, on_edges = _face_polygon(Z, on_plane, normal)
        faces.append(Face(tuple(polygon), normal, tuple(on_plane)))
        on_some_edge.update(on_edges)

    vertex_set = set()
    for f in faces:
        vertex_set.update(f.vertices)
    vertices = sorted(vertex_set)

    # A point of a facet plane is a vertex of P iff it is a vertex of that face.
    for f in faces:
        _check(set(f.points_on_plane) & vertex_set == set(f.vertices),
               f"vertex/face incidence inconsistent for face {f.vertices}")

    on_boundary = set()
    for f in faces:
        on_boundary.update(f.points_on_plane)
    redundant = {}
    for i in range(n):
        if i not in vertex_set:
            redundant[i] = ("edge" if i in on_some_edge else
                            "face" if i in on_boundary else "interior")

    # Edges: every directed polygon edge a->b must occur exactly once, and its
    # reverse b->a must occur in the neighbouring face (closed 2-manifold).
    directed: Dict[Tuple[int, int], int] = {}
    for fi, f in enumerate(faces):
        k = len(f.vertices)
        for t in range(k):
            a, b = f.vertices[t], f.vertices[(t + 1) % k]
            _check((a, b) not in directed, f"directed edge {a}->{b} occurs in two faces")
            directed[(a, b)] = fi
    edge_faces: Dict[Tuple[int, int], Tuple[int, int]] = {}
    for (a, b), fi in directed.items():
        _check((b, a) in directed, f"edge {a}-{b} belongs to only one face")
        if a < b:
            edge_faces[(a, b)] = (fi, directed[(b, a)])
    edges = sorted(edge_faces)

    neighbours: Dict[int, List[int]] = {v: [] for v in vertices}
    for a, b in edges:
        neighbours[a].append(b)
        neighbours[b].append(a)
    for v in vertices:
        neighbours[v].sort()
        _check(len(neighbours[v]) >= 3, f"vertex {v} has degree {len(neighbours[v])} < 3")

    V, E, F = len(vertices), len(edges), len(faces)
    _check(V - E + F == 2, f"Euler characteristic V-E+F = {V - E + F} != 2")

    # The (twice) vector area of each face polygon is a positive multiple of
    # its outward normal (orientation check), and the vector areas of all
    # faces sum to zero (the faces close up to a surface).
    total = (0, 0, 0)
    for f in faces:
        a = Z[f.vertices[0]]
        area = (0, 0, 0)
        for t in range(1, len(f.vertices) - 1):
            area = vadd(area, cross(vsub(Z[f.vertices[t]], a), vsub(Z[f.vertices[t + 1]], a)))
        _check(cross(area, f.normal) == (0, 0, 0) and dot(area, f.normal) > 0,
               f"face {f.vertices} is not oriented counter-clockwise seen from outside")
        total = vadd(total, area)
    _check(total == (0, 0, 0), "face vector areas do not sum to zero")

    return Hull(n, L, vertices, faces, edges, edge_faces, neighbours, redundant)


# =============================================================================
# Volume and centroid of the solid
# =============================================================================

def volume_and_centroid(points: Sequence[Point], hull: Hull) -> Tuple[Fraction, Point]:
    """Exact volume and centre of mass of the solid conv(points).

    Fan-triangulate every face (a, b, c ordered counter-clockwise seen from
    outside) and decompose P into the tetrahedra (O, a, b, c), O = average of
    the vertices (an interior point).  det(a-O, b-O, c-O) = 6 vol(O,a,b,c) > 0
    for each of them (checked), and the centroid of a tetrahedron is the
    average of its four vertices.
    """
    vs = [points[v] for v in hull.vertices]
    nv = len(vs)
    O = tuple(sum(p[k] for p in vs) / nv for k in range(3))
    six_vol = Fraction(0)
    moment = (Fraction(0), Fraction(0), Fraction(0))
    for f in hull.faces:
        a = points[f.vertices[0]]
        for t in range(1, len(f.vertices) - 1):
            b, c = points[f.vertices[t]], points[f.vertices[t + 1]]
            d = det3(vsub(a, O), vsub(b, O), vsub(c, O))
            _check(d > 0, "non-positive tetrahedron in the volume decomposition")
            six_vol += d
            moment = vadd(moment, vscale(d, vadd(vadd(O, a), vadd(b, c))))
    volume = six_vol / 6
    centroid = vscale(1 / (4 * six_vol), moment)
    return volume, centroid  # type: ignore[return-value]


# =============================================================================
# Classification of faces, edges and vertices
# =============================================================================

def _decide(margins: Iterable[Fraction]) -> str:
    """'eq' if all margins > 0, 'none' if some margin < 0, else 'degenerate'."""
    margins = list(margins)
    if all(m > 0 for m in margins):
        return "eq"
    if any(m < 0 for m in margins):
        return "none"
    return "degenerate"


@dataclass
class FaceResult:
    index: int
    vertices: Tuple[int, ...]
    normal: Tuple[int, int, int]
    status: str                        # STABLE | NONEQ | DEGENERATE
    foot: Point                        # orthogonal projection of c onto the face plane
    edge_margins: List[Fraction]       # >0: foot strictly on the inner side of polygon edge k
    witness_edge: Optional[Tuple[int, int]]   # an edge whose line separates the foot (NONEQ)
    margin_float: float                # min signed distance foot -> edge lines (display only)


@dataclass
class EdgeResult:
    vertices: Tuple[int, int]
    faces: Tuple[int, int]
    status: str                        # SADDLE | NONEQ | DEGENERATE
    t: Fraction                        # foot q = a + t (b - a)
    foot: Point
    max_value: Fraction                # max over w != a, b of (w - q).(q - c); must be < 0
    witness: int                       # a vertex attaining max_value
    reason: str
    height_float: float                # max_value / |q - c| (display only)
    margin_float: float                # min(t|ab|, (1-t)|ab|, -height)    (display only)


@dataclass
class VertexResult:
    index: int
    status: str                        # UNSTABLE | NONEQ | DEGENERATE
    max_value: Fraction                # max over w != v of (w - v).(v - c)
    witness: int                       # a vertex attaining it
    adjacent_max_value: Fraction       # the same maximum over the neighbours of v only
    adjacent_witness: int
    height_float: float                # max_value / |v - c| (display only)

    @property
    def margin_float(self) -> float:   # display only; > 0 iff unstable
        return -self.height_float


def _argmax(values: Dict[int, Fraction]) -> int:
    """Index with the largest value; ties broken by the smallest index."""
    return max(values, key=lambda w: (values[w], -w))


def classify_faces(points: Sequence[Point], hull: Hull, c: Point) -> List[FaceResult]:
    """Stable equilibria.

    With outward normal n and polygon edges a -> b (counter-clockwise seen from
    outside), the foot q = c - ((c-a).n / n.n) n of c lies strictly on the
    inner side of edge ab iff det(n, b - a, q - a) > 0.  Since q - c is
    parallel to n, det(n, b-a, q-a) = det(n, b-a, c-a): the projection need
    not even be formed (it is formed anyway, for the report and as a check).
    """
    results = []
    for fi, f in enumerate(hull.faces):
        n = f.normal
        nn = dot(n, n)
        p0 = points[f.vertices[0]]
        h = dot(n, vsub(c, p0))
        _check(h < 0, f"centroid not strictly inside the plane of face {fi}")
        foot = vsub(c, vscale(h / nn, n))
        _check(dot(n, vsub(foot, p0)) == 0, "foot point is not on the face plane")

        k = len(f.vertices)
        margins: List[Fraction] = []
        dists: List[float] = []
        for t in range(k):
            a, b = points[f.vertices[t]], points[f.vertices[(t + 1) % k]]
            m = det3(n, vsub(b, a), vsub(c, a))
            _check(m == det3(n, vsub(b, a), vsub(foot, a)), "projection invariance violated")
            margins.append(m)
            dists.append(float(m) / (math.sqrt(nn) * _norm_float(vsub(b, a))))
        decision = _decide(margins)
        status = {"eq": STABLE, "none": NONEQ, "degenerate": DEGENERATE}[decision]
        witness_edge = None
        if status == NONEQ:
            t = next(t for t, m in enumerate(margins) if m < 0)
            witness_edge = (f.vertices[t], f.vertices[(t + 1) % k])
        results.append(FaceResult(fi, f.vertices, n, status, foot, margins, witness_edge, min(dists)))
    return results


def classify_edges(points: Sequence[Point], hull: Hull, c: Point) -> List[EdgeResult]:
    """Saddle equilibria.

    Edge [a, b], d = b - a.  The foot of c on the line ab is q = a + t d with
    t = (c - a).d / d.d, and u = q - c is orthogonal to d.  [a, b] carries a
    saddle iff 0 < t < 1 and the plane through q with normal u supports P and
    meets it exactly in [a, b], i.e. (w - q).u < 0 for every vertex w != a, b
    (a and b themselves lie on that plane because (a - q) and (b - q) are
    parallel to d).  The margins are therefore t, 1 - t and -(w - q).u.

    Self-check (a lemma): P lies in the dihedral wedge bounded by the planes of
    the two faces incident to [a, b], so the sign of max_w (w - q).u over all
    vertices equals its sign over the vertices of these two faces.
    """
    results = []
    for (a, b) in hull.edges:
        pa, pb = points[a], points[b]
        d = vsub(pb, pa)
        t = dot(vsub(c, pa), d) / dot(d, d)
        q = vadd(pa, vscale(t, d))
        u = vsub(q, c)
        _check(dot(u, d) == 0, "foot of perpendicular is not orthogonal")
        _check(u != (0, 0, 0), "centroid lies on an edge line")
        _check(dot(vsub(pa, q), u) == 0 and dot(vsub(pb, q), u) == 0, "edge not in its own plane")

        values = {w: dot(vsub(points[w], q), u) for w in hull.vertices if w not in (a, b)}
        witness = _argmax(values)
        max_value = values[witness]

        f1, f2 = hull.edge_faces[(a, b)]
        local = set(hull.faces[f1].vertices) | set(hull.faces[f2].vertices)
        local_max = max(values[w] for w in local - {a, b})
        _check(sign(local_max) == sign(max_value),
               f"edge {a}-{b}: local and global support tests disagree")

        decision = _decide([t, 1 - t] + [-x for x in values.values()])
        status = {"eq": SADDLE, "none": NONEQ, "degenerate": DEGENERATE}[decision]
        # Human-readable reason, listing a strictly violated condition first.
        if t < 0 or t > 1:
            reason = f"foot of c outside the edge (t = {float(t):.6g})"
        elif max_value > 0:
            reason = f"plane through the foot cuts P (vertex {witness} strictly beyond it)"
        elif t == 0 or t == 1:
            reason = "foot of c is exactly an endpoint of the edge"
        elif max_value == 0:
            reason = f"plane through the foot also contains vertex {witness}"
        else:
            reason = "saddle"
        height = float(max_value) / _norm_float(u)
        length = _norm_float(d)
        margin = min(float(t) * length, float(1 - t) * length, -height)
        results.append(EdgeResult((a, b), (f1, f2), status, t, q, max_value, witness, reason,
                                  height, margin))
    return results


def classify_vertices(points: Sequence[Point], hull: Hull, c: Point) -> List[VertexResult]:
    """Unstable equilibria.

    Vertex v is unstable iff (w - v).(v - c) < 0 for every vertex w != v.  A
    vertex w with (w - v).(v - c) > 0 "shadows" v.

    Self-check (a lemma): P is contained in the cone v + cone{w - v : w
    adjacent to v}, so the maximum of (w - v).(v - c) over all vertices has
    the same sign as its maximum over the neighbours of v.  In particular a
    non-unstable vertex is always shadowed by an ADJACENT vertex.
    """
    results = []
    for v in hull.vertices:
        p = points[v]
        u = vsub(p, c)
        values = {w: dot(vsub(points[w], p), u) for w in hull.vertices if w != v}
        witness = _argmax(values)
        adjacent = {w: values[w] for w in hull.neighbours[v]}
        adjacent_witness = _argmax(adjacent)
        _check(sign(values[witness]) == sign(adjacent[adjacent_witness]),
               f"vertex {v}: adjacent-shadowing lemma violated")
        decision = _decide(-x for x in values.values())
        status = {"eq": UNSTABLE, "none": NONEQ, "degenerate": DEGENERATE}[decision]
        results.append(VertexResult(v, status, values[witness], witness,
                                    adjacent[adjacent_witness], adjacent_witness,
                                    float(values[witness]) / _norm_float(u)))
    return results


# =============================================================================
# Top-level verification
# =============================================================================

@dataclass
class VerificationResult:
    name: str
    points: List[Point]
    fingerprint: str
    hull: Hull
    volume: Fraction
    centroid: Point
    faces: List[FaceResult]
    edges: List[EdgeResult]
    vertices: List[VertexResult]
    degeneracies: List[str]

    @property
    def stable_faces(self) -> List[FaceResult]:
        return [f for f in self.faces if f.status == STABLE]

    @property
    def saddle_edges(self) -> List[EdgeResult]:
        return [e for e in self.edges if e.status == SADDLE]

    @property
    def unstable_vertices(self) -> List[int]:
        return [v.index for v in self.vertices if v.status == UNSTABLE]

    @property
    def counts(self) -> Tuple[int, int, int]:
        """(S, H, U)."""
        return len(self.stable_faces), len(self.saddle_edges), len(self.unstable_vertices)

    @property
    def poincare_hopf_ok(self) -> bool:
        S, H, U = self.counts
        return S - H + U == 2

    @property
    def has_degeneracy(self) -> bool:
        return bool(self.degeneracies)

    @property
    def certified(self) -> bool:
        """No degenerate candidate and S - H + U = 2 (all self-checks passed)."""
        return not self.has_degeneracy and self.poincare_hopf_ok

    @property
    def mono_unstable(self) -> bool:
        return self.certified and self.counts[2] == 1

    @property
    def mono_stable(self) -> bool:
        return self.certified and self.counts[0] == 1

    def smallest_margin(self) -> Tuple[float, str]:
        """Informational robustness indicator (floating point, NOT used in any decision).

        Every element carries a signed float margin in length units, of the
        same sign as its exact decision (> 0 equilibrium, < 0 not): for a face
        the smallest signed distance of the foot of c to the face's edge lines;
        for an edge the smallest of the distances of the foot to the endpoints
        (signed) and the depth of the other vertices below the plane through
        the foot; for a vertex the depth of the other vertices below the plane
        through it.  Returns the smallest |margin| and the element attaining it.
        """
        items = ([(abs(f.margin_float), f"face {f.index} {list(f.vertices)}") for f in self.faces]
                 + [(abs(e.margin_float), f"edge {e.vertices[0]}-{e.vertices[1]}") for e in self.edges]
                 + [(abs(v.margin_float), f"vertex {v.index}") for v in self.vertices])
        return min(items, key=lambda item: item[0])

    def to_json(self) -> Dict[str, Any]:
        """Machine-readable report.  Exact rationals are given as strings 'p/q'."""
        fr = str
        vec = lambda p: [fr(x) for x in p]            # noqa: E731
        vecf = lambda p: [float(x) for x in p]        # noqa: E731
        S, H, U = self.counts
        hull = self.hull
        return {
            "name": self.name,
            "input_sha256": self.fingerprint,
            "num_input_points": hull.num_points,
            "common_denominator": hull.common_denominator,
            "V": len(hull.vertices), "E": len(hull.edges), "F": len(hull.faces),
            "euler_characteristic": len(hull.vertices) - len(hull.edges) + len(hull.faces),
            "vertices": hull.vertices,
            "redundant_points": [{"index": i, "location": loc} for i, loc in sorted(hull.redundant.items())],
            "volume": fr(self.volume), "volume_float": float(self.volume),
            "centroid": vec(self.centroid), "centroid_float": vecf(self.centroid),
            "faces": [{
                "index": f.index, "vertices": list(f.vertices), "outward_normal": list(f.normal),
                "status": f.status, "foot": vec(f.foot), "foot_float": vecf(f.foot),
                "edge_margins": [fr(m) for m in f.edge_margins],
                "witness_edge": list(f.witness_edge) if f.witness_edge else None,
                "min_distance_to_edge_lines_float": f.margin_float,
            } for f in self.faces],
            "edges": [{
                "vertices": list(e.vertices), "faces": list(e.faces), "status": e.status,
                "t": fr(e.t), "t_float": float(e.t), "foot": vec(e.foot),
                "max_value": fr(e.max_value), "witness": e.witness,
                "max_height_float": e.height_float, "margin_float": e.margin_float,
                "reason": e.reason,
            } for e in self.edges],
            "vertex_classification": [{
                "index": v.index, "status": v.status,
                "max_value": fr(v.max_value), "witness": v.witness,
                "adjacent_max_value": fr(v.adjacent_max_value),
                "adjacent_witness": v.adjacent_witness,
                "max_height_float": v.height_float,
            } for v in self.vertices],
            "stable_faces": [list(f.vertices) for f in self.stable_faces],
            "saddle_edges": [list(e.vertices) for e in self.saddle_edges],
            "unstable_vertices": self.unstable_vertices,
            "S": S, "H": H, "U": U,
            "poincare_hopf_ok": self.poincare_hopf_ok,
            "smallest_margin_float": {"value": self.smallest_margin()[0],
                                      "element": self.smallest_margin()[1]},
            "degeneracies": self.degeneracies,
            "has_degeneracy": self.has_degeneracy,
            "certified": self.certified,
            "mono_unstable": self.mono_unstable,
            "mono_stable": self.mono_stable,
        }


def verify(points: Sequence[Sequence[Any]], name: str = "") -> VerificationResult:
    """Run the complete exact verification on a list of points.

    `points` may contain Fractions, ints or strings (see parse_coordinate).
    Raises InputError for invalid input and SelfCheckError if an internal
    consistency check fails.
    """
    pts = parse_points(points)
    hull = convex_hull(pts)
    volume, c = volume_and_centroid(pts, hull)
    faces = classify_faces(pts, hull, c)
    edges = classify_edges(pts, hull, c)
    verts = classify_vertices(pts, hull, c)

    degeneracies = []
    for f in faces:
        if f.status == DEGENERATE:
            degeneracies.append(f"face {f.index} {list(f.vertices)}: projection of c lies on "
                                f"the boundary of the face")
    for e in edges:
        if e.status == DEGENERATE:
            degeneracies.append(f"edge {e.vertices[0]}-{e.vertices[1]}: {e.reason}")
    for v in verts:
        if v.status == DEGENERATE:
            degeneracies.append(f"vertex {v.index}: the plane through it perpendicular to v - c "
                                f"also contains vertex {v.witness}")
    return VerificationResult(name, pts, input_fingerprint(pts), hull, volume, c,
                              faces, edges, verts, degeneracies)


# =============================================================================
# Human-readable report
# =============================================================================

def _frac_str(x: Fraction, limit: int = 0) -> str:
    s = str(x)
    if limit and len(s) > limit:
        return (f"<{len(str(abs(x.numerator)))}-digit / {len(str(x.denominator))}-digit "
                f"fraction, see JSON report>")
    return s


def format_report(res: VerificationResult, verbose: bool = False) -> str:
    hull = res.hull
    V, E, F = len(hull.vertices), len(hull.edges), len(hull.faces)
    S, H, U = res.counts
    out: List[str] = []
    w = out.append
    bar = "=" * 78
    w(bar)
    w(f"EXACT EQUILIBRIUM CERTIFICATE: {res.name}")
    w(bar)
    w(f"input points           : {hull.num_points}")
    w(f"input sha256           : {res.fingerprint}")
    w(f"common denominator     : {hull.common_denominator}")
    w(f"convex hull            : V = {V}, E = {E}, F = {F};  V - E + F = {V - E + F}  [OK]")
    sizes: Dict[int, int] = {}
    for f in hull.faces:
        sizes[len(f.vertices)] = sizes.get(len(f.vertices), 0) + 1
    w("face sizes             : " + ", ".join(f"{cnt} x {k}-gon" for k, cnt in sorted(sizes.items())))
    if hull.redundant:
        w(f"REDUNDANT input points : {len(hull.redundant)} (not vertices of the hull; ignored)")
        for i, loc in sorted(hull.redundant.items()):
            where = {"interior": "in the interior of P", "face": "in the relative interior of a face",
                     "edge": "in the relative interior of an edge"}[loc]
            w(f"    point {i}: {where}")
    else:
        w("redundant input points : none (every input point is a vertex)")
    w(f"volume                 : {_frac_str(res.volume, 200)}")
    w(f"                       ~ {float(res.volume):.15g}")
    w("centroid (exact)       :")
    for axis, x in zip("xyz", res.centroid):
        w(f"    {axis} = {_frac_str(x, 400)}")
    w("centroid (float)       : (" + ", ".join(f"{float(x):.15g}" for x in res.centroid) + ")")
    w("")

    if verbose:
        w(f"FACES (F = {F})   margin = min signed distance of the foot of c to the edge lines")
        for f in res.faces:
            extra = f"  outside edge {f.witness_edge}" if f.witness_edge else ""
            w(f"  face {f.index:3d} {str(list(f.vertices)):40s} {f.status:16s} "
              f"margin {f.margin_float:+.3e}{extra}")
        w("")
        w(f"EDGES (E = {E})   t = position of the foot on the edge; "
          f"h = max height of another vertex above the plane")
        for e in res.edges:
            w(f"  edge {e.vertices[0]:3d}-{e.vertices[1]:<3d} {e.status:16s} t = {float(e.t):+.6f}  "
              f"h = {e.height_float:+.3e}  ({e.reason})")
        w("")

    w(f"VERTICES (V = {V})   h = max over w of (w - v).(v - c)/|v - c|  "
      f"(height of the highest other vertex above the plane through v perpendicular to v - c)")
    for v in res.vertices:
        if v.status == UNSTABLE:
            w(f"  vertex {v.index:3d}  UNSTABLE         h = {v.height_float:+.3e} "
              f"(closest: vertex {v.witness})")
        else:
            adj = "adjacent" if v.witness in hull.neighbours[v.index] else "non-adjacent"
            kind = "touched" if v.status == DEGENERATE else "shadowed"
            w(f"  vertex {v.index:3d}  {v.status:16s} h = {v.height_float:+.3e}  "
              f"{kind} by vertex {v.witness} ({adj}); adjacent witness: vertex {v.adjacent_witness}")
    w("")

    w("SUMMARY")
    w(f"  stable faces      S = {S}: " + "; ".join(
        f"face {f.index} {list(f.vertices)}" for f in res.stable_faces))
    w(f"  saddle edges      H = {H}: " + ", ".join(
        f"{e.vertices[0]}-{e.vertices[1]}" for e in res.saddle_edges))
    w(f"  unstable vertices U = {U}: " + ", ".join(str(i) for i in res.unstable_vertices))
    w(f"  (S, H, U) = ({S}, {H}, {U});  S - H + U = {S - H + U}  "
      f"[Poincare-Hopf {'OK' if res.poincare_hopf_ok else 'FAILED'}]")
    if res.degeneracies:
        w(f"  DEGENERACIES FOUND: {len(res.degeneracies)}")
        for d in res.degeneracies:
            w(f"    - {d}")
    else:
        w("  degeneracies: none (every decision was a strict inequality)")
    value, element = res.smallest_margin()
    w(f"  smallest decisive margin (float, informational only): {value:.3e} at {element}")
    w(f"  CERTIFIED (exact, nondegenerate): {res.certified}")
    w(f"  MONO-UNSTABLE: {res.mono_unstable}")
    w(f"  MONO-STABLE:   {res.mono_stable}")
    if not res.certified:
        w("  (verdicts are True only for a certified, nondegenerate polyhedron)")
    w(bar)
    return "\n".join(out)


# =============================================================================
# Float -> rational conversion
# =============================================================================

def rationalize(points: Iterable[Sequence[float]], denom: int) -> List[Point]:
    """Round floating-point points to rationals with common denominator `denom`.

    Each coordinate x becomes round(x * denom) / denom, where x is taken at its
    exact binary value and the rounding is round-half-to-even.  Consequently
    exact zeros stay zero, and x and -x are mapped to exact negatives of each
    other, so exact coplanarities with coordinate planes and mirror symmetries
    through them are preserved.
    """
    denom = int(denom)
    if denom <= 0:
        raise ValueError("denominator must be a positive integer")
    out: List[Point] = []
    for i, p in enumerate(points):
        row = [float(x) for x in p]
        if len(row) != 3:
            raise ValueError(f"point {i} does not have 3 coordinates")
        if not all(math.isfinite(x) for x in row):
            raise ValueError(f"point {i} has a non-finite coordinate")
        out.append(tuple(Fraction(round(Fraction(x) * denom), denom) for x in row))  # type: ignore[misc]
    return out


def write_polyhedron_json(path: str, points: Sequence[Sequence[Any]], name: str = "") -> None:
    """Write points in the input format of this program."""
    pts = parse_points(points)
    rows = ",\n".join("  " + json.dumps([str(x) for x in p]) for p in pts)   # one vertex per line
    with open(path, "w", encoding="utf-8") as fh:
        fh.write('{\n "name": ' + json.dumps(name) + ',\n "vertices": [\n' + rows + "\n ]\n}\n")


def npy_to_json(npy_path: str, denom: int, out_path: str, name: Optional[str] = None) -> List[Point]:
    """Convert a (V, 3) float array stored with numpy.save into the JSON input format."""
    import numpy as np  # lazy: only this conversion path needs numpy

    arr = np.load(npy_path)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise InputError(f"{npy_path}: expected an array of shape (V, 3), got {arr.shape}")
    pts = rationalize(arr.tolist(), denom)
    if name is None:
        name = f"{os.path.splitext(os.path.basename(npy_path))[0]} (denominator {denom})"
    write_polyhedron_json(out_path, pts, name)
    return pts


# =============================================================================
# Command-line interface
# =============================================================================

EXIT_CERTIFIED = 0        # nondegenerate, Poincare-Hopf OK
EXIT_NOT_CERTIFIED = 1    # degeneracy found or Poincare-Hopf failed
EXIT_INPUT_ERROR = 2      # unreadable / invalid input
EXIT_SELF_CHECK = 3       # internal consistency check failed (bug)


def _parse_denominator(text: str) -> int:
    try:
        value = Fraction(text)
    except (ValueError, ZeroDivisionError):
        raise argparse.ArgumentTypeError(f"invalid denominator {text!r}")
    if value.denominator != 1 or value <= 0:
        raise argparse.ArgumentTypeError(f"denominator must be a positive integer, got {text!r}")
    return int(value)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Exact-arithmetic certificate of the equilibria (stable faces, saddle edges, "
                    "unstable vertices) of a homogeneous convex polyhedron with rational vertices.",
        epilog="Exit codes: 0 certified (nondegenerate and S-H+U=2), 1 not certified, "
               "2 input error, 3 internal self-check failure.")
    parser.add_argument("input", nargs="?", help="polyhedron JSON file")
    parser.add_argument("--report", metavar="FILE", help="also write a machine-readable JSON report")
    parser.add_argument("--verbose", action="store_true", help="list every face and edge")
    parser.add_argument("--from-npy", metavar="NPY",
                        help="convert a (V,3) float array to the JSON input format and exit")
    parser.add_argument("--denom", type=_parse_denominator, help="common denominator for --from-npy")
    parser.add_argument("--out", metavar="FILE", help="output JSON file for --from-npy")
    parser.add_argument("--name", help="polyhedron name to store with --from-npy")
    args = parser.parse_args(argv)

    if args.from_npy:
        if args.denom is None or args.out is None:
            parser.error("--from-npy requires --denom and --out")
        try:
            pts = npy_to_json(args.from_npy, args.denom, args.out, args.name)
        except (OSError, ValueError) as exc:
            print(f"INPUT ERROR: {exc}", file=sys.stderr)
            return EXIT_INPUT_ERROR
        print(f"wrote {len(pts)} rational points with denominator {args.denom} to {args.out}")
        return 0

    if not args.input:
        parser.error("an input JSON file (or --from-npy) is required")
    try:
        name, pts = load_polyhedron_json(args.input)
        result = verify(pts, name)
    except InputError as exc:
        print(f"INPUT ERROR: {exc}", file=sys.stderr)
        return EXIT_INPUT_ERROR
    except SelfCheckError as exc:
        print(f"SELF-CHECK FAILED (no result is reported): {exc}", file=sys.stderr)
        return EXIT_SELF_CHECK

    print(format_report(result, verbose=args.verbose))
    if args.report:
        with open(args.report, "w", encoding="utf-8") as fh:
            json.dump(result.to_json(), fh, indent=1)
            fh.write("\n")
        print(f"JSON report written to {args.report}")
    return EXIT_CERTIFIED if result.certified else EXIT_NOT_CERTIFIED


if __name__ == "__main__":
    sys.exit(main())
