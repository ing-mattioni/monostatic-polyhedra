#!/usr/bin/env python3
"""
Tests for exact_verify.py.   Run with:   python verify/test_exact_verify.py

The exact verifier itself uses only the standard library.  These tests also
use numpy and scipy, to cross-check the exact results against a completely
independent floating-point implementation (Qhull convex hull + numpy
arithmetic) -- see `float_reference` below.
"""

import io
import json
import math
import os
import random
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout, redirect_stderr
from fractions import Fraction

import numpy as np
from scipy.spatial import ConvexHull

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exact_verify as ev  # noqa: E402

EXPLORE = os.path.join(os.path.dirname(HERE), "explore")
SCRIPT = os.path.join(HERE, "exact_verify.py")


# =============================================================================
# Independent floating-point reference implementation
# =============================================================================

def float_reference(points):
    """Classify equilibria in floating point, independently of exact_verify.

    The hull comes from scipy/Qhull (triangulated; triangles with equal plane
    equations are merged into polygonal faces), the centroid from a
    tetrahedral decomposition in numpy.  Returns a dict with the volume,
    centroid, the sets of stable faces (frozensets of vertex indices), saddle
    edges (frozensets) and unstable vertices, and `min_margin`, the smallest
    |signed margin| of any decision (so callers can tell whether the float
    answer can be trusted).
    """
    P = np.array([[float(x) for x in p] for p in points])
    hull = ConvexHull(P)
    verts = sorted({int(i) for i in hull.simplices.ravel()})

    O = P[verts].mean(axis=0)
    six, mom = 0.0, np.zeros(3)
    for s in hull.simplices:
        a, b, c = P[s]
        d = abs(np.linalg.det(np.array([a - O, b - O, c - O])))
        six += d
        mom += d * (O + a + b + c)
    centroid, volume = mom / (4 * six), six / 6

    groups = []                                    # coplanar triangles -> one face
    for s, eq in zip(hull.simplices, hull.equations):
        for g in groups:
            if np.max(np.abs(g["eq"] - eq)) < 1e-9:
                g["tris"].append(s)
                break
        else:
            groups.append({"eq": eq, "tris": [s]})
    tri_edges = {}
    for gi, g in enumerate(groups):
        for s in g["tris"]:
            for k in range(3):
                i, j, o = int(s[k]), int(s[(k + 1) % 3]), int(s[(k + 2) % 3])
                tri_edges.setdefault(frozenset((i, j)), []).append((gi, i, j, o))
    poly_edges = [e for e, lst in tri_edges.items() if lst[0][0] != lst[1][0]]

    margins = []
    stable = set()
    for gi, g in enumerate(groups):
        n, off = g["eq"][:3], g["eq"][3]
        q = centroid - (n @ centroid + off) * n
        dists = []
        for e, lst in tri_edges.items():
            for (gj, i, j, o) in lst:
                if gj == gi and e in poly_edges:
                    side = np.cross(P[j] - P[i], P[o] - P[i]) @ n
                    val = np.cross(P[j] - P[i], q - P[i]) @ n
                    dists.append(math.copysign(1, side) * val / np.linalg.norm(P[j] - P[i]))
        m = min(dists)
        margins.append(m)
        if m > 0:
            stable.add(frozenset(int(i) for s in g["tris"] for i in s))

    saddles = set()
    for e in poly_edges:
        i, j = sorted(e)
        d = P[j] - P[i]
        t = (centroid - P[i]) @ d / (d @ d)
        q = P[i] + t * d
        u = q - centroid
        h = max((P[w] - q) @ u for w in verts if w not in (i, j)) / np.linalg.norm(u)
        m = min(t * np.linalg.norm(d), (1 - t) * np.linalg.norm(d), -h)
        margins.append(m)
        if m > 0:
            saddles.add(e)

    unstable = set()
    for v in verts:
        u = P[v] - centroid
        h = max((P[w] - P[v]) @ u for w in verts if w != v) / np.linalg.norm(u)
        margins.append(-h)
        if h < 0:
            unstable.add(v)

    return {"volume": volume, "qhull_volume": hull.volume, "centroid": centroid,
            "stable": stable, "saddles": saddles, "unstable": unstable,
            "counts": (len(stable), len(saddles), len(unstable)),
            "min_margin": min(abs(m) for m in margins)}


# =============================================================================
# Helpers
# =============================================================================

def F(*xs):
    return tuple(Fraction(x) for x in xs)


def cube(a=1):
    return [F(x, y, z) for x in (-a, a) for y in (-a, a) for z in (-a, a)]


def box(lx, ly, lz):
    return [F(x, y, z) for x in (0, lx) for y in (0, ly) for z in (0, lz)]


def prism(triangle, h):
    return [F(x, y, 0) for x, y in triangle] + [F(x, y, h) for x, y in triangle]


def random_ellipsoid_points(n, seed, axes=(1, 1, 1), denom=10**6):
    rng = random.Random(seed)
    pts = []
    for _ in range(n):
        g = [rng.gauss(0, 1) for _ in range(3)]
        r = math.sqrt(sum(x * x for x in g))
        pts.append([axes[k] * g[k] / r for k in range(3)])
    return ev.rationalize(pts, denom)


def exact_sets(res):
    return ({frozenset(f.vertices) for f in res.stable_faces},
            {frozenset(e.vertices) for e in res.saddle_edges},
            set(res.unstable_vertices))


def rational_rotation(a, b, c, d):
    """Rotation matrix with rational entries from the integer quaternion (a,b,c,d)."""
    n = a * a + b * b + c * c + d * d
    R = [[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
         [2*(b*c + a*d), a*a - b*b + c*c - d*d, 2*(c*d - a*b)],
         [2*(b*d - a*c), 2*(c*d + a*b), a*a - b*b - c*c + d*d]]
    return [[Fraction(x, n) for x in row] for row in R]


def run_main(argv):
    """Call ev.main(argv) capturing stdout/stderr; return (exit code, stdout, stderr)."""
    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        code = ev.main(argv)
    return code, out.getvalue(), err.getvalue()


# =============================================================================
# Tests
# =============================================================================

class CrossCheckMixin:
    def assert_matches_float(self, res, points=None, require_robust=True):
        """Exact result agrees with the independent float computation."""
        ref = float_reference(points if points is not None else res.points)
        # volume and centroid against scipy/Qhull
        self.assertAlmostEqual(float(res.volume), ref["qhull_volume"],
                               delta=1e-11 * max(1.0, ref["qhull_volume"]))
        scale = max(abs(float(x)) for p in res.points for x in p)
        for k in range(3):
            self.assertAlmostEqual(float(res.centroid[k]), ref["centroid"][k], delta=1e-11 * scale)
        if ref["min_margin"] < 1e-8:
            self.assertFalse(require_robust, "float reference is not robust for this input")
            return ref
        self.assertEqual(exact_sets(res), (ref["stable"], ref["saddles"], ref["unstable"]))
        return ref


class TestClassicalSolids(CrossCheckMixin, unittest.TestCase):

    def check(self, points, counts, VEF, volume=None, centroid=None):
        res = ev.verify(points, "test")
        hull = res.hull
        self.assertEqual((len(hull.vertices), len(hull.edges), len(hull.faces)), VEF)
        self.assertEqual(res.counts, counts)
        self.assertFalse(res.has_degeneracy)
        self.assertTrue(res.poincare_hopf_ok)
        self.assertTrue(res.certified)
        if volume is not None:
            self.assertEqual(res.volume, Fraction(volume))
        if centroid is not None:
            self.assertEqual(res.centroid, F(*centroid))
        self.assert_matches_float(res)
        return res

    def test_cube(self):
        res = self.check(cube(), (6, 12, 8), (8, 12, 6), 8, (0, 0, 0))
        self.assertTrue(all(len(f.vertices) == 4 for f in res.hull.faces))
        self.assertFalse(res.mono_stable or res.mono_unstable)

    def test_tetrahedron(self):
        pts = [F(1, 1, 1), F(1, -1, -1), F(-1, 1, -1), F(-1, -1, 1)]
        self.check(pts, (4, 6, 4), (4, 6, 4), Fraction(8, 3), (0, 0, 0))

    def test_octahedron(self):
        pts = [F(*(s * (k == i) for k in range(3))) for i in range(3) for s in (1, -1)]
        self.check(pts, (8, 12, 6), (6, 12, 8), Fraction(4, 3), (0, 0, 0))

    def test_box_1x2x3(self):
        self.check(box(1, 2, 3), (6, 12, 8), (8, 12, 6), 6, ("1/2", 1, "3/2"))

    def test_box_rational_shifted(self):
        pts = [tuple(x + Fraction(1, 7) for x in p) for p in box(Fraction(1, 3), 2, Fraction(5, 2))]
        self.check(pts, (6, 12, 8), (8, 12, 6), Fraction(5, 3),
                   (Fraction(1, 6) + Fraction(1, 7), 1 + Fraction(1, 7), Fraction(5, 4) + Fraction(1, 7)))


class TestPyramidsAndPrisms(CrossCheckMixin, unittest.TestCase):

    def test_right_square_pyramids(self):
        # For a right pyramid over the square [-1,1]^2 one checks by hand that
        # all 5 faces carry stable points, all 5 vertices are unstable and all
        # 8 edges carry saddles, for EVERY height h: (S,H,U) = (5,8,5).
        for h in (Fraction(1, 100), Fraction(1, 10), 1, 10, 1000):
            pts = [F(-1, -1, 0), F(1, -1, 0), F(1, 1, 0), F(-1, 1, 0), F(0, 0, h)]
            res = ev.verify(pts)
            self.assertEqual(res.counts, (5, 8, 5), f"h = {h}")
            self.assertTrue(res.certified)
            self.assertEqual(res.volume, Fraction(4) * h / 3)
            self.assertEqual(res.centroid, F(0, 0, Fraction(h) / 4))
            self.assert_matches_float(res)

    def test_oblique_pyramids(self):
        # Moving the apex sideways (tall/skinny vs flat) changes the counts;
        # S - H + U = 2 must hold whenever the polyhedron is nondegenerate.
        seen = set()
        base = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        for h in (Fraction(1, 10), 1, 10):
            for ax in (0, Fraction(1, 2), 1, Fraction(3, 2), 3, 10):
                for ay in (Fraction(1, 3), Fraction(-7, 5)):
                    pts = [F(x, y, 0) for x, y in base] + [F(ax, ay, h)]
                    res = ev.verify(pts)
                    if res.has_degeneracy:
                        # an exact degeneracy must show up as a ~zero float margin
                        self.assertLess(float_reference(pts)["min_margin"], 1e-9)
                        continue
                    self.assertTrue(res.poincare_hopf_ok)
                    self.assert_matches_float(res)
                    seen.add(res.counts)
        self.assertGreater(len(seen), 2, f"counts did not vary: {seen}")

    def test_flat_vs_tall_skinny_pyramid_over_triangle(self):
        tri = [(0, 0), (5, 0), (1, 2)]
        counts = {}
        for h in (Fraction(1, 50), 50):
            pts = [F(x, y, 0) for x, y in tri] + [F(4, 1, h)]
            res = ev.verify(pts)
            self.assertTrue(res.certified)
            self.assert_matches_float(res)
            counts[h] = res.counts
        for c in counts.values():
            self.assertEqual(c[0] - c[1] + c[2], 2)

    def test_right_triangular_prism(self):
        # Right prism over the right triangle (0,0),(4,0),(0,3): every face is
        # stable (the centroid of the triangle projects inside each side),
        # every vertex unstable, every edge a saddle: (5, 9, 6).
        pts = prism([(0, 0), (4, 0), (0, 3)], 2)
        res = ev.verify(pts)
        self.assertEqual(res.counts, (5, 9, 6))
        self.assertTrue(res.certified)
        self.assertEqual(sorted(len(f.vertices) for f in res.hull.faces), [3, 3, 4, 4, 4])
        self.assertEqual(res.volume, 12)
        self.assertEqual(res.centroid, F(Fraction(4, 3), 1, 1))
        self.assert_matches_float(res)

    def test_obtuse_triangular_prism(self):
        # The foot of the triangle's centroid on side (6,0)-(5,1) is outside the
        # side, so that rectangular face carries no equilibrium.
        pts = prism([(0, 0), (6, 0), (5, 1)], 3)
        res = ev.verify(pts)
        self.assertTrue(res.certified)
        quad = next(f for f in res.faces if set(f.vertices) == {1, 2, 4, 5})
        self.assertEqual(quad.status, ev.NONEQ)
        self.assertEqual(res.counts[0], 4)
        self.assert_matches_float(res)

    def test_degenerate_prism(self):
        # Triangle B=(0,0), C=(1,0), A=(-1,1): its centroid (0,1/3) projects
        # onto the line BC exactly at B.  In the prism this gives exact
        # degeneracies: the face y=0 (foot of c on its edge BB'), the edges
        # BC and B'C' (foot at an endpoint), the edge BB' (its plane contains
        # the whole face y=0) and the vertices B, B' (orthoplane contains C, C').
        pts = prism([(0, 0), (1, 0), (-1, 1)], 1)       # indices B=0 C=1 A=2 B'=3 C'=4 A'=5
        res = ev.verify(pts)
        self.assertTrue(res.has_degeneracy)
        self.assertFalse(res.certified)
        self.assertFalse(res.mono_stable or res.mono_unstable)
        degenerate_faces = {frozenset(f.vertices) for f in res.faces if f.status == ev.DEGENERATE}
        degenerate_edges = {e.vertices for e in res.edges if e.status == ev.DEGENERATE}
        degenerate_verts = {v.index for v in res.vertices if v.status == ev.DEGENERATE}
        self.assertEqual(degenerate_faces, {frozenset({0, 1, 3, 4})})
        self.assertEqual(degenerate_edges, {(0, 1), (3, 4), (0, 3)})
        self.assertEqual(degenerate_verts, {0, 3})
        reasons = {e.vertices: e.reason for e in res.edges}
        self.assertIn("endpoint", reasons[(0, 1)])
        self.assertIn("also contains vertex", reasons[(0, 3)])
        self.assertEqual(res.centroid, F(0, Fraction(1, 3), Fraction(1, 2)))
        self.assertLess(float_reference(pts)["min_margin"], 1e-12)


class TestHullRobustness(unittest.TestCase):

    def test_redundant_points(self):
        pts = cube() + [F(0, 0, 0), F(1, Fraction(1, 2), Fraction(1, 3)), F(1, 1, 0)]
        res = ev.verify(pts)
        self.assertEqual(res.hull.redundant, {8: "interior", 9: "face", 10: "edge"})
        self.assertEqual(res.hull.vertices, list(range(8)))
        self.assertEqual(res.counts, (6, 12, 8))
        self.assertTrue(res.certified)

    def test_large_coplanar_faces_and_collinear_points(self):
        # A prism over a regular-ish 12-gon with extra points on its edges.
        poly = [(10, 0), (9, 5), (5, 9), (0, 10), (-5, 9), (-9, 5), (-10, 0),
                (-9, -5), (-5, -9), (0, -10), (5, -9), (9, -5)]
        pts = prism(poly, 3) + [F(Fraction(19, 2), Fraction(5, 2), 0), F(10, 0, Fraction(3, 2))]
        res = ev.verify(pts)
        sizes = sorted(len(f.vertices) for f in res.hull.faces)
        self.assertEqual(sizes, [4] * 12 + [12, 12])
        self.assertEqual(res.hull.redundant, {24: "edge", 25: "edge"})
        self.assertTrue(res.poincare_hopf_ok)

    def test_input_errors(self):
        bad_inputs = {
            "too few": [F(0, 0, 0), F(1, 0, 0), F(0, 1, 0)],
            "duplicate": cube() + [F(1, 1, 1)],
            "coplanar": [F(0, 0, 0), F(1, 0, 0), F(0, 1, 0), F(1, 1, 0), F("1/2", "1/3", 0)],
            "collinear": [F(k, 2 * k, 3 * k) for k in range(5)],
            "bad string": [("a", 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
            "float": [(0.5, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
            "short row": [(0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
        }
        for label, pts in bad_inputs.items():
            with self.assertRaises(ev.InputError, msg=label):
                ev.verify(pts)

    def test_similarity_invariance(self):
        # Rotating (by a rational rotation), scaling and translating must map
        # the classification to itself and transform volume/centroid exactly.
        pts = random_ellipsoid_points(20, seed=7, axes=(1, 2, 3))
        R = rational_rotation(1, 2, 3, 4)
        for i in range(3):
            for j in range(3):
                self.assertEqual(sum(R[i][k] * R[j][k] for k in range(3)), int(i == j))
        lam, shift = Fraction(7, 3), F("1/2", -5, "2/9")
        def g(p):
            return tuple(lam * sum(R[i][k] * p[k] for k in range(3)) + shift[i] for i in range(3))
        r1, r2 = ev.verify(pts), ev.verify([g(p) for p in pts])
        self.assertEqual(exact_sets(r1), exact_sets(r2))
        self.assertEqual(r2.volume, lam ** 3 * r1.volume)
        self.assertEqual(r2.centroid, g(r1.centroid))
        self.assertTrue(r1.certified and r2.certified)


class TestRandomPolyhedra(CrossCheckMixin, unittest.TestCase):

    def test_random_ellipsoid_polyhedra(self):
        cases = [(12, 1, (1, 1, 1)), (25, 2, (1, 2, 3)), (40, 3, (3, 1, 1)), (30, 4, (1, 1, 4))]
        for n, seed, axes in cases:
            pts = random_ellipsoid_points(n, seed, axes)
            res = ev.verify(pts)
            self.assertEqual(len(res.hull.vertices), n)
            self.assertTrue(res.certified, f"seed {seed}")
            self.assertTrue(res.poincare_hopf_ok)
            self.assert_matches_float(res)

    def test_random_integer_clouds(self):
        # Points of a coarse lattice: many interior (redundant) points, and
        # coplanar / collinear configurations on the hull are likely.
        for seed in range(4):
            rng = random.Random(100 + seed)
            pts = sorted({F(*(rng.randint(-12, 12) for _ in range(3))) for _ in range(45)})
            res = ev.verify(pts)
            ref = float_reference(pts)
            self.assertAlmostEqual(float(res.volume), ref["qhull_volume"], delta=1e-9)
            for k in range(3):
                self.assertAlmostEqual(float(res.centroid[k]), ref["centroid"][k], delta=1e-9)
            self.assertTrue(any(loc == "interior" for loc in res.hull.redundant.values()))
            if res.has_degeneracy:
                self.assertLess(ref["min_margin"], 1e-9)
            else:
                self.assertTrue(res.poincare_hopf_ok)
                if ref["min_margin"] > 1e-8:
                    self.assertEqual(res.counts, ref["counts"])
                    self.assertEqual(set(res.unstable_vertices), ref["unstable"])


class TestRationalizeAndCLI(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def path(self, name):
        return os.path.join(self.dir, name)

    def test_rationalize(self):
        pts = [[0.0, -0.0, 0.123456789], [-0.123456789, 1.5, -2.25], [0.1234565, -0.1234565, 1e-12]]
        r = ev.rationalize(pts, 10**6)
        self.assertEqual(r[0][0], 0)
        self.assertEqual(r[0][1], 0)
        self.assertEqual(r[0][2], Fraction(123457, 10**6))
        self.assertEqual(r[1][0], -r[0][2])                 # mirror symmetry preserved
        self.assertEqual(r[2][0], -r[2][1])                 # also for round-half cases
        self.assertEqual(r[2][2], 0)
        for p in r:
            for x in p:
                self.assertEqual(10**6 % x.denominator, 0)
        with self.assertRaises(ValueError):
            ev.rationalize([[float("nan"), 0, 0]], 10)

    def test_json_roundtrip_and_decimal_literals(self):
        pts = ev.rationalize(np.random.default_rng(0).normal(size=(10, 3)).tolist(), 10**9)
        ev.write_polyhedron_json(self.path("p.json"), pts, "roundtrip")
        name, back = ev.load_polyhedron_json(self.path("p.json"))
        self.assertEqual((name, back), ("roundtrip", pts))
        # JSON decimal literals are read exactly (0.1 -> 1/10, not a binary float)
        with open(self.path("d.json"), "w") as fh:
            fh.write('{"name": "dec", "vertices": [[0.1, 0, 0], [0, 0.1, 0], [0, 0, 0.1], [0, 0, 0]]}')
        _, pts = ev.load_polyhedron_json(self.path("d.json"))
        self.assertEqual(pts[0][0], Fraction(1, 10))

    def test_cli_exit_codes_and_report(self):
        ev.write_polyhedron_json(self.path("cube.json"), cube(), "cube")
        code, out, _ = run_main([self.path("cube.json"), "--report", self.path("r.json")])
        self.assertEqual(code, 0)
        self.assertIn("(S, H, U) = (6, 12, 8)", out)
        with open(self.path("r.json")) as fh:
            rep = json.load(fh)
        self.assertEqual((rep["S"], rep["H"], rep["U"], rep["certified"]), (6, 12, 8, True))
        self.assertEqual(rep["volume"], "8")

        ev.write_polyhedron_json(self.path("deg.json"), prism([(0, 0), (1, 0), (-1, 1)], 1))
        self.assertEqual(run_main([self.path("deg.json")])[0], 1)

        with open(self.path("dup.json"), "w") as fh:
            json.dump({"vertices": [[0, 0, 0], [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]}, fh)
        code, _, err = run_main([self.path("dup.json")])
        self.assertEqual(code, 2)
        self.assertIn("duplicate", err)

    def test_cli_from_npy_and_subprocess(self):
        np.save(self.path("t.npy"), np.array([[1.0, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / 3)
        code, out, _ = run_main(["--from-npy", self.path("t.npy"), "--denom", "1e6",
                                 "--out", self.path("t.json")])
        self.assertEqual(code, 0)
        _, pts = ev.load_polyhedron_json(self.path("t.json"))
        self.assertEqual(pts[0], F("333333/1000000", "333333/1000000", "333333/1000000"))
        # the real command line, in a fresh interpreter
        proc = subprocess.run([sys.executable, SCRIPT, self.path("t.json")],
                              capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("(S, H, U) = (4, 6, 4)", proc.stdout)


class TestSeedPolyhedra(CrossCheckMixin, unittest.TestCase):
    """The flat pyramids found numerically in explore/ (float analysis: mono-unstable)."""

    def check_seed(self, filename, n_base, denoms):
        path = os.path.join(EXPLORE, filename)
        if not os.path.exists(path):
            self.skipTest(f"{path} not found")
        arr = np.load(path)
        results = {}
        for denom in denoms:
            pts = ev.rationalize(arr.tolist(), denom)
            res = ev.verify(pts, f"{filename} / {denom}")
            S, H, U = res.counts
            print(f"\n    {filename} denom={denom:.0e}: (S,H,U) = ({S},{H},{U}), "
                  f"unstable vertices {res.unstable_vertices}, "
                  f"saddle edges {[e.vertices for e in res.saddle_edges]}, "
                  f"stable faces {[f.vertices if len(f.vertices) < 5 else f'{len(f.vertices)}-gon' for f in res.stable_faces]}, "
                  f"degeneracies {len(res.degeneracies)}, "
                  f"smallest margin {res.smallest_margin()[0]:.2e} ({res.smallest_margin()[1]})",
                  end="")
            # the planar base must be ONE face with all n_base base points as vertices
            base = [f for f in res.hull.faces if len(f.vertices) > 3]
            self.assertEqual(len(base), 1)
            self.assertEqual(sorted(base[0].vertices), list(range(n_base)))
            self.assertEqual(len(res.hull.vertices), n_base + 1)
            self.assertEqual(res.hull.redundant, {})
            # exact symmetry facts: mirror symmetry x -> -x, apex height 1/100 above y = 0
            self.assertEqual(res.centroid[0], 0)
            self.assertEqual(res.centroid[1], Fraction(1, 400))
            self.assertFalse(res.has_degeneracy)
            self.assertTrue(res.poincare_hopf_ok)
            self.assertTrue(res.certified)
            self.assertEqual(res.unstable_vertices, [0])
            self.assertTrue(res.mono_unstable)
            self.assert_matches_float(res)
            results[denom] = res.counts
        self.assertEqual(len(set(results.values())), 1, f"counts depend on the denominator: {results}")

    def test_seed_V18(self):
        self.check_seed("seed_V18.npy", 17, [10**6, 10**9])

    def test_seed_V16(self):
        self.check_seed("seed_V16.npy", 15, [10**8, 10**9])


if __name__ == "__main__":
    unittest.main(verbosity=2)
