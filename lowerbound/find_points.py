"""
Step 5: numerical search for actual solutions of the systems that remain
uncertified, and geometric analysis of the solutions found.

Unknowns: r_i in R^3 for i != v, with r_v = (1,0,0).  We maximise the minimum
slack t (local optimisation, scipy SLSQP, many starts) in several modes:

  "relaxed"  exactly the system (4.3) of the paper, for the case's face F:
                 |r_i|^2 - r_i.r_j(i) <= -t            (i != v)
                 1 + sum_{k in F} r_{k,1} <= -t
  "convex"   relaxed + convexity with the combinatorics of G: for every face
             (a,b,c) of G (consistently oriented) and every other vertex k,
                 det(r_b - r_a, r_c - r_a, r_k - r_a) >= t,
             i.e. the points are in strictly convex position and the boundary
             complex of their hull is exactly G (a GEOMETRIC solution: the
             paper's option (4)).
  "hull3d"   (exploratory strengthening, per (G, v, tree), independent of F)
             shadow inequalities + "the origin is a convex combination of the
             centroids c_T of ALL tetrahedra conv(v, F'), F' not containing v"
             (all three coordinates; this is necessary since the true centroid
             is the volume-weighted average of the c_T).  Implemented with
             extra variables lam_T >= 0, sum lam_T = 1, sum lam_T c_T = 0
             (bilinear, so still a quadratic system).
  "hull3d_convex"  hull3d + the convexity constraints above.
  "full"     (exploratory, sample only) convexity + the TRUE centroid of the
             solid at the origin, i.e. an actual mono-unstable polyhedron with
             this graph, unstable vertex v and shadowing tree.  Because flat
             (degenerate) polyhedra drive every slack to 0, convexity is
             imposed as a hard margin: normalised determinants
             det/(|u||w||s|) >= DELTA_FULL, and t is the minimum normalised
             shadow slack (r_i.r_j - |r_i|^2)/|r_j|^2.

A mode is "solved" for a system if some start reaches t > 1e-7, re-evaluated
independently of the optimiser's status (all inequalities then hold strictly).
Local search: failure to find a point is evidence, not proof.

For every point set found we compute with scipy.spatial.ConvexHull: convex
position, whether the hull graph equals G (labelled), the true centroid of the
solid, the unstable vertices w.r.t. the true centroid (mono-unstable?) and
whether the origin lies in the convex hull of the tetrahedra centroids.

Usage:  python find_points.py [V] [--starts 24] [--workers 4] [--full-sample 40]
        python find_points.py retry [V]   (second pass, more starts, for the
                                           systems where "convex" failed)
Reads uncertified_V<V>.json (from strengthen.py); writes V<V>_feasible_points.json.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import Counter, defaultdict

import numpy as np
from scipy.optimize import linprog, minimize
from scipy.spatial import ConvexHull

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import enumerate_cases as ec  # noqa: E402

SOLVED_TOL = 1e-7
DELTA_FULL = 0.02
MODES_PER_SYSTEM = ("relaxed", "convex")
MODES_PER_TREE = ("hull3d", "hull3d_convex")


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def _cross(a, b):
    """Row-wise cross product (much faster than np.cross for small arrays)."""
    return np.stack([a[:, 1] * b[:, 2] - a[:, 2] * b[:, 1],
                     a[:, 2] * b[:, 0] - a[:, 0] * b[:, 2],
                     a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]], axis=1)


def oriented_faces(G):
    """Faces of the triangulation G, all oriented consistently (as traversed
    by networkx's planar embedding)."""
    import networkx as nx
    ok, emb = nx.check_planarity(G)
    assert ok
    seen = set()
    out = []
    for u, w in emb.edges():
        if (u, w) in seen:
            continue
        out.append(tuple(emb.traverse_face(u, w, mark_half_edges=seen)))
    return out


def solid_centroid_from_faces(R, ofaces):
    """Centroid and volume of the polyhedron with vertex positions R and
    consistently oriented triangular boundary `ofaces` (signed volumes)."""
    F = np.asarray(ofaces)
    Ra, Rb, Rc = R[F[:, 0]], R[F[:, 1]], R[F[:, 2]]
    d = np.einsum("ij,ij->i", Ra, _cross(Rb, Rc)) / 6.0
    vol = d.sum()
    mom = (d[:, None] * (Ra + Rb + Rc)).sum(0) / 4.0
    return mom / vol, abs(vol)


def origin_in_hull(P):
    """Is 0 in conv(rows of P)?  LP feasibility; returns (bool, lambda)."""
    m = len(P)
    A_eq = np.vstack([P.T, np.ones((1, m))])
    b_eq = np.array([0.0, 0.0, 0.0, 1.0])
    r = linprog(np.zeros(m), A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * m,
                method="highs")
    return r.status == 0, (r.x.tolist() if r.status == 0 else None)


def hull_analysis(R, G, v, decomposition):
    """Convex-hull based analysis of the point set R (rows = vertices)."""
    V = len(R)
    out = {}
    cT = np.array([(R[v] + R[a] + R[b] + R[c]) / 4 for a, b, c in decomposition])
    out["origin_in_conv_tetra_centroids"] = origin_in_hull(cT)[0]
    try:
        H = ConvexHull(R)
    except Exception as e:                      # coplanar etc.
        rank = np.linalg.matrix_rank(R - R.mean(0), tol=1e-9)
        out.update(degenerate=True, affine_rank=int(rank), error=str(e)[:80])
        return out
    hv = sorted(int(i) for i in H.vertices)
    out["degenerate"] = False
    out["hull_vertices"] = hv
    out["convex_position"] = len(hv) == V
    interior_pt = R[hv].mean(0)
    vol, mom = 0.0, np.zeros(3)
    edges = set()
    for s in H.simplices:
        a, b, c = R[s]
        d = abs(np.linalg.det(np.array([a - interior_pt, b - interior_pt,
                                        c - interior_pt]))) / 6.0
        vol += d
        mom += d * (a + b + c + interior_pt) / 4.0
        for x, y in ((s[0], s[1]), (s[1], s[2]), (s[0], s[2])):
            edges.add(tuple(sorted((int(x), int(y)))))
    cen = mom / vol
    out["volume"] = float(vol)
    out["true_centroid"] = cen.tolist()
    out["hull_graph_equals_G"] = bool(out["convex_position"] and
                                      edges == {tuple(sorted(e)) for e in G.edges()})
    # equilibria w.r.t. the true centroid (vertex i unstable iff
    # (r_j - r_i).(r_i - c) < 0 for all other hull vertices j)
    unstable = [i for i in hv
                if all((R[j] - R[i]) @ (R[i] - cen) < 0 for j in hv if j != i)]
    out["unstable_vertices_true_centroid"] = unstable
    out["mono_unstable"] = len(unstable) == 1
    out["farthest_from_true_centroid"] = int(np.argmax(np.linalg.norm(R - cen, axis=1)))
    out["v_unstable_wrt_true_centroid"] = v in unstable
    # and w.r.t. the origin (the centroid assumed by the system)
    out["unstable_vertices_origin"] = [
        i for i in hv if all((R[j] - R[i]) @ R[i] < 0 for j in hv if j != i)]
    return out


# ---------------------------------------------------------------------------
# The optimisation problems
# ---------------------------------------------------------------------------

class Problem:
    """Constraint functions (in the form  g(x) >= 0) and Jacobians.
    Variable vector x = (coordinates of r_i, i != v  [3(V-1)],
                         lam_T  [only in hull3d modes],  t)."""

    def __init__(self, V, v, parent, face, ofaces, mode):
        self.V, self.v, self.parent, self.face, self.mode = V, v, parent, face, mode
        self.others = [i for i in range(V) if i != v]
        self.idx = {i: k for k, i in enumerate(self.others)}
        self.n = 3 * (V - 1)
        self.ofaces = ofaces
        self.decomp = [tuple(sorted(f)) for f in ofaces if v not in f]
        self.nl = len(self.decomp) if mode.startswith("hull3d") else 0
        self.nx = self.n + self.nl + 1
        self.use_tree = True
        self.use_com = mode in ("relaxed", "convex")
        self.use_conv = mode in ("convex", "hull3d_convex", "full")
        # convexity triples (oriented face, other vertex).  The orientation
        # sign is fixed: reflection in the plane z = 0 fixes r_v and all
        # shadow / centroid constraints and flips every determinant, so the
        # opposite sign only yields the mirror images.
        self.ctrip = [(f, k) for f in ofaces for k in range(V) if k not in f]
        self._ct = tuple(np.array(a) for a in zip(*[(f[0], f[1], f[2], k)
                                                    for f, k in self.ctrip]))
        self._oth = np.array(self.others)
        self._par = np.array([parent[i] for i in self.others])
        self._col = np.full(V, -10 ** 6)
        for i in self.others:
            self._col[i] = 3 * self.idx[i]

    # -- helpers --
    def R(self, x):
        R = np.zeros((self.V, 3))
        R[self.v] = (1.0, 0.0, 0.0)
        R[self.others] = x[:self.n].reshape(-1, 3)
        return R

    def lam(self, x):
        return x[self.n:self.n + self.nl]

    def _scatter(self, J, rows, verts, grads):
        """J[rows, coordinates of verts] += grads (v is fixed: skipped)."""
        for c in range(3):
            cols = self._col[verts] + c
            ok = cols >= 0
            np.add.at(J, (rows[ok], cols[ok]), grads[ok, c])

    def _conv_parts(self, R):
        """Orientation determinants det(r_b-r_a, r_c-r_a, r_k-r_a) and their
        gradients w.r.t. r_a, r_b, r_c, r_k (vectorised over all triples)."""
        A, B, Cc, K = self._ct
        u, w, s = R[B] - R[A], R[Cc] - R[A], R[K] - R[A]
        uw = _cross(u, w)
        dets = np.einsum("ij,ij->i", uw, s)
        gb, gc, gk = _cross(w, s), _cross(s, u), uw
        return dets, (-(gb + gc + gk), gb, gc, gk), (u, w, s)

    # -- inequality constraints g(x) >= 0 --
    def ineq(self, x):
        R, t = self.R(x), x[-1]
        vals = []
        Ri, P = R[self._oth], R[self._par]
        slack = np.einsum("ij,ij->i", Ri, P) - np.einsum("ij,ij->i", Ri, Ri)
        if self.mode == "full":       # normalised shadow slack, hard convexity margin
            vals.append(slack / np.einsum("ij,ij->i", P, P) - t)
            dets, _, (u, w, s) = self._conv_parts(R)
            nrm = (np.linalg.norm(u, axis=1) * np.linalg.norm(w, axis=1)
                   * np.linalg.norm(s, axis=1))
            vals.append(dets / np.maximum(nrm, 1e-12) - DELTA_FULL)
            return np.concatenate(vals)
        vals.append(slack - t)
        if self.use_com:
            vals.append([-(1.0 + R[list(self.face), 0].sum()) - t])
        if self.use_conv:
            vals.append(self._conv_parts(R)[0] - t)
        return np.concatenate(vals)

    def ineq_jac(self, x):
        if self.mode == "full":       # finite differences are fine here
            return None
        R = self.R(x)
        m = len(self.others)
        rows = m + (1 if self.use_com else 0) + (len(self.ctrip) if self.use_conv else 0)
        J = np.zeros((rows, self.nx))
        J[:, -1] = -1.0
        rr = np.arange(m)
        Ri, P = R[self._oth], R[self._par]
        self._scatter(J, rr, self._oth, -(2 * Ri - P))
        self._scatter(J, rr, self._par, Ri)
        row = m
        if self.use_com:
            for k in self.face:
                if k != self.v:
                    J[row, self._col[k]] -= 1.0
            row += 1
        if self.use_conv:
            _, grads, _ = self._conv_parts(R)
            rr = row + np.arange(len(self.ctrip))
            for verts, g in zip(self._ct, grads):
                self._scatter(J, rr, verts, g)
        return J

    # -- equality constraints --
    def eq(self, x):
        R = self.R(x)
        if self.mode == "full":
            return solid_centroid_from_faces(R, self.ofaces)[0]
        if self.nl:
            lam = self.lam(x)
            cT = np.array([(R[self.v] + R[a] + R[b] + R[c]) / 4
                           for a, b, c in self.decomp])
            return np.append(lam @ cT, lam.sum() - 1.0)
        return None

    def eq_jac(self, x):
        R = self.R(x)
        if self.mode == "full":
            return self._centroid_jac(R)
        lam = self.lam(x)
        J = np.zeros((4, self.nx))
        for T, (a, b, c) in enumerate(self.decomp):
            cT = (R[self.v] + R[a] + R[b] + R[c]) / 4
            J[:3, self.n + T] = cT
            for k in (a, b, c):
                if k != self.v:
                    J[:3, self._col[k]:self._col[k] + 3] += np.eye(3) * lam[T] / 4
        J[3, self.n:self.n + self.nl] = 1.0
        return J

    def _centroid_jac(self, R):
        """Jacobian of the solid centroid c = mom/vol (vol = sum d_f,
        mom = sum d_f (r_a+r_b+r_c)/4, d_f = r_a.(r_b x r_c)/6)."""
        F = np.array(self.ofaces)
        Ra, Rb, Rc = R[F[:, 0]], R[F[:, 1]], R[F[:, 2]]
        d = np.einsum("ij,ij->i", Ra, _cross(Rb, Rc)) / 6.0
        S4 = (Ra + Rb + Rc) / 4.0
        vol = d.sum()
        c = (d[:, None] * S4).sum(0) / vol
        gd = (_cross(Rb, Rc) / 6.0, _cross(Rc, Ra) / 6.0, _cross(Ra, Rb) / 6.0)
        J = np.zeros((3, self.nx))
        for slot in range(3):
            for f in range(len(F)):
                p = F[f, slot]
                if p == self.v:
                    continue
                k = self._col[p]
                J[:, k:k + 3] += (np.outer(S4[f] - c, gd[slot][f])
                                  + np.eye(3) * d[f] / 4.0) / vol
        return J

    def achieved_slack(self, x):
        """min slack at x, recomputed from scratch; equality violations count
        as negative slack."""
        y = x.copy()
        y[-1] = 0.0
        g = self.ineq(y)
        m = len(self.others)
        if self.mode == "full":
            t = float(np.min(g[:m]))                 # normalised shadow slack
            margin = float(np.min(g[m:]))            # convexity margin (hard)
            if margin < -1e-9:
                t = min(t, margin)
        else:
            t = float(np.min(g))
        e = self.eq(y)
        if e is not None and np.max(np.abs(e)) > 1e-9:
            t = min(t, -float(np.max(np.abs(e))))
        if self.nl and np.min(self.lam(y)) < -1e-12:
            t = min(t, float(np.min(self.lam(y))))
        return t


# ---------------------------------------------------------------------------
# Starting points
# ---------------------------------------------------------------------------

def tree_sample(V, v, parent, rng):
    """Random configuration satisfying all shadow inequalities: each r_i is
    uniform in the Thales ball of its parent (ball with diameter [0, r_j])."""
    R = np.zeros((V, 3))
    R[v] = (1.0, 0.0, 0.0)
    done = {v}
    while len(done) < V:
        for i in range(V):
            if i not in done and parent[i] in done:
                c = R[parent[i]] / 2
                u = rng.normal(size=3)
                u *= rng.random() ** (1 / 3) / np.linalg.norm(u)
                R[i] = c + np.linalg.norm(c) * u
                done.add(i)
    return R


def spiral_sample(V, v, parent, rng):
    """Every child rotated by a random angle (mostly in one rotational sense)
    about a random axis close to e_3, and shrunk by ~cos(angle): a good start
    for the spiralling solutions."""
    R = np.zeros((V, 3))
    R[v] = (1.0, 0.0, 0.0)
    done = {v}
    k = np.array([0.0, 0.0, 1.0]) + 0.3 * rng.normal(size=3)
    k /= np.linalg.norm(k)
    while len(done) < V:
        for i in range(V):
            if i not in done and parent[i] in done:
                p = R[parent[i]]
                th = rng.uniform(0.2, 0.9) * rng.choice([1, 1, 1, -1])
                rot = (p * np.cos(th) + np.cross(k, p) * np.sin(th)
                       + k * (k @ p) * (1 - np.cos(th)))
                R[i] = 0.97 * np.cos(th) * rot
                done.add(i)
    return R


def convex_realization(prob, rng, max_tries=5000):
    """Random convex realization of G with r_v = (1,0,0), oriented like
    prob.ofaces.  Rejection sampling: the hull of V random points on the unit
    sphere is a triangulation; keep it if isomorphic to G, relabel by a random
    isomorphism, map r_v to (1,0,0) by a Householder reflection and reflect in
    z = 0 if the orientation is wrong.  Returns the R array or None."""
    import networkx as nx
    from networkx.algorithms.isomorphism import GraphMatcher
    V = prob.V
    G = nx.Graph()
    G.add_nodes_from(range(V))
    for a, b, c in prob.ofaces:
        G.add_edges_from(((a, b), (b, c), (a, c)))
    e1 = np.array([1.0, 0.0, 0.0])
    for _ in range(max_tries):
        P = rng.normal(size=(V, 3))
        P /= np.linalg.norm(P, axis=1)[:, None]
        try:
            H = ConvexHull(P)
        except Exception:
            continue
        if len(H.vertices) < V:
            continue
        Hg = nx.Graph()
        Hg.add_nodes_from(range(V))
        for t in H.simplices:
            Hg.add_edges_from(((t[0], t[1]), (t[1], t[2]), (t[0], t[2])))
        isos = list(GraphMatcher(G, Hg).isomorphisms_iter())
        if not isos:
            continue
        phi = isos[rng.integers(len(isos))]
        R = np.array([P[phi[i]] for i in range(V)])
        w = R[prob.v] - e1
        if np.linalg.norm(w) > 1e-12:
            R = R @ (np.eye(3) - 2 * np.outer(w, w) / (w @ w)).T
        dets = prob._conv_parts(R)[0]
        if np.all(dets < 0):
            R[:, 2] *= -1
        elif not np.all(dets > 0):
            continue
        return R
    return None


def _neg_t(x):
    return -x[-1]


def _neg_t_jac(x):
    g = np.zeros_like(x)
    g[-1] = -1.0
    return g


def solve_mode(prob, starts, rng, R_init=None, maxiter=400):
    """Maximise t from several starts; returns (best t, best x)."""
    Rs = []
    if R_init is not None:
        Rs.append(R_init)
        for _ in range(max(1, starts // 4)):
            Rs.append(R_init + np.vstack([np.zeros((1, 3)) if i == prob.v else
                                          rng.normal(scale=0.05, size=(1, 3))
                                          for i in range(prob.V)]))
    for s in range(starts):
        if prob.use_conv and s % 2 == 0:
            R = convex_realization(prob, rng)
            if R is not None:
                Rs.append(R)
                continue
        Rs.append(spiral_sample(prob.V, prob.v, prob.parent, rng) if s % 2 == 0
                  else tree_sample(prob.V, prob.v, prob.parent, rng))
    cons = [dict(type="ineq", fun=prob.ineq)]
    if prob.mode != "full":
        cons[0]["jac"] = prob.ineq_jac
    if prob.mode == "full" or prob.nl:
        cons.append(dict(type="eq", fun=prob.eq, jac=prob.eq_jac))
    bounds = [(-1.5, 1.5)] * prob.n + [(0.0, 1.0)] * prob.nl + [(-2.0, 1.0)]
    best = (-np.inf, None)
    for R0 in Rs:
        x0 = np.concatenate([R0[prob.others].ravel(),
                             np.full(prob.nl, 1.0 / max(prob.nl, 1)), [0.0]])
        x0[-1] = min(0.0, prob.achieved_slack(x0))
        try:
            r = minimize(_neg_t, x0, jac=_neg_t_jac, constraints=cons,
                         bounds=bounds, method="SLSQP",
                         options=dict(maxiter=maxiter, ftol=1e-12))
        except Exception:
            continue
        t = prob.achieved_slack(r.x)
        if t > best[0]:
            best = (t, r.x)
        if best[0] > 1e-3:
            break
    return best


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

_GDS_CACHE: dict = {}


def _graph(V, gid):
    if V not in _GDS_CACHE:
        _GDS_CACHE[V] = ec.graph_data(V, with_trees=False)
    gd = _GDS_CACHE[V][gid]
    return gd.graph, gd


def run_modes(args):
    """Run `modes` (in order, each warm-started from the previous mode's
    solution when there is one) on one (case, ...) task."""
    case, V, starts, modes, seed = args
    G, gd = _graph(V, case["graph"])
    of = oriented_faces(G)
    rng = np.random.default_rng(seed)
    out = dict(graph=case["graph"], root=case["root"], parent=case["parent"],
               face=case.get("face"))
    R_prev = None
    for mode in modes:
        prob = Problem(V, case["root"], case["parent"], case.get("face"), of, mode)
        t, x = solve_mode(prob, starts if mode != "full" else max(4, starts // 4),
                          rng, R_init=R_prev, maxiter=300 if mode == "full" else 400)
        rec = dict(max_min_slack=t, solved=bool(t > SOLVED_TOL))
        if x is not None:
            R = prob.R(x)
            rec["points"] = R.tolist()
            if prob.use_com:
                rec["tetra_x_sum"] = float(R[case["root"], 0]
                                           + sum(R[k, 0] for k in case["face"]))
            if prob.nl:
                rec["lambda"] = prob.lam(x).tolist()
            rec["analysis"] = hull_analysis(R, G, case["root"], prob.decomp)
            if rec["solved"]:
                R_prev = R
        out[mode] = rec
    return out


def summarize_mode(res, mode):
    tried = [r for r in res if mode in r]
    sol = [r for r in tried if r[mode]["solved"]]
    print(f"\nmode {mode}: solution found for {len(sol)} of {len(tried)} tried")
    ts = [r[mode]["max_min_slack"] for r in tried]
    if tried:
        print(f"   max-min slack: min {min(ts):.3g}, median {np.median(ts):.3g}, "
              f"max {max(ts):.3g}")
    if not sol:
        return
    an = [r[mode]["analysis"] for r in sol]
    nd = [a for a in an if not a.get("degenerate", False)]
    print(f"   degenerate (affinely dependent) point sets: {len(an) - len(nd)}")
    print("   convex position:", sum(a["convex_position"] for a in nd),
          "| hull graph == G:", sum(a["hull_graph_equals_G"] for a in nd),
          "| mono-unstable w.r.t. true centroid:", sum(a["mono_unstable"] for a in nd),
          "| v unstable w.r.t. true centroid:",
          sum(a["v_unstable_wrt_true_centroid"] for a in nd))
    print("   #unstable vertices w.r.t. true centroid:",
          dict(sorted(Counter(len(a["unstable_vertices_true_centroid"]) for a in nd).items())))
    print("   #unstable vertices w.r.t. the origin:",
          dict(sorted(Counter(len(a["unstable_vertices_origin"]) for a in nd).items())))
    print("   origin in conv(tetrahedra centroids):",
          sum(a["origin_in_conv_tetra_centroids"] for a in an), "of", len(an))
    if nd:
        cx = np.array([a["true_centroid"] for a in nd])
        vol = np.array([a["volume"] for a in nd])
        print("   true centroid: x in [%.3f, %.3f], |c| in [%.3f, %.3f]; volume in [%.2e, %.2e]"
              % (cx[:, 0].min(), cx[:, 0].max(), np.linalg.norm(cx, axis=1).min(),
                 np.linalg.norm(cx, axis=1).max(), vol.min(), vol.max()))
    by = defaultdict(int)
    for r in sol:
        by[(r["graph"], r["root"])] += 1
    print("   solved, by (graph, root):", dict(sorted(by.items())))


def _retry_task(args):
    """Second pass for a system whose "convex" mode failed: warm starts from
    the relaxed solution lifted out of its plane by random heights, plus more
    random starts, with a different seed."""
    rec, V, starts, seed = args
    G, gd = _graph(V, rec["graph"])
    of = oriented_faces(G)
    rng = np.random.default_rng(seed)
    prob = Problem(V, rec["root"], rec["parent"], rec["face"], of, "convex")
    R0 = np.array(rec["relaxed"]["points"])
    # unit normal of the (nearly planar) relaxed solution
    nrm = np.linalg.svd(R0 - R0.mean(0))[2][-1]
    best = (-np.inf, None)
    for k in range(starts // 2):
        h = rng.normal(scale=rng.choice([0.01, 0.03, 0.1]), size=V)
        h[rec["root"]] = 0.0
        t, x = solve_mode(prob, 1, rng, R_init=R0 + np.outer(h, nrm))
        if t > best[0]:
            best = (t, x)
        if best[0] > 1e-3:
            break
    if best[0] <= SOLVED_TOL:
        t, x = solve_mode(prob, starts // 2, rng)
        if t > best[0]:
            best = (t, x)
    t, x = best
    out = dict(max_min_slack=t, solved=bool(t > SOLVED_TOL), retried=True)
    if x is not None:
        R = prob.R(x)
        out["points"] = R.tolist()
        out["tetra_x_sum"] = float(R[rec["root"], 0] + sum(R[k, 0] for k in rec["face"]))
        out["analysis"] = hull_analysis(R, G, rec["root"], prob.decomp)
    return out


def retry_convex(V, starts, workers):
    """Re-run mode "convex" with more starts for the systems where the first
    pass found no geometric solution; updates V<V>_feasible_points.json."""
    path = os.path.join(HERE, f"V{V}_feasible_points.json")
    with open(path) as f:
        data = json.load(f)
    ps = data["per_system"]
    todo = [k for k, r in enumerate(ps) if not r["convex"]["solved"]]
    print(f"retrying mode 'convex' for {len(todo)} systems with {starts} starts")
    tasks = [(ps[k], V, starts, 777000 + k) for k in todo]
    import multiprocessing as mp
    with mp.get_context("spawn").Pool(min(workers, 4)) as pool:
        res = list(pool.imap(_retry_task, tasks, chunksize=2))
    for k, r in zip(todo, res):
        if r["max_min_slack"] > ps[k]["convex"]["max_min_slack"]:
            ps[k]["convex"] = r
        else:
            ps[k]["convex"]["retried"] = True
    print(f"newly solved: {sum(r['solved'] for r in res)} of {len(todo)}")
    with open(path, "w") as f:
        json.dump(data, f, indent=1)
    summarize_mode(ps, "convex")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "retry":
        V = int(sys.argv[2]) if len(sys.argv) > 2 else 7
        retry_convex(V, starts=64, workers=4)
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("V", nargs="?", type=int, default=7)
    ap.add_argument("--starts", type=int, default=24)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--full-sample", type=int, default=40,
                    help="number of (graph, root, tree) triples, evenly spaced, "
                         "also run in mode 'full'")
    ap.add_argument("--level", default="bc", choices=["base", "b", "bc"],
                    help="take the systems still uncertified at this level")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    with open(os.path.join(HERE, f"uncertified_V{a.V}.json")) as f:
        unc = json.load(f)
    key = {"base": "cert_base", "b": "cert_b", "bc": "cert_bc"}[a.level]
    cases = [c for c in unc if not c[key]]
    if a.limit:
        cases = cases[:a.limit]
    # distinct (graph, root, tree) triples for the per-tree modes
    trees = {}
    for c in cases:
        trees.setdefault((c["graph"], c["root"], tuple(c["parent"])),
                         dict(graph=c["graph"], root=c["root"], parent=c["parent"]))
    trees = list(trees.values())
    print(f"{len(cases)} systems uncertified at level {a.level}; "
          f"{len(trees)} distinct (graph, root, tree)")
    nfs = min(a.full_sample, len(trees))
    full_idx = set(np.linspace(0, len(trees) - 1, nfs).round().astype(int)) if nfs else set()
    tasks = [(c, a.V, a.starts, list(MODES_PER_SYSTEM), 1000 + k)
             for k, c in enumerate(cases)]
    tasks += [(c, a.V, a.starts,
               list(MODES_PER_TREE) + (["full"] if k in full_idx else []),
               500000 + k) for k, c in enumerate(trees)]
    t0 = time.time()
    if a.workers > 1:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(min(a.workers, 4)) as pool:
            res = list(pool.imap(run_modes, tasks, chunksize=2))
    else:
        res = [run_modes(t) for t in tasks]
    print(f"done in {time.time() - t0:.0f}s")
    per_system, per_tree = res[:len(cases)], res[len(cases):]
    with open(os.path.join(HERE, f"V{a.V}_feasible_points.json"), "w") as f:
        json.dump(dict(per_system=per_system, per_tree=per_tree,
                       notes=__doc__), f, indent=1)
    for mode in MODES_PER_SYSTEM:
        summarize_mode(per_system, mode)
    for mode in list(MODES_PER_TREE) + ["full"]:
        summarize_mode(per_tree, mode)


if __name__ == "__main__":
    main()
