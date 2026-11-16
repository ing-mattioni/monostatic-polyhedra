"""
core.py -- trusted kernel of the pentahedron proof: rooted trees, rigorous evaluation of the
chart data on a box (ball arithmetic, python-flint `arb`) and the rigorous certificate check.

A CERTIFICATE for (chart, box X, rooted tree (s, pi)) is an affine multiplier function
        y_i(x) = y0_i + sum_k Y1[i][k] (x_k - c_k)        (c = centre of X; i runs over the 4 faces != s)
such that  y_i(x) >= 0  on X  and
        N(x) = sum_i y_i(x) * sym( nu_i(x) eta_{i -> pi(i)}(x)^T )       is NEGATIVE DEFINITE for all x in X.
(See shapes.py for nu, eta and the README for why this refutes the tree on every shape in X.)

Rigorous check (function `check`):
  * y_i(X) is enclosed by the ball  y0_i + sum_k Y1[i][k] * [-h_k, h_k]  and must be >= 0;
  * entrywise mean value theorem on the convex box:  N(x) = N(c) + sum_k (x_k - c_k) dN/dx_k (xi),
    dN/dx_k (x) = sum_i ( Y1[i][k] C_i(x) + y_i(x) dC_i/dx_k (x) ),   C_i = sym(nu_i eta^T);
    a ball enclosure G_k of dN/dx_k over X is computed with ball arithmetic; with D_k = mid(G_k),
    E_k = rad(G_k):   N(x) = L(delta) + Rem,   L(delta) = N(c) + sum_k delta_k D_k  (delta = x - c),
    |Rem_ab| <= sum_k h_k E_k,ab =: Ebar_ab,  hence  ||Rem||_2 <= ||Ebar||_F =: r;
  * L is affine in delta, so -L(delta) - r I > 0 on the box iff it holds at the 8 corners (convexity of
    the positive definite cone); positive definiteness is checked by Sylvester's criterion with
    ball arithmetic (a comparison of balls is True only if it holds for all points of the balls).
All inputs (box corners, multipliers) are binary floating point numbers, i.e. exact dyadic rationals.
"""
import itertools
from flint import arb

from shapes import get_chart, unpack, ADJ_PAIRS, FACES

ADJ = {i: [j for j in range(5) if j != i and (min(i, j), max(i, j)) in ADJ_PAIRS] for i in range(5)}
DIRECTED = [(i, j) for i in range(5) for j in ADJ[i]]          # 18 directed edges
DIR_INDEX = {e: k for k, e in enumerate(DIRECTED)}
SYM_IDX = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]     # storage order of symmetric 3x3 matrices


def rooted_trees():
    """All rooted spanning trees of the face adjacency graph K5 minus {A,B}:
    list of (root s, parent tuple) with parent[s] = -1 and parent[i] adjacent to i."""
    out = []
    for s in range(5):
        others = [i for i in range(5) if i != s]
        for choice in itertools.product(*[ADJ[i] for i in others]):
            par = [-1] * 5
            for i, j in zip(others, choice):
                par[i] = j
            ok = True
            for i in others:
                x, steps = i, 0
                while x != s and steps <= 5:
                    x = par[x]
                    steps += 1
                if x != s:
                    ok = False
                    break
            if ok:
                out.append((s, tuple(par)))
    return out


TREES = rooted_trees()
if len(TREES) != 375:          # explicit (not an assert: must survive python -O)
    raise RuntimeError("expected 375 rooted spanning trees")


def tree_edges(tid):
    s, par = TREES[tid]
    return [DIR_INDEX[(i, par[i])] for i in range(5) if i != s]


def _sym_outer(a, b):
    """sym(a b^T) in SYM_IDX order"""
    return [(a[p] * b[q] + a[q] * b[p]) / 2 for (p, q) in SYM_IDX]


def _edge_matrices(vals):
    """vals = chart evaluator output -> for each directed edge: (C, [dC/dx0, dC/dx1, dC/dx2])"""
    nu, eta = unpack(vals)
    out = []
    for (i, j) in DIRECTED:
        key, sg = ((i, j), 1) if i < j else ((j, i), -1)
        n, dn = nu[i]
        e, de = eta[key]
        e = [sg * t for t in e]
        de = [[sg * t for t in row] for row in de]
        C = _sym_outer(n, e)
        dC = []
        for k in range(3):
            a = _sym_outer(dn[k], e)
            b = _sym_outer(n, de[k])
            dC.append([a[m] + b[m] for m in range(6)])
        out.append((C, dC))
    return out


class BoxData:
    """Rigorous data of a chart on the box [lo, hi] (floats)."""

    def __init__(self, chart_name, lo, hi):
        ch = get_chart(chart_name)
        self.lo, self.hi = tuple(lo), tuple(hi)
        self.c = [(arb(l) + arb(h)) / 2 for l, h in zip(lo, hi)]         # exact (dyadic)
        self.h = [(arb(h) - arb(l)) / 2 for l, h in zip(lo, hi)]
        X = [c + h * arb(0, 1) for c, h in zip(self.c, self.h)]         # balls enclosing the box
        self.at_c = _edge_matrices(ch["fn"](*self.c))
        self.at_X = _edge_matrices(ch["fn"](*X))
        self.unit = [h * arb(0, 1) for h in self.h]                      # x_k - c_k for x in the box


def _posdef(M):
    """M in SYM_IDX order (arb entries).  True only if certainly positive definite."""
    xx, xy, xz, yy, yz, zz = M
    if not (xx > 0):
        return False
    m2 = xx * yy - xy * xy
    if not (m2 > 0):
        return False
    det = xx * (yy * zz - yz * yz) - xy * (xy * zz - yz * xz) + xz * (xy * yz - yy * xz)
    return bool(det > 0)


def check(bd, tid, y0, Y1):
    """Rigorous check of a certificate (see module docstring).  y0: 4 floats, Y1: 4x3 floats."""
    edges = tree_edges(tid)
    y0a = [arb(t) for t in y0]
    Y1a = [[arb(t) for t in row] for row in Y1]
    yX = []
    for i in range(4):
        b = y0a[i] + Y1a[i][0] * bd.unit[0] + Y1a[i][1] * bd.unit[1] + Y1a[i][2] * bd.unit[2]
        if not (b.lower() >= 0):
            return False
        yX.append(b)
    Nc = [arb(0)] * 6
    for i in range(4):
        C = bd.at_c[edges[i]][0]
        Nc = [Nc[m] + y0a[i] * C[m] for m in range(6)]
    D, Ebar = [], [arb(0)] * 6
    for k in range(3):
        G = [arb(0)] * 6
        for i in range(4):
            CX, dCX = bd.at_X[edges[i]]
            G = [G[m] + Y1a[i][k] * CX[m] + yX[i] * dCX[k][m] for m in range(6)]
        D.append([g.mid() for g in G])
        Ebar = [Ebar[m] + bd.h[k] * G[m].rad() for m in range(6)]
    r2 = Ebar[0] ** 2 + Ebar[3] ** 2 + Ebar[5] ** 2 + 2 * (Ebar[1] ** 2 + Ebar[2] ** 2 + Ebar[4] ** 2)
    r = r2.upper().sqrt().upper()          # rigorous upper bound of the Frobenius norm
    for sg in itertools.product((-1, 1), repeat=3):
        M = [-(Nc[m] + sg[0] * bd.h[0] * D[0][m] + sg[1] * bd.h[1] * D[1][m] + sg[2] * bd.h[2] * D[2][m])
             for m in range(6)]
        M[0] = M[0] - r
        M[3] = M[3] - r
        M[5] = M[5] - r
        if not _posdef(M):
            return False
    return True
