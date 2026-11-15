"""
propose.py -- UNTRUSTED helper: propose affine multipliers (y0, Y1) for a (box, tree) by solving a
small semidefinite programme in floating point.  Whatever it returns is checked rigorously by core.check.

SDP (sample points: centre and the 8 corners of the box, xi = (x - c)/h in {0} u {-1,1}^3):
    maximise lam   s.t.   - sum_i yh_i(xi_s) Chat_i(x_s)  >=  lam * I   for all samples s,
                          yh_i(xi) = yh0_i + Yh1[i] . xi  >= 0 at the corners,   sum_i yh0_i = 1,
where Chat_i = C_i / scale_i (scale_i = Frobenius norm of C_i at the centre).
"""
import itertools
import numpy as np
import cvxpy as cp

from shapes import get_chart
from core import _edge_matrices, tree_edges, SYM_IDX

CORNERS = [np.array(s, float) for s in itertools.product((-1, 1), repeat=3)]
SAMPLES = [np.zeros(3)] + CORNERS


def _full(sym6):
    M = np.zeros((3, 3))
    for m, (p, q) in enumerate(SYM_IDX):
        M[p, q] = M[q, p] = sym6[m]
    return M


def edge_mats_float(chart_name, x):
    ch = get_chart(chart_name)
    em = _edge_matrices([float(t) for t in ch["fn"](*[float(t) for t in x])])
    return [_full(C) for C, dC in em]


class Proposer:
    def __init__(self):
        self.P = [[cp.Parameter((3, 3), symmetric=True) for _ in range(4)] for _ in SAMPLES]
        self.y0 = cp.Variable(4)
        self.Y1 = cp.Variable((4, 3))
        self.lam = cp.Variable()
        cons = [cp.sum(self.y0) == 1]
        for s, xi in enumerate(SAMPLES):
            expr = 0
            for i in range(4):
                expr = expr + (self.y0[i] + self.Y1[i, :] @ xi) * self.P[s][i]
            cons.append(-expr - self.lam * np.eye(3) >> 0)
        for xi in CORNERS:
            cons.append(self.y0 + self.Y1 @ xi >= 1e-6)      # strictly positive: robust to rounding
        self.prob = cp.Problem(cp.Maximize(self.lam), cons)

    def propose(self, chart_name, lo, hi, tid, sample_mats=None):
        """-> (lam, y0, Y1) in true (unscaled) units, or None"""
        lo, hi = np.asarray(lo, float), np.asarray(hi, float)
        c, h = (lo + hi) / 2, (hi - lo) / 2
        edges = tree_edges(tid)
        if sample_mats is None:
            sample_mats = [edge_mats_float(chart_name, c + h * xi) for xi in SAMPLES]
        scale = [max(np.linalg.norm(sample_mats[0][e]), 1e-300) for e in edges]
        for s in range(len(SAMPLES)):
            for i, e in enumerate(edges):
                self.P[s][i].value = sample_mats[s][e] / scale[i]
        try:
            self.prob.solve(solver="CLARABEL")
        except Exception:
            return None
        if self.lam.value is None or self.lam.value <= 0:
            return None
        y0 = np.array(self.y0.value) / np.array(scale)
        Y1 = np.array(self.Y1.value) / np.array(scale)[:, None] / h[None, :]
        return float(self.lam.value), [float(t) for t in y0], [[float(t) for t in row] for row in Y1]
