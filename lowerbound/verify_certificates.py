"""
Independent exact verification of the infeasibility certificates stored in a
results CSV (results_V5.csv, results_V6.csv, results_V7.csv) written by
certify.py.

For each row marked certified, the multipliers y are read as exact decimal
rationals (exactly the digits stored in the file), the reduced homogenized
matrix  M(y) = sum_k y_k M_k  (size V: constant + first coordinates of the
V-1 non-root vertices) is assembled from scratch -- independently of
certify.py -- and checked to be positive definite in exact integer arithmetic
(fraction-free Gaussian elimination: all leading principal minors > 0).
Since all y_k are checked to be >= 0 as well, Theorem 8 of arXiv:2401.17906
together with the decoupling identity

    sum_k y_k Qhat_k  ~  blockdiag(M(y), A(y), A(y)),  A(y) = trailing block of M(y)

(simultaneous permutation of rows/columns; see certify.py) proves the
infeasibility of the corresponding system.

The system of a row (root v, parent array p, tetrahedron (v, a, b, c)):
    |r_i|^2 - r_i . r_{p[i]} <= 0   for i != v, in increasing order of i,
    1 + r_{a,1} + r_{b,1} + r_{c,1} <= 0,
with r_v = (1, 0, 0); y lists the multipliers in this order.

Usage:  python verify_certificates.py results_V7.csv
"""

import csv
import sys
from fractions import Fraction
from math import lcm


def reduced_matrix_times2(V, v, parent, face):
    """Return the list of integer matrices 2*M_k (k = shadow inequalities in
    increasing vertex order, then the centroid inequality)."""
    others = [i for i in range(V) if i != v]
    pos = {i: 1 + k for k, i in enumerate(others)}       # index in M
    mats = []
    for i in others:
        M = [[0] * V for _ in range(V)]
        a = pos[i]
        M[a][a] += 2                                      # |r_i|^2 (x-part)
        j = parent[i]
        if j == v:                                        # - r_i . e_1 = - x_i
            M[0][a] -= 1
            M[a][0] -= 1
        else:                                             # - r_i . r_j
            b = pos[j]
            M[a][b] -= 1
            M[b][a] -= 1
        mats.append(M)
    M = [[0] * V for _ in range(V)]
    M[0][0] += 2                                          # constant 1 (= r_{v,1})
    for k in face:                                        # + x_k
        M[0][pos[k]] += 1
        M[pos[k]][0] += 1
    mats.append(M)
    return mats


def is_pd(N):
    """Exact positive definiteness of a symmetric integer matrix (Bareiss)."""
    N = [row[:] for row in N]
    d, prev = len(N), 1
    for k in range(d):
        if N[k][k] <= 0:
            return False
        for i in range(k + 1, d):
            for j in range(k + 1, d):
                N[i][j] = (N[i][j] * N[k][k] - N[i][k] * N[k][j]) // prev
        prev = N[k][k]
    return True


def verify_row(r):
    parent = [int(t) for t in r["parent"].split()]
    V = len(parent)
    tet = [int(t) for t in r["tetrahedron"].split()]
    v, face = tet[0], tet[1:]
    assert parent[v] == -1 and v not in face
    y = [Fraction(t) for t in r["y"].split()]
    if min(y) < 0:
        return False
    den = 1
    for t in y:
        den = lcm(den, t.denominator)
    Y = [int(t * den) for t in y]
    mats = reduced_matrix_times2(V, v, parent, face)
    assert len(mats) == len(Y)
    N = [[sum(Y[k] * mats[k][i][j] for k in range(len(Y))) for j in range(V)]
         for i in range(V)]
    return is_pd(N)


def main(path):
    n = n_cert = ok = 0
    bad = []
    with open(path) as f:
        for r in csv.DictReader(f):
            n += 1
            if r["certified"] != "1":
                continue
            n_cert += 1
            if verify_row(r):
                ok += 1
            else:
                bad.append((r["graph"], r["root"], r["parent"], r["tetrahedron"]))
    print(f"{path}: {n} rows, {n_cert} marked certified, {ok} verified exactly "
          f"from the stored decimal multipliers, {len(bad)} failed")
    for b in bad[:20]:
        print("   FAILED:", b)
    return len(bad) == 0


if __name__ == "__main__":
    sys.exit(0 if all(main(p) for p in sys.argv[1:]) else 1)
