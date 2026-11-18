#!/usr/bin/env python3
"""Independent exact check of the fourteen-vertex construction and 0<t<=1 family.

Python 3.10+, standard library only. Rebuilds the triangular hull from every
point triple. Face feet and edge feet are checked against their entire support
conditions using integer numerators, independently of exact_verify.py.

    python verify/check_c.py
    python verify/check_c.py --report examples/V14_family_report.json

The parameter certificate proves U=1 on the whole interval, not constant S/H.
"""
import argparse
from itertools import combinations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def signs(values):
    return ''.join('+' if v > 0 else '-' if v < 0 else '0' for v in values)


def classify(values):
    if all(v > 0 for v in values):
        return 'equilibrium'
    require(any(v < 0 for v in values), 'degenerate equilibrium')
    return 'non-equilibrium'


def check(points):
    n = len(points)
    faces = []
    for i, j, k in combinations(range(n), 3):
        normal = cross(sub(points[j], points[i]), sub(points[k], points[i]))
        if normal == (0, 0, 0):
            continue
        values = [dot(normal, sub(points[l], points[i])) for l in range(n) if l not in (i, j, k)]
        if all(v < 0 for v in values):
            faces.append((i, j, k))
        elif all(v > 0 for v in values):
            faces.append((i, k, j))
        elif all(v <= 0 for v in values) or all(v >= 0 for v in values):
            raise ValueError('nontriangular facet or non-full-dimensional input')
    require(n == 14 and len(faces) == 24, 'unexpected hull counts')
    directed = [(f[i], f[(i+1) % 3]) for f in faces for i in range(3)]
    require(len(set(directed)) == 72, 'repeated oriented edge')
    require(all((j, i) in directed for i, j in directed), 'boundary is not closed')
    edges = sorted({tuple(sorted(e)) for e in directed})
    require(len(edges) == 36 and set(sum((list(f) for f in faces), [])) == set(range(n)), 'missing vertex')
    D, C = 0, [0, 0, 0]
    for f in faces:
        a, b, c = (points[i] for i in f)
        d = dot(a, cross(b, c))
        D += d
        C = [C[l] + d*(a[l]+b[l]+c[l]) for l in range(3)]
    require(D == 90438831 and C == [-948565765, 351349436, 24043586], 'mass identity mismatch')
    W = 4*D
    centered = [tuple(W*p[l]-C[l] for l in range(3)) for p in points]
    vertex_rows, unstable = [], []
    for i in range(n):
        margins = [-dot(sub(points[j], points[i]), centered[i]) for j in range(n) if j != i]
        state = classify(margins)
        if state == 'equilibrium':
            unstable.append(i)
        vertex_rows.append({'vertex': i, 'state': state, 'support_signs': signs(margins)})
    face_rows, stable = [], []
    for f in faces:
        a, b, c = (points[i] for i in f)
        normal = cross(sub(b, a), sub(c, a))
        norm2 = dot(normal, normal)
        h = W*dot(normal, a)-dot(normal, C)
        require(h > 0, 'centroid not strictly interior')
        den = W*norm2
        foot = tuple(norm2*C[l]+h*normal[l] for l in range(3))
        values = []
        for i, j in zip(f, f[1:]+f[:1]):
            along = sub(points[j], points[i])
            to_foot = tuple(foot[l]-den*points[i][l] for l in range(3))
            values.append(dot(normal, cross(along, to_foot)))
        state = classify(values)
        if state == 'equilibrium':
            stable.append(sorted(f))
        face_rows.append({'face': list(f), 'edge_signs': signs(values), 'state': state,
                          'integer_margins': [str(v) for v in values]})
    edge_rows, saddles = [], []
    for i, j in edges:
        a = points[i]
        e = sub(points[j], a)
        L = dot(e, e)
        den = W*L
        tnum = -dot(e, centered[i])
        radial = tuple(L*centered[i][l]+tnum*e[l] for l in range(3))
        support = [-dot(sub(points[k], a), radial) for k in range(n) if k not in (i, j)]
        values = [tnum, den-tnum] + support
        state = classify(values)
        if state == 'equilibrium':
            saddles.append([i, j])
        edge_rows.append({'edge': [i, j], 'segment_signs': signs(values[:2]),
                         'support_signs': signs(support), 'state': state,
                         'integer_margins': [str(v) for v in values]})
    require(unstable == [0], 'unexpected unstable vertices')
    require(sorted(stable) == [[0, 5, 7], [0, 6, 13]], 'unexpected stable faces')
    require(saddles == [[6, 7]], 'unexpected saddle edges')
    parents = [None, 0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13, 0]
    family = []
    root_family = []
    for i, j in [(i, parents[i]) for i in range(1, n)] + [(0, j) for j in range(1, n)]:
        delta = sub(points[j], points[i])
        A = delta[0]*centered[i][0]+delta[1]*centered[i][1]
        B = delta[2]*centered[i][2]
        require((A < 0 and A+B < 0) if i == 0 else (A > 0 and A+B > 0), 'family endpoint sign')
        if i != 0:
            require(tuple(sorted((i, j))) in edges, 'witness is not a true edge')
        row = {'i': i, 'j': j, 'A': str(A), 'B': str(B), 'A_plus_B': str(A+B)}
        (root_family if i == 0 else family).append(row)
    return {'schema': 1, 'V': n, 'E': len(edges), 'F': len(faces), 'S': len(stable),
            'H': len(saddles), 'U': len(unstable), 'D': str(D), 'C': [str(x) for x in C],
            'centroid_denominator': str(W), 'faces': [list(f) for f in faces],
            'stable_faces': sorted(stable), 'saddle_edges': saddles, 'unstable_vertices': unstable,
            'vertex_checks': vertex_rows, 'face_checks': face_rows, 'edge_checks': edge_rows,
            'family_interval': '0 < t <= 1', 'family_witnesses': family, 'family_root': root_family}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    raw = json.loads((ROOT/'examples/V14_integer.json').read_text(encoding='utf-8'))
    points = [tuple(int(c) for c in row) for row in raw['vertices']]
    result = check(points)
    if args.report:
        args.report.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print('PASS: rebuilt hull V=14 E=36 F=24; exact S=2 H=1 U=1; no degenerate equilibria.')
    print('PASS: every vertical contraction 0<t<=1 has the same face lattice and U=1.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
