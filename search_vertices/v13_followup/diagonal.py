"""Exact search over every positive diagonal stretch after deleting one vertex.

Diagonal metrics up to scale form a triangle. For each vertex other than a
proposed root, retain the union of half planes admitting a nonnegative escape
dot product. Including equality also tests the boundaries between chambers.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "verify"))
import exact_verify as ev


def clip(poly, row):
    out = []
    prev = poly[-1]
    val = lambda p: row[0]*p[0] + row[1]*p[1] + row[2]*(1-p[0]-p[1])
    a = val(prev)
    for cur in poly:
        b = val(cur)
        if a*b < 0:
            t = a/(a-b)
            out.append(tuple(x+t*(y-x) for x,y in zip(prev,cur)))
        if b >= 0:
            out.append(cur)
        prev, a = cur, b
    return list(dict.fromkeys(out))


def inside(poly):
    # The barycentre of a convex set of nonnegative triples is strictly
    # positive iff no coordinate vanishes throughout that set.
    return poly and any(x>0 for x,y in poly) and any(y>0 for x,y in poly) and any(x+y<1 for x,y in poly)


def solve(rows, poly):
    if not rows:
        n = len(poly)
        return tuple(sum(p[k] for p in poly)/n for k in (0,1))
    options = []
    for i, choices in enumerate(rows):
        pieces = []
        for row in choices:
            piece = clip(poly, row)
            if piece == poly:
                return solve(rows[:i]+rows[i+1:], poly)
            if inside(piece):
                pieces.append(piece)
        if not pieces:
            return None
        options.append((len(pieces), i, pieces))
    _, i, pieces = min(options, key=lambda item:item[0])
    for piece in pieces:
        found = solve(rows[:i]+rows[i+1:], piece)
        if found is not None:
            return found
    return None


def test(points):
    hull = ev.convex_hull(points)
    volume, c = ev.volume_and_centroid(points, hull)
    rows = [[tuple((p[k]-c[k])*(q[k]-p[k]) for k in range(3))
             for j,q in enumerate(points) if i != j] for i,p in enumerate(points)]
    feasible = []
    for root in range(len(points)):
        xy = solve(rows[:root]+rows[root+1:], [(Q(0),Q(0)),(Q(1),Q(0)),(Q(0),Q(1))])
        if xy is not None:
            weights = (*xy, 1-sum(xy))
            maxima = [max(sum(a*b for a,b in zip(w,row)) for row in choices)
                      for choices in rows for w in [weights]]
            feasible.append(dict(root=root, weights=list(map(str,weights)), zero=maxima.count(0),
                                 unstable=sum(v<0 for v in maxima)))
    return dict(vertices=len(hull.vertices), volume=str(volume), possible_roots=feasible)


if __name__ == "__main__":
    _, points = ev.load_polyhedron_json(str(ROOT/"examples/V14_integer.json"))
    # Positive control: restricting the admissible root to 0 must retain
    # at least one feasible metric for the undeleted exact example.
    control = test(points)
    if 0 not in [r["root"] for r in control["possible_roots"]]:
        raise ValueError("Known positive control failed")
    results = [test(points[:i]+points[i+1:]) for i in range(len(points))]
    (HERE/"diagonal_results.json").write_text(json.dumps(dict(control=control, deletions=results),indent=2)+"\n")
    print([(i,len(r["possible_roots"])) for i,r in enumerate(results)])
