"""Cross-check added-point centroids and exactly rank the angular escape margins."""
from fractions import Fraction as Q
from itertools import combinations
import json
import math
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from replacement import check_box,facets,ROOT,HERE,ev,dot,sub,cross


def check(ok,detail):
    if not ok:
        raise ValueError(detail)


def score(path):
    _,p=ev.load_polyhedron_json(str(path))
    res=ev.verify(p)
    c=res.centroid
    centered=[sub(v,c) for v in p]
    radii=[dot(g,g) for g in centered]
    root=max(range(len(p)),key=radii.__getitem__)
    margins=[]
    distances=[]
    for i,v in enumerate(p):
        if i==root:
            continue
        row=[]
        for j,w in enumerate(p):
            if i==j:
                continue
            d=sub(w,v)
            norm2=dot(d,d)
            num=dot(centered[i],d)
            row.append((num*abs(num)/(radii[i]*norm2),j))
            distances.append(norm2)
        margin,witness=max(row)
        margins.append((margin,i,witness))
    s,vertex,witness=min(margins)
    return dict(file=path.name, exact_V=len(res.hull.vertices), counts=list(res.counts),
                unstable=res.unstable_vertices,nondegenerate=res.certified,
                farthest=root,worst=vertex,best_direction=witness,signed_squared_margin=str(s),
                angular_margin=math.copysign(math.sqrt(float(abs(s))),s),
                separation_over_diameter=math.sqrt(float(min(distances)/max(distances))))


def check_nonaffine_move():
    source = HERE.parent / "v13_followup/spatial_0_144_best.json"
    _, initial = ev.load_polyhedron_json(str(source))
    _, final = ev.load_polyhedron_json(str(HERE / "moves_best.json"))
    fixed = [i for i, (a, b) in enumerate(zip(initial, final)) if a == b]
    moved = [i for i, (a, b) in enumerate(zip(initial, final)) if a != b]
    check(len(initial) == len(final) == 13 and moved, "Vertex correspondence")
    for indices in combinations(fixed, 4):
        a, b, c, d = (initial[i] for i in indices)
        det = dot(sub(b, a), cross(sub(c, a), sub(d, a)))
        if det:
            return dict(source=source.name, fixed=fixed, moved=moved,
                        affine_basis=list(indices), determinant=str(det))
    raise ValueError("No full-dimensional fixed affine frame")


def cross_checks():
    p=[tuple(map(int,r)) for r in json.loads((ROOT/"examples/V14_integer.json").read_text())["vertices"]]
    check(max(dot(z,z) for z in p)<=1024**2,"Exterior radius premise")
    tests=0
    # These 14 positive controls insert the missing original point and must
    # reproduce the known mono-unstable polyhedron at its exact centroid.
    for k in range(14):
        base=p[:k]+p[k+1:]
        D,C,fs=facets(base)
        t=p[k]
        check(check_box(base,D,C,fs,[],1,t,t) is None,"Positive control excluded")
    for a,b in combinations(range(14),2):
        base=[v for k,v in enumerate(p) if k not in (a,b)]
        D,C,fs=facets(base)
        for t in (p[a],tuple((x+y)//2 for x,y in zip(p[a],p[b])),
                  (45056-147*a, -164*b,31011),(-1666,2901,31*(a-b))):
            if t in base:
                continue
            delta=[max(0,dot(n,t)-o) for n,o,s in fs]
            full=D+sum(delta)
            add_c=tuple(Q(C[k]+sum(w*(s[k]+t[k]) for w,(_,_,s) in zip(delta,fs)),4*full)
                        for k in range(3))
            pts=[tuple(map(Q,v)) for v in base+[t]]
            h=ev.convex_hull(pts)
            vol,c=ev.volume_and_centroid(pts,h)
            check(6*vol==full and c==add_c,"Mismatch in volume or centroid")
            # Also test a box of nonzero width containing the independently
            # classified sample (an inconclusive box is allowed).
            state=check_box(base,D,C,fs,[],1,tuple(v-2 for v in t),tuple(v+3 for v in t))
            if state=="excluded":
                vs=ev.classify_vertices(pts,h,c)
                check(sum(v.status==ev.UNSTABLE for v in vs)>=2,"Incorrect exclusion")
            if state=="inside12":
                check(len(h.vertices)<=12,"Incorrect redundancy")
            tests+=1
    return tests


if __name__=="__main__":
    paths=[HERE.parent/"v13_followup/focused_144_192_best.json",HERE/"continuation_count.json",
           HERE/"moves_best.json",HERE/"moves_sep0.021_best.json",HERE/"reroute_best.json",HERE/"lifting_best.json"]
    out=dict(positive_controls=14,centroids_cross_checked=cross_checks(),
             nonaffine_move=check_nonaffine_move(),candidates=[score(p) for p in paths])
    (HERE/"audit.json").write_text(json.dumps(out,indent=1)+"\n")
    print("PASS:",out["positive_controls"],"positive controls;",out["centroids_cross_checked"],"exact centroid cross-checks")
    for r in out["candidates"]:
        print(r["file"],r["counts"],r["angular_margin"],r["separation_over_diameter"])
