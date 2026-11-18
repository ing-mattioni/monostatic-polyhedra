"""Replay every affine deletion certificate using only Python integer arithmetic."""
import gzip
from itertools import combinations
import json
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])


def geometry(p):
    facets = {}
    for i,j,k in combinations(range(len(p)),3):
        n = cross(sub(p[j],p[i]),sub(p[k],p[i]))
        side = [dot(n,sub(v,p[i])) for v in p]
        if max(side)==0 and min(side)<0:
            facets[tuple(h for h,v in enumerate(side) if v==0)] = n
        if min(side)==0 and max(side)>0:
            facets[tuple(h for h,v in enumerate(side) if v==0)] = tuple(-x for x in n)
    neighbours = [set() for _ in p]
    D, C = 0, [0,0,0]
    for face,normal in facets.items():
        # Orient the boundary of each supporting polygon by testing its
        # other points against every proposed directed edge.
        successor = {}
        for i in face:
            candidates = [j for j in face if i!=j and all(
                dot(normal,cross(sub(p[j],p[i]),sub(p[k],p[i])))>0
                for k in face if k not in (i,j))]
            require(len(candidates)==1,"Bad face boundary")
            j = successor[i] = candidates[0]
            neighbours[i].add(j)
            neighbours[j].add(i)
        order = [face[0]]
        for _ in face[1:]:
            order.append(successor[order[-1]])
        require(set(order)==set(face) and successor[order[-1]]==order[0],"Open face")
        a = p[order[0]]
        for j,k in zip(order[1:],order[2:]):
            b,c = p[j],p[k]
            d = dot(a,cross(b,c))
            D += d
            for l in range(3):
                C[l] += d*(a[l]+b[l]+c[l])
    require(D>0 and all(len(n)>=3 for n in neighbours),"Bad hull or volume")
    return [sorted(n) for n in neighbours],[tuple(4*D*z-y for z,y in zip(v,C)) for v in p]


def metric_row(g,d):
    r = [u*v for u,v in zip(g,d)]+[g[i]*d[j]+g[j]*d[i] for i,j in ((0,1),(0,2),(1,2))]
    factor = gcd(*r)
    return tuple(v//factor for v in r)


def replay(points, proof):
    dropped,frame = proof["deleted"],proof["frame"]
    require(isinstance(dropped,int) and 0<=dropped<14 and isinstance(frame,int),"Bad input chart")
    # Invertible triangular change of coordinates, determinant 100.
    points = [(x+frame*y,y+100*z if frame else y,100*z)
              for i,(x,y,z) in enumerate(points) if i!=dropped]
    neighbours,g = geometry(points)
    rows = {(i,j):metric_row(g[i],sub(points[j],points[i]))
            for i,others in enumerate(neighbours) for j in others}
    conflicts = []
    for terms in proof["terms"]:
        v = [0]*6
        keys = set()
        for i,j,w in terms:
            require((i,j) in rows and type(w) is int and w>0,"Invalid dual term")
            keys.add((i,j))
            for k,z in enumerate(rows[i,j]):
                v[k] -= w*z
        # Twice the NEGATIVE sum of the symmetric witness matrices. Its
        # nonzero positive semidefiniteness contradicts a positive metric.
        a,b,c = [2*z for z in v[:3]]
        d,e,f = v[3:]
        minors = (a,b,c,a*b-d*d,a*c-e*e,b*c-f*f,
                  a*b*c+2*d*e*f-a*f*f-b*e*e-c*d*d)
        require(a+b+c>0 and min(minors)>=0,"Dual sum fails principal-minor test")
        conflicts.append(keys)
    for label,stream in proof["trees"].items():
        root = int(label)
        require(0<=root<13,"Invalid exempt vertex")
        it = iter(stream)
        def verify(path):
            token = next(it,None)
            require(type(token) is int,"Incomplete tree")
            if token==-1:
                for v in path:
                    seen = set()
                    while v in path:
                        if v in seen:
                            return
                        seen.add(v)
                        v = path[v]
                raise ValueError("Missing directed cycle")
            if token < -1:
                k = -token-2
                require(k<len(conflicts) and conflicts[k]<=set(path.items()),"Wrong leaf certificate")
                return
            require(0<=token<13 and token!=root and token not in path,"Invalid split")
            for j in neighbours[token]:
                verify({**path,token:j})
        verify({})
        require(next(it,None) is None,"Surplus tree data")
    return {(dropped,int(k)) for k in proof["trees"]}


def main():
    raw = json.loads((ROOT/"examples/V14_integer.json").read_text())
    points = [tuple(map(int,v)) for v in raw["vertices"]]
    require(len(points)==14,"Unexpected source")
    covered = set()
    for file in sorted(HERE.glob("affine_*.proof.json.gz")):
        with gzip.open(file,"rt",encoding="utf-8") as source:
            covered.update(replay(points,json.load(source)))
    expected = {(i,j) for i in range(14) for j in range(13)}
    require(covered==expected,f"Missing cases: {sorted(expected-covered)}")
    print("PASS: 182/182 deletion/root cases. Every affine image has at least two strictly unstable vertices.")


if __name__=="__main__":
    main()
