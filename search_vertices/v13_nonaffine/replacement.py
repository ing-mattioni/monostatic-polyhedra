"""Exact box search: retain twelve original vertices and introduce one freely positioned point."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/"verify"))
import exact_verify as ev


def dot(x,y):
    return sum(i*j for i,j in zip(x,y))


def sub(x,y):
    return tuple(a-b for a,b in zip(x,y))


def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])


def triangles(p):
    h=ev.convex_hull([tuple(map(Fraction,v)) for v in p])
    out=[]
    for f in h.faces:
        a=p[f.vertices[0]]
        for j,k in zip(f.vertices[1:],f.vertices[2:]):
            b,c=p[j],p[k]
            out.append((a,b,c))
    return out


def facets(p):
    D=0
    C=[0,0,0]
    faces=[]
    for a,b,c in triangles(p):
        n=cross(sub(b,a),sub(c,a))
        d=dot(a,cross(b,c))
        offset=dot(n,a)
        if max(dot(n,v)-offset for v in p)>0:
            raise ValueError("Wrong facet orientation")
        s=tuple(sum(z) for z in zip(a,b,c))
        D+=d
        for l in range(3):
            C[l]+=d*s[l]
        faces.append((n,offset,s))
    if D<=0:
        raise ValueError("Nonpositive starting volume")
    return D,C,faces


def interval(n,lo,hi):
    return (sum(a*(l if a>=0 else h) for a,l,h in zip(n,lo,hi)),
            sum(a*(h if a>=0 else l) for a,l,h in zip(n,lo,hi)))


def multiply(a,b,c,d):
    products=(a*c,a*d,b*c,b*d)
    return min(products),max(products)


def check_box(p,D,C,faces,outer,M,lo,hi):
    # Coordinates use denominator M. D and C are the original integer
    # determinant sum and first moment, before adding any visible pyramids.
    if any(interval(n,lo,hi)[0]>offset*M for n,offset,_ in outer):
        return "outside"
    dl=dh=D*M
    # Center the extra moment about c0 to reduce interval dependence.
    cl=[0]*3
    ch=[0]*3
    sl=[None]*3
    sh=[None]*3
    visible=False
    for n,offset,s in faces:
        a,b=interval(n,lo,hi)
        a,b=max(0,a-offset*M),max(0,b-offset*M)
        if b==0:
            continue
        visible=True
        dl+=a
        dh+=b
        for k in range(3):
            sl[k]=s[k] if sl[k] is None else min(sl[k],s[k])
            sh[k]=s[k] if sh[k] is None else max(sh[k],s[k])
            constant=M*(D*s[k]-C[k])
            l,h=multiply(a,b,constant+D*lo[k],constant+D*hi[k])
            cl[k]+=l
            ch[k]+=h
    if not visible:
        return "inside12"
    for k in range(3):
        # Scaled centroid = M*C0/(4*D0) + extra/(4*D0*augmented_D).
        ld=dl if cl[k]<0 else dh
        hd=dl if ch[k]>0 else dh
        cl[k]=(M*C[k]*ld+cl[k])//(4*D*ld)
        ch[k]=-((-M*C[k]*hd-ch[k])//(4*D*hd))
        # Independent enclosure: 4*c = C/D + eta*(t+sbar-C/D),
        # eta=1-D/D_aug, with sbar a convex average of visible facet sums.
        l=M*(D*sl[k]-C[k])+D*lo[k]
        h=M*(D*sh[k]-C[k])+D*hi[k]
        ld=dh if l<0 else dl
        hd=dh if h>0 else dl
        cl[k]=max(cl[k],(M*C[k]*ld+l*(ld-D*M))//(4*D*ld))
        ch[k]=min(ch[k],-((-M*C[k]*hd-h*(hd-D*M))//(4*D*hd)))
    # An arc i->j is possible only if the escape product could be >=0.
    # A closed group contains a strictly unstable point (its farthest).
    # Two closed sink components therefore exclude a unique unstable point.
    n=len(p)+1
    reach=[1<<i for i in range(n)]
    for i,v in enumerate(p):
        gl=sub(v,ch)
        gh=sub(v,cl)
        for j,w in enumerate(p):
            if i!=j and interval(sub(w,v),gl,gh)[1]>=0:
                reach[i]|=1<<j
        if sum(multiply(lo[k]-v[k],hi[k]-v[k],gl[k],gh[k])[1] for k in range(3))>=0:
            reach[i]|=1<<(n-1)
        if sum(multiply(v[k]-hi[k],v[k]-lo[k],lo[k]-ch[k],hi[k]-cl[k])[1] for k in range(3))>=0:
            reach[-1]|=1<<i
    for k in range(n):
        bit=1<<k
        for j in range(n):
            if reach[j]&bit:
                reach[j]|=reach[k]
    common=reach[0]
    for mask in reach[1:]:
        common &= mask
    return "excluded" if not common else None


def experiment(job):
    pair,depth,budget,cube=job
    vertices=[tuple(map(int,p)) for p in json.loads((ROOT/"examples/V14_integer.json").read_text())["vertices"]]
    p=[v for i,v in enumerate(vertices) if i not in pair]
    D,C,faces=facets(p)
    outer=facets(vertices)[2]
    M=2**(depth+4)
    ps=[tuple(M*z for z in v) for v in p]
    lo0,hi0=[tuple(fn(a)*M for a in zip(*vertices)) for fn in (min,max)]
    if cube:
        lo0,hi0=(-cube*M,)*3,(cube*M,)*3
        outer=[]
    stack=[(lo0,hi0,0)]
    counts=dict(outside=0,inside12=0,excluded=0,unresolved=0)
    nodes=0
    first=None
    while stack:
        lo,hi,level=stack.pop()
        nodes+=1
        state=check_box(ps,D,C,faces,outer,M,lo,hi)
        if state:
            counts[state]+=1
        elif level==depth or nodes>=budget:
            counts["unresolved"]+=1
            if first is None:
                first=[list(lo),list(hi),M]
        else:
            # Bisect along the longest axis relative to its original range.
            if cube:
                # Fine layers near the original sheet, coarse far from it.
                widths=(1400,1500,25)
                axis=max(range(3),key=lambda k:Fraction(hi[k]-lo[k],widths[k]*M+max(0,lo[k],-hi[k])))
            else:
                axis=max(range(3),key=lambda k:Fraction(hi[k]-lo[k],hi0[k]-lo0[k]))
            mid=(hi[axis]+lo[axis])//2
            a,b=list(lo),list(hi)
            a[axis]=b[axis]=mid
            stack.append((tuple(a),hi,level+1))
            stack.append((lo,tuple(b),level+1))
    return dict(removed=pair,nodes=nodes,leaves=counts,first_unresolved=first)


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--depth",type=int,default=24)
    ap.add_argument("--budget",type=int,default=3000)
    ap.add_argument("--only",type=int,nargs=2)
    ap.add_argument("--cube",type=int,default=0)
    ap.add_argument("--certify",action="store_true",help="Require complete coverage and the analytic exterior cutoff")
    ap.add_argument("--workers",type=int,default=6)
    a=ap.parse_args()
    pairs=[tuple(a.only)] if a.only else list(combinations(range(14),2))
    if a.certify:
        vertices=json.loads((ROOT/"examples/V14_integer.json").read_text())["vertices"]
        if len(vertices)!=14 or a.only or a.cube<41*1024 or any(
                sum(int(z)**2 for z in v)>1024**2 for v in vertices):
            raise ValueError("Domain or coordinate bound does not cover the stated family")
    out=[]
    with ProcessPoolExecutor(a.workers) as pool:
        for r in pool.map(experiment,[(p,a.depth,a.budget,a.cube) for p in pairs]):
            out.append(r)
    stem=f"replacement_{a.depth}_{a.budget}"+(f"_cube{a.cube}" if a.cube else "")+(f"_{pairs[0][0]}_{pairs[0][1]}" if a.only else "")
    (HERE/(stem+".json")).write_text(json.dumps(out,indent=1)+"\n")
    print("Completed",sum(r["leaves"]["unresolved"]==0 for r in out),"of",len(out),"domains;",
          sum(r["nodes"] for r in out),"boxes")
    if a.certify and (len(out)!=91 or any(r["leaves"]["unresolved"] for r in out)):
        raise ValueError("Incomplete certificate")
