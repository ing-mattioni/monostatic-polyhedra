"""Optimize all heights by a linear program for a fixed polygon and two hull triangulations.

The volume is fixed to one and the three affine height freedoms are removed.
The resulting horizontal witnesses suffice for a thin body when all are positive.
"""
from concurrent.futures import ProcessPoolExecutor
import argparse
import json
import os
from pathlib import Path
import sys
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[key]="1"
import numpy as np
from scipy.optimize import linprog, minimize
from scipy.spatial import ConvexHull

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from continuation import read, ev, evaluate, SOURCES


def outline(angles, left):
    arms=[]
    for sign,steps in ((-1,angles[:left]),(1,angles[left:])):
        phi=np.cumsum(steps)
        radius=np.cumprod(np.cos(steps))
        arms.append(radius[:,None]*np.c_[sign*np.sin(phi),np.cos(phi)])
    return np.vstack([[0.,1.],arms[0],arms[1][::-1]])


def triangulations(q,h):
    hull=ConvexHull(np.c_[q,h])
    signs=hull.equations[:,2]
    if min(abs(signs))<1e-10 or len(ConvexHull(q).vertices)!=len(q):
        return None
    return [hull.simplices[signs*side>0] for side in (1,-1)]


def coefficients(q,triangles):
    n=len(q)
    volume=np.zeros(n)
    moment=np.zeros((2,n))
    creases=[]
    for side,tri in zip((1,-1),triangles):
        tri_q=q[tri]
        cross=np.linalg.det(tri_q[:,1:]-tri_q[:,:1])
        area=abs(cross)/2
        if min(area)<1e-8:
            return None
        common={}
        for f,coords,A in zip(tri,tri_q,area):
            volume[f]+=side*A/3
            moment[:,f]+=side*A/12*(coords+coords.sum(0)).T
            for i in range(3):
                edge=tuple(sorted(np.delete(f,i)))
                if edge not in common:
                    common[edge]=f[i]
                    continue
                j=common[edge]
                bary=np.linalg.solve(np.c_[coords,np.ones(3)].T,np.r_[q[j],1.])
                row=np.zeros(n)
                row[j]=side
                row[f]=-side*bary
                creases.append(row)
    return volume,moment,np.array(creases)


def solve(q,triangles,parent,root,focus=None):
    n=len(q)
    edges=[i for i in range(n) if i!=root]
    coeff=coefficients(q,triangles)
    if coeff is None:
        return None
    volume,moment,creases=coeff
    vec=q[parent[edges]]-q[edges]
    vec/=np.linalg.norm(vec,axis=1)[:,None]
    vec/=np.maximum(np.linalg.norm(q[edges],axis=1),.05)[:,None]
    active=np.ones(n-1) if focus is None else (np.asarray(edges)==focus).astype(float)
    rhs=np.r_[np.zeros(len(creases)),np.sum(vec*q[edges],axis=1)-.0004*(1-active)]
    constraints=np.vstack([np.c_[creases,np.zeros(len(creases))],np.c_[vec@moment,active]])
    # Select a well-conditioned triple to define the zero plane.
    gauge=[0,n//3,2*n//3]
    bounds=[(0.,0.) if i in gauge else (None,None) for i in range(n+1)]
    sol=linprog(np.r_[np.zeros(n),-1.],A_ub=constraints,b_ub=rhs,
                A_eq=[np.r_[volume,0.]],b_eq=[1.],bounds=bounds,method="highs")
    if not sol.success:
        return None
    h=sol.x[:-1]
    return sol.x[-1],h,moment@h


def job(seed):
    rng=np.random.default_rng(25000+seed)
    left=[6,5,7,4,8,3,9][seed%7]
    raw=np.r_[np.full(left,2.60/left),np.full(12-left,2.60/(12-left))]
    angles=raw*rng.uniform(.85,1.15,12)
    q=outline(angles,left)
    par=np.r_[0,np.arange(left),np.arange(left+2,13),0]
    best=None
    programs=0
    # Each trial releases all heights, beginning at a different pair of upper
    # and lower triangulations. Iterate near the preceding optimal lift.
    for k in range(40):
        h=rng.normal(0,1,13)
        if k%3 and best is not None:
            h=best[1]/max(abs(best[1]))+h*rng.choice([.025,.12,.4])
        tri=triangulations(q,h)
        if tri is None:
            continue
        r=solve(q,tri,par,0)
        programs+=1
        if r is not None and (best is None or r[0]>best[0]):
            best=(*r,tri)
    if best is None:
        return dict(seed=seed,programs=programs)
    # Move all twelve planar steps independently. Maintain the chosen
    # triangulations while recomputing the optimum of the height LP.
    tri=best[3]
    def objective(a):
        nonlocal programs
        qq=outline(a,left)
        if max(a[:left].sum(),a[left:].sum())>3.12:
            return .3+sum(abs(a-angles))
        turn=np.roll(qq,-1,axis=0)-qq
        nxt=np.roll(turn,-1,axis=0)
        if np.any(turn[:,0]*nxt[:,1]-turn[:,1]*nxt[:,0]<=1e-6):
            return .3+sum(abs(a-angles))
        programs+=1
        r=solve(qq,tri,par,0)
        return .3 if r is None else -r[0]
    fit=minimize(objective,angles,method="SLSQP",bounds=[(.10,1.32)]*12,
                 constraints=[dict(type="ineq",fun=lambda a:np.array([3.10-sum(a[:left]),3.10-sum(a[left:])]))],
                 options=dict(maxiter=24,ftol=2e-8,eps=1e-5))
    if fit.fun < -best[0]:
        qq=outline(fit.x,left)
        r=solve(qq,tri,par,0)
        if r is not None:
            q,best=qq,(*r,tri)
    # Try switching the triangulations at the optimized polygon, and then
    # protect 11 witnesses while concentrating the deficit on one vertex.
    for k in range(20):
        tri=triangulations(q,best[1]/max(abs(best[1]))+rng.normal(0,.10,13))
        if tri is None:
            continue
        r=solve(q,tri,par,0)
        programs+=1
        if r is not None and r[0]>best[0]:
            best=(*r,tri)
    choices=[best[:3]]
    for i in range(1,13):
        extra=solve(q,best[3],par,0,i)
        programs+=1
        if extra is not None:
            choices.append(extra)
    records=[]
    for lp,h,c in choices:
        h=h/max(abs(h))
        x=np.c_[q,.002*h]
        val=evaluate(x,np.ones(3))
        if val is not None:
            records.append((val,x))
    if not records:
        return dict(seed=seed,programs=programs)
    (u,margin),x=min(records,key=lambda item:(item[0][0],-item[0][1]))
    pts=ev.rationalize(x,10**10)
    cert=ev.verify(pts)
    return dict(seed=seed,left=left,programs=programs,score=float(margin),flat_lp=float(best[0]),
                counts=list(cert.counts),verts=len(cert.hull.vertices),certified=cert.certified,xyz=x.tolist())


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--count",type=int,default=84)
    a=ap.parse_args()
    # Positive control: select the known upper and lower triangulations,
    # using exact-input coordinates in their published order.
    p=read(SOURCES[0])/1000
    parents=np.array([0,0,1,2,3,4,5,8,9,10,11,12,13,0])
    control=solve(p[:,:2],triangulations(p[:,:2],p[:,2]),parents,0)
    if control is None or control[0]<=0:
        raise RuntimeError("Failed V14 control")
    print("Known 14-vertex control",control[0],flush=True)
    records=[]
    best=(99,2.)
    with ProcessPoolExecutor(6) as pool:
        for r in pool.map(job,range(a.count)):
            if "xyz" in r:
                rank=(r["counts"][2],-r["score"])
                if r["certified"] and r["verts"]==13 and rank<best:
                    best=rank
                    pts=ev.rationalize(r["xyz"],10**10)
                    ev.write_polyhedron_json(str(HERE/"lifting_best.json"),pts,"LP lifting search")
                    print(r["seed"],r["counts"],r["score"],r["flat_lp"],flush=True)
                r.pop("xyz")
            records.append(r)
    (HERE/"lifting_runs.json").write_text(json.dumps(dict(control=control[0],runs=records),indent=1)+"\n")
