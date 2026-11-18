"""Free replacement and coordinate searches with the exact geometric formula for an added point."""
import os
for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[k]="1"
from concurrent.futures import ProcessPoolExecutor
import argparse
from itertools import combinations
import json
from pathlib import Path
import sys
import numpy as np
from numba import njit
from scipy.spatial import ConvexHull
from scipy.stats import qmc

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from search import cma_run
import core
from continuation import SOURCES,read,ev


def base_parameters(x):
    hull=ConvexHull(x)
    tri=hull.simplices
    a,b,c=x[tri].transpose(1,0,2)
    n=np.cross(b-a,c-a)
    n*=np.sign((n*hull.equations[:,:3]).sum(1))[:,None]
    offset=(n*a).sum(1)
    s=a+b+c
    return offset.sum(),offset@s,n,offset,s,tri


@njit(cache=True)
def scan(trials,base,D,C,n,off,s,tri,axes,focus,min_sep):
    out=np.empty(len(trials))
    x=np.empty((len(base)+1,3))
    x[:-1]=base*axes
    for k in range(len(trials)):
        t=trials[k]
        if np.max(np.abs(t))>6.:
            out[k]=-5.
            continue
        d=D
        cc=C.copy()
        hidden=np.ones(len(base),dtype=np.bool_)
        for f in range(len(n)):
            excess=n[f]@t-off[f]
            if excess>0:
                d+=excess
                cc+=excess*(t+s[f])
            else:
                for idx in tri[f]:
                    hidden[idx]=False
        if d<=D*(1+1.e-10) or np.any(hidden):
            out[k]=-4.
            continue
        cc=cc*axes/(4*d)
        x[-1]=t*axes
        m,parent,root=core.margins_c(x,cc)
        m[root]=2.
        pair,diam=core.shape_stats(x)
        if pair<min_sep:
            out[k]=-4.
            continue
        if focus:
            sorted_m=np.sort(m)
            # Protect the eleven other escape witnesses while improving
            # the second unstable vertex. Reward staying in the U=2 basin.
            out[k]=sorted_m[0]+10.*min(0.,sorted_m[1]-.0001)
        else:
            out[k]=m.min()
    return out


def vary(x,initial,rng,focus,min_sep):
    mean=x.mean(0)
    u,s,v=np.linalg.svd(x-mean,full_matrices=False)
    axes=s/s.max()
    base=u*np.sqrt(len(x))
    scale=s/np.sqrt(len(x))
    params=base_parameters(base)
    def objective(y):
        return scan(y,base,*params,axes,focus,min_sep)
    y0=(initial-mean)@v.T/scale
    unit=qmc.Sobol(3,scramble=True,seed=int(rng.integers(2**30))).random_base2(10)
    starts=np.vstack([y0,(unit-.5)*3.5,(unit-.5)*10., y0+(unit-.5)*.30])
    values=objective(starts)
    best=values.max()
    start=starts[values.argmax()]
    for scale_cma in (.30,.04):
        val,y,_=cma_run(objective,start,scale_cma,2200,seed=int(rng.integers(2**30)),stall=65)
        if val>best:
            start,best=y,val
    candidate=(start*scale)@v+mean if best>-3. else initial.copy()
    return candidate


def census(x):
    h=ConvexHull(x)
    if len(h.vertices)!=13:
        return (99,-3.)
    D,C,*_=base_parameters(x)
    m,_,r=core.margins_c(x,C/(4*D))
    u=int((m<0).sum())
    m[r]=2.
    return u,float(m.min())


def job(args):
    seed,min_sep=args
    rng=np.random.default_rng(221000+seed)
    pairs=list(combinations(range(14),2))
    if seed<len(pairs):
        i,j=pairs[seed]
        pts=read(SOURCES[0])/1000
        base=np.delete(pts,[i,j],axis=0)
        tip=vary(base,(pts[i]+pts[j])/2,rng,False,min_sep)
        x=np.vstack([base,tip])
        family="pair"
    else:
        paths=[HERE/"continuation_count.json",HERE.parent/"v13_followup/focused_144_192_best.json",
               HERE.parent/"v13_followup/spatial_0_144_best.json"]
        x=read(paths[(seed-91)%3])
        family="sweeps"
    initial=census(x)
    prev=initial
    prev_sep=core.shape_stats(x)[0]
    original=x.copy()
    for k,i in enumerate(rng.permutation(13) if family=="pair" else np.r_[rng.permutation(13),rng.permutation(13)]):
        proposal=x.copy()
        proposal[i]=vary(np.delete(x,i,axis=0),x[i],rng,True,min_sep)
        if np.linalg.norm(proposal[i]-original.mean(0))>3*max(np.linalg.norm(original-original.mean(0),axis=1)):
            continue
        new=census(proposal)
        sep=core.shape_stats(proposal)[0]
        repair=(min_sep>.012 and prev_sep<min_sep and sep>=min_sep and new[0]<99)
        if repair or ((min_sep<=.012 or sep>=min_sep) and (new[0],-new[1])<(prev[0],-prev[1])):
            x,prev=proposal,new
            prev_sep=sep
        if prev[0]==1:
            break
    pts=ev.rationalize(x,10**10)
    cert=ev.verify(pts)
    return dict(seed=seed,family=family,initial=initial,score=prev[1],counts=list(cert.counts),
                vertices=len(cert.hull.vertices),certified=cert.certified,
                separation=prev_sep,changed=int((np.linalg.norm(x-original,axis=1)>1e-9).sum()),xyz=x.tolist())


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--min-sep",type=float,default=.012)
    args=ap.parse_args()
    tag="moves" if args.min_sep==.012 else f"moves_sep{args.min_sep:g}"
    dummy=read(SOURCES[0])[:12].copy()
    scan(np.zeros((1,3)),dummy,*base_parameters(dummy),np.ones(3),True,args.min_sep)
    out=[]
    rank=(100,3.)
    with ProcessPoolExecutor(6) as pool:
        for r in pool.map(job,[(i,args.min_sep) for i in range(139)]):
            rr=(r["counts"][2],-r["score"])
            if r["certified"] and r["vertices"]==13 and r["separation"]>=args.min_sep-1e-10 and rr<rank:
                rank=rr
                pts=ev.rationalize(r["xyz"],10**10)
                ev.write_polyhedron_json(str(HERE/(tag+"_best.json")),pts,"Independent vertex relocations")
                (HERE/(tag+"_best_report.json")).write_text(json.dumps(ev.verify(pts).to_json(),indent=1)+"\n")
                print(r["seed"],r["counts"],r["score"],r["changed"],flush=True)
            r.pop("xyz")
            out.append(r)
    (HERE/(tag+"_runs.json")).write_text(json.dumps(out,indent=1)+"\n")
