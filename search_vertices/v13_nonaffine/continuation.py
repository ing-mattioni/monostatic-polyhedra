"""Retract an auxiliary vertex while independently repositioning the other thirteen."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import json
import os
from pathlib import Path
import sys

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
from scipy.spatial import ConvexHull,QhullError

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0,str(HERE.parent/"v13_followup"))
from spatial import mass, cosines, optimize, ev

SOURCES = [ROOT/"examples/V14_integer.json", HERE.parent/"certified/V14_resh_d6.json",
           HERE.parent/"certified/V14_spiral76_d7.json"]


def read(path):
    return np.array([[float(Fraction(k)) for k in p] for p in
                     json.loads(path.read_text())["vertices"]])


def chart(points):
    u,s,v = np.linalg.svd(points-points.mean(0),full_matrices=False)
    return u*np.sqrt(len(points)),s/s.max()


def augmented_mass(weights,offset,t):
    if t==0:
        return mass
    def geometry(x):
        tip = weights@x+t*offset
        c,dc,hull = mass(np.vstack([x,tip]))
        return c,dc[:-1]+weights[:,None,None]*dc[-1],hull
    return geometry


def tree(x,axes,c):
    radii = np.linalg.norm((x-c)*axes,axis=1)
    root = int(radii.argmax())
    angular = cosines(x*axes,c*axes)
    angular[radii[:,None]>=radii[None,:]] = -3.
    return root,angular.argmax(1)


def evaluate(x,axes):
    c,_,h = mass(x)
    if len(h.vertices)!=len(x):
        return None
    p = x*axes
    m = cosines(p,c*axes).max(1)
    far = int(np.linalg.norm((x-c)*axes,axis=1).argmax())
    return int((m<0).sum()),float(np.min(np.delete(m,far)))


def keep(record,x,axes,stage):
    value = evaluate(x,axes)
    if value is None:
        return
    u,score = value
    if (u,-score) < (record.get("U",100),-record.get("score",-3)):
        record.update(U=u,score=score,xyz=(x*axes).tolist(),stage=stage)


def run(seed):
    rng = np.random.default_rng(18000+seed)
    source,removed,mode = seed//14%3,seed%14,seed//42
    full,axes = chart(read(SOURCES[source]))
    x = np.delete(full,removed,axis=0)
    tip = full[removed]
    # Retract towards an interior point near the discarded vertex.
    dist = np.linalg.norm((x-tip)*np.maximum(axes,.12),axis=1)
    w = .08*np.ones(13)/13
    close = np.argsort(dist)[:3]
    w[close] += .92*(1/dist[close])/(1/dist[close]).sum()
    offset = tip-w@x
    rec = dict(seed=seed,source=SOURCES[source].name,removed=removed,mode=mode)
    keep(rec,x,axes,"start")
    for k,t in enumerate((1.,.75,.48,.25,.09,0.,0.)):
        geom = augmented_mass(w,offset,t)
        if mode in (2,3) and k in (1,3):
            # Bend near the missing cap, alternating the sign of the transverse
            # displacements so neighbouring layers can exchange order.
            x = x.copy()
            local = np.argsort(dist)[:5]
            x[local] += rng.normal(0,[.06,.06,.4],size=(5,3))
        c = geom(x)[0]
        root,parent = tree(x,axes,c)
        angle = cosines(x*axes,c*axes).max(1)
        angle[root] = 2.
        worst = int(angle.argmin()) if mode%2 else None
        trial = optimize(x,axes,root,parent,65,worst,geometry=geom)
        keep(rec,trial,axes,f"{k}:{t}")
        if rec.get("U")==1 and rec["score"]>0:
            break
        if len(ConvexHull(trial).vertices)==13:
            x = trial
        else:
            for scale in (.75,.5,.25):
                mix = x+scale*(trial-x)
                if len(ConvexHull(mix).vertices)==13:
                    x=mix
                    break
    return rec


def save_exact(r,tag):
    pts = ev.rationalize(r["xyz"],10**10)
    result = ev.verify(pts)
    r.update(exact_V=len(result.hull.vertices), exact_counts=list(result.counts),
             exact_certified=result.certified)
    if tag:
        ev.write_polyhedron_json(str(HERE/(tag+".json")),pts,tag)
        (HERE/(tag+"_report.json")).write_text(json.dumps(result.to_json(),indent=1)+"\n")
    return result


if __name__=="__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--count",type=int,default=168)
    ap.add_argument("--workers",type=int,default=6)
    a=ap.parse_args()
    out=[]
    ranks={}
    with ProcessPoolExecutor(a.workers) as pool:
        for r in pool.map(run,range(a.count)):
            result=save_exact(r,None)
            if r["exact_V"]==13 and result.certified:
                for tag,rank in (("continuation_count",(-result.counts[2],r["score"])),
                                 ("continuation_margin",(r["score"],))):
                    if tag not in ranks or rank>ranks[tag]:
                        ranks[tag]=rank
                        save_exact(r,tag)
                        print(tag,r["seed"],r["stage"],result.counts,r["score"],flush=True)
            r.pop("xyz")
            out.append(r)
    (HERE/"continuation_runs.json").write_text(json.dumps(out,indent=1)+"\n")
