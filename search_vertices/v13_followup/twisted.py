"""Two arms exchanged by a half turn, with free spatial twists and a tight radial recurrence."""
import os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key] = "1"
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys
import numpy as np
from numba import njit

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import core
from search import cma_run


@njit(cache=True)
def make(p, parents, t):
    n = len(parents)
    x = np.zeros((2*n+1,3))
    x[0,1] = 1.
    ph = np.zeros(n+1)
    rad = np.ones(n+1)
    penalty = 0.
    for i in range(1,n+1):
        j = parents[i-1]
        ph[i] = ph[j]+p[i-1]
        theta = 0. if i==1 else p[n+i-2]
        u = np.array([np.sin(ph[i])*np.cos(t*theta),np.cos(ph[i]),np.sin(ph[i])*np.sin(t*theta)/t])
        rad[i] = x[j,0]*u[0]+x[j,1]*u[1]+t*t*x[j,2]*u[2]
        x[i] = rad[i]*u
        x[i+n] = x[i]
        x[i+n,0] *= -1
        x[i+n,2] *= -1
        penalty += max(0.,.08-p[i-1])+max(0.,p[i-1]-1.3)+max(0.,.015-rad[i])
        penalty += max(0.,.001-(x[j,1]-x[i,1]))
    return x, penalty


@njit(cache=True)
def evaluate(candidates,parents,t):
    result = np.empty(len(candidates))
    for i in range(len(candidates)):
        x,penalty = make(candidates[i],parents,t)
        if penalty>0:
            result[i] = -1-penalty
        else:
            # The centroid component along the symmetry axis is invariant
            # under rescaling the transverse axis. Condition before joggling.
            height = np.abs(x[:,2]).max()
            if height < 1e-13:
                result[i] = -1.
                continue
            x[:,2] /= height
            volume,c,count,vertices,ok = core.hull_centroid(x,core.JOG)
            result[i] = -c[1] if (ok and vertices.all() and volume>1e-7 and
                                   x[:,1].min()<c[1]<1 and abs(c[0])+abs(c[2])<1e-7) else -1.
    return result


def job(seed):
    rng = np.random.default_rng(seed)
    parents = np.arange(6)
    if seed%4==2:
        parents[-1] = 4
    if seed%4==3:
        parents[-2] = 3
    t = [1e-4,.05,.25,.75][(seed//4)%4]
    x0 = np.r_[rng.uniform(.36,.50,6),rng.uniform(-1,1,5)]
    best,p,ne = cma_run(lambda z:evaluate(z,parents,t),x0,.16,38000,seed=1000+seed,stall=280)
    x,_ = make(p,parents,t)
    x[:,2] *= t
    return dict(seed=seed,scale=t,parents=parents.tolist(),score=best,params=p.tolist(),evals=ne,xyz=x.tolist())


if __name__=="__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--count",type=int,default=48)
    a = ap.parse_args()
    # Compile once before starting child processes.
    evaluate(np.array([[.44]*6+[.1,-.5,.6,-.5,.2]]),np.arange(6),.1)
    with ProcessPoolExecutor(4) as pool:
        runs = list(pool.map(job,range(a.count)))
    best = max(runs,key=lambda r:r["score"])
    sys.path.insert(0,str(HERE.parents[1]/"verify"))
    import exact_verify as ev
    from fractions import Fraction
    for r in runs:
        # Decimal reconstruction retains thin coordinates and the half-turn.
        pts = [[Fraction(format(z,".16g")) for z in row] for row in r["xyz"]]
        result = ev.verify(pts)
        r["exact_centroid_y"] = str(result.centroid[1])
        if abs(float(result.centroid[1])+r["score"])>1e-6 and r["score"]>-0.99:
            raise ValueError(f"Seed {r['seed']} failed the independent centroid check")
        r["exact_U"] = result.counts[2]
        r["exact_V"] = len(result.hull.vertices)
        if result.mono_unstable and r["exact_V"]==13:
            ev.write_polyhedron_json(str(HERE/f"twisted_certified_{r['seed']}.json"),pts)
    (HERE/"twisted_best.json").write_text(json.dumps(best,indent=1)+"\n")
    (HERE/"twisted_runs.json").write_text(json.dumps([{k:v for k,v in r.items() if k not in ("params","xyz")} for r in runs],indent=1)+"\n")
    print(best["score"],"at seed",best["seed"],"total evaluations",sum(r["evals"] for r in runs))
