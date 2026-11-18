"""Search non-affine projective images of all fourteen single-vertex deletions."""
import os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key] = "1"
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys
import numpy as np
from numba import njit
from scipy.spatial import ConvexHull

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(HERE.parent),str(ROOT/"verify")]
from search import cma_run
import core
import exact_verify as ev


@njit(cache=True)
def shape(x,p):
    direction = p[:3].copy()
    levels = x@direction
    # Keep all denominators positive. The map is invertible wherever
    # 1+h.x>0. Parameterization covers interior directions on both sides.
    levels /= 1.0-levels.min()
    den = 1.0+levels
    y = x/den.reshape((-1,1))
    # Triangular factor of the affine metric, up to overall scale.
    a = np.exp(max(-6.,min(3.,p[3])))
    b = np.exp(max(-9.,min(3.,p[4])))
    z = y.copy()
    z[:,0] = y[:,0]+p[5]*y[:,1]+p[6]*y[:,2]
    z[:,1] = a*(y[:,1]+p[7]*y[:,2])
    z[:,2] = b*y[:,2]
    return y,z,a,b


@njit(cache=True)
def evaluate(batch,x,tri):
    out = np.empty(len(batch))
    for k in range(len(batch)):
        p = batch[k]
        y,z,a,b = shape(x,p)
        # Origin is interior and remains fixed. Integrate the unscaled
        # projective image then send its volume centroid by the linear map.
        total = 0.
        c = np.zeros(3)
        for face in tri:
            i,j,l = face
            v = abs(y[i]@np.cross(y[j],y[l]))
            total += v
            c += v*(z[i]+z[j]+z[l])
        c /= 4*total
        value,root = core.margins_pen(z,c,0.03*core.shape_stats(z)[1],.2)
        sep,diam = core.shape_stats(z)
        out[k] = value - 2*max(0.,.025-sep) - max(0.,np.max(np.abs(p[:3]))-35)*.003
    return out


def inputs():
    raw = json.loads((ROOT/"examples/V14_integer.json").read_text())
    return np.asarray(raw["vertices"],dtype=float)


def seed_run(seed):
    dropped = seed%14
    rng = np.random.default_rng(16000+seed)
    x = np.delete(inputs(),dropped,axis=0)
    x -= x.mean(0)
    x /= np.max(np.abs(x),axis=0)
    tri = ConvexHull(x).simplices
    p = np.r_[rng.normal(0,.7,3),rng.uniform(-.6,.6),rng.uniform(-4,.3),rng.normal(0,.25,3)]
    f = lambda batch:evaluate(batch,x,tri)
    _,p,ne = cma_run(f,p,.5,26000,seed=seed+1,stall=230)
    _,z,_,_ = shape(x,p)
    pts = ev.rationalize(z,10**10)
    cert = ev.verify(pts)
    record = dict(seed=seed,deleted=dropped,params=p.tolist(),evaluations=ne,
                  vertices=len(cert.hull.vertices),counts=cert.counts,certified=cert.certified)
    if cert.mono_unstable:
        ev.write_polyhedron_json(str(HERE/f"projective_hit_{seed}.json"),pts)
    return record


if __name__=="__main__":
    pool_x = inputs()[:13].copy()
    pool_x -= pool_x.mean(0)
    evaluate(np.zeros((1,8)),pool_x,ConvexHull(pool_x).simplices)
    with ProcessPoolExecutor(6) as pool:
        results = list(pool.map(seed_run,range(84)))
    (HERE/"projective_runs.json").write_text(json.dumps(results,indent=1)+"\n")
    print("Runs",len(results),"evaluations",sum(r["evaluations"] for r in results),
          "minimum U",min(r["counts"][2] for r in results))
