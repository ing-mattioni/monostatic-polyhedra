"""Seeded search with free 3D vertices, prescribed escape trees and a volume-centroid Jacobian."""
import os
for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[key] = "1"

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import ConvexHull

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/"verify"))
import exact_verify as ev


def mass(x):
    hull = ConvexHull(x)
    tri = hull.simplices.copy()
    a,b,c = x[tri].transpose(1,0,2)
    normals = np.cross(b-a,c-a)
    sign = np.sign(np.einsum("ij,ij->i",normals,hull.equations[:,:3]))
    grad = np.stack([np.cross(b,c), np.cross(c,a), np.cross(a,b)],axis=1)*sign[:,None,None]
    det = np.einsum("ij,ij->i",a,grad[:,0,:])
    total = det.sum()
    means = (a+b+c)/4
    center = det@means/total
    dc = np.zeros((len(x),3,3))
    for k in range(3):
        deriv = (means-center)[:,:,None]*grad[:,k,None,:] + det[:,None,None]*np.eye(3)/4
        np.add.at(dc,tri[:,k],deriv/total)
    return center, dc, hull


def cosines(x, c):
    g = x-c
    d = x[None,:,:]-x[:,None,:]
    den = np.linalg.norm(g,axis=1)[:,None]*np.linalg.norm(d,axis=2)
    return np.divide(np.einsum("ik,ijk->ij",g,d),den,out=np.full(den.shape,-2.),where=den>0)


def optimize(x0, axes, root, parent, iterations, bottleneck=None, geometry=mass):
    idx = np.delete(np.arange(len(x0)), root)
    target = parent[idx]
    n = len(x0)
    cache = [None,None]
    def evaluate(z):
        if np.array_equal(z,cache[0]):
            return cache[1]
        x = z[:-1].reshape(n,3)
        c,dc,_ = geometry(x)
        g = (x[idx]-c)*axes
        d = (x[target]-x[idx])*axes
        rg = np.maximum(np.linalg.norm(g,axis=1),1e-12)
        rd = np.maximum(np.linalg.norm(d,axis=1),1e-12)
        angular = (d*g).sum(1)/(rg*rd)
        v = (d/(rd*rg)[:,None]-angular[:,None]*g/(rg**2)[:,None])*axes
        w = (g/(rd*rg)[:,None]-angular[:,None]*d/(rd**2)[:,None])*axes
        jac = -np.einsum("ik,jkl->ijl",v,dc)
        ii = np.arange(n-1)
        jac[ii,idx] += v-w
        jac[ii,target] += w
        # Prevent coincident escape edges. A positive cosine on a numerically
        # collapsed edge is meaningless. Bounds in lifted coordinates keep
        # the mass calculation well-conditioned even for flat physical shapes.
        sep = np.zeros((n-1,n,3))
        sep[ii,target] = d*axes/rd[:,None]
        sep[ii,idx] = -sep[ii,target]
        active = np.ones(n-1) if bottleneck is None else (idx==bottleneck).astype(float)
        vals = np.r_[angular-z[-1]*active-.0003*(1-active),rd-.045]
        J = np.c_[np.vstack((jac.reshape(n-1,-1),sep.reshape(n-1,-1))), np.r_[-active,np.zeros(n-1)]]
        cache[:] = [z.copy(), (vals,J)]
        return vals,J
    z = np.r_[x0.ravel(),-1.]
    ft = np.zeros_like(z)
    ft[-1] = -1.
    res = minimize(lambda z:-z[-1], z, jac=lambda z:ft, method="SLSQP",
                   constraints=[dict(type="ineq",fun=lambda z:evaluate(z)[0],jac=lambda z:evaluate(z)[1])],
                   bounds=[(v-.43,v+.43) for v in x0.ravel()]+[(-1.,1.)],
                   options=dict(ftol=2e-9,maxiter=iterations))
    return res.x[:-1].reshape(n,3)


def generate(seed):
    rng = np.random.default_rng(seed)
    mode = seed%6
    if mode in (0,1):
        data = json.loads((ROOT/"examples/V14_integer.json").read_text())
        x = np.delete(np.array(data["vertices"],float),seed//6%14,axis=0)
    elif mode == 2:
        files = sorted((ROOT/"search_vertices/sols").glob("V13_*.npy"))
        files = [f for f in files if not f.stem.endswith("_flat")]
        x = np.load(files[rng.integers(len(files))])
    else:
        # A branching tight-spiral tree in space. Children can change the
        # plane of their parent's spiral or branch from an internal vertex.
        narm = 2 if mode == 3 else 3 if mode == 4 else 1
        lengths = rng.multinomial(12-narm,np.ones(narm)/narm)+1
        x = [np.array([0.,0.,1.])]
        parents = [0]
        for arm, length in enumerate(lengths):
            start = 0 if arm==0 or rng.random()<.5 else rng.integers(len(x))
            u = x[start]/np.linalg.norm(x[start])
            direction = np.array([np.cos(arm*2*np.pi/narm),np.sin(arm*2*np.pi/narm),0.])
            v = direction-direction.dot(u)*u
            v /= np.linalg.norm(v)
            p = start
            for k in range(length):
                u = x[p]/np.linalg.norm(x[p])
                twist = rng.normal(0,.24)
                v = np.cos(twist)*v+np.sin(twist)*np.cross(u,v)
                ang = rng.uniform(.7,1.3)*min(.85,3/max(length,1))
                child = np.linalg.norm(x[p])*np.cos(ang)*(np.cos(ang)*u+np.sin(ang)*v)
                v = -np.sin(ang)*u+np.cos(ang)*v
                parents.append(p)
                x.append(child)
                p = len(x)-1
        x = np.asarray(x)
    # SVD only changes the frame. Assign a physical thickness independently
    # while retaining all three free coordinates in the optimization.
    u,s,_ = np.linalg.svd(x-x.mean(0),full_matrices=False)
    x = u*np.sqrt(len(x))/2
    thick = [0.012,.35, .12, .06, .22, .45][mode]*rng.uniform(.6,1.4)
    axes = np.array([1.,np.clip(s[1]/s[0],.30,1.),thick])
    x += rng.normal(0,.022 if mode<3 else .012,x.shape)
    if mode<3:
        c = mass(x)[0]
        ang = cosines(x*axes,c*axes)
        root = int(np.linalg.norm((x-c)*axes,axis=1).argmax())
        # Force an acyclic first tree by requiring strictly increasing radius.
        radius = np.linalg.norm((x-c)*axes,axis=1)
        ang[radius[None,:]<=radius[:,None]] = -3.
        parents = ang.argmax(1)
    else:
        root = 0
    return x, axes, root, np.array(parents)


def run(job):
    seed,focus = job
    x,axes,root,parents = generate(seed)
    best = None
    for _ in range(3):
        other = [i for i in range(len(x)) if i!=root]
        bottleneck = other[(seed//6)%len(other)] if focus else None
        x = optimize(x,axes,root,parents,85,bottleneck)
        c,_,hull = mass(x)
        physical = x*axes
        ang = cosines(physical,c*axes)
        # Inspect only actual extreme points for the recorded score.
        vertices = hull.vertices
        av = ang[np.ix_(vertices, vertices)].max(1)
        root = vertices[np.linalg.norm((x[vertices]-c)*axes,axis=1).argmax()]
        score = float(av[vertices!=root].min())
        unstable = int((av<0).sum())
        rank = (-unstable if focus else 0,score)
        if len(vertices)==13 and (best is None or rank>best[3]):
            best = score,physical.tolist(),unstable,rank
        if best is not None and best[0]>0:
            break
        parents = ang.argmax(1)
    if best is None:
        return dict(seed=seed)
    return dict(seed=seed, score=best[0], xyz=best[1], unstable=best[2])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--start",type=int,default=0)
    parser.add_argument("--count",type=int,default=120)
    parser.add_argument("--workers",type=int,default=6)
    parser.add_argument("--focus",action="store_true",help="Protect eleven vertices and optimize one escape")
    args = parser.parse_args()
    runs = []
    tag = ("focused" if args.focus else "spatial")+f"_{args.start}_{args.count}"
    best_rank = (-100,-3.)
    with ProcessPoolExecutor(args.workers) as pool:
        for r in pool.map(run,[(i,args.focus) for i in range(args.start,args.start+args.count)]):
            rank = (-r.get("unstable",99) if args.focus else 0,r.get("score",-3.))
            if "xyz" in r and rank>best_rank:
                best_rank = rank
                path = HERE/(tag+"_best.json")
                pts = ev.rationalize(r["xyz"],10**9)
                ev.write_polyhedron_json(str(path),pts,f"13-vertex search, seed {r['seed']}")
                exact = ev.verify(pts).to_json()
                print(r["seed"],r["score"],"exact V,S,H,U",*[exact[k] for k in ("V","S","H","U")],flush=True)
                (HERE/(tag+"_exact.json")).write_text(json.dumps(exact,indent=1)+"\n")
            r.pop("xyz",None)
            runs.append(r)
    (HERE/(tag+"_runs.json")).write_text(json.dumps(runs,indent=1)+"\n")
