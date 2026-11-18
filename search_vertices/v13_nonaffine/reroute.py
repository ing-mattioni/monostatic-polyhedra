"""Refine two-unstable seeds after changing the escape tree and the shape together."""
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
from continuation import chart,read,mass,cosines,optimize,ev,np,HERE
from scipy.spatial import ConvexHull, QhullError


def valid(x,axes):
    physical=x*axes
    d=np.linalg.norm(physical[:,None]-physical[None,:],axis=2)
    d[np.diag_indices(len(x))]=np.inf
    return d.min()>=.012*np.max(d[np.isfinite(d)]) and len(ConvexHull(x).vertices)==13


def witness_tree(m,root,rng,noise):
    # A best local witness can point back along a cycle. Start from the
    # proposed root and attach one child at a time, allowing radii to reorder.
    n=len(m)
    parent=np.full(n,root)
    included=np.zeros(n,dtype=bool)
    included[root]=True
    priorities=m+rng.normal(0,noise,(n,n))
    for _ in range(n-1):
        available=np.where((~included[:,None]) & included[None,:],priorities,-np.inf)
        i,j=np.unravel_index(available.argmax(),available.shape)
        parent[i]=j
        included[i]=True
    return parent


def run(seed):
    rng=np.random.default_rng(180001+seed)
    names=[HERE/"moves_best.json",HERE/"continuation_count.json",
           HERE.parent/"v13_followup/focused_144_192_best.json"]
    x,axes=chart(read(names[seed%3]))
    base=x.copy()
    # Independently bend a local cluster, then allow every vertex to move.
    center=rng.integers(13)
    loc=np.argsort(np.linalg.norm((x-x[center])*axes,axis=1))[:4]
    x[loc]+=rng.normal(0,[.085,.085,.45],(4,3))
    if not valid(x,axes):
        x=base.copy()
    axes=axes.copy()
    axes[2]*=(1.,3.,10.,.3)[seed//3%4]
    best=None
    for k in range(4):
        c=mass(x)[0]
        m=cosines(x*axes,c*axes)
        weakest=np.argsort(m.max(1))
        root=weakest[int(seed%7==0 and k==0)]
        parent=witness_tree(m,root,rng,(.0,.025,.075,.15)[seed//12%4])
        other=next(int(i) for i in weakest if i!=root)
        try:
            y=optimize(x,axes,root,parent,100,other if k!=1 or seed%5 else None)
            for alpha in (1.,.85,.55,.3):
                z=(1-alpha)*x+alpha*y
                if not valid(z,axes):
                    continue
                c=mass(z)[0]
                scores=cosines(z*axes,c*axes).max(1)
                far=np.argmax(np.linalg.norm((z-c)*axes,axis=1))
                rank=(-int((scores<0).sum()),float(np.min(np.delete(scores,far))))
                if best is None or rank>best[0]:
                    best=rank,z.copy()
                x=z
                break
        except QhullError:
            x=base.copy()
        if best and best[0][0]==-1:
            break
    if best is None:
        return dict(seed=seed)
    rank,z=best
    p=ev.rationalize(z*axes,10**10)
    c=ev.verify(p)
    return dict(seed=seed,counts=list(c.counts),vertices=len(c.hull.vertices),certified=c.certified,
                score=rank[1],xyz=(z*axes).tolist())


if __name__=="__main__":
    best=(100,3.)
    out=[]
    with ProcessPoolExecutor(6) as pool:
        for r in pool.map(run,range(144)):
            if "xyz" in r:
                rank=r["counts"][2],-r["score"]
                if r["certified"] and r["vertices"]==13 and rank<best:
                    best=rank
                    pts=ev.rationalize(r["xyz"],10**10)
                    ev.write_polyhedron_json(str(HERE/"reroute_best.json"),pts,"Rerouted escape tree")
                    (HERE/"reroute_best_report.json").write_text(json.dumps(ev.verify(pts).to_json(),indent=1)+"\n")
                    print(r["seed"],rank,flush=True)
                r.pop("xyz")
            out.append(r)
    (HERE/"reroute_runs.json").write_text(json.dumps(out,indent=1)+"\n")
