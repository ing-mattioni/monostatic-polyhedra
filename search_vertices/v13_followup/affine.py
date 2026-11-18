"""Search all positive-definite metrics on a fixed seed using witness disjunctions.

An excluded LP branch is accepted only after exact rational reconstruction
of a Farkas certificate. Unresolved branches are reported explicitly.
"""
import os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key] = "1"
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
from pathlib import Path
import argparse
import itertools
import gzip
import json
import math
import sys
import numpy as np
from scipy.optimize import linprog

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0,str(ROOT/"verify"))
import exact_verify as ev


def coefficients(a,b):
    return tuple([a[i]*b[i] for i in range(3)]+
                 [a[i]*b[j]+a[j]*b[i] for i,j in ((0,1),(0,2),(1,2))])


def primitive(row):
    den = math.lcm(*(Q(v).denominator for v in row))
    vals = [int(v*den) for v in row]
    g = math.gcd(*vals)
    return tuple(v//g for v in vals)


def reconstruct(cols, guess):
    keep = np.flatnonzero(guess>1e-10)
    # A basic feasible solution has at most seven nonzero coefficients.
    if not 1<=len(keep)<=7:
        return False
    a = [[Q(cols[j][i]) for j in keep]+[Q(int(i==6))] for i in range(7)]
    for k in range(len(keep)):
        pivot = next((j for j in range(k,7) if a[j][k]),None)
        if pivot is None:
            return False
        a[k],a[pivot] = a[pivot],a[k]
        a[k] = [v/a[k][k] for v in a[k]]
        for i in range(7):
            if i!=k and a[i][k]:
                a[i] = [x-a[i][k]*y for x,y in zip(a[i],a[k])]
    coeffs = [a[i][-1] for i in range(len(keep))]
    if all(v>=0 for v in coeffs) and all(
            sum(v*cols[j][i] for j,v in zip(keep,coeffs))==int(i==6) for i in range(7)):
        return {int(j):v for j,v in zip(keep, coeffs) if v}
    return False


def search(points, budget, roots):
    hull = ev.convex_hull(points)
    _,c = ev.volume_and_centroid(points,hull)
    rows = [[primitive(coefficients(ev.vsub(points[i],c),ev.vsub(points[j],points[i])))
             for j in hull.neighbours[i]] for i in range(len(points))]
    directions = [v for v in itertools.product((-1,0,1),repeat=3)
                  if any(v) and next(k for k in v if k)!=-1]
    cuts = list(dict.fromkeys(primitive(coefficients(v,v)) for v in directions))
    totals = dict(nodes=0, exact_leaves=0, cycle_leaves=0, unresolved=0, feasible=[], excluded_roots=[])
    # Share closed conjunctions across different choices of exempt root.
    excluded = {}
    proofs, trees, trace = [], {}, []
    def visit(root,chosen,depth):
        key = frozenset(chosen)
        for conflict in excluded:
            if conflict <= key:
                trace.append(-2-excluded[conflict])
                return True
        if chosen:
            neighbour = {i:hull.neighbours[i][j] for i,j in chosen}
            seen = set()
            v = chosen[-1][0]
            while v in neighbour:
                if v in seen:
                    # A nonnegative escape to a distinct vertex strictly
                    # increases squared distance in a positive metric.
                    totals["cycle_leaves"] += 1
                    trace.append(-1)
                    return True
                seen.add(v)
                v = neighbour[v]
        totals["nodes"] += 1
        if totals["nodes"]>budget:
            totals["unresolved"] += 1
            return False
        for attempt in range(20):
            all_rows = [rows[i][j] for i,j in chosen] + cuts
            M = np.array(all_rows,float)
            M /= np.max(abs(M),axis=1)[:,None]
            # Maximize a common margin on rank-one positive tests. The
            # chosen escape halfspaces need only be nonnegative.
            k = len(chosen)
            res = linprog([0]*6+[-1],A_ub=np.c_[-M,np.r_[np.zeros(k),np.ones(len(cuts))]],
                          b_ub=np.zeros(len(M)),A_eq=[[1,1,1,0,0,0,0]],b_eq=[1],
                          bounds=(None,None),method="highs")
            if res.status==2 or (res.success and res.x[-1]<1e-9):
                # sum lambda*A + sum mu*vv' = 0, sum mu=1. All weights
                # nonnegative. A positive-definite metric is impossible.
                eq = np.vstack([M.T,np.r_[np.zeros(k),np.ones(len(cuts))]])
                dual = linprog(np.zeros(len(M)),A_eq=eq,b_eq=[0]*6+[1],method="highs")
                if dual.success:
                    cols = [(*row, int(n>=k)*max(map(abs,row))) for n,row in enumerate(all_rows)]
                    # Row normalization rescales each nonnegative multiplier.
                    weights = reconstruct(cols,dual.x)
                    if weights:
                        totals["exact_leaves"] += 1
                        terms = [(chosen[l],w) for l,w in weights.items() if l<k]
                        ints = primitive([w for _,w in terms])
                        conflict = frozenset(ij for ij,_ in terms)
                        index = len(proofs)
                        excluded[conflict] = index
                        proofs.append([[i,hull.neighbours[i][j],w] for ((i,j),_),w in zip(terms,ints)])
                        trace.append(-2-index)
                        return True
            if not res.success:
                break
            m = res.x[:6]
            mat = np.array([[m[0],m[3],m[4]],[m[3],m[1],m[5]],[m[4],m[5],m[2]]])
            vals,vecs = np.linalg.eigh(mat)
            # Branch on an unstable vertex first, even for a relaxed metric.
            # Every possible witness neighbour is covered independently.
            missing = [i for i in range(len(points)) if i!=root and all(i!=h for h,_ in chosen)]
            if missing:
                def deficit(i):
                    A = np.array(rows[i],float)
                    return (A@m/np.max(abs(A),axis=1)).max()
                worst = min(missing,key=deficit)
                if deficit(worst)<-1e-7 or (vals[0]>1e-9 and deficit(worst)<1e-7):
                    trace.append(worst)
                    ok = True
                    for j in range(len(rows[worst])):
                        if not visit(root,chosen+[(worst,j)],depth+1):
                            ok = False
                    return ok
            if vals[0]>1e-9:
                totals["feasible"].append(dict(root=root,metric=m.tolist()))
                return False
            scale = 32*2**min(attempt,12)
            v = tuple(int(x) for x in np.rint(scale*vecs[:,0]))
            cut = primitive(coefficients(v,v))
            if len(cuts)>=55:
                break
            if cut not in cuts:
                cuts.append(cut)
        totals["unresolved"] += 1
        return False
    for root in roots:
        trace.clear()
        if visit(root,[],0):
            totals["excluded_roots"].append(root)
            trees[root] = trace.copy()
    return totals,dict(terms=proofs,trees=trees)


def deletion(args):
    i,budget,frame,roots = args
    _,pts = ev.load_polyhedron_json(str(ROOT/"examples/V14_integer.json"))
    # Fixed invertible change of units for conditioning, absorbable in G.
    pts = [(p[0],p[1],100*p[2]) for p in pts]
    if i>=0:
        pts.pop(i)
    if frame:
        pts = [(x+frame*y,y+z,z) for x,y,z in pts]
    r,proof = search(pts,budget,range(len(pts)) if roots is None else roots)
    proof.update(deleted=i,frame=frame)
    suffix = f"_frame{frame}" if frame else ""
    with gzip.open(HERE/f"affine_{i}{suffix}.proof.json.gz","wt",encoding="utf-8") as f:
        json.dump(proof,f,separators=(",",":"))
    print(i,r["nodes"],r["exact_leaves"],r["unresolved"],len(r["feasible"]),flush=True)
    return r


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--indices",type=int,nargs="+",default=list(range(14)))
    parser.add_argument("--budget",type=int,default=15000)
    parser.add_argument("--workers",type=int,default=6)
    parser.add_argument("--roots",type=int,nargs="+")
    parser.add_argument("--frame",type=int,default=0)
    args = parser.parse_args()
    with ProcessPoolExecutor(args.workers) as pool:
        suffix = f"_frame{args.frame}" if args.frame else ""
        for i,r in zip(args.indices,pool.map(deletion,[(i,args.budget,args.frame,args.roots) for i in args.indices])):
            (HERE/f"affine_{i}{suffix}.json").write_text(json.dumps(r,indent=1)+"\n")
