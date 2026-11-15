import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
sys.path.insert(0, 'explore')
from penta_scan import TREES, best_margin, FACES
def prism_covectors(s):
    s = np.asarray(s, float)
    B = np.array([[0,0],[1,0],[0,1]], float); A = 0.5
    S1 = s.sum(); 
    V = A*S1/3
    xbar = A*(S1*B[:,0].sum() + (s*B[:,0]).sum())/12/V
    ybar = A*(S1*B[:,1].sum() + (s*B[:,1]).sum())/12/V
    zbar = (A/2)*((s*s).sum() + s[0]*s[1]+s[0]*s[2]+s[1]*s[2])/6/V
    c = np.array([xbar, ybar, zbar])
    # planes n.x <= d  ->  covector l = n/(d - n.c)
    planes = []
    # Q1: face containing lateral edges 2,3 (opposite vertex 1=B[0]=(0,0)): x + y <= 1
    planes.append((np.array([1.,1.,0.]), 1.0))
    # Q2: opposite vertex 2=(1,0): x >= 0
    planes.append((np.array([-1.,0.,0.]), 0.0))
    # Q3: opposite vertex 3=(0,1): y >= 0
    planes.append((np.array([0.,-1.,0.]), 0.0))
    # A (top): z <= s1 + (s2-s1)x + (s3-s1)y
    planes.append((np.array([-(s[1]-s[0]), -(s[2]-s[0]), 1.]), s[0]))
    # B (bottom): -z <= 0
    planes.append((np.array([0.,0.,-1.]), 0.0))
    L = np.array([n/(d - n@c) for n, d in planes])
    return L
for s in [(1,1,1),(1,.5,.2),(1,.9,.1),(1,.5,.01),(1,.5,1e-4),(1,.01,.01),(1,1e-3,1e-4),(1,.999,.998),(1,1,.01),(1,.3,.3),(1,.7,.4),(1,.2,.05)]:
    L = prism_covectors(s)
    vals = [best_margin(L, tr) for tr in TREES]
    k = int(np.argmax(vals))
    print(f"prism s={s}: max SDP margin over trees = {vals[k]:+.3e}   root {FACES[TREES[k][0]]} {dict((FACES[i],FACES[j]) for i,j in TREES[k][1].items())}", flush=True)
