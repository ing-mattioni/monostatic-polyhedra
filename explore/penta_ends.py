import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
sys.path.insert(0, 'explore')
from penta_scan import covectors, TREES, best_margin, FACES
paths = {
 "(e,e,e)":        lambda e: (e, e, e),
 "(e,2e,.5)":      lambda e: (e, 2*e, .5),
 "(e,.5,.5)":      lambda e: (e, .5, .5),
 "(e,e^2,.5)":     lambda e: (e, e*e, .5),
 "(1-e,.5,.5)":    lambda e: (1-e, .5, .5),
 "(1-e,1-2e,.5)":  lambda e: (1-e, 1-2*e, .5),
 "(1-e,1-2e,1-3e)":lambda e: (1-e, 1-2*e, 1-3*e),
 "(1-e,1-e,1-e)":  lambda e: (1-e, 1-e, 1-e),
 "(e,.5,1-e)":     lambda e: (e, .5, 1-e),
 "(e,1-e,1-2e)":   lambda e: (e, 1-e, 1-2*e),
 "(e,2e,1-e)":     lambda e: (e, 2*e, 1-e),
 "(1-e,1-e^2,.5)": lambda e: (1-e, 1-e*e, .5),
}
for name, f in paths.items():
    row = []
    for e in (1e-1, 1e-2, 1e-3, 1e-4):
        L, g = covectors(f(e))
        vals = [best_margin(L, tr) for tr in TREES]
        row.append(max(vals))
    print(f"{name:18s} " + "  ".join(f"{v:+.2e}" for v in row), flush=True)
