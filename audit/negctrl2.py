import sys, gzip, json, random
from recheck import traverse, check_cert
sys.path.insert(0, r"C:\Projects\Claude4\penta")
from core import BoxData, check
d = json.load(gzip.open("snap_certs_Z.json.gz", "rt"))
recs, st = traverse(d)
rng = random.Random(11)
tab = {}
for r in rng.sample(recs, 300):
    lo, hi, depth, (tid, y0, Y1), idx = r
    bd = BoxData("Z", [float(t) for t in lo], [float(t) for t in hi])
    for trial in range(3):
        t2 = rng.randrange(375)
        mine = check_cert("Z", lo, hi, t2, y0, Y1, maxsub=0)[0]
        theirs = check(bd, t2, y0, Y1)
        tab[(mine, theirs)] = tab.get((mine, theirs), 0) + 1
    # perturbed multipliers: shrink y0 margin
    for fac in (0.0, 0.5, 2.0):
        Y1p = [[t * (1 + fac * rng.uniform(-1, 1)) for t in row] for row in Y1]
        mine = check_cert("Z", lo, hi, tid, y0, Y1p, maxsub=0)[0]
        theirs = check(bd, tid, y0, Y1p)
        tab[("pert", mine, theirs)] = tab.get(("pert", mine, theirs), 0) + 1
print("(mine, core.check) counts for random other trees:", {k: v for k, v in tab.items() if len(k) == 2})
print("perturbed multipliers ('pert', mine, core):", {k: v for k, v in tab.items() if len(k) == 3})
