import sys, gzip, json, random
sys.argv = ["x"]
from recheck import traverse, check_cert
d = json.load(gzip.open("snap_certs_Z.json.gz", "rt"))
recs, st = traverse(d)
rng = random.Random(7)
res = {"wrong_tree": 0, "scaled_Y1_x50": 0, "neg_y0": 0, "box_enlarged": 0, "original_ok": 0}
n = 200
for r in rng.sample(recs, n):
    lo, hi, depth, (tid, y0, Y1), idx = r
    res["original_ok"] += check_cert("Z", lo, hi, tid, y0, Y1, maxsub=0)[0]
    res["wrong_tree"] += check_cert("Z", lo, hi, (tid + 1) % 375, y0, Y1, maxsub=0)[0]
    res["scaled_Y1_x50"] += check_cert("Z", lo, hi, tid, y0, [[50 * t for t in row] for row in Y1], maxsub=0)[0]
    res["neg_y0"] += check_cert("Z", lo, hi, tid, [-t for t in y0], Y1, maxsub=0)[0]
    # enlarge the box 4x around its centre (clipped to [0,1]); same affine multiplier
    c = [(a + b) / 2 for a, b in zip(lo, hi)]; h = [(b - a) / 2 for a, b in zip(lo, hi)]
    lo2 = [max(0, c[k] - 4 * h[k]) for k in range(3)]; hi2 = [min(1, c[k] + 4 * h[k]) for k in range(3)]
    c2 = [(a + b) / 2 for a, b in zip(lo2, hi2)]
    y0b = [y0[i] + sum(Y1[i][k] * float(c2[k] - c[k]) for k in range(3)) for i in range(4)]
    res["box_enlarged"] += check_cert("Z", lo2, hi2, tid, y0b, Y1, maxsub=0)[0]
print(f"out of {n}: number that PASS the exact check:", res)
