"""(Run from the repository root; needs the local, git-ignored lb_joint/uncertified_V9.json.)
Generate hybrid certificates for a few uncertified V=9 trees and test verify_joint.check_hybrid_certificate,
including tampered variants that must be rejected."""
import os, sys, json, random, copy
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, 'lowerbound'); sys.path.insert(0, 'lb_joint')
import enumerate_cases as ec, hybrid_cert as hc
import verify_joint as vj
V = 9
unc = json.load(open('lb_joint/uncertified_V9.json'))
rng = random.Random(21)
sample = rng.sample(unc, 5)
gds = ec.graph_data(V, with_trees=False)
ok = bad_caught = 0; tamper_total = 0
for c in sample:
    gd = gds[c["graph"]]; v = c["root"]; parent = c["parent"]
    rec, z = hc.certify_tree(V, v, parent, gd.faces)
    assert rec is not None, "hybrid failed"
    E, F, nbrs = vj.check_sphere(V, [list(e) for e in gd.edges], [list(f) for f in gd.faces])
    vj.check_hybrid_certificate(V, F, nbrs, v, tuple(parent), rec); ok += 1
    tampers = []
    t = copy.deepcopy(rec); t["pairs"] = t["pairs"][:-1]; tampers.append(("drop a pair", t))
    t = copy.deepcopy(rec); t["single"] = t["single"][:-1]; tampers.append(("drop a good triangle", t))
    t = copy.deepcopy(rec); t["mu"] = [0] * V; t["mu2"] = [0] * V; tampers.append(("mu = 0", t))
    t = copy.deepcopy(rec); t["q"] = v; tampers.append(("q = v", t))
    t = copy.deepcopy(rec); t["single"][0][1] = [-1 if i != v else 0 for i in range(V)]; tampers.append(("negative y", t))
    t = copy.deepcopy(rec); t["pairs"].append(copy.deepcopy(t["pairs"][0])); tampers.append(("duplicate pair", t))
    t = copy.deepcopy(rec); f0 = t["single"][0][0]; t["pairs"].append([f0, t["pairs"][0][1], t["pairs"][0][2]]); tampers.append(("good also bad", t))
    t = copy.deepcopy(rec); t["single"][0][1] = [0] * V; tampers.append(("y = 0 (matrix not PD)", t))
    for name, tt in tampers:
        tamper_total += 1
        try:
            vj.check_hybrid_certificate(V, F, nbrs, v, tuple(parent), tt)
            print("  NOT CAUGHT:", name)
        except vj.Fail as e:
            bad_caught += 1
print(f"valid hybrid certificates accepted: {ok}/{len(sample)}; tampered variants rejected: {bad_caught}/{tamper_total}")
