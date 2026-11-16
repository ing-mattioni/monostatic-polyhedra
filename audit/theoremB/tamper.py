#!/usr/bin/env python3
"""
tamper.py -- negative controls for check_theoremB.py.

Part 1 (end-to-end, V = 6): tampered copies of lb_joint/certs_V6.jsonl.gz
  are written to audit/theoremB/tamper/ and the checker is run on each as a
  separate process; each run must exit with status != 0 and print the
  expected reason.  One untampered copy must pass.
Part 2 (record level, V = 8, certG of Prop. 4.3): single records of
  lb_joint/exp1_V8_certG.jsonl.gz are tampered and passed to the checker's
  record test check_record (same code path as the full run).
Part 3 (sharpness): for one Prop. 4.2 record (V = 6) and one Prop. 4.3
  record (V = 8) the smallest integer perturbation delta that destroys
  positive definiteness is located exactly; delta - 1 must be accepted and
  delta rejected.
Part 4 (end-to-end, V = 8): a tampered copy of the certG file with several
  tampered records and one dropped record; the full checker must report
  exactly these.
"""
import copy
import gzip
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import check_theoremB as C  # noqa: E402

PY = sys.executable
TD = os.path.join(HERE, "tamper")
os.makedirs(TD, exist_ok=True)


def read(fn):
    with gzip.open(fn, "rt") as f:
        return [json.loads(l) for l in f if l.strip()]


def write(fn, recs):
    with gzip.open(fn, "wt") as f:
        for r in recs:
            f.write(json.dumps(r) + "\n")


def run(V, files):
    p = subprocess.run([PY, "-I", os.path.join(HERE, "check_theoremB.py"), str(V)] + files,
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


results = []


def report(name, ok, detail):
    results.append((name, ok))
    print("%-58s %s   %s" % (name, "REJECTED as expected" if ok else "*** NOT DETECTED ***", detail),
          flush=True)


# ---------------------------------------------------------------- part 1
print("Part 1: end-to-end, V=6 (certs_V6.jsonl.gz)", flush=True)
base = read(os.path.join(ROOT, "lb_joint", "certs_V6.jsonl.gz"))
certs = [i for i, r in enumerate(base) if r["kind"] == "cert"]
gidx = [i for i, r in enumerate(base) if r["kind"] == "graph"]


def first_cert(pred=lambda r: True):
    return next(i for i in certs if pred(base[i]))


def t_untampered(recs):
    return recs


def t_drop_tree(recs):
    del recs[certs[-1]]
    return recs


def t_neg_y(recs):
    r = recs[certs[0]]
    y = r["Y"][0][1]
    k = next(j for j, t in enumerate(y) if t > 0)
    y[k] = -y[k]
    return recs


def t_neg_mu(recs):
    r = recs[certs[0]]
    r["mu"] = [-t for t in r["mu"]]
    return recs


def t_missing_block(recs):
    r = recs[first_cert(lambda r: len(r["Y"]) > 1)]
    r["Y"].pop()
    return recs


def t_dup_block(recs):
    r = recs[first_cert(lambda r: len(r["Y"]) > 1)]
    r["Y"][-1] = copy.deepcopy(r["Y"][0])
    return recs


def t_block_with_apex(recs):
    r = recs[certs[0]]
    g = next(x for x in recs if x["kind"] == "graph" and x["gid"] == r["gid"])
    F = next(f for f in g["faces"] if r["root"] in f)
    r["Y"][0][0] = F
    return recs


def t_yroot(recs):
    r = recs[certs[0]]
    r["Y"][0][1][r["root"]] = 1
    return recs


def t_float(recs):
    r = recs[certs[0]]
    r["mu"][0] = r["mu"][0] + 0.5
    return recs


def t_bad_parent(recs):
    # point a vertex to a non-neighbour that is the root (stays acyclic)
    for i in certs:
        r = recs[i]
        g = next(x for x in recs if x["kind"] == "graph" and x["gid"] == r["gid"])
        E = {tuple(sorted(e)) for e in g["edges"]}
        v = r["root"]
        for a in range(6):
            if a != v and tuple(sorted((a, v))) not in E and r["parent"][a] != -1:
                r["parent"][a] = v
                return recs
    raise RuntimeError("no non-neighbour found")


def t_cycle_parent(recs):
    r = recs[certs[0]]
    v = r["root"]
    g = next(x for x in recs if x["kind"] == "graph" and x["gid"] == r["gid"])
    E = {tuple(sorted(e)) for e in g["edges"]}
    for (a, b) in E:
        if v not in (a, b):
            r["parent"][a] = b
            r["parent"][b] = a
            return recs


def t_flip_sphere(recs):
    # an edge flip in sphere gid 0: still a sphere, but now isomorphic to gid 1
    # or not in the list -> listed != generated / trees do not match
    g = recs[gidx[0]]
    faces = [tuple(sorted(f)) for f in g["faces"]]
    E = {tuple(sorted(e)) for e in g["edges"]}
    for (a, c) in sorted(E):
        two = [f for f in faces if a in f and c in f]
        b = [x for x in two[0] if x not in (a, c)][0]
        d = [x for x in two[1] if x not in (a, c)][0]
        if tuple(sorted((b, d))) not in E:
            nf = [f for f in faces if f not in two] + [tuple(sorted((a, b, d))), tuple(sorted((c, b, d)))]
            g["faces"] = [list(f) for f in nf]
            g["edges"] = [list(e) for e in sorted((E - {(a, c)}) | {tuple(sorted((b, d)))})]
            return recs


def t_not_sphere(recs):
    g = recs[gidx[0]]
    g["faces"] = g["faces"][:-1]
    return recs


def t_drop_root(recs):
    g = recs[gidx[0]]
    u = g["roots"].pop()
    g["aut"].pop(str(u), None)
    return recs


def t_bad_aut(recs):
    g = recs[gidx[0]]
    u = next(k for k in g["aut"])
    g["aut"][u] = list(range(6))[::-1]
    return recs


tests = [
    ("untampered copy (must PASS)", t_untampered, None),
    ("dropped tree record", t_drop_tree, "required trees without valid certificate"),
    ("negated multiplier y_i", t_neg_y, "negative multiplier y"),
    ("negated direction vector mu", t_neg_mu, "not positive definite"),
    ("missing triangle block", t_missing_block, "do not cover all triangles"),
    ("duplicated block replacing another", t_dup_block, "has two blocks"),
    ("block triangle containing the apex", t_block_with_apex, "contains the apex"),
    ("y_root != 0", t_yroot, "y_root != 0"),
    ("non-integer entry in mu", t_float, "not an integer vector"),
    ("parent not adjacent", t_bad_parent, "not an edge of G"),
    ("parent map with a 2-cycle", t_cycle_parent, "cycle"),
    ("sphere changed by an edge flip", t_flip_sphere, "FAIL"),
    ("sphere with a triangle removed", t_not_sphere, "not a simplicial 2-sphere"),
    ("root removed from root list", t_drop_root, "orbits without listed root"),
    ("wrong stated automorphism", t_bad_aut, "stated automorphism"),
]
for name, fn, expect in tests:
    recs = fn(copy.deepcopy(base))
    path = os.path.join(TD, "V6_%s.jsonl.gz" % fn.__name__)
    write(path, recs)
    code, out = run(6, [path])
    if expect is None:
        ok = (code == 0 and "RESULT V=6: PASS" in out)
        results.append((name, ok))
        print("%-58s %s" % (name, "PASS as expected" if ok else "*** UNEXPECTED FAIL ***"), flush=True)
        continue
    ok = code != 0 and expect in out
    line = next((l for l in out.splitlines() if l.startswith("FAIL")), "")
    report(name, ok, line[:110])

# ---------------------------------------------------------------- part 2
print("\nPart 2: record level, V=8, Prop. 4.3 records (exp1_V8_certG.jsonl.gz)", flush=True)
graphs = {}
with gzip.open(os.path.join(ROOT, "lb_joint", "certs_V8.jsonl.gz"), "rt") as f:
    for l in f:
        r = json.loads(l)
        if r["kind"] == "graph":
            graphs[r["gid"]] = C.make_graph(8, r)
        elif r["kind"] == "cert":
            break
C._init(8, graphs)
with gzip.open(os.path.join(ROOT, "lb_joint", "exp1_V8_certG.jsonl.gz"), "rt") as f:
    G8 = []
    for l in f:
        r = json.loads(l)
        if r["kind"] == "certG":
            G8.append(r)
        if len(G8) >= 50:
            break
r0 = G8[0]
k, key, err, n = C.check_record(r0)
print("untampered certG record: error = %r, PD tests = %d" % (err, n), flush=True)
results.append(("untampered certG record accepted", err is None))


def tg(name, mut, expect):
    r = copy.deepcopy(r0)
    mut(r)
    _, _, err, _ = C.check_record(r)
    report(name, err is not None and expect in err, "error: %s" % err)


def m_wsum(r):
    r["groups"][0]["W"][0][1] += 1
    r["groups"][0]["W"][1][0] += 1


def m_wasym(r):
    r["groups"][0]["W"][0][1] += 1
    r["groups"][1]["W"][0][1] -= 1


def m_drop_block(r):
    r["groups"][1]["blocks"].pop()


def m_neg_y(r):
    y = r["groups"][1]["blocks"][0][1]
    j = next(i for i, t in enumerate(y) if t > 0)
    y[j] = -y[j]


def m_s1(r):
    r["groups"][0]["blocks"][0][2] = [[1, 2, 5]]


def m_s1flag(r):
    r["s1"] = True


def m_dupq(r):
    r["groups"][1]["q"] = r["groups"][0]["q"]


def m_apex_tri(r):
    q = r["groups"][1]["q"]
    gr = graphs[r["gid"]]
    F = next(t for t in gr["tris"] if q in t)
    r["groups"][1]["blocks"][0][0] = list(F)


def m_drop_group(r):
    r["groups"].pop()


def m_yroot(r):
    r["groups"][1]["blocks"][0][1][r["root"]] = 3


def m_mu(r):
    r["groups"][1]["mu"] = [-t for t in r["groups"][1]["mu"]]


tg("W_q not summing to zero (one symmetric entry +1)", m_wsum, "sum of W_q is not 0")
tg("W_q not symmetric (sum still zero)", m_wasym, "W not symmetric")
tg("missing triangle block (second apex)", m_drop_block, "do not cover all triangles")
tg("negated multiplier y_i (second apex)", m_neg_y, "negative multiplier y")
tg("y_root != 0 in a block of the second apex", m_yroot, "y_root != 0")
tg("negated mu^(q)", m_mu, "not positive definite")
tg("S1 term injected", m_s1, "S1 terms present")
tg("s1 flag set", m_s1flag, "s1 flag set")
tg("two groups with the same apex", m_dupq, "apexes not distinct")
tg("block triangle containing its apex", m_apex_tri, "contains the apex")
tg("second group dropped (W_v alone != 0)", m_drop_group, "sum of W_q is not 0")

# ---------------------------------------------------------------- part 3
print("\nPart 3: sharpness of the exact positive-definiteness test", flush=True)


def smallest_breaking(ok_at):
    """ok_at(d) monotone: True for small d, False for large d.  Returns least
    d >= 1 with ok_at(d) False (exact integer search)."""
    hi = 1
    while ok_at(hi):
        hi *= 2
    lo = hi // 2   # ok_at(lo) True (or lo = 0)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ok_at(mid):
            lo = mid
        else:
            hi = mid
    return hi


# Prop. 4.2 record, V = 6: decrease mu_a (a in no block-triangle? no: any a)
C._init(6, {r["gid"]: C.make_graph(6, r) for r in base if r["kind"] == "graph"})
rc = base[certs[0]]
a = rc["root"]


def ok42(d):
    r = copy.deepcopy(rc)
    r["mu"][a] -= d
    return C.check_record(r)[2] is None


d42 = smallest_breaking(ok42)
r = copy.deepcopy(rc)
r["mu"][a] -= d42 - 1
e1 = C.check_record(r)[2]
r = copy.deepcopy(rc)
r["mu"][a] -= d42
e2 = C.check_record(r)[2]
print("Prop 4.2 record (V=6, gid %d, root %d): mu_root - %d accepted (%r); mu_root - %d rejected (%r)"
      % (rc["gid"], rc["root"], d42 - 1, e1, d42, e2), flush=True)
results.append(("sharp PD threshold, Prop. 4.2", e1 is None and e2 is not None))

# Prop. 4.3 record, V = 8: move delta from W_q[b][b] of group 0 to group 1
C._init(8, graphs)
b = 3


def ok43(d):
    r = copy.deepcopy(r0)
    r["groups"][0]["W"][b][b] -= d
    r["groups"][1]["W"][b][b] += d
    return C.check_record(r)[2] is None


d43 = smallest_breaking(ok43)
r = copy.deepcopy(r0)
r["groups"][0]["W"][b][b] -= d43 - 1
r["groups"][1]["W"][b][b] += d43 - 1
e1 = C.check_record(r)[2]
r = copy.deepcopy(r0)
r["groups"][0]["W"][b][b] -= d43
r["groups"][1]["W"][b][b] += d43
e2 = C.check_record(r)[2]
print("Prop 4.3 record (V=8, gid %d, root %d): W_q[%d][%d] shifted by %d accepted (%r); by %d rejected (%r)"
      % (r0["gid"], r0["root"], b, b, d43 - 1, e1, d43, e2), flush=True)
results.append(("sharp PD threshold, Prop. 4.3 (sum of W still 0)", e1 is None and e2 is not None))

# ---------------------------------------------------------------- part 4
if "--skip-v8" not in sys.argv:
    print("\nPart 4: end-to-end V=8 with a tampered certG file (full run, several minutes)", flush=True)
    allG = read(os.path.join(ROOT, "lb_joint", "exp1_V8_certG.jsonl.gz"))
    idx = [i for i, r in enumerate(allG) if r["kind"] == "certG"]
    muts = [m_wsum, m_drop_block, m_neg_y, m_mu, m_wasym]
    tampered_keys = []
    for m, i in zip(muts, idx[100:105]):
        m(allG[i])
        tampered_keys.append((allG[i]["gid"], allG[i]["root"], tuple(allG[i]["parent"])))
    dropped = allG[idx[200]]
    dropped_key = (dropped["gid"], dropped["root"], tuple(dropped["parent"]))
    del allG[idx[200]]
    path = os.path.join(TD, "V8_exp1_tampered.jsonl.gz")
    write(path, allG)
    code, out = run(8, [os.path.join(ROOT, "lb_joint", "certs_V8.jsonl.gz"), path])
    flagged = [l for l in out.splitlines() if l.startswith("FAIL")]
    for l in flagged:
        print("   " + l[:150])
    summary = [l for l in out.splitlines() if l.startswith("(D)") or l.startswith("RESULT")]
    for l in summary:
        print("   " + l)
    ok = (code != 0 and "5 invalid certificate records" in out
          and "6 required trees without valid certificate" in out)
    report("V=8: 5 tampered + 1 dropped certG records", ok, "exit status %d" % code)

print("\nSUMMARY: %d controls, %d behaved as expected" % (len(results), sum(ok for _, ok in results)))
sys.exit(0 if all(ok for _, ok in results) else 1)
