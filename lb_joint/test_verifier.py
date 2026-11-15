"""Negative tests for verify_joint.py: tampered certificate files must be
rejected (exit status 1).  Uses certs_V6.jsonl.gz; temporary files go to the
system temp directory."""
import gzip
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "certs_V6.jsonl.gz")
V = 6


def load():
    with gzip.open(SRC, "rt") as f:
        return [json.loads(l) for l in f]


def run(recs, label, expect_ok=False):
    fd, path = tempfile.mkstemp(suffix=".jsonl.gz")
    os.close(fd)
    with gzip.open(path, "wt") as f:
        for r in recs:
            f.write(json.dumps(r) + "\n")
    p = subprocess.run([sys.executable, os.path.join(HERE, "verify_joint.py"), str(V), path],
                       capture_output=True, text=True)
    os.remove(path)
    last = [l for l in p.stdout.splitlines() if l.strip()][-1]
    ok = (p.returncode == 0) == expect_ok
    print(f"{'PASS' if ok else '*** UNEXPECTED ***'}  {label}: exit {p.returncode}: {last[:110]}")
    return ok


def main():
    base = load()
    certs = [i for i, r in enumerate(base) if r["kind"] == "cert"]
    graphs = [i for i, r in enumerate(base) if r["kind"] == "graph"]
    results = [run(base, "untampered file", expect_ok=True)]

    r = json.loads(json.dumps(base)); del r[certs[17]]
    results.append(run(r, "one certificate deleted"))

    r = json.loads(json.dumps(base)); c = r[certs[5]]; c["mu"] = [0] * V
    results.append(run(r, "mu set to zero"))

    r = json.loads(json.dumps(base)); c = r[certs[9]]
    y = c["Y"][0][1]; k = next(i for i, t in enumerate(y) if t > 0); y[k] = -y[k]
    results.append(run(r, "negative multiplier"))

    r = json.loads(json.dumps(base)); c = r[certs[3]]; c["Y"].pop()
    results.append(run(r, "one tetrahedron missing"))

    r = json.loads(json.dumps(base)); c = r[certs[3]]; c["Y"][1] = c["Y"][0]
    results.append(run(r, "tetrahedron listed twice"))

    r = json.loads(json.dumps(base)); c = r[certs[11]]
    c["Y"] = [[f, [0] * V] for f, y in c["Y"]]       # M_T = sym(mu 1_T^T) has rank <= 2
    results.append(run(r, "all shadow multipliers zero"))

    # scaling a certificate by a positive integer keeps it valid
    r = json.loads(json.dumps(base)); c = r[certs[11]]
    c["Y"] = [[f, [7 * t for t in y]] for f, y in c["Y"]]; c["mu"] = [7 * t for t in c["mu"]]
    results.append(run(r, "certificate scaled by 7 (valid)", expect_ok=True))

    # replace graph 1 (and its certificates) by an exact copy of graph 0 with
    # its certificates: everything else is valid, so only the
    # non-isomorphism check (i) can reject it
    r = json.loads(json.dumps(base))
    g0 = r[graphs[0]]
    keep = [x for x in r if not (x["kind"] in ("graph", "cert") and x.get("gid") == 1)]
    g1 = dict(g0); g1["gid"] = 1
    extra = [dict(x, gid=1) for x in r if x["kind"] == "cert" and x["gid"] == 0]
    r = keep[:1] + [g0, g1] + [x for x in keep[1:] if x["kind"] == "cert"] + extra
    results.append(run(r, "isomorphic duplicate graph (with its certificates)"))

    # bogus automorphism: compose with a transposition that breaks the edge set
    r = json.loads(json.dumps(base)); g = r[graphs[0]]
    E = {frozenset(e) for e in g["edges"]}
    u = next(iter(g["aut"])); p = g["aut"][u]
    for a in range(V):
        for b in range(a + 1, V):
            if int(u) in (a, b):
                continue
            q = list(p); q[a], q[b] = q[b], q[a]
            if {frozenset(q[x] for x in e) for e in E} != E:
                break
        else:
            continue
        break
    g["aut"][u] = q
    results.append(run(r, "corrupted automorphism"))

    # drop a face (no longer a sphere)
    r = json.loads(json.dumps(base)); r[graphs[1]]["faces"].pop()
    results.append(run(r, "triangle removed from a sphere"))

    # wrong parent (not a tree of G / not enumerated) -> missing coverage
    r = json.loads(json.dumps(base)); c = r[certs[0]]
    other = r[certs[1]]; c["parent"] = other["parent"]
    results.append(run(r, "certificate attached to another tree"))

    print("all negative tests behave as expected" if all(results) else "SOME TESTS FAILED")
    return all(results)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
