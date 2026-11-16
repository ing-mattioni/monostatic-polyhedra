#!/usr/bin/env python3
"""
depth_claim.py -- checks the (non-load-bearing) remark of paper Section 4.3:
"For V=8, the trees not covered by Proposition 4.2 are exactly those in which
some triangle avoiding v has all three vertices at depth >= 5 in the tree."
Only the combinatorial half can be tested without an SDP solver: every tree
certified by Prop. 4.3 (certG) should have such a triangle, and no tree that
has one should carry a Prop. 4.2 certificate (cert).
"""
import gzip
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def depth(parent, v):
    d = [0] * len(parent)
    for i in range(len(parent)):
        x, k = i, 0
        while x != v:
            x = parent[x]
            k += 1
        d[i] = k
    return d


faces = {}
cnt = {"cert": [0, 0], "certG": [0, 0]}
for fn in ("certs_V8.jsonl.gz", "exp1_V8_certG.jsonl.gz"):
    with gzip.open(os.path.join(ROOT, "lb_joint", fn), "rt") as f:
        for l in f:
            r = json.loads(l)
            if r["kind"] == "graph":
                faces[r["gid"]] = [tuple(t) for t in r["faces"]]
            if r["kind"] not in cnt:
                continue
            v = r["root"]
            d = depth(r["parent"], v)
            deep = any(v not in t and min(d[a] for a in t) >= 5 for t in faces[r["gid"]])
            cnt[r["kind"]][1 if deep else 0] += 1
for k, (no, yes) in cnt.items():
    print("%-6s trees: with a deep (>=5) triangle avoiding v: %d, without: %d" % (k, yes, no))
ok = cnt["cert"][1] == 0 and cnt["certG"][0] == 0
print("remark consistent with the certificate files:", ok)
sys.exit(0 if ok else 1)
