# Monostatic polyhedra: new bounds for the Conway–Guy problems

A homogeneous convex polyhedron is **mono-stable** if it can rest on exactly one face and
**mono-unstable** if it can balance on exactly one vertex. How small can such polyhedra be?
(Conway & Guy, 1966/69; all cases open.)

| quantity | before | now (this repository) |
|---|---|---|
| min. number of vertices of a mono-unstable polyhedron | 7 ≤ V ≤ 18 | **9 ≤ V ≤ 14** |
| min. number of faces of a mono-stable polyhedron | 5 ≤ F ≤ 14 | **6 ≤ F ≤ 14** |

* **Theorem 1.** Every homogeneous convex polyhedron with at most 8 vertices has at least two unstable
  vertices. (Exact integer certificates: 24 795 cases for 7 vertices, 432 041 for 8.)
* **Theorem 2.** There is a homogeneous mono-unstable polyhedron with 14 vertices; one with integer
  coordinates of absolute value ≤ 1000 is in `examples/V14_integer.json`. (Exact rational verification.)
* **Theorem 3.** Every homogeneous convex pentahedron has at least two stable faces — the first
  improvement of the lower bound in Conway and Guy's original problem since 1969. (Rigorous
  ball-arithmetic cover of the 3-dimensional space of affine shapes, 375 rolling patterns.)

![Top view of the 14-vertex mono-unstable polyhedron: a thin plate whose outline is a double spiral around the centre of mass](examples/V14_plate.svg)

The mathematics (statements, proofs of all reductions, discussion, references) is in
[`paper/paper.md`](paper/paper.md); the proof of Theorem 3 is detailed in
[`penta/README.md`](penta/README.md).

**Status / caveats.** These are computer-assisted results produced in an AI-driven research session
(Claude, with the repository owner); they have **not been peer reviewed**. Each theorem is checked by a
small independent verifier (below). Theorem 1 uses one external fact (the number of triangulations of
the sphere with ≤ 8 vertices, OEIS A000109); Theorem 3 relies on `python-flint`/Arb ball arithmetic and
on `sympy` for polynomial algebra. Not found: a 13-vertex mono-unstable or a 13-face mono-stable
polyhedron (limited searches; no conclusion). The problems themselves remain open: the exact minima are
still unknown.

## Verify

```bash
python -m venv .venv && . .venv/bin/activate       # Windows: .venv\Scripts\activate
pip install sympy python-flint                     # only needed for Theorem 3

# Theorem 1 (standard library only, exact integers)
python lb_joint/verify_joint.py 7 lb_joint/certs_V7.jsonl.gz                                   # ~5 s
python lb_joint/verify_joint.py 8 lb_joint/certs_V8.jsonl.gz lb_joint/exp1_V8_certG.jsonl.gz   # ~3 min

# Theorem 2 (standard library only, exact rationals)
python verify/exact_verify.py examples/V14_integer.json

# Theorem 3
python penta/verify.py
```

## Layout

| path | content |
|---|---|
| `paper/paper.md` | write-up |
| `examples/` | the 14-vertex polyhedron of Theorem 2 with its exact certificate report |
| `verify/` | `exact_verify.py`: exact-rational classifier of all equilibria of a polyhedron (+ tests, certified 16/18-vertex examples) |
| `lb_joint/` | Theorem 1: certificate generator, standalone verifier, certificates for V = 4…8 |
| `penta/` | Theorem 3: chart derivation, rigorous kernel, search, verifier, Monte-Carlo sanity check |
| `search_vertices/` | search for small mono-unstable polyhedra; certified 14- and 15-vertex bodies |
| `search_faces/` | search for small mono-stable polyhedra; certified 14/15-face bodies with improved margins; Reshetov's bodies re-verified |
| `lowerbound/` | replication of the published SDP method (Bozóki–Domokos–Papp–Regős) and analysis of where it stops |
| `explore/` | exploratory scripts (flat pyramids, ideal-frustum bound, prototypes) |
