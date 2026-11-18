# Review guide: fourteen vertices

The publication now has one main theorem and one contraction corollary.
Former A/B claims are not part of this paper or deposit.

| Claim or reduction | Paper | Exact implementation / evidence |
|---|---|---|
| Uniform volume density; distinction from equal point masses | Section 1 | Centroid is integrated over volume, never averaged over vertices |
| 14 exposed vertices, 36 edges, 24 facets | Section 3.1, coordinate table and facet list | V14_short_check.py steps 1–2; independently reconstructed in both other checkers |
| Complete closed outward surface | Section 3.1, facet-adjacency argument | 264 strict support tests, opposite edge incidences, all points present |
| Volume and centroid | Section 3.2 | Signed determinant sums in all three checkers |
| Vertex equilibrium criterion | Lemma 1 | Tangent-cone argument; strict dot products |
| Thirteen non-equilibrium witnesses; vertex 0 unstable | Section 3.3, Table 3 | V14_short_check.py step 4; check_c.py checks all vertex pairs |
| Two stable faces | Section 3.4, Appendix A | Exact projected feet; 24 cyclic sign triples |
| One saddle edge | Section 3.4, Appendix A | Exact segment feet and 12 global support signs for each of 36 true edges |
| No degenerate equilibria | Section 3.4 | Every accepted feature has strict signs; every rejected feature has a strict violation |
| Same face lattice under positive contraction | Section 4 | Invertible linear map diag(1,1,t) |
| U = 1 for every 0 < t <= 1 | Corollary 1, Table 3, Section 4 | Positive affine function of t squared at both endpoints for 13 witnesses; direct root signs also checked |
| Numerical tables match supplied exact data | Verification section | verify/check_manuscript.py |
| Citation scope and limits | Sections 1 and 5 | Bibliography includes public 16-vertex example and relevant research archive |

## Trust boundary

For the shortest proof of existence, inspect examples/V14_short_check.py
(about 100 lines), its hard-coded points and triangles, and the supporting
mathematical argument. This script proves U = 1, not the full S/H census.

For the complete theorem and corollary, inspect verify/check_c.py. It imports
only the standard library and does not import the other two checkers.
It reconstructs the hull from all triples, integrates signed tetrahedra,
and evaluates explicit integer predicates. Its precomputed reports are outputs,
not mathematical assumptions.

verify/exact_verify.py is a separate, more general implementation, using rational
arithmetic and polygonal hull logic. It also handles examples beyond this paper.
Its optional NumPy conversion mode is not used here.

All checking paths ultimately trust Python's integer arithmetic, the interpreter
and the usual software/hardware execution environment. They are not formal
proof-assistant certificates. Multiple implementations and passing tests do not
replace review of the reductions.

reproduce.py, verify/check_manuscript.py and paper/latex/make_ledger.py orchestrate
integrity and presentation consistency; no optimizer, plot library or TeX
engine is trusted for the existence theorem. Hashes check the distributed bytes,
not novelty or correctness by themselves.

## Run and inspect

```sh
python reproduce.py
python examples/V14_short_check.py
python verify/exact_verify.py examples/V14_integer.json --verbose
python verify/check_c.py
```

The exact face/edge numerators and all contraction coefficients are in
examples/V14_family_report.json. Appendix A prints their signs. S = 2 and H = 1
are proved only for t = 1; the interval statement concerns U = 1.

The retained short checker's labels “Theorem C” and “Lemma 2.1” are historical:
they now correspond to Theorem 1 and Lemma 1. Its trusted source is unchanged.
