# Guide for reviewers

This file is for a mathematician who wants to check the proofs in [`paper/paper.md`](paper/paper.md). The
computations are re-run by `python reproduce.py` (see the README). What a human has to check is the
**hand-proved reductions** that turn each theorem into the finite computation, and the claim that the
verifiers check what the paper says they check. Below, every such step is listed with the place where it
is proved and the place where it is implemented. Line references are to the trusted files only.

Recommended order:
1. Theorem C (§1, half an hour). It is the easiest entry point.
2. Theorem B (§2).
3. Theorem A (§3).

## 1. Theorem C (14 vertices)

| step | paper | code |
|---|---|---|
| hull facets and their completeness (strict support plus closed oriented surface) | §5.3 items 1–2 | `examples/V14_short_check.py`, checks (1)–(2) |
| volume and centroid from the oriented facets (divergence theorem) | §5.3 item 3 | check (3) |
| vertex equilibria via Lemma 2.1, with one witness per vertex | §2, §5.3 item 4 | check (4) |

The script has about 100 lines, uses only integers and computes no hull itself. The witness table is
printed in §5.3 of the paper.

## 2. Theorem B (at least 9 vertices)

| # | step | paper | what to check | code |
|---|---|---|---|---|
| B1 | shadowing lemma | Lemma 2.1 | the tangent-cone argument; weak inequality $\ge 0$ for "not unstable"; $\lvert r\rvert$ strictly increases | — |
| B2 | boundary triangulation | §4.1 (1) | triangulating faces by diagonals gives a simplicial 2-sphere on the same vertices, whose graph contains all edges of $P$ (hence the tree) | — |
| B3 | cone identity | Lemma 4.1 | tiling from a vertex $q$; flat tetrahedra allowed; some $T^*$ has positive volume | — |
| B4 | certificates | Props. 4.2, 4.3 | the direction of the inequalities; $\mu$ shared by all $T$; $\sum W_q=0$; the verifier stores $2M_T$ (integer) | `verify_joint.py`: `check_certificate`, `check_general_certificate` |
| B5 | exact positive definiteness | §4.3 | Bareiss without pivoting: the $k$-th pivot is the $k$-th leading principal minor, and every division is exact (checked) | `is_positive_definite` |
| B6 | all spheres | §4.2 | (i) the checks of a simplicial 2-sphere; (ii) the exhaustive labelled enumeration and its necessary conditions (min degree ≥ 3, each edge in ≥ 2 triangles); (iii) Wagner's flip theorem for maximal planar graphs; flips that would create a double edge are excluded | `check_sphere`, `enumerate_spheres.py`, `flip_closure.py` |
| B7 | all roots | §4.2 | the listed permutations are verified automorphisms (edges *and* triangles) mapping each vertex to a root | `main` (graph records) |
| B8 | all trees | §4.2 | exhaustive enumeration of acyclic parent maps, cross-checked with Kirchhoff's determinant | `rooted_spanning_trees`, `kirchhoff` |
| B9 | relabelling | proof of Theorem B | (4.1) and (4.2) depend only on the labelled combinatorics, so isomorphisms carry certificates | — |

Record kinds used for $V\le 8$: `cert` (Prop. 4.2) for 448,187 trees and `certG` with two apexes and no
steepest-neighbour terms (Prop. 4.3) for 10,207 trees with $V=8$. The kinds `certH` and the `s1` terms are
implemented and documented in the verifier, but they are used only in the (unfinished) $V=9$ experiments.

## 3. Theorem A (pentahedra)

| # | step | paper | what to check | code |
|---|---|---|---|---|
| A1 | affine invariance, metric $W$ | §3.1, (3.1) | the covectors transform as $A^{-\mathsf T}\ell$; every $W\succ0$ arises | — |
| A2 | stability uses adjacent faces only | Lemma 2.2 | faces meeting face $i$ only at a vertex (pyramids) impose no condition | — |
| A3 | rolling is acyclic, giving a rooted tree | Lemma 3.1 | weak inequality; the face with the largest $\ell^{\mathsf T}W\ell$ is stable | — |
| A4 | certificate | Prop. 3.2 | $\operatorname{tr}(NW)<0$ against $\sum y_i\nu_i^{\mathsf T}W\eta_i\ge0$; the sign of $\eta$ for edges $i>j$ | `core._edge_matrices` (`sg = -1`) |
| A5 | combinatorial types | Lemmas 3.3, 3.4 | Euler argument; concurrency or parallelism of the lateral lines; $0<t_i<1$ with a single A/B labelling; pyramid as degenerate prism, with adjacency a subgraph of $K_5-AB$ | `core.rooted_trees` (375 trees, checked) |
| A6 | **coverage by the two charts** | Lemma 3.5 | every pentahedron lands in $D_Z\cap\{uvw\le1-\rho_0\}$ or in $D_P$, after sorting; the map from chart Z to chart P for $t_1>1-\rho_0$; pyramids at $w=1$ and $q=0$; parallelograms at $\rho=q=0$ | `verify.py` (grids, skip-leaf condition) |
| A7 | **positivity of $\nu$, $\eta$** | Lemma 3.6, App. B | $\Delta_k=\Lambda(d_k-n_k\cdot g)>0$; the only division is by the monomial $u$; scaling is by positive integers; the explicit determinant table proves the signs in the three-tetrahedron decomposition; at pyramid boundaries the middle tetrahedron is flat and the moment formulas extend by continuity | `shapes._tet_data`, `_finish`, `_monomial_content`, `_integer_scale` |
| A8 | box check | Lemma 3.7 | mean value theorem entrywise; Frobenius bounds spectral norm; convexity over the 8 corners; Arb semantics (balls contain exact values, comparisons are "certainly true") | `core.BoxData`, `core.check`, `core._posdef` |
| A9 | the cover | §3.4 | grid exactness; children computed by the verifier; each of the 375 trees certified on every root-to-leaf path; skip leaves only where $\mathrm{lo}_u\mathrm{lo}_v\mathrm{lo}_w>1-\rho_0$ exactly | `verify.verify_cell`, `verify.main` |
| A10 | generated code | §3.3 | the `lambdify` output contains no divisions and no float literals (checked by the audit) | `shapes.build_evaluator` |

Chart Z can be checked entirely by hand from Appendix B of the paper. For chart P, `audit/derive_check.py`
re-derives the proportionality factors from the vertices, using an independent centroid computation.
`penta/test_shapes.py` cross-checks numerically.

## 4. What does not need to be trusted

* The certificate generators and searches: `lb_joint/certify_joint.py`, `general_cert.py`,
  `hybrid_cert.py`, `penta/search.py`, `propose.py`, `compact.py`, and `search_*`.
* Every SDP solver, and every floating-point computation except the dyadic box corners and multipliers
  of Theorem A. Those are exact binary rationals and enter the ball arithmetic exactly.

A wrong certificate can only make a verifier fail. It cannot make one pass.

## 5. Trusted computing base, and known limitations

* The trusted base consists of the Python interpreter, `sympy` (exact polynomial algebra: expansion,
  differentiation, `lambdify` to a division-free polynomial evaluator), Arb through `python-flint`
  (ball arithmetic), and the trusted files named above.
* Both computer-assisted theorems were re-checked by independent implementations written by separate AI
  agents: Theorem A in `audit/` (exact rationals, all 1,339,090 certificates) and Theorem B in
  `audit/theoremB/` (all 458,394 cases, with different algorithms for every part). See Appendix A.3 of
  the paper. Neither agent found a mathematical error. Their suggestions for the exposition are
  incorporated.
* All verification so far has been done by machines and by AI agents. No human expert has yet checked
  the reductions in §§2–4 of the paper. That check is what this guide is for.
