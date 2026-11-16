# Independent audit of the pentahedron proof (Theorem 3)

This is the report of an independent referee pass carried out by a separate AI agent (Claude Opus) that
was asked to try to break the proof in `penta/`: mathematics, code, numerics, verifier. Its scripts and
logs are in this directory; it shares no code with `penta/core.py` / `penta/shapes.py`. The report is
reproduced as delivered; a section at the end lists what was changed in response.

**Summary.** No mathematical error and no counterexample were found. All 575,441 chart-Z certificates
verify in the independent exact checker, in both the original and the compacted file; samples of chart P
also all pass. The one real defect was in `verify.py`: under `python -O` or `PYTHONOPTIMIZE=1` it accepted
an essentially empty certificate set (checks written as `assert`).

## (a) Mathematics (`penta/README.md` §1–§3)

No error and no gap affecting correctness was found after checking every listed point.

1. **Foot and stability criterion.** The foot `x_i = W l_i/(l_iᵀWl_i)` is right. Face *i* stable ⇔
   `l_iᵀW(l_j−l_i) < 0` for all *j ≠ i*. A foot on the face boundary counts as "not stable", which matches
   the `≥ 0` used in the certificate.
2. **"Not stable ⇒ an adjacent face j".** Correct, including pyramids. Inside the plane of face *i* the
   irredundant constraints are its edges; faces meeting *i* only at a vertex (opposite lateral faces of a
   pyramid, A and B at the apex) give redundant constraints.
3. **Acyclicity and the tree.** `l_jᵀWl_j > l_iᵀWl_i + 2l_iᵀW(l_j−l_i) ≥ l_iᵀWl_i`. So the face with the
   largest value is always stable, and with exactly one stable face *s* the choice of π gives a spanning
   tree rooted at *s*.
4. **Certificate direction.** `<N,W> < 0` for W ⪰ 0, W ≠ 0, against `Σ y_i ν_iᵀWη_i ≥ 0`. Correct.
5. **Signs.** Edges i→j with i > j use −η_(j,i), which is right. The only monomial divided out is `u`, for
   η(Q1,A) in chart Z; u > 0 at genuine shapes. Chart P divides out no monomial.
6. **Classification.** Euler's formula leaves only the triangular prism and the quadrilateral pyramid. The
   prism's lateral lines are concurrent or parallel. The README did not prove the concurrency step and
   used two facts silently: O lies outside each lateral segment, and the same triangle is nearer to O on
   all three lines. These justify 0 < t_i ≤ 1 with a single A/B labelling. *Add a short lemma.*
7. **Coverage by the two charts.** No hole found. Sorting and swapping A/B are automorphisms of K5−{AB},
   so the set of 375 rooted trees is invariant. Toric coordinates u = t1/t2, v = t2/t3, w = t3, all in
   (0,1]. w = 1 gives exactly the pyramids whose base has a pair of non-parallel opposite sides; their
   adjacency graph is a subgraph. Chart P with ρ = 0 gives the parallel-edge prisms; ρ = 0, q = 0 gives the
   pyramids over trapezoids and parallelograms. A chart-Z shape with uvw > 1−ρ₀ maps to ρ = 1−t1,
   p = (1−t2)/(1−t1), q = (1−t3)/(1−t2), all in range. Skip leaves are sound. The Z→P affine equivalence
   was checked exactly at 36 random overlap points.
8. **Positivity of ν and η** (only asserted in the README). `derive_check.py` derives `l_k` from the
   vertices with its own divergence-theorem centroid and proves exactly, in sympy, that each ν_k ∥ l_k and
   η_ij ∥ (l_j − l_i) in both charts; the ratios factor into a positive constant times centroid-to-plane
   factors over a volume factor, plus 1/u. Each factor is nonzero at genuine shapes and the genuine region
   is connected, so the sign is constant; it is positive (also checked exactly at 760 random points).
   *The README should include this argument.*
9. **Numerics in `core.check`** (`arb_semantics.py`, python-flint 0.9.0): `arb(float)` is exact;
   comparisons are "certainly true" (overlapping balls and NaN compare False); `rad()` returns an upper
   bound for the radius; `lower()`/`upper()` round outward. The mean-value / Frobenius / 8-corner argument
   is sound and the Sylvester formula is correct. The generated lambdify sources contain no `/`, no float
   literals and no function calls; the only powers are 2, 3, 4. Harmless oddity: `arb(inf) > 0` is True,
   but `verify.py` rejects non-finite floats.

## (b) Independent exact re-check (`recheck.py`)

Method: rebuilds every node's box exactly from the bisection trees; Taylor-expands the polynomials exactly
on [−1,1]³; bounds the degree ≥ 2 remainder by absolute coefficients and the spectral norm by the maximal
row sum; checks the 8 corners with exact Sylvester and y ≥ 0 exactly; would subdivide on failure (never
needed). The polynomials enter only as data dumped from `shapes.py` and are verified by (a)8.

| run | certificates | result |
|---|---|---|
| Chart Z, original file | all 575,441 (30,232 nodes) | 0 failures |
| Chart Z, compacted file (the final proof) | all 575,441 | 0 failures |
| Chart Z, stratified sample | 20,540 (20,000 uniform, 190 at depth ≥ 12, 350 on the boundary) | 0 failures |
| Chart P, original file, sample | 6,556 (199 deep, 357 boundary) | 0 failures |
| Chart P, compacted file, sample | 5,555 (196 deep, 359 boundary) | 0 failures |

Coverage: an independent traversal confirms the full grid, exact skip leaves, and all 375 trees on every
root-to-leaf path; counts match `penta/verify.log` exactly. Pointwise end-to-end: 460 exact random shapes
on the original files and 300 on the compacted files (pyramids, prisms, Z→P overlap shapes): at every
shape all 375 trees are certified and every certificate on the path gives N(x) ≺ 0 exactly. Negative
controls: the checker rejects negated y0 200/200 and most certificates with Y1 × 50 or a box enlarged 4×;
in 1,800 perturbed or wrong-tree comparisons `core.check` never accepted anything the independent checker
rejected. *Not done by the auditor: a full exact run over all chart-P certificates (sampled only).*

## (c) Falsification (`falsify.py`, `falsify2.py`)

Random sampling: 1,123,143 valid convex pentahedra (frusta, prisms, pyramids; singular values in e^±4, edge
ratios down to 1e-5); none has fewer than 2 stable faces; smallest second-best margin ≈ 3e-5.
Unconstrained optimisation (400 evolution-strategy chains × 3,000 steps per family, then Nelder–Mead): no
run reached ≤ 0; the value 0 is approached only in the flat degenerate limit (vol/diam³ ≈ 1e-8). With
vol/diam³ ≥ 1e-3 the minimiser always sits on that constraint, with best second-best margins 3.0e-4…4.9e-4.

## (d) `verify.py` audit

The logic is sound: grid check, recomputed exact midpoints, per-path pending sets, Z-only skip leaves with
an exact Fraction test, integer tree ids in [0,375), finite-float multipliers, rejection of FAILED nodes.
Tampered files run through the real `verify.py`:

| case | caught? |
|---|---|
| Missing cells (one legal Z skip cell, ρ₀ = 15/16, P empty) | yes |
| Same file with `python -O` / `PYTHONOPTIMIZE=1` | **no — "ALL CHECKS PASSED"** |
| Illegal skip leaf | yes |
| NaN multipliers | yes |
| Only one child | yes |
| Non-dict node | crash with exit 1 (fail-safe) |
| ρ₀ mismatch between files | yes, but only through a bare `assert` |

**Must fix:** the grid, ρ₀, chart-field and grid-exactness checks in `verify.py`, `len(TREES) == 375` in
`core.py` and the identity checks in `shapes.py` are bare `assert`s, which Python strips under -O.
This does not invalidate the current certificates. `compact.py` cannot affect soundness.

## (e) Remarks on the write-up

1. "First improvement … since 1969" should be hedged ("to our knowledge").
2. §5 of the paper is an outline; add the concurrency/labelling lemma, the pyramid adjacency remark and the
   ν/η positivity argument to `penta/README.md`.
3. State that the multipliers were rounded and that `verify.py` passed on the final files.
4. The Monte-Carlo sanity run had no stored log.
5. "One SDP per tree decides the question" overstates it; say "reduces to".
6. Paper §3 called OEIS A000109 the one external fact while the README said the list is re-derived.

---

## Follow-up by the main session (2026-11-16)

* **`assert` → explicit checks.** `penta/verify.py`, `penta/core.py`, `penta/shapes.py` (and
  `lb_joint/enumerate_spheres.py`) no longer use `assert` for any check. The full verification was re-run
  with `python -O`: all checks pass (`penta/verify.log`). The tamper suite was re-run against the fixed
  verifier: every tampered file is rejected, also under `-O` (`tamper_after_fix.log`).
* **Chart P, all certificates.** The auditor's independent checker was run on *all* 763,649 chart-P
  certificates of the final file: see `recheck_compactP_all.log` for the outcome.
* **Write-up.** The lemma on the prism type, the pyramid description and the positivity argument were
  added to `penta/README.md`; the statements criticised in (e) were corrected; `penta/sanity.log` stores
  the Monte-Carlo run; the completeness of the sphere lists in Theorem 1 is verified by
  `lb_joint/enumerate_spheres.py` (no external count).
