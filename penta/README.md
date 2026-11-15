# No homogeneous convex pentahedron is mono-stable

**Theorem C.** Let *P* be a convex polyhedron with five faces and uniform density, and let *c* be its
centre of mass. Then at least two faces of *P* contain the orthogonal projection of *c* onto their plane
in their relative interior (i.e. *P* can rest on at least two faces). Consequently a homogeneous
mono-stable polyhedron has at least **6** faces. (The previous lower bound, 5, is Conway and Guy's
theorem that no homogeneous tetrahedron is mono-stable, 1969.)

The proof is computer-assisted: a finite cover of the 3-dimensional space of affine shapes of
pentahedra, on each piece of which every possible "rolling pattern" is refuted by a 3×3 matrix
certificate checked in rigorous ball arithmetic.

## 1. Reduction: affine shape × metric

Put the origin at the centroid and write `P = {x : l_i·x ≤ 1, i = 1..5}` (covectors `l_i`).
The centroid is affinely covariant, so instead of moving *P* we fix an affine normal form of *P* and let
the Euclidean structure be an arbitrary inner product `<x,y> = xᵀ M y`; put `W = M⁻¹` (positive definite).

* The foot of the centroid on the plane of face *i* is `x_i = W l_i / (l_iᵀ W l_i)`.
* Face *i* is *stable* iff `l_j·x_i < 1` for every other face *j*. If it is not stable, there is an
  **adjacent** face *j* (sharing an edge with *i*) with `l_iᵀ W (l_j − l_i) ≥ 0`
  (the constraints of non-adjacent faces are redundant in the plane of face *i*).
* If `l_iᵀ W (l_j − l_i) ≥ 0` and `l_i ≠ l_j` then `l_jᵀ W l_j > l_iᵀ W l_i` (expand
  `(l_j − l_i)ᵀ W (l_j − l_i) > 0`). So the relation "*i* rolls to *j*" has no cycles, and the face with the
  largest `l_iᵀ W l_i` (the face nearest to *c*) is stable.

Hence, if *P* had at most one stable face *s*, choosing for every other face *i* an adjacent face `π(i)`
with `l_iᵀ W (l_π(i) − l_i) ≥ 0` gives a spanning tree of the face-adjacency graph rooted at *s*.

**Certificate.** Fix a shape and a rooted spanning tree `(s, π)`. Let `ν_i` be a positive multiple of `l_i`
and `η_{i→j}` a positive multiple of `l_j − l_i`. If there are numbers `y_i ≥ 0` (*i ≠ s*) such that

    N = Σ_{i≠s} y_i · sym(ν_i η_{i→π(i)}ᵀ)        (sym(abᵀ) = (abᵀ + baᵀ)/2)

is **negative definite**, the tree is impossible for that shape: `<N, W> < 0` for every `W ⪰ 0, W ≠ 0`,
whereas `<N, W> = Σ y_i ν_iᵀ W η_{i→π(i)} ≥ 0`.

## 2. Shapes

A convex polyhedron with five faces is combinatorially a triangular prism (two triangles *A*, *B*, three
quadrilaterals *Q1, Q2, Q3*; 6 vertices) or a quadrilateral pyramid (5 vertices) — Euler's formula.
In the first case the three lateral edges lie on three pairwise coplanar lines, which are concurrent
(apex *O*) or parallel; the pyramid is the limit in which one lateral edge has length 0.
Face adjacency in the prism type: all pairs except {*A*, *B*}; this graph has 75 spanning trees, i.e.
**375 rooted trees**; the adjacency graph of a pyramid is a subgraph, so its rooted trees are among them.
All 375 are refuted at every shape. Lateral edges are labelled so that the parameters below are sorted
(the set of 375 trees is invariant under relabelling). `Q_i` is the lateral face not containing edge *i*.

* **Chart Z (cone frame).** `O = 0`, `B_i = e_i`, `A_i = t_i e_i`, `0 < t_1 ≤ t_2 ≤ t_3 ≤ 1`
  (`t_3 = 1`: pyramid), with toric coordinates `t_3 = w, t_2 = wv, t_1 = wvu`, `(u,v,w) ∈ [0,1]³`.
  Every prism-type pentahedron with concurrent lateral edges, and every pyramid whose base has a pair of
  non-parallel opposite sides, is affinely equivalent to exactly such a shape.
* **Chart P (prism frame).** `B_1 = (0,0,0), B_2 = (1,0,0), B_3 = (0,1,0)`, `G = (1/3,1/3)`,
  `σ = ρ(1, p, pq)`, `s = (1, p, pq)`, `A_i = ((1−σ_i)(B_i)_{xy} + σ_i G, s_i)`,
  `(ρ,p,q) ∈ [0,ρ₀]×[0,1]²`. For `ρ > 0` this is the affine image of the cone frame with `t_i = 1 − σ_i`
  (apex `(G, 1/ρ)`); `ρ = 0` gives the prisms with parallel lateral edges (and, with `q = 0`, the pyramids
  over trapezoids and parallelograms).

Chart Z degenerates at the corner `t → (1,1,1)` (where chart P is used); the toric coordinates make all
direction data continuous on the rest of the closed cube. A shape of chart Z with `t_1 = uvw > 1 − ρ₀`
is a shape of chart P with `ρ = 1 − t_1 < ρ₀`. So it suffices to refute all trees on
`{(u,v,w) ∈ [0,1]³ : uvw ≤ 1 − ρ₀}` in chart Z and on the whole box in chart P.

`shapes.py` derives, exactly (sympy), for each chart polynomial vectors `ν_f` and `η_{(i,j)}`: if plane *k*
is `n_k·x ≤ d_k` in the chart frame and the centroid is `g = mom/vol`, then with
`Δ_k = d_k·vol − n_k·mom > 0` (the centroid is interior) `l_k = vol·n_k/Δ_k` and
`l_j − l_i = vol·(Δ_i n_j − Δ_j n_i)/(Δ_iΔ_j)`; so `ν_k = n_k` and
`η_{(i,j)} = (Δ_i n_j − Δ_j n_i)/monomial` (the monomial content in the chart variables is divided out so
that the vectors do not vanish identically on the boundary of the parameter box; identity checks are
executed in exact arithmetic), scaled to integer coefficients. `test_shapes.py` cross-checks them
numerically against covectors computed directly from the vertices.

## 3. Rigorous cover

For a box *X* (centre *c*, half-widths *h*) of a chart and a rooted tree, a certificate is an **affine**
multiplier `y_i(x) = y0_i + Σ_k Y1[i][k](x_k − c_k)` with `y_i ≥ 0` on *X* and `N(x) ≺ 0` on *X*, where
`N(x) = Σ y_i(x) C_i(x)`, `C_i = sym(ν_i η_{i→π(i)}ᵀ)`. `core.check` verifies this rigorously:

1. `y_i(X) ≥ 0` by ball arithmetic;
2. entrywise mean value theorem: `N(x) = N(c) + Σ_k (x_k − c_k) ∂_k N(ξ)`, with
   `∂_k N = Σ_i (Y1[i][k] C_i + y_i ∂_k C_i)`; a ball enclosure `G_k` of `∂_k N` over *X* is computed;
   with `D_k = mid G_k`, `E_k = rad G_k`: `N(x) = L(δ) + Rem`, `L(δ) = N(c) + Σ δ_k D_k`, `δ = x − c`,
   `|Rem_ab| ≤ Σ_k h_k (E_k)_ab`, hence `‖Rem‖₂ ≤ r :=` Frobenius norm of that bound;
3. `−L(δ) − r·I ≻ 0` is checked at the 8 corners `δ = (±h_1, ±h_2, ±h_3)` by Sylvester's criterion in
   ball arithmetic; `L` is affine, so it then holds on the whole box, whence `N(x) ≺ 0` on *X*.

Ball arithmetic is `python-flint`'s `arb` (Arb library): every operation returns a ball containing the
exact result; a comparison of balls is true only if it holds for all points. Box corners and multipliers
are binary floating-point numbers, i.e. exact dyadic rationals; the chart polynomials have integer
coefficients and are evaluated without divisions.

`search.py` (untrusted) builds the cover by bisection; multipliers come from a small SDP (`propose.py`).
`verify.py` (trusted, independent of the search code) re-traverses the stored bisection trees, recomputes
every box itself, and re-checks every certificate; skip leaves of chart Z are accepted only if
`lo_u·lo_v·lo_w > 1 − ρ₀` holds exactly.

**The cover in this repository** (`ρ₀ = 1/2`):

| chart | top cells | boxes (nodes) | certificates | skip leaves | max. bisection depth |
|---|---|---|---|---|---|
| Z | 8×8×8 on `[0,1]³` | 32 136 | 575 441 | 870 | 13 |
| P | 4×8×8 on `[0,½]×[0,1]²` | 37 120 | 763 649 | – | 12 |

The search reported no failed (box, tree) pair. `compact.py` (untrusted) rounds the multipliers to 6
significant digits where the rounded certificate still verifies; `verify.log` is the output of the final
run of `verify.py` on the files shipped here.

## 4. Running it

    python shapes.py          # print the chart polynomials (derivation + identity checks)
    python test_shapes.py     # numerical cross-check of the chart data against direct geometry
    python search.py Z ; python search.py P      # (re)generate certs_Z.json.gz, certs_P.json.gz
    python verify.py          # THE PROOF CHECK: exit status 0 iff all certificates verify
    python sanity.py 2000000  # independent Monte-Carlo sanity check (random pentahedra)

## 5. What the computer check does and does not cover

Trusted: `shapes.py`, `core.py`, `verify.py`, the Python interpreter, sympy (polynomial algebra) and
python-flint/Arb (ball arithmetic). The mathematical reductions of Sections 1–2 are proved by hand above
and are not machine-checked. The result has not been independently refereed.
