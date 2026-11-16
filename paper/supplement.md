# Supplement to "Monostatic polyhedra: no mono-stable pentahedron, and $9\le V_U\le14$"

This supplement collects experiments and side results that are not needed for Theorems A–C. Unless
stated otherwise they are numerical observations, not proofs. Every polyhedron that is called
*certified* below was checked in exact rational arithmetic by `verify/exact_verify.py`.

## S1. The vertex lower bound at $V=9$

Proving $V_U\ge 10$ requires excluding $V=9$. The relevant counts are:

* 50 simplicial 2-spheres with 9 vertices; the list is complete by the flip-closure check,
  `lb_joint/flip_closure.py`, with 667 flips;
* 9,915,284 (sphere, root, tree) cases.

**Stage 1** (`lb_joint/run_V9_stage1.log`). The single-cone certificate of Proposition 4.2 excludes
9,163,771 cases (92.4%). The certificate file is 1.46 GB and is not part of the repository.

**Stage 2.** For the remaining 751,513 trees three further sound ingredients were tried. All of them are
implemented, with their soundness arguments, in `lb_joint/verify_joint.py`.

1. *Coupled cones* (record kind `certG`; Proposition 4.3).
2. *Hybrid certificates* (record kind `certH`). The triangles avoiding $v$ are split into "good" and
   "bad" ones. A good $T$ needs a positive definite block of Proposition 4.2. A bad $T$ needs, for every
   tetrahedron $T'$ of the cone decomposition from a second apex $q$, a positive definite block that
   combines both directions $\mu$ and $\mu'$. This is about 4 times faster than `certG`, because no
   coupling matrix is needed.
3. *Steepest-neighbour inequalities* (S1). Assign to each vertex $i\ne v$ its steepest shadowing
   neighbour $p(i)\in\arg\max_{k\sim i} r_k\cdot r_i$. This is again an acyclic shadowing tree, and it
   satisfies $r_i\cdot r_k\le r_i\cdot r_{p(i)}$ for all neighbours $k$. Every configuration has such a
   tree, so these inequalities may be added for each enumerated tree.

**Obstruction** (`lb_joint/v9_obstruction.py`, `.log`, `.json`). We sampled 3,000 of the stage-2 trees
uniformly (seed 2026):

| outcome | trees |
|---|---|
| excluded by a hybrid certificate | 2,911 |
| excluded only with all nine cone decompositions + S1 | 22 |
| **not excluded by any certificate available here** | **67** |

For the 67 survivors the SDP relaxation is feasible even with the cone identities from all nine vertices
and the S1 inequalities. No certificate of these types exists for them.

Extrapolated, about 17,000 trees (0.17% of all cases) survive. All survivors in the sample have maximal
tree depth 6 or 7. A typical one is sphere 0, root 0, with the path
$0\leftarrow7\leftarrow6\leftarrow5\leftarrow4\leftarrow3\leftarrow2\leftarrow1$ and the leaf
$0\leftarrow8$. The 294 sampled trees of depth 8 were all excluded.

The relaxation never uses the actual volume weights of the tetrahedra (only their nonnegativity) or the
convexity of the configuration beyond the shadowing inequalities. These are the natural places to look
for stronger valid constraints.

## S2. Thin plates in more detail

### S2.1 Flat pyramids: 16 vertices

Take the outline polygon $B$ in a plane and one apex lifted slightly above a point $A'\in B$. The centre
of mass of the pyramid projects to $(3g+A')/4$, where $g$ is the centroid of $B$. Put $A'$ at the bottom
of a double tight spiral and optimise the spiral angles. This works for outlines with 15 vertices (arm
splits 7+7, 8+6, 9+5, 10+4), giving certified mono-unstable polyhedra with 16 vertices
(`verify/certificates/seed_V16*`).

For 14-gons the numerical optimum stays infeasible for every arm split. The previous record, 18 vertices
[DKLRV20], is a flat pyramid over Conway and Guy's 17-gon; certified versions are in
`verify/certificates/seed_V18*`.

### S2.2 Plates with heights on the outline

The construction of Theorem C (§5.2 of the paper) needs no apex. CMA-ES on the planar surrogate (the
thin-plate limit of §5.1) and then in full 3D finds 14-vertex examples for arm splits 7+6, 8+5 and 9+4.
A version with a larger safety margin is `search_vertices/certified/V14_spiral76_d7.json`. A certified
15-vertex body is `search_vertices/certified/V15_insert_d7.json`.

### S2.3 Thirteen vertices: negative evidence for thin plates

* *Search.* 328 optimisation runs with about $3\cdot10^7$ evaluations covered all arm splits and many
  height patterns. They ended at clean minimax points with all margins equal and negative: the best
  angular margin was $-7.9\cdot10^{-3}$, against $+1.07\cdot10^{-2}$ at 14 vertices.
* *Flat limit.* In three dimensions the optimiser always drives the thickness to zero: the flat limit is
  locally optimal.
* *Ideal thickness bound.* For a given outline $K$, the thickness function that pushes the weighted
  centroid furthest in a given direction is a frustum: a plateau below a level, with the minimal concave
  extension above it. Even with this ideal thickness, which would cost extra vertices, symmetric tight
  spirals need a 13-gon outline. In units with $|T-O|=1$ the ideal centroid height is $+0.0285$ for the
  best 11-gon and $-0.0034$ for the best 13-gon; it must be negative.

None of this is a proof, and genuinely three-dimensional (non-plate) designs were explored only lightly.

## S3. Polarity between the two record bodies

**Lemma (polarity).** Let $0\in\operatorname{int}P$, $P=\{x:\ell_i\cdot x\le1\}$, and let
$P^\circ=\operatorname{conv}\{\ell_i\}$ be the polar body. Face $i$ of $P$ is stable with respect to the
origin iff $\ell_i$ is an unstable vertex of $P^\circ$ with respect to the origin. Both conditions say
$\ell_j\cdot\ell_i<|\ell_i|^2$ for all $j\ne i$ (Lemma 2.2 of the paper).

So stable faces of $P$ correspond to unstable vertices of $P^\circ$, but with respect to the centre of
mass of $P$, which in general is not the centre of mass of $P^\circ$. This is why the face problem and the
vertex problem are different problems.

Reshetov's 14-face mono-stable polyhedron is a long thin spindle whose 14 faces are all "lateral" and
slightly tilted. Its polar about the centre of mass is a thin plate with 14 vertices, mono-unstable with
respect to the origin. The centre of mass of that plate misses the origin by only 0.7% of the diameter,
and a local optimisation turns it into a genuine mono-unstable polyhedron
(`search_vertices/certified/V14_resh_d6.json`).

In the limits "infinitely long spindle" and "infinitely thin plate" both problems become planar problems
about a double spiral. They differ only in which point must lie below the spiral centre: the centroid of
the plate, or the Santaló point of the plate.

## S4. Faces: Reshetov's bodies and improved margins

Reshetov's four published bodies, with 14–17 faces and integer plane data [Res14], were re-verified
exactly (`search_faces/certified/reshetov1*.json`).

We measure the safety margin as in [Res14]: the distance of the centroid projections from the face
boundaries, divided by $\operatorname{vol}^{1/3}$. On this measure the published bodies have margins
$1.4\cdot10^{-4}$ (14 faces) and $6.4\cdot10^{-4}$ (15 faces); his reported best values are
$1.2\cdot10^{-3}$ and $4.2\cdot10^{-3}$.

Local optimisation reaches $2.2\cdot10^{-3}$ (14 faces) and $5.5\cdot10^{-3}$ (15 faces) in floating
point. The bodies `search_faces/certified/improved14.json` and `improved15.json` are rounded to 4-digit
integer plane data and are exactly certified mono-stable, with margins $1.2\cdot10^{-3}$ and
$4.6\cdot10^{-3}$.

About 100 optimisation runs aimed at 13 faces failed. This is a small effort compared with Reshetov's
$10^{10}$ evaluations and proves nothing.

## S5. Simplices

For $d$-simplices, the conditions of §3.1 of the paper involve only the Gram matrix of the scaled facet
normals, which has zero row sums. "Is there a monostatic $d$-simplex with a given rolling pattern" then
reduces to one small semidefinite feasibility problem.

All 719 rooted trees on 10 nodes are infeasible numerically, consistent with the absence of a monostatic
9-simplex. This case was left open by Dawson, Finbow and Mak. It is not new: it is equivalent to
[PRDB23, Cor. 6], which states that mono-unstable 0-skeleta need 11 vertices in every dimension.

## S6. Exploration code

| directory | content |
|---|---|
| `lowerbound/` | replication of the SDP method of [BDPR24]; where it stops at $V=7$ |
| `lb_joint/` | certificate generators (`certify_joint.py`, `general_cert.py`, `hybrid_cert.py`, `stage2_hybrid.py`), the trusted verifier and the V = 9 experiments |
| `penta/` | `search.py`, `propose.py`, `compact.py` (untrusted, they produce the cover); `sanity.py` (Monte-Carlo test) |
| `search_vertices/`, `search_faces/`, `explore/` | optimisation searches and prototypes |
| `audit/` | the independent referee pass on Theorem A (scripts, logs, `REPORT.md`) |
