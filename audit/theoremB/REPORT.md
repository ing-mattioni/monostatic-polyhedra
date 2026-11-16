# Referee report on Theorem B (at most 8 vertices implies at least two unstable vertices)

*Independent referee pass by a separate AI agent (Claude Opus). The agent was instructed to try to break
the proof, and it neither read nor imported the repository's verifier code. The report is reproduced as
delivered; the agent's runtime environment did not let it write files outside its scripts and logs, so
the main session saved this text. The changes made in response are listed at the end.*

**Scope.** `paper/paper.md` §2 (Lemma 2.1) and §4.1–4.3; the docstring of `lb_joint/verify_joint.py`,
read for the file format only; the certificate files `lb_joint/certs_V4.jsonl.gz` … `certs_V8.jsonl.gz`
and `lb_joint/exp1_V8_certG.jsonl.gz`. I tried to break the proof, and I did not read, import or copy the
code of `verify_joint.py` or of any other repository file.

**Files in `audit/theoremB/`:**

| file | content |
|---|---|
| `check_theoremB.py` | the independent exact checker (standard library only, 788 lines) |
| `recheck.log` | its output for V = 4..8, with the SHA-256 of the inputs and of the checker |
| `tamper.py`, `negative_controls.log` | negative controls (the tampered files are regenerated in `tamper/`) |
| `sanity_sample.py`, `sanity.log` | floating-point sampling and an adversarial search (numpy/scipy) |
| `depth_claim.py`, `depth_claim.log` | a check of a side remark in §4.3 |

## (a) Mathematics

I checked every hand-proved step and found **no error**. Three places are argued too briefly (A2, A3, A4);
each can be fixed in a sentence or two, and none affects correctness.

**A1. Lemma 2.1 is correct.**
* The centre of mass c is interior, so no $r_i$ is 0.
* $P$ lies in the tangent cone $r_i+\operatorname{cone}\{r_j-r_i : j\sim i\}$. So if $(r_j-r_i)\cdot r_i<0$
  for all neighbours $j$, then $(x-r_i)\cdot r_i<0$ for every $x\in P\setminus\{r_i\}$.
* Conversely, take a neighbour $j$ with $(r_j-r_i)\cdot r_i\ge 0$. If the value is 0, the plane contains
  $r_j$. If it is positive, the open segment $(0,r_j)$ lies in the interior of $P$ and crosses the plane.
* $|r_j|^2-|r_i|^2=|r_j-r_i|^2+2(r_j-r_i)\cdot r_i$ gives $|r_j|>|r_i|$.
* Strict and weak inequalities are used consistently: "unstable" is strict, so "not unstable" means
  weakly shadowed ($\ge0$), which is exactly (4.1). The strict definition makes Theorem B the stronger
  statement.
* Adjacency means adjacency in $P$, and the edges of $P$ are edges of $G$.
* $|r|$ strictly increases along $p$, so $p$ has no cycles and its only sink is $v$: it is a spanning tree
  oriented towards $v$.

**A2. Small gap: the case $U=0$ is never excluded explicitly.**
* The proof assumes "exactly one unstable vertex", but the negation of "at least two" is "at most one".
* $U=0$ is excluded by Lemma 2.1: the vertex with the largest $|r_i|$ has no shadowing neighbour.
* Suggestion: add "In particular the vertex farthest from c is unstable, so U ≥ 1."

**A3. Small gap: §4.1 step 1 asserts without proof that the triangulated boundary is a simplicial
complex.**
* A diagonal of a face $F$ is never an edge of $P$. An edge of $P$ through two vertices of $F$ meets $F$
  in a face of $F$ that contains both, so the two vertices are adjacent in $F$.
* A diagonal of $F$ is never a diagonal of another face $F'$, and no triangle appears twice:
  $F\cap F'$ is empty, a vertex or an edge, so it has at most two vertices, and if it has two they are
  adjacent.
* So the result is a simplicial 2-sphere with $3V-6$ edges and $2V-4$ triangles. Its graph contains all
  edges of $P$, hence the tree.

**A4. Lemma 4.1 is correct, but its proof is terse.** A precise version:
* For $x\in\operatorname{int}P$, the ray from $q$ through $x$ leaves $P$ at a unique boundary point $y$,
  and the open segment $(q,y)$ is interior.
* If $y$ lay in a triangle or face containing $q$, the segment $[q,y]$ would lie on $\partial P$.
  So $x\in\operatorname{conv}(q,F)$ for some $F\ni y$ lying in a face not containing $q$: the tetrahedra
  cover $P$.
* The interiors are disjoint: if $x$ is interior to $\operatorname{conv}(q,F_1)$ and to
  $\operatorname{conv}(q,F_2)$, then $y$ lies in the relative interiors of both $F_1$ and $F_2$, so
  $F_1=F_2$.
* $\operatorname{conv}(q,F)$ is flat exactly when $F$ lies in a face of $P$ containing $q$. This covers
  non-simplicial faces and coplanar vertices.
* Additivity of first moments gives $\sum_T\operatorname{vol}(T)\sum_{k\in T}r_k=4\operatorname{vol}(P)\,c=0$,
  and $\sum_T\operatorname{vol}(T)=\operatorname{vol}(P)>0$, so some $T$ has positive volume.
* Suggested rewording: "$F\not\ni q$ unless $x$ is on the boundary" should read "$y$ lies in no triangle
  or face containing $q$ unless $x\in\partial P$".

**A5. Proposition 4.2 is correct.**
* With $S_{i\to j}=E_{ii}-\tfrac12(E_{ij}+E_{ji})$, $\langle S_{i\to j},R\rangle=R_{ii}-R_{ij}$, so the factor
  $\tfrac12$ is right.
* $\langle\operatorname{sym}(\mu1_T^{\mathsf T}),R\rangle=w\cdot s_T$ with $s_T=\sum_{k\in T}r_k$ ($T$
  includes the apex).
* $\mu$ must be shared by all $T$, because only for a single vector $w$ does $\sum_T\operatorname{vol}(T)\,w\cdot s_T$
  vanish. With a different direction for each $T$ the weighted sum carries no sign information.
* Some $T^*$ with $\operatorname{vol}(T^*)>0$ has $w\cdot s_{T^*}\le0$. With $y\ge0$ and (4.1) this gives
  $\langle M_{T^*},R\rangle\le0$, contradicting $M\succ0$, $R\succeq0$, $R\ne0$.
* $y_v=0$ is needed because $v$ has no parent and hence no inequality.
* The files store $2M_T$ as an integer matrix; positive definiteness is unaffected by this scaling.

**A6. Proposition 4.3 is correct.**
* For each apex $q$, choose $T_q$ by the argument of A5. The sum of the chosen matrices is positive
  definite, and its pairing with $R$ is at most $\langle\sum_qW_q,R\rangle=0$.
* $\sum W_q=0$ must hold exactly; the $W_q$ need not be semidefinite.
* The files store $W_q+2(\sum yS+\operatorname{sym}(\mu1_T^{\mathsf T}))$. This is Prop. 4.3 with $W_q/2$
  in place of $W_q$, and the sum is still 0.
* In every block, $y_v=0$ refers to the root $v$.
* $Q$ need not contain $v$, but in all 10,207 records $Q=\{v,q\}$.

**A7. The isomorphism invariance is correct.**
* A bijection carrying triangles to triangles, followed by a listed automorphism, conjugates $R$ by a
  permutation matrix and maps $p$ to $\psi p\psi^{-1}$.
* It carries over the shadowing inequalities (4.1), the cone identities (4.2) for every apex, and the
  triangle sets.
* The isomorphism has to be one of complexes, because the certificates depend on the triangles. For
  $V\ge4$ this coincides with graph isomorphism (Whitney), and my checker works with complexes
  directly.

**A8. I found no polyhedron that the case analysis could miss.**
* $V<4$ is impossible.
* Non-simplicial faces and coplanar vertices are handled by A3 and by the flat tetrahedra of A4.
* Equality in the shadowing condition is allowed, and no general position is assumed anywhere.
* Every triangulation is one of the listed spheres (completeness, below). Every tree of $G$ is
  enumerated, including trees that use diagonals.

**A9. Documentation inconsistencies in the trusted code (not errors in the proof).**
* `verify_joint.py` also accepts "S1" terms and "certH" records, whose soundness is argued only in its
  docstring. I checked both arguments and believe them sound. Neither is used for $V\le8$: every `certG`
  record has `s1 = false` and no S1 terms, and there are no `certH` records. My checker rejects S1 terms
  outright.
* The `verify_joint.py` docstring says that sphere completeness relies on the OEIS count, whereas the
  paper (§4.2) proves completeness with `enumerate_spheres.py` and with `flip_closure.py` (via Wagner's
  theorem). The docstring should be brought in line. (My own completeness proof uses neither the OEIS
  count nor Wagner's theorem.)

**A10. Side claims that agree with my computations.**
* The labelled sphere counts $\sum V!/|\mathrm{Aut}|$ are 1, 10, 195, 5,712 and 223,440, as stated.
* There are 139,636 (tree, tetrahedron) systems at $V=7$, as stated in §1.3.
* The total number of cases is 458,394, as stated.
* The depth-$\ge5$ remark of §4.3 agrees with the files in the half that can be tested without an SDP
  solver: all 10,207 Prop. 4.3 trees have such a triangle, and none of the 421,834 Prop. 4.2 trees do.

## (b) Independent re-check

Exact arithmetic only, no floating point in any decision, 16 worker processes.

1. **Each listed complex is a sphere.** The triangles are distinct; the edge list equals the set of
   triangle sides; each edge lies in exactly two triangles; every vertex link is a single cycle through
   all neighbours; the complex is connected; and $V-E+F=2$. By the classification of surfaces this
   identifies the 2-sphere.
2. **The list is complete, by my own generation.** Lemma (proved in the checker's docstring): every
   simplicial 2-sphere with $n\ge5$ vertices arises from one with $n-1$ vertices by coning a fan disk. The
   disk is one triangle, two triangles sharing an edge, or three consecutive triangles around a vertex.
   The proof:
   * since $\sum\deg=6n-12$, some vertex has degree 3, 4 or 5;
   * after its star is deleted, the chords of the link polygon in the complementary disk cannot cross,
     so a fan of non-edges always fills the hole;
   * coning that fan disk gives back the original sphere.

   All extensions starting from the tetrahedron were generated, each re-validated as a sphere and
   reduced to an exact canonical form: the least relabelled triangle list over all relabellings that
   sort the vertices by (degree, multiset of neighbour degrees). The listed spheres are pairwise distinct
   and, as a set, equal to the generated classes. The class counts are 1, 1, 2, 5 and 14; their agreement
   with OEIS A000109 is incidental and not relied on.
3. **Roots.** The full automorphism group is computed by brute force. Every orbit contains a listed root,
   and every stated automorphism is valid and maps its vertex to a listed root.
4. **Trees.** Spanning trees are enumerated as acyclic $(V-1)$-subsets of the edges, by union–find. Their
   counts equal the matrix-tree determinant, computed by exact rational elimination. Each tree is
   oriented to every listed root.
5. **Certificates.** Every record, not just one per tree, must pass all of these checks:
   * the entries are integers (booleans and floats are rejected);
   * the parent array is a tree of $G$ oriented to the root, with every parent adjacent;
   * $y\ge0$ and $y_{\text{root}}=0$;
   * there is exactly one block per triangle avoiding the apex;
   * for Prop. 4.3: the apexes are distinct, the $W_q$ are symmetric with $\sum W_q=0$ exactly, and there
     are no S1 terms;
   * every block matrix is positive definite, tested by **exact $LDL^{\mathsf T}$ over the rationals**
     (Schur-complement recursion), not by Bareiss elimination.

**Results.** All runs exit with status 0; $V=8$ takes about 230–250 s.

| V | spheres listed / generated | (sphere, root) | trees | Kirchhoff agrees | Prop. 4.2 valid | Prop. 4.3 valid | invalid | PD tests | missing |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 1/1 | 1 | 16 | yes | 16 | – | 0 | 16 | 0 |
| 5 | 1/1 | 2 | 150 | yes | 150 | – | 0 | 375 | 0 |
| 6 | 2/2 | 4 | 1,392 | yes | 1,392 | – | 0 | 5,568 | 0 |
| 7 | 5/5 | 16 | 24,795 | yes | 24,795 | – | 0 | 139,636 | 0 |
| 8 | 14/14 | 63 | 432,041 | yes | 421,834 | 10,207 | 0 | 3,313,358 | 0 |

## (c) Negative controls

`negative_controls.log`: all **30 of 30 controls behaved as expected.**

**Part 1: tampered copies of the V=6 file, each run end to end in its own process.**
* The untampered copy passes.
* Each of these 14 tampers is rejected with the correct reason:
  * certificate tampers: a dropped tree, negated $y$, negated $\mu$, a missing triangle block, a duplicated
    block, a block containing the apex, $y_{\text{root}}\ne0$, a non-integer entry, a non-adjacent parent,
    a parent 2-cycle;
  * sphere and root tampers: an edge-flipped sphere (reported as isomorphic to the other listed sphere),
    a sphere with a triangle removed, a root removed from the list, a wrong stated automorphism.

**Part 2: tampered Prop. 4.3 records at V = 8, through the checker's record test.**
* The untampered record is accepted.
* All 11 tampers are rejected: $W_q$ not summing to 0; $W_q$ asymmetric with the sum still 0; a missing
  block; negated $y$; $y_{\text{root}}\ne0$; negated $\mu^{(q)}$; an injected S1 term; the s1 flag set; a
  duplicated apex; a block triangle containing its apex; a dropped group.

**Part 3: sharpness of the positive-definiteness test.**
* For a Prop. 4.2 record, $\mu_{\text{root}}-26$ is accepted and $\mu_{\text{root}}-27$ is rejected.
* For a Prop. 4.3 record, moving 94 from $W_v[3][3]$ to $W_q[3][3]$ (sum still 0) is accepted, and moving
  95 is rejected.

**Part 4: end to end at V = 8.** The Prop. 4.3 file had 5 tampered records and 1 dropped. The checker
reported exactly the 5 invalid records and the 6 trees left without a certificate, and exited with
status 1.

## (d) Sanity sampling (floating point, not part of the proof)

The test uses the definition directly: vertex $i$ is unstable iff $(r_j-r_i)\cdot(r_i-c)<0$ for every
other vertex $j$.
* **Samples:** 2,400 hulls with exactly 5–8 vertices, 100 per (family, V). The families are uniform in a
  cube, Gaussian, points on a sphere, elongated (aspect up to 100), nearly flat (thickness down to
  $10^{-3}$), and points on needle- or plate-shaped ellipsoids.
* **Unstable vertices:** none of the 2,400 has fewer than two. The minimum per cell is 2 for the
  elongated, flat and ellipsoid families, and 3–4 for the others.
* **Lemma 2.1:** the neighbour test and the definition disagree 0 times.
* **Cone identity:** the relative error is at most about $5\cdot10^{-15}$.
* **Closest sample:** the smallest second-unstable margin was $2.0\cdot10^{-4}$.
* **Adversarial hill-climb:** 12 runs × 400 steps per V maximised the second-smallest margin; a score
  $\ge0$ would mean a mono-unstable example. The best score was $-6.8\cdot10^{-5}$, and scores approach 0
  only in degenerate limits.

## (e) Verdict

**Mathematics.** The reduction is correct: Lemma 2.1, §4.1, Lemma 4.1, Propositions 4.2 and 4.3, and the
isomorphism argument. I found no error and no missing case. There are three small expository gaps
(A2: $U\ge1$; A3: the triangulation is a simplicial complex; A4: a precise proof of Lemma 4.1) and two
documentation inconsistencies in `verify_joint.py` (A9).

**Computation.** The independent checker shares no code with the repository. It uses different methods
for each part:
* sphere completeness: generation by coning fan disks;
* root coverage: orbits recomputed from the full automorphism group;
* trees: enumeration of edge subsets;
* positive definiteness: rational $LDL^{\mathsf T}$.

It confirms every one of the 458,394 cases for V = 4..8 with zero invalid records, and the negative
controls show that it is not vacuous.

**Theorem B is, in my opinion, correct.** It depends only on classical facts used in the hand proofs:
the tangent cone of a polytope, the classification of surfaces, and the Jordan curve theorem in a disk.

---

## Changes made in response (main session, 2026-11-16)

* **A2.** The paper now states after Lemma 2.1 that the vertex farthest from $c$ is unstable, so
  $U\ge1$. The proof of Theorem B starts from "at most one unstable vertex".
* **A3.** §4.1 now proves that the triangulated boundary is a simplicial 2-sphere.
* **A4.** The proof of Lemma 4.1 was replaced by the precise version above.
* **A9.** The docstring and messages of `verify_joint.py` now say that completeness of the sphere lists
  is proved by `enumerate_spheres.py` and `flip_closure.py`. The OEIS count is kept only as an
  additional cross-check. This changes only comments and printed text, no check; `CHECKSUMS.sha256` was
  regenerated.
