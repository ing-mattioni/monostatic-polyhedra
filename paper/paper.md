# Monostatic polyhedra: no mono-stable pentahedron, and $9\le V_U\le 14$

*Claudio Mattioni (claudio.mattioni@gmail.com)* &nbsp; *Preprint, November 2026; not yet peer reviewed.*
*Code, certificates and logs: archived at [doi:10.5281/zenodo.24286668](https://doi.org/10.5281/zenodo.24286668)
(version v1.0; see Appendix A).*

**Abstract.** A homogeneous convex polyhedron is *mono-stable* if it can rest on only one of its faces,
and *mono-unstable* if it can be balanced on only one of its vertices. No homogeneous tetrahedron is
mono-stable (Conway, Goldberg and Guy, 1969), and Goldberg and Guy (1969) asked whether every homogeneous
pentahedron has at least two stable faces. We prove that it does, so a homogeneous mono-stable polyhedron
has at least 6 faces; the best known example, due to Reshetov, has 14. For the vertex version we prove that
every homogeneous convex polyhedron with at most 8 vertices has at least two unstable vertices, which
raises the lower bound of Bozóki, Domokos, Papp and Regős from 7 to 9. We also give a mono-unstable polyhedron
with 14 vertices and integer coordinates; the previous record had 18. The two lower bounds are
computer-assisted.
* *Pentahedra.* We separate the affine shape of the polyhedron from the metric. Over a three-parameter
  space of affine normal forms, each of the 375 possible "rolling patterns" is excluded by a $3\times3$
  linear matrix inequality. These inequalities are verified on a finite cover with rigorous ball
  arithmetic: 1,339,090 certificates.
* *Vertices.* We strengthen the Gram-matrix method of Bozóki et al. in two ways: the centroid identity is
  tested in a direction that depends on the configuration, and the identities of several cone
  decompositions are coupled. Integer certificates, checked in exact arithmetic, exclude all 458,394
  combinatorial cases with at most 8 vertices.

The 14-vertex example is a thin plate shaped as a double Conway spiral, and its verification fits on one
page.

---

## 1. Introduction

### 1.1 The problems

Let $P\subset\mathbb R^3$ be a convex polyhedron of uniform density, with centre of mass $c$.

* A face $F$ is a **stable equilibrium** if the orthogonal projection of $c$ onto the plane of $F$ lies
  in the relative interior of $F$: then $P$ can rest on $F$ on a horizontal table.
* A vertex $v$ is an **unstable equilibrium** if the plane through $v$ orthogonal to $v-c$ meets $P$
  only in $v$: then $P$ can be balanced on $v$.

In the generic (nondegenerate) case the numbers $S$, $H$, $U$ of stable, saddle-type and unstable
equilibria satisfy $S-H+U=2$. $P$ is **mono-stable** if $S=1$ and **mono-unstable** if $U=1$.

Conway and Guy [CG66] asked to show that every homogeneous tetrahedron has at least two stable faces, and
to construct a mono-stable polyhedron. The published solution [CGG69] settled the first claim (see also
[Daw85]) and gave Guy's 19-face mono-stable polyhedron. In the same note Goldberg and Guy asked, among other questions:

> *Can it be proved that each uniform pentahedron has at least two stable faces?*

(This is the first problem of Goldberg and Guy in Lángi's list [Lán21, slide 4].)
Shephard [She68] asked for the mono-stable polyhedron with the fewest faces. The examples were later
reduced to 18 faces (Bezdek [Bez11]) and to 14 faces (Reshetov [Res14], found by a large computer search).
So the minimum number $F_S$ of faces of a homogeneous mono-stable polyhedron satisfied $5\le F_S\le 14$.

For the vertex version, a homogeneous mono-unstable polyhedron with 18 vertices was constructed by
Domokos, Kovács, Lángi, Regős and Varga [DKLRV20]. Bozóki, Domokos, Kovács and Regős [BDKR22] and
Papp, Regős, Domokos and Bozóki [PRDB23] settled the version with point masses at the vertices (the
minimum is 11 vertices). For homogeneous polyhedra, Bozóki, Domokos, Papp and Regős [BDPR24] proved
$V_U\ge7$ with a semidefinite relaxation of the Gram matrix of the vertices. So $7\le V_U\le18$.

Related recent work includes Lángi's solution of several problems of Conway on monostable polyhedra and
his existence proof of a polyhedral Gömböc [Lán22], the first explicit mono-monostatic polyhedron
(56,946 faces, verified in exact rational arithmetic) by Schettini Gherardini [Sch26], and a physical
(necessarily inhomogeneous) monostable tetrahedron [ADD25].

### 1.2 Results

**Theorem A (pentahedra).** Every homogeneous convex polyhedron with five faces has at least two stable
faces. Hence $F_S\ge 6$.

**Theorem B (lower bound for vertices).** Every homogeneous convex polyhedron with at most eight vertices
has at least two unstable vertices. Hence $V_U\ge 9$.

**Theorem C (upper bound for vertices).** The convex hull of the 14 points

    (  -4, 1000,  -1)  (-477,  638,  -3)  (-633,  344,  -2)  (-657,   62,  -4)  (-577, -182,   5)
    (-432, -355, -11)  (-242, -454,  13)  ( 181, -490, -12)  ( 359, -426,  12)  ( 518, -295,  -6)
    ( 630, -103,   5)  ( 670,  128,   3)  ( 626,  374,   2)  ( 474,  636,   2)

is a homogeneous mono-unstable polyhedron with 14 vertices; its only unstable vertex is the first point.
Hence $V_U\le 14$.

| quantity | previously known | now |
|---|---|---|
| $F_S$, fewest faces of a mono-stable polyhedron | $5\le F_S\le 14$ [CGG69, Res14] | $6\le F_S\le 14$ |
| $V_U$, fewest vertices of a mono-unstable polyhedron | $7\le V_U\le 18$ [BDPR24, DKLRV20] | $9\le V_U\le 14$ |

To our knowledge Theorem A is the first improvement since 1969 of the lower bound in the original
Conway–Guy problem, and it answers the question of Goldberg and Guy quoted above. Theorem B also settles the
case $V=7$, for which [BDPR24] conjectured non-existence.

### 1.3 Methods

**Theorem A: affine shape times metric (§3).** Centres of mass are affinely covariant. So instead of
deforming the polyhedron we fix an affine normal form and let the inner product of $\mathbb R^3$ vary.
Pentahedra have a three-dimensional space of affine normal forms. For a fixed normal form and a fixed
rooted spanning tree of the face-adjacency graph (a "rolling pattern"), the condition that every face
other than the root rolls to its parent is a linear matrix inequality in the $3\times3$ inverse Gram matrix
of the metric. A dual certificate excludes it. We cover the shape space by boxes, and on each box we
certify all 375 rooted trees with multipliers that depend affinely on the shape, using a mean-value
enclosure in ball arithmetic. Earlier impossibility proofs concern simplices. All simplices are affinely
equivalent, so for them the same reduction leaves no shape parameter at all.

**Theorem B: stronger Gram-matrix certificates (§4).** If $v$ is the only unstable vertex, every other
vertex is "shadowed" by an adjacent vertex: this gives a rooted spanning tree, and with it linear
inequalities on the Gram matrix $R\succeq0$ of the vertex vectors [BDPR24]. The centre of mass enters
through the cone decomposition of $P$ from a vertex,
$\sum_T\operatorname{vol}(T)\sum_{k\in T}r_k=0$. [BDPR24] use one scalar consequence of this vector
identity, its pairing with $r_v$, and test each tetrahedron separately. We add two pieces of information:

1. The identity is paired with a direction $w=\sum_k\mu_kr_k$ that depends linearly on the
   configuration, with one coefficient vector $\mu$ for all tetrahedra of the decomposition. Finding $\mu$
   and the multipliers is a single convex semidefinite programme per tree.
2. The same polyhedron is tiled from every vertex, with different unknown volume weights. The
   identities of two cone decompositions are coupled by matrices $W_q$ with $\sum_qW_q=0$.

At $V=7$ our numerical reproduction (`lowerbound/`) found feasible solutions, to numerical tolerance,
for 1,426 of the 139,636 (tree, tetrahedron) systems of [BDPR24]. With ingredient 1, every case with
$V\le7$ and 97.6% of the trees with $V=8$ are excluded; ingredient 2 excludes the remaining 10,207 trees
with $V=8$. At $V=9$, the sampled searches in §6.1 leave some trees uncertified and suggest that stronger
geometric constraints may be useful. They do not establish feasibility of the relaxation or
non-existence of certificates for those trees.

**Theorem C: thin plates (§5).** As a polyhedron is flattened, the shadowing conditions become planar
conditions with respect to a weighted centroid of the outline, while the heights still shape the mass
distribution. This turns the construction into a planar problem with a free weight. We solve it with a
double "Conway spiral" whose outline vertices carry small alternating heights.

### 1.4 Nature of the proofs

The proofs of Theorems A and B consist of hand-proved reductions (§§2–4, written out in full here) and
finite computations. The computations are checked by small verifiers that do not import the code that
produced the certificates:
* Theorem B: about 550 lines of standard-library Python, exact integer arithmetic.
* Theorem A: about 540 lines, using `sympy` for exact polynomial algebra and Arb ball arithmetic through
  `python-flint`.

Theorem C is checked by a 100-line integer-only script (§5.3). Appendix A lists the trusted code, the
commands, the SHA-256 hashes of the certificate files and the running times. Independent
implementations that share no code with ours re-checked all 1,339,090 certificates of Theorem A in exact
rational arithmetic, and all 458,394 cases of Theorem B with different algorithms. Neither found a failure
(Appendix A.3).

The research, the code and the first draft of this manuscript were produced with substantial assistance
of an AI system (Claude, Anthropic); see Appendix A.4.

---

## 2. Preliminaries

All polyhedra are convex, three-dimensional and of uniform density. Equilibria are taken with respect to
the centre of mass, and coordinates are chosen with the centre of mass at the origin unless stated
otherwise.

**Lemma 2.1 (shadowing [BDPR24]).** Let $r_1,\dots,r_V$ be the vertices of $P$, with the centre of mass at
the origin. Vertex $r_i$ is an unstable equilibrium if and only if $(r_j-r_i)\cdot r_i<0$ for every vertex
$r_j$ adjacent to $r_i$. If $(r_j-r_i)\cdot r_i\ge0$ ("$j$ shadows $i$"), then $|r_j|>|r_i|$.

*Proof.* $P$ is contained in its tangent cone $r_i+\operatorname{cone}\{r_j-r_i:j\sim i\}$ at $r_i$.

* If $(r_j-r_i)\cdot r_i<0$ for all neighbours $j$, then $(x-r_i)\cdot r_i<0$ for every
  $x\in P\setminus\{r_i\}$. So the plane through $r_i$ orthogonal to $r_i$ meets $P$ only in $r_i$.
* If some neighbour $j$ has $(r_j-r_i)\cdot r_i\ge 0$, there are two cases. If the value is $0$, the plane
  contains the point $r_j\neq r_i$ of $P$. If it is positive, the plane separates $r_j$ from the interior
  point $0$, so it meets the interior of $P$.

Finally,
$0<|r_j-r_i|^2=|r_j|^2-|r_i|^2-2(r_j-r_i)\cdot r_i\le|r_j|^2-|r_i|^2$. $\square$

In particular, the vertex farthest from $c$ has no shadowing neighbour, so it is unstable and $U\ge1$.
Consequently, if $v$ is the only unstable vertex, choosing for every $i\neq v$ an adjacent shadowing vertex
$p(i)$ gives a map along which $|r|$ strictly increases. So $p$ has no cycles: it is a spanning tree of the
graph of $P$, oriented towards $v$.

**Lemma 2.2 (faces).** Let $P=\{x:\ell_i\cdot x\le 1,\ i=1,\dots,F\}$ with the centre of mass at the
origin, each $\ell_i$ defining a face. The foot of the perpendicular from the origin onto the plane of
face $i$ is $\ell_i/|\ell_i|^2$. Face $i$ is stable if and only if $\ell_j\cdot\ell_i<|\ell_i|^2$ for every
face $j$ that shares an edge with face $i$.

*Proof.* Inside the plane of face $i$, the face is cut out by the half-planes $\ell_j\cdot x\le1$ of the
faces $j$ adjacent along its edges. A point of the plane is in the relative interior of the face if and
only if it satisfies these inequalities strictly. Faces meeting face $i$ in a vertex, or not at all,
impose no further condition. $\square$

---

## 3. No mono-stable pentahedron (Theorem A)

### 3.1 Affine shape times metric

Let $P_0=\{x:\ell_i\cdot x\le1\}$ be a polyhedron with centre of mass at the origin, and let $A$ be an
invertible linear map. Then $AP_0=\{x:(A^{-\mathsf T}\ell_i)\cdot x\le1\}$, and its centre of mass is
again the origin, because the centre of mass of a uniform body is affinely covariant. Put
$W=A^{-1}A^{-\mathsf T}\succ0$. By Lemma 2.2, face $i$ of $AP_0$ is stable iff
$(A^{-\mathsf T}\ell_j)\cdot(A^{-\mathsf T}\ell_i)<|A^{-\mathsf T}\ell_i|^2$ for all adjacent $j$, that is,
iff

$$\ell_i^{\mathsf T}W(\ell_j-\ell_i)<0\quad\text{for all faces } j \text{ adjacent to } i. \tag{3.1}$$

Every $W\succ0$ arises in this way (take $A=W^{-1/2}$), and every polyhedron affinely equivalent to $P_0$
is a translate of some $AP_0$. Call face $i$ *$W$-stable* if (3.1) holds. Then:

> *Theorem A holds iff for every convex pentahedron $P_0$ in a fixed set of affine normal forms (§3.2)
> and every $W\succ0$, at least two faces are $W$-stable.*

**Lemma 3.1 (rolling).** Fix $P_0$ and $W\succ0$, and write $\phi_i=\ell_i^{\mathsf T}W\ell_i$. If face
$i$ is not $W$-stable, it has an adjacent face $j$ with $\ell_i^{\mathsf T}W(\ell_j-\ell_i)\ge0$, and then
$\phi_j>\phi_i$. Consequently the face with the largest $\phi_i$ is $W$-stable. If only one face $s$ is
$W$-stable, choosing such a $j=\pi(i)$ for each $i\ne s$ defines a spanning tree of the face-adjacency
graph rooted at $s$.

*Proof.* The first claim is (3.1). Expanding the positive quantity
$(\ell_j-\ell_i)^{\mathsf T}W(\ell_j-\ell_i)$ gives
$\phi_j=\phi_i+2\ell_i^{\mathsf T}W(\ell_j-\ell_i)+(\ell_j-\ell_i)^{\mathsf T}W(\ell_j-\ell_i)>\phi_i$.
So $\pi$ strictly increases $\phi$: it has no cycles, and every orbit ends at the only face without a
parent, $s$. $\square$

**Proposition 3.2 (certificate).** Fix $P_0$ and a rooted spanning tree $(s,\pi)$. Let $\nu_i=\alpha_i\ell_i$
and $\eta_i=\beta_i(\ell_{\pi(i)}-\ell_i)$ with $\alpha_i,\beta_i>0$. Suppose there are $y_i\ge0$
($i\ne s$) such that

$$N=\sum_{i\ne s}y_i\,\operatorname{sym}(\nu_i\eta_i^{\mathsf T})\ \prec\ 0,\qquad
\operatorname{sym}(ab^{\mathsf T})=\tfrac12(ab^{\mathsf T}+ba^{\mathsf T}).$$

Then no $W\succ0$ realises the tree, i.e. has
$\ell_i^{\mathsf T}W(\ell_{\pi(i)}-\ell_i)\ge0$ for all $i\ne s$.

*Proof.* For $W\succ0$ and $N\prec0$ we have $\operatorname{tr}(NW)<0$. But
$\operatorname{tr}(NW)=\sum_iy_i\,\nu_i^{\mathsf T}W\eta_i=\sum_iy_i\alpha_i\beta_i\,\ell_i^{\mathsf T}W(\ell_{\pi(i)}-\ell_i)\ge0$
if the tree is realised. $\square$

So Theorem A follows if, for every normal form and every rooted spanning tree of its face-adjacency
graph, such a certificate exists.

### 3.2 Combinatorics and affine normal forms

**Lemma 3.3.** A convex polyhedron with five faces is combinatorially either a *triangular prism* (two
triangles $A$, $B$, three quadrilaterals $Q_1,Q_2,Q_3$, six vertices) or a *quadrilateral pyramid* (five
vertices).

*Proof.* Euler's formula gives $E=V+3$. Since $2E\ge3V$, $V\le6$. Since $2E\ge3F=15$, $V\ge5$. A face has at
most four neighbours, hence at most four sides. With $V=6$ ($E=9$, $2E=18$) the side counts sum to 18, so
they are $3,3,4,4,4$. Each quadrilateral is then adjacent to all four other faces, so the quadrilaterals
pairwise share an edge (three edges) and the remaining six edges join a quadrilateral to a triangle. The
two triangles are therefore disjoint, which is the triangular prism. With $V=5$ ($E=8$, $2E=16$)
they are $3,3,3,3,4$, which is the pyramid. $\square$

**Lemma 3.4 (prism type).** Let $P$ be of prism type with triangles $A=A_1A_2A_3$, $B=B_1B_2B_3$ and lateral
edges $A_iB_i$ on lines $L_i$. Then $L_1,L_2,L_3$ are either concurrent in a point $O$ or parallel.

* In the concurrent case, $O$ lies on none of the segments $[A_i,B_i]$, and after possibly exchanging the
  names $A$, $B$ we have $A_i=O+t_i(B_i-O)$ with $0<t_i<1$ for all $i$.
* In the parallel case, $A_i=B_i+s_id$ with a common vector $d$ and all $s_i>0$.

*Proof.* $L_i$ and $L_j$ lie in the plane of the lateral face that contains both edges, so the three lines
are pairwise coplanar; they do not all lie in one plane.

* If $L_1\cap L_2=\{O\}$, suppose $L_3$ met $L_1$ and $L_2$ in two different points, or met one of them and
  were parallel to the other. Then $L_3$ would lie in the plane of $L_1,L_2$. So $L_3\ni O$.
* If $L_1\parallel L_2$, the same argument shows that $L_3$ is parallel to both.

In a lateral face, $L_i$ and $L_j$ carry opposite sides $A_iB_i$, $A_jB_j$ of a convex quadrilateral. The
lines through two opposite sides of a convex quadrilateral meet, if at all, outside both sides, so
$O\notin[A_i,B_i]$. Moreover $A_i$ is nearer to $O$ than $B_i$ iff $A_j$ is nearer to $O$ than $B_j$;
otherwise the sides $A_iA_j$ and $B_iB_j$ would cross. The parallel case is analogous. $\square$

A quadrilateral pyramid with apex $X$ has the same description with one lateral edge of length zero.
Choose two opposite sides of the base as $A_1B_1$, $A_2B_2$ and put $A_3=B_3=X$. Then $A=A_1A_2X$,
$B=B_1B_2X$, $Q_1=A_2B_2X$, $Q_2=A_1B_1X$, and $Q_3$ is the base. Here $Q_1,Q_2$ meet only in $X$, as do
$A,B$.

*Face-adjacency graph.* In the prism type every two faces share an edge except $A$ and $B$; the graph is
$K_5$ minus an edge. It has $3\cdot5^2=75$ spanning trees and therefore $5\cdot75=375$ rooted spanning trees.
The adjacency graph of a pyramid is a subgraph, so its rooted spanning trees are among these 375. Renaming
the lateral edges permutes $Q_1,Q_2,Q_3$, and exchanging $A$ and $B$ is an automorphism of $K_5-AB$; both
leave the set of 375 trees invariant. We therefore always certify **all 375 trees**, and labellings may
be chosen freely.

**The two charts.** Fix $\rho_0=\tfrac12$. Throughout, $Q_i$ denotes the lateral face that does not contain
lateral edge $i$.

* **Chart Z (cone frame).** $O=0$, $B_i=e_i$, $A_i=t_ie_i$ with $0<t_1\le t_2\le t_3\le1$; $t_3=1$ is a
  pyramid. We use the toric coordinates $t_3=w$, $t_2=wv$, $t_1=wvu$, with
  $(u,v,w)\in D_Z=\{(u,v,w)\in[0,1]^3: uvw\le1-\rho_0\}$.
* **Chart P (prism frame).** $B_1=(0,0,0)$, $B_2=(1,0,0)$, $B_3=(0,1,0)$, $G=(\tfrac13,\tfrac13)$,
  $\sigma=\rho\,(1,p,pq)$, $s=(1,p,pq)$ and $A_i=\big((1-\sigma_i)B_i+\sigma_iG,\ s_i\big)$, with
  $(\rho,p,q)\in D_P=[0,\rho_0]\times[0,1]^2$.
  * For $\rho>0$ the lateral lines meet in the apex $(G,1/\rho)$ and $A_i=O+(1-\sigma_i)(B_i-O)$: this is an
    affine image of the shape $t=(1-\rho,\,1-\rho p,\,1-\rho pq)$ of chart Z.
  * For $\rho=0$ it is the prism with parallel lateral edges and heights $s$.

**Lemma 3.5 (coverage).** Every convex pentahedron is affinely equivalent, after a suitable labelling of
its faces, to the shape of a point of $D_Z$ with $u,v,w>0$, or to the shape of a point of $D_P$ with $p>0$.

*Proof.* The four cases below cover every pentahedron. In each, relabel the lateral edges so that the
parameters are sorted, and use the free choice of an affine frame.

* *Prism type, concurrent lines.* By Lemma 3.4, the affine map with $O\mapsto0$ and $B_i\mapsto e_i$ gives
  $0<t_1\le t_2\le t_3<1$. The vectors $B_i-O$ are independent because the lines are not coplanar.
  * If $t_1\le1-\rho_0$: $u=t_1/t_2$, $v=t_2/t_3$, $w=t_3$ is a point of $D_Z$.
  * Otherwise $\rho=1-t_1\in(0,\rho_0)$, $p=(1-t_2)/(1-t_1)$ and $q=(1-t_3)/(1-t_2)$ lie in $(0,1]$, and the
    chart-P shape $(\rho,p,q)$ is affinely equivalent to it.
* *Prism type, parallel lines.* Map $B_1,B_2,B_3$ to the standard triangle and $d$ to a vertical vector,
  scaled so that the largest $s_i$ equals 1. With $s=(1,p,pq)$ sorted, this is $(0,p,q)\in D_P$.
* *Pyramid whose base has two non-parallel opposite sides.* Take them as the lateral edges $A_1B_1$,
  $A_2B_2$; their lines meet in a point $O$ of the base plane, outside both sides. The map $O\mapsto0$,
  $B_1\mapsto e_1$, $B_2\mapsto e_2$, $X\mapsto e_3$ gives $t_3=1$ and $0<t_1\le t_2<1$, by the convexity
  argument of Lemma 3.4 applied to the base. Then proceed as in the concurrent case; in chart P this
  gives $q=0$.
* *Pyramid over a parallelogram.* Take two parallel opposite sides as $A_1B_1$, $A_2B_2$. This is the
  parallel case with $s_3=0$, i.e. $(0,p,0)\in D_P$. $\square$

The closed sets $D_Z$, $D_P$ also contain degenerate points that are not pentahedra (for example $u=0$
in $D_Z$, or $p=0$ in $D_P$). The certificates below hold on the closed sets, but they are *interpreted*
only at genuine pentahedra.

**The data $\nu$ and $\eta$.** In a chart, face $k$ is $n_k\cdot x\le d_k$, with $n_k$, $d_k$ polynomial in
the chart parameters, and the centre of mass is $g=m/\Lambda$, with $\Lambda=24\operatorname{vol}(P)$ and
$m=24\operatorname{vol}(P)\,g$ also polynomial (Appendix B). Put

$$\Delta_k=d_k\Lambda-n_k\cdot m=\Lambda\,(d_k-n_k\cdot g).$$

In coordinates centred at $g$ the plane of face $k$ is $\ell_k\cdot y=1$ with
$\ell_k=n_k/(d_k-n_k\cdot g)=\Lambda n_k/\Delta_k$. Therefore

$$n_k=\frac{\Delta_k}{\Lambda}\,\ell_k,\qquad
\Delta_in_j-\Delta_jn_i=\frac{\Delta_i\Delta_j}{\Lambda}\,(\ell_j-\ell_i).$$

We take $\nu_k=n_k$, and $\eta_{(i,j)}=(\Delta_in_j-\Delta_jn_i)/\omega_{ij}$ for adjacent $i<j$ (and
$\eta_{(j,i)}=-\eta_{(i,j)}$). Here $\omega_{ij}$ is the largest monomial in the chart variables dividing all
three components. It is $u$ for the pair $(Q_1,A)$ of chart Z and $1$ otherwise. Dividing by it keeps the
vector from vanishing identically on a boundary face of the parameter box. Both vectors are finally
multiplied by positive integers to clear denominators.

**Lemma 3.6 (positivity).** At every genuine pentahedron of either chart, $\nu_k$ is a positive multiple
of $\ell_k$ and $\eta_{(i,j)}$ is a positive multiple of $\ell_j-\ell_i$.

*Proof.* $\Lambda>0$. Also $\Delta_k=\Lambda(d_k-n_k\cdot g)>0$, because the centre of mass is an
interior point and $n_k\cdot x\le d_k$ on $P$. The monomial $u$ is positive at genuine shapes. $\square$

The polynomials are derived exactly with `sympy` (`penta/shapes.py`). Chart Z is short enough to check by
hand (Appendix B). The audit re-derived the proportionality factors of both charts from the vertices, with
an independent centroid computation (Appendix A.3).

### 3.3 Rigorous certificates on boxes

For a box $X=\hat x+[-h_1,h_1]\times[-h_2,h_2]\times[-h_3,h_3]$ of a chart and a rooted tree $(s,\pi)$, a
certificate consists of **affine** multipliers

$$y_i(x)=y_i^0+\sum_kY_{ik}(x_k-\hat x_k)\qquad(i\ne s).$$

They must satisfy $y_i\ge0$ on $X$ and $N(x)=\sum_iy_i(x)C_i(x)\prec0$ on $X$, where
$C_i=\operatorname{sym}(\nu_i\eta_{i\to\pi(i)}^{\mathsf T})$. By Lemma 3.6 and Proposition 3.2 the tree is
then excluded at every genuine shape in $X$.

**Lemma 3.7 (box check).** Let $G_k$ be an interval matrix containing $\partial N/\partial x_k$ on all of
$X$. Write $D_k=\operatorname{mid}G_k$ and $E_k=\operatorname{rad}G_k$ (entrywise), and put
$\bar E=\sum_kh_kE_k$ and $r=\|\bar E\|_F$. Suppose that, at each of the 8 corners
$\delta\in\{\pm h_1\}\times\{\pm h_2\}\times\{\pm h_3\}$,

$$-\Big(N(\hat x)+\sum_k\delta_kD_k\Big)-r\,I\succ0 .$$

Then $N(x)\prec0$ for all $x\in X$.

*Proof.* By the mean value theorem applied to each entry on the convex box,
$N(x)=N(\hat x)+\sum_k\delta_k\,\partial_kN(\xi)$, with $\delta=x-\hat x$ and a point $\xi$ that depends on the
entry. Hence $N(x)=L(\delta)+\mathrm{Rem}$, where $L(\delta)=N(\hat x)+\sum_k\delta_kD_k$ and
$|\mathrm{Rem}_{ab}|\le\bar E_{ab}$, so $\|\mathrm{Rem}\|_2\le\|\bar E\|_F=r$. The map
$\delta\mapsto-L(\delta)-rI$ is affine and positive definite at the corners, hence on the whole box, by
convexity of the positive definite cone. So $N(x)\preceq L(\delta)+rI\prec0$. $\square$

Note that $\partial_kN=\sum_i\big(Y_{ik}C_i+y_i\,\partial_kC_i\big)$. The check is carried out in Arb ball
arithmetic [Joh17] (`penta/core.py`):
* every operation returns a ball that contains the exact result, and a comparison is accepted only if it
  holds for all points of the balls;
* $y_i\ge0$ on $X$ is checked on the ball $y_i^0+\sum_kY_{ik}[-h_k,h_k]$;
* positive definiteness at the corners is checked by Sylvester's criterion;
* box corners and multipliers are binary floating-point numbers, i.e. exact dyadic rationals;
* the chart polynomials have integer coefficients and are evaluated without divisions.

### 3.4 The cover

The certificates are stored on bisection trees over the cells of a grid: $8\times8\times8$ cells on
$[0,1]^3$ for chart Z, and $4\times8\times8$ cells on $[0,\rho_0]\times[0,1]^2$ for chart P. A certificate
stored at a node is valid on that node's box. The verifier `penta/verify.py` performs these checks:

1. Both files use the same dyadic $\rho_0\in(0,1)$, and the top cells are exactly the cells of the stated
   grid.
2. In every cell the bisection tree is traversed. The boxes of the children are computed by the verifier
   (bisection of the stated axis at the exact midpoint), so the leaves of a cell partition it.
3. Along every root-to-leaf path, each of the 375 rooted trees receives a certificate that passes
   Lemma 3.7 on the box of the node where it is stored.
4. The only exception is a *skip leaf* of chart Z. It is accepted only if
   $\mathrm{lo}_u\,\mathrm{lo}_v\,\mathrm{lo}_w>1-\rho_0$ holds exactly, so the whole box lies outside $D_Z$.

So for every point of $D_Z$ and of $D_P$ and every rooted tree, some certificate holds on a box containing
that point.

| chart | grid | boxes | certificates | skip leaves | max. depth |
|---|---|---|---|---|---|
| Z | $8\times8\times8$ on $[0,1]^3$ | 32,136 | 575,441 | 870 | 13 |
| P | $4\times8\times8$ on $[0,\frac12]\times[0,1]^2$ | 37,120 | 763,649 | – | 12 |

The cover was produced by an untrusted search (`penta/search.py`) that bisects a box until a small
semidefinite programme yields multipliers passing Lemma 3.7. The multipliers were then rounded to six
significant digits wherever the rounded certificate still passes (`penta/compact.py`).

*Proof of Theorem A.* Let $P$ be a homogeneous convex pentahedron and suppose it had at most one stable
face.
1. By Lemma 3.5, $P$ is affinely equivalent to a genuine shape $x$ in $D_Z$ or $D_P$, with some labelling
   of the faces.
2. By §3.1 there is $W\succ0$ for which this shape has at most one $W$-stable face.
3. By Lemma 3.1 there is a rooted spanning tree $(s,\pi)$ of $K_5-AB$ realised by $W$.
4. By §3.4 some box $X\ni x$ carries a certificate for $(s,\pi)$. By Lemmas 3.6 and 3.7 it gives $y(x)\ge0$
   with $N(x)\prec0$.
5. Proposition 3.2 then says that $W$ cannot realise $(s,\pi)$, a contradiction. $\square$

*Checks beyond the proof.*
* `penta/verify.py` passes on the final certificate files, also under `python -O` (it contains no
  `assert` statements).
* An independent exact-rational checker re-verified every certificate (Appendix A.3).
* Random sampling of more than $3\cdot10^6$ pentahedra, with two independent codes, and
  evolution-strategy optimisation found no pentahedron with fewer than two stable faces. The second-best
  margin approaches 0 only in the flat degenerate limit.

---

## 4. At least nine vertices (Theorem B)

### 4.1 The Gram-matrix reduction

Let $P$ have $V$ vertices $r_1,\dots,r_V$, centre of mass $0$, and exactly one unstable vertex $v$.

1. Triangulate each face of $P$ by diagonals, adding no vertices. The boundary becomes a simplicial
   2-sphere $G$ on the $V$ vertices: a maximal planar graph with $3V-6$ edges and $2V-4$ triangles. It is
   a simplicial complex because two faces $F\ne F'$ of $P$ meet in the empty set, a vertex or an edge.
   * A diagonal of $F$ is never an edge of $P$: an edge of $P$ through two vertices of $F$ meets $F$ in a
     face of $F$ containing both, so they are adjacent in $F$.
   * No diagonal or triangle occurs in two faces.
2. By Lemma 2.1 there is a spanning tree $p$, oriented towards $v$, that uses edges of $P$ (hence of $G$)
   and satisfies $(r_{p(i)}-r_i)\cdot r_i\ge0$ for $i\ne v$.
3. Let $R=(r_a\cdot r_b)_{a,b}$ be the Gram matrix. Then $R\succeq0$, $R\ne0$, and with
   $S_{i\to j}=E_{ii}-\operatorname{sym}(E_{ij})$ (where $E_{ij}$ are the matrix units)

$$\langle S_{i\to p(i)},R\rangle=|r_i|^2-r_i\cdot r_{p(i)}\le0\qquad(i\ne v). \tag{4.1}$$

**Lemma 4.1 (cone identities).** For every vertex $q$, the tetrahedra $T=\operatorname{conv}(\{q\}\cup F)$,
with $F$ a triangle of $G$ not containing $q$, have disjoint interiors and fill $P$ (some may be flat).
Hence

$$\sum_{T}\operatorname{vol}(T)\sum_{k\in T}r_k=4\operatorname{vol}(P)\,c=0. \tag{4.2}$$

*Proof.*
* *Cover.* For $x\in\operatorname{int}P$, the ray from $q$ through $x$ leaves $P$ at a unique boundary
  point $y$, and the open segment $(q,y)$ is interior. If $y$ lay in a face containing $q$, the segment
  $[q,y]$ would lie in that face, so on $\partial P$. Hence $x\in\operatorname{conv}(q,F)$ for a triangle
  $F\ni y$ of a face not containing $q$.
* *Disjoint interiors.* If $x$ is interior to both $\operatorname{conv}(q,F_1)$ and
  $\operatorname{conv}(q,F_2)$, then $y$ is in the relative interior of $F_1$ and of $F_2$, so $F_1=F_2$.
* *Flat tetrahedra.* $\operatorname{conv}(q,F)$ is flat exactly when $F$ lies in a face of $P$ containing
  $q$.

The identity follows from the additivity of first moments, since the centroid of a tetrahedron is the
mean of its vertices. $\square$

Write $1_T\in\{0,1\}^V$ for the indicator vector of the four vertices of $T$.

**Proposition 4.2 (one cone, joint direction).** Fix $(G,v,p)$. Suppose there are $\mu\in\mathbb R^V$ and,
for each tetrahedron $T$ of the cone decomposition from $v$, multipliers $y^T_i\ge0$ ($i\ne v$) with

$$M_T=\sum_{i\ne v}y^T_iS_{i\to p(i)}+\operatorname{sym}(\mu1_T^{\mathsf T})\ \succ\ 0\qquad\text{for all }T.$$

Then no polyhedron realises $(G,v,p)$.

*Proof.* Put $w=\sum_k\mu_kr_k$. Pairing (4.2) with $w$ gives
$\sum_T\operatorname{vol}(T)\langle\operatorname{sym}(\mu1_T^{\mathsf T}),R\rangle=0$. Since the volumes
are nonnegative and not all zero, some $T^*$ with $\operatorname{vol}(T^*)>0$ has
$\langle\operatorname{sym}(\mu1_{T^*}^{\mathsf T}),R\rangle\le0$. With (4.1) this gives
$\langle M_{T^*},R\rangle\le0$, which contradicts $M_{T^*}\succ0$, $R\succeq0$, $R\neq0$. $\square$

**Proposition 4.3 (several cones).** Fix $(G,v,p)$ and a set $Q$ of vertices. Suppose there are symmetric
matrices $W_q$ ($q\in Q$) with $\sum_{q\in Q}W_q=0$, vectors $\mu^{(q)}$, and multipliers $y^{q,T}\ge0$,
such that for every $q\in Q$ and every tetrahedron $T$ of the cone decomposition from $q$

$$W_q+\sum_{i\ne v}y^{q,T}_iS_{i\to p(i)}+\operatorname{sym}(\mu^{(q)}1_T^{\mathsf T})\ \succ\ 0.$$

Then no polyhedron realises $(G,v,p)$.

*Proof.* For each $q$ choose $T_q$ as in the previous proof, now using (4.2) for the apex $q$. Add the $|Q|$
matrices for $(q,T_q)$: the $W_q$ cancel, the sum is positive definite, and its pairing with $R$ is
$\le0$. $\square$

For fixed $(G,v,p)$, finding the data of Proposition 4.2 is one convex semidefinite programme, because
$\mu$ is shared by all blocks. Proposition 4.3 with fixed $Q$ is one as well.

**What is new compared with [BDPR24].** In our notation, the method of [BDPR24] is Proposition 4.2 with
$\mu=e_v$, i.e. $w=r_v$, decided tetrahedron by tetrahedron.
* Allowing a general $\mu$ uses the full vector identity (4.2). It says that the centre of mass lies at
  the origin in every direction, not only along $r_v$. A direction is never fixed in advance; it is a
  linear function of the unknown configuration, the same for all tetrahedra, chosen per tree.
* Proposition 4.3 adds a second piece of geometric information. The same body is tiled from different
  vertices with different unknown weights, and both identities must hold simultaneously.

Neither ingredient uses the actual tetrahedron volumes or the convexity of the configuration beyond the
shadowing inequalities. These are natural places to strengthen the method (§6.1).

### 4.2 Finitely many cases

The cases are triples (sphere $G$, root $v$, rooted spanning tree $p$).

* **Spheres.** For $V=4,\dots,8$ there are 1, 1, 2, 5 and 14 simplicial 2-spheres up to isomorphism. The
  verifier checks that each listed complex is a simplicial 2-sphere and that the listed ones are pairwise
  non-isomorphic. Completeness is proved twice, without relying on the literature count (OEIS A000109):
  * `lb_joint/enumerate_spheres.py` enumerates all labelled simplicial 2-spheres on $\{0,\dots,V-1\}$
    exhaustively: all edge sets with $3V-6$ edges, and all exact double covers by triangles. It compares
    their number with $\sum_kV!/|\operatorname{Aut}(G_k)|$ over the listed spheres. The counts are 1, 10,
    195, 5,712 and 223,440 labelled spheres for $V=4,\dots,8$, in agreement in every case.
  * `lb_joint/flip_closure.py` checks that the list is closed under edge flips up to isomorphism (159 flips
    for $V=8$). By Wagner's theorem [Wag36] any two triangulations of the sphere with $V$ vertices are
    connected by flips, so the list is complete.
* **Roots.** Every vertex is either a listed root or is mapped to one by a listed automorphism, which is
  verified to preserve edges and triangles.
* **Trees.** For every (sphere, root), all maps $p$ from the other vertices to neighbours without directed
  cycles are enumerated exhaustively. Their number is compared with Kirchhoff's matrix-tree determinant.

### 4.3 Verification and proof of Theorem B

Every certificate has integer entries. Every matrix that must be positive definite is checked exactly, by
fraction-free Gaussian elimination: all leading principal minors must be positive. The verifier
`lb_joint/verify_joint.py` uses only the Python standard library and exact integers.

| $V$ | spheres | (sphere, root) | rooted trees | certificates |
|---|---|---|---|---|
| 4 | 1 | 1 | 16 | Prop. 4.2 |
| 5 | 1 | 2 | 150 | Prop. 4.2 |
| 6 | 2 | 4 | 1,392 | Prop. 4.2 |
| 7 | 5 | 16 | 24,795 | Prop. 4.2 |
| 8 | 14 | 63 | 432,041 | Prop. 4.2 for 421,834; Prop. 4.3 with $Q=\{v,q\}$ for 10,207 |

For $V=8$, the trees assigned a two-cone certificate in the stored files are exactly those in which some
triangle avoiding $v$ has all three vertices at depth $\ge5$ in the tree. Each is excluded using a second
apex $q$. This observation about the stored certificates does not assert that a one-cone certificate
is impossible for any of these trees.

*Proof of Theorem B.* Suppose $P$ has $V\le8$ vertices and at most one unstable vertex. Then $V\ge4$, and
by the remark after Lemma 2.1 there is exactly one unstable vertex $v$.
1. §4.1 produces a triangulated boundary $G$ and a shadowing tree $p$.
2. By §4.2, an isomorphism maps $(G,v,p)$ onto a listed sphere, a listed root and an enumerated tree. The
   isomorphism preserves (4.1) and (4.2), because these depend only on the labelling.
3. That case carries a verified certificate of Proposition 4.2 or 4.3, which is a contradiction. $\square$

---

## 5. A mono-unstable polyhedron with 14 vertices (Theorem C)

### 5.1 The flat-limit principle

Let $x_i(\varepsilon)=(p_i,\varepsilon h_i)$ with $p_i\in\mathbb R^2$, $h_i\in\mathbb R$, and let
$P_\varepsilon=\operatorname{conv}\{x_i(\varepsilon)\}$. Then $P_\varepsilon$ is a linear image of $P_1$,
so its centre of mass is $(c',\varepsilon c_h)$ with $c'$ independent of $\varepsilon$.

* $c'$ is the centroid of the outline $K=\operatorname{conv}\{p_i\}$ weighted by the *thickness function*
  of $P_1$: the upper concave envelope of the heights minus the lower convex envelope.
* $(x_j-x_i)\cdot(x_i-c)=(p_j-p_i)\cdot(p_i-c')+\varepsilon^2(h_j-h_i)(h_i-c_h)$.

Suppose every vertex but one is strictly shadowed, in the plane and with respect to $c'$, by a neighbour
along the outline. Then $P_\varepsilon$ is mono-unstable for all small $\varepsilon>0$, because the
vertex farthest from $c$ is always unstable. The heights do not enter the shadowing conditions in the
limit; they only shape the mass distribution.

*Tight Conway spirals.* In the plane, a chain $T=P_0,P_1,\dots$ about a point $O$, with polar angles
increasing by $\varphi_k$ and radii $r_k=r_{k-1}\cos\varphi_k$, has $\angle OP_kP_{k-1}=90^\circ$. So $P_k$ is
on the verge of being shadowed by $P_{k-1}$. Moving the reference point from $O$ by a small $\delta$ with
$\delta\cdot(P_{k-1}-P_k)<0$ for all $k$ makes every shadowing strict.

Two such arms on either side of $T$ form a convex polygon that is mono-unstable with respect to every
reference point slightly "below" $O$. The polygon's own centroid is never such a point: a homogeneous
planar convex body has at least four equilibria with respect to its centroid, two of them unstable
[DPR94]. This is why the third dimension is needed: its job
is to pull the weighted centroid $c'$ down.

### 5.2 The construction

Give the outline vertices themselves small heights:
* close to zero at $T$ and on the upper parts of both arms;
* alternating in sign, with growing amplitude, along the last three or four vertices of each arm.

The plate is then thick near the bottom of the spiral (locally it resembles an antiprism) and thin near
$T$, without any extra vertex. Numerical optimisation (CMA-ES on the planar surrogate, then in full 3D)
finds 14-vertex examples robustly, for arm sizes $7+6$, $8+5$ and $9+4$. The polyhedron of Theorem C is
such a plate with arms $7+6$, rounded to small integers.

![Top view of the polyhedron of Theorem C](../examples/V14_plate.svg)

### 5.3 Verification

The convex hull has 24 triangular faces. Each is given below as an outward-oriented triple of point
indices (0-based, in the order of Theorem C):

    (1,0,5) (1,6,0) (5,0,7) (6,13,0) (0,9,7) (0,13,9) (3,2,1) (2,4,1) (3,1,5) (4,6,1) (3,4,2) (3,5,4)
    (4,5,6) (5,7,6) (6,7,8) (6,8,13) (7,9,8) (8,9,10) (8,10,11) (8,11,12) (8,12,13) (9,11,10) (9,12,11) (13,12,9)

The check, in integers:

1. *Facets.* For each triple $abc$, every other point $r_k$ satisfies
   $((b-a)\times(c-a))\cdot(r_k-a)<0$. So each triangle is a facet of the hull, and no other point lies on
   its plane.
2. *All facets.* Every edge of a triangle lies in exactly two triangles, with opposite orientations, and
   every point is a corner. Each edge of the hull lies in exactly two facets, and the facets are connected
   through edges, so these are all the facets. The polyhedron has $V=14$, $E=36$, $F=24$.
3. *Volume and centre of mass.* With outward orientation, $D=\sum\det(a,b,c)=6\operatorname{vol}=90\,438\,831$
   and $C=\sum\det(a,b,c)(a+b+c)=24\operatorname{vol}\cdot c=(-948\,565\,765,\ 351\,349\,436,\ 24\,043\,586)$.
   Therefore $g_i=4Dr_i-C=4D(r_i-c)$ is an integer vector.
4. *Equilibria.* By Lemma 2.1, vertex $i$ is not unstable iff some neighbour $j$ has
   $s(i,j)=(r_j-r_i)\cdot g_i\ge0$. The witnesses are listed below. For vertex 0, all six neighbours
   $1,5,6,7,9,13$ have $s<0$; the largest value is $s(0,1)=-130\,591\,697\,728\,185$.

| $i$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| witness $j$ | 0 | 1 | 2 | 3 | 4 | 5 | 8 |
| $s(i,j)$ | 2249229145437 | 909298793102 | 541515043628 | 456118206758 | 180771149147 | 46753300478 | 351981830882 |

| $i$ | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|
| witness $j$ | 9 | 10 | 11 | 12 | 13 | 0 |
| $s(i,j)$ | 488326853183 | 512415268962 | 462201185684 | 597174613290 | 789702550352 | 1200773226408 |

So vertex 0 is the only unstable vertex. `examples/V14_short_check.py` (standard library, integers only)
performs exactly steps 1–4. The general exact classifier `verify/exact_verify.py` computes the hull itself
and classifies all equilibria: $(S,H,U)=(2,1,1)$ with no degeneracies, consistent with $S-H+U=2$. This
proves Theorem C. $\square$

---

## 6. Limits of the methods and open problems

### 6.1 Where the Gram-matrix method stops: $V=9$

There are 50 simplicial 2-spheres with 9 vertices; the list is complete by the flip-closure check. There
are 9,915,284 (sphere, root, tree) cases.

* Proposition 4.2 settles 9,163,771 of them (92.4%).
* For the remaining 751,513 we tried Proposition 4.3, and a cheaper sound variant that uses the second
  cone only for the "bad" tetrahedra of the first. We also tried the steepest-neighbour inequalities
  $r_i\cdot r_k\le r_i\cdot r_{p(i)}$, which are valid once each vertex is assigned its steepest shadowing
  neighbour (see `lb_joint/verify_joint.py`).

On a random sample of 3,000 of these trees, 2,911 were excluded by the variant, a further 22 were excluded
using all nine cone decompositions plus the steepest-neighbour inequalities, and **67 remained
uncertified by the tested searches**. The last outcome records failure to obtain an exact certificate;
it does not prove that the relaxation is feasible or that no certificate of these types exists.
Extrapolating the observed failure rate suggests about 17,000 cases (0.17% of the total) would remain
uncertified by the same procedure; this is a sample-based estimate, not an exhaustive count. All 67
sampled survivors have maximal tree depth 6 or 7, while the 294 sampled trees of depth 8 were all
excluded. Relations between the tetrahedron weights and the geometry, or further convexity constraints,
are possible directions for extending the method. Details are in the supplement.

### 6.2 Searches for 13 vertices and 13 faces

We did not find a mono-unstable polyhedron with 13 vertices.
* A search over thin plates made 328 optimisation runs and about $3\cdot10^7$ evaluations. It ended at
  minimax points with all margins equal and negative.
* An exploratory frustum thickness model helps explain the centroid shift in spiral constructions,
  but supplies no certified lower bound on the number of vertices.

About 100 optimisation runs aimed at a mono-stable polyhedron with 13 faces also failed. Compared with
Reshetov's $10^{10}$ evaluations this is a small effort. Neither search proves anything. Both, together
with a polarity between Reshetov's 14-face spindle and our 14-vertex plates, are described in the
supplement.

### 6.3 Open problems

1. Close $9\le V_U\le14$. For the lower bound, resolve the uncertified $V=9$ cases (§6.1), possibly with
   stronger geometric constraints. For the upper bound,
   is there a mono-unstable polyhedron that is not a thin plate?
2. Close $6\le F_S\le14$. Theorem A uses a three-dimensional shape space; for hexahedra it is
   six-dimensional, with more combinatorial types. Is there a Gram-type relaxation for faces, as there is
   for vertices?
3. Prove that thin plates need 14 vertices. This is a planar variational problem.

---

## Appendix A. Verification and reproducibility

### A.1 Trusted code and what it checks

| claim | trusted code | dependencies | checks |
|---|---|---|---|
| Theorem A | `penta/verify.py`, `penta/core.py`, `penta/shapes.py` (≈540 lines) | `sympy` 1.14, `python-flint` 0.9 (Arb) | §3.4, items 1–4 |
| Theorem B | `lb_joint/verify_joint.py` (≈550 lines) | standard library | §4.2–4.3 |
| Theorem B, sphere lists | `lb_joint/enumerate_spheres.py`, `lb_joint/flip_closure.py` | standard library | §4.2 |
| Theorem C | `examples/V14_short_check.py` (≈100 lines); `verify/exact_verify.py` | standard library | §5.3 |

The certificate generators, searches and optimisers (`lb_joint/certify_joint.py`, `penta/search.py`,
`search_*/`, ...) are not trusted. None of the verifiers imports them, and none contains an `assert` used
as a check. Each verifier exits with status 0 only if every check passes. Tampered certificate files
(dropped cells or trees, illegal skip leaves, NaN or negated multipliers, changed matrices) are rejected
(`lb_joint/test_verifier.py`, `lb_joint/test_certH.py`, `audit/tamper_after_fix.log`).

### A.2 Reproduction

All files are archived at [doi:10.5281/zenodo.24286668](https://doi.org/10.5281/zenodo.24286668)
(version v1.0). In the unpacked archive:

    python -m venv venv-verify
    venv-verify/bin/python -m pip install -r requirements-verify.txt   # sympy, mpmath, python-flint (pinned)
    venv-verify/bin/python reproduce.py                                # all checks; --quick: seconds

`reproduce.py` checks the SHA-256 hashes below, then runs every verifier as a separate process.

| file | SHA-256 |
|---|---|
| `lb_joint/certs_V4.jsonl.gz` … `certs_V8.jsonl.gz`, `exp1_V8_certG.jsonl.gz` | see `CHECKSUMS.sha256` |
| `penta/certs_Z.json.gz` (33.8 MB, 575,441 certificates) | `3fda788e324a303b798668cdbc6ed08f35e3e432ff2b37c2663ff95c3605df78` |
| `penta/certs_P.json.gz` (42.9 MB, 763,649 certificates) | `796f057b6d428c25086664b0e0ee7e084a0044b8a72b62ff296b689d10028a85` |
| `lb_joint/certs_V8.jsonl.gz` (48.4 MB) | `089be2691aa339e6d11e061e07a4a06820a7c62236fecdffe65b6ad0789c5aa2` |
| `lb_joint/exp1_V8_certG.jsonl.gz` (5.8 MB) | `05851381f6b07dcecd9c6d6f2d4a82627686ad4667b4bbffe22943d8e0de0bc6` |

Running times on a 16-thread desktop are:

* Theorem C: under 1 s;
* Theorem B: about 3.5 min for $V=8$, and seconds for $V\le7$;
* the sphere enumeration: 30 s for $V=8$;
* Theorem A: about 2 min with 14 processes.

A complete run (about 7 min) from a fresh clone of the repository, in a fresh virtual environment, is
recorded in `reproduce_clean_run.log`.

### A.3 Independent re-checks

**Theorem A.** A separate referee pass was done by an independent AI agent (a separate Claude instance,
`audit/REPORT.md`). It was asked to break the proof and shares no code with `penta/`. In summary:

* It re-derived the face covectors from the vertices, using its own divergence-theorem centroid, and
  proved exactly that $\nu_k\parallel\ell_k$ and $\eta_{ij}\parallel(\ell_j-\ell_i)$ with positive
  factors, in both charts.
* It re-verified **every** certificate of both charts with an exact-rational implementation: exact
  Taylor expansion of the polynomials on each box, a remainder bound by absolute coefficients, and exact
  Sylvester tests at the corners. All 1,339,090 passed, and no subdivision was needed.
* It confirmed the coverage logic (grid, skip leaves, all 375 trees on every path) and the Z→P
  equivalence at random overlap points.
* It tried to falsify the theorem with $1.1\cdot10^6$ random pentahedra and with optimisation.

It found no mathematical error and no counterexample. Its one code finding was fixed: checks written as
`assert` are disabled under `python -O`.

**Theorem B.** A second independent agent (`audit/theoremB/REPORT.md`) reviewed §2 and §4 and wrote its own
exact checker from the mathematical statements, without reading `verify_joint.py`'s code. The checker
is 788 lines, uses only the standard library, and works differently in each part:
* it proves sphere completeness by generating all triangulations through coning fan disks, with
  canonical forms;
* it recomputes the automorphism orbits;
* it enumerates trees as edge subsets;
* it tests positive definiteness by rational $LDL^{\mathsf T}$.

It confirmed all 458,394 cases with zero invalid records, and 30 negative controls behaved as expected.
It found no error and no missing case. Its three expository suggestions are incorporated in Lemma 2.1,
§4.1 and Lemma 4.1.

### A.4 Use of AI

The problem choice, the methods, the code, the computations and the first draft of this text were
produced in an AI-assisted research session (Claude, Anthropic) directed by the author. The
audits of A.3 were performed by separate AI instances. No human expert has yet checked the proofs. The
hand-written reductions (§§2–4 and Appendix B) are the parts that most need human review;
`REVIEW_GUIDE.md` lists them with pointers to the code.

## Appendix B. Chart Z explicitly

With $t=(uvw,\,vw,\,w)$ and $\tau=t_1t_2t_3=uv^2w^3$:

$$P=\{x:\ x_1,x_2,x_3\ge0,\ \ x_1+x_2+x_3\le1,\ \ x_1/t_1+x_2/t_2+x_3/t_3\ge1\},$$

* $\operatorname{vol}(P)=(1-\tau)/6$;
* $g_i=(1-\tau t_i)/(4(1-\tau))$, the big simplex minus the small one;
* $\Lambda=4(1-\tau)$ and $m_i=1-\tau t_i$.

The faces and the data are:

| face $k$ | $n_k\cdot x\le d_k$ | $\Delta_k$ |
|---|---|---|
| $Q_i$ ($i=1,2,3$) | $-x_i\le0$ | $1-\tau t_i$ |
| $A$ | $-(x_1+u\,x_2+uv\,x_3)\le-uvw$ | $1+u+uv-4uvw+\tau t_1$ |
| $B$ | $x_1+x_2+x_3\le1$ | $1-4\tau+\tau(t_1+t_2+t_3)$ |

Here $\tau t_1=u^2v^3w^4$, $\tau t_2=uv^3w^4$ and $\tau t_3=uv^2w^4$. Then $\nu_k=n_k$ and
$\eta_{(i,j)}=\Delta_in_j-\Delta_jn_i$, except that $\eta_{(Q_1,A)}$ is divided by $u$. For example
$\eta_{(Q_1,Q_2)}=(1-\tau t_2,\ -(1-\tau t_1),\ 0)$. `python penta/shapes.py` prints all of them, and
those of chart P, whose vertices and planes are given in §3.2 and in `penta/shapes.py`.

The orientation choices in `shapes._tet_data` can also be checked explicitly. In either chart the
polyhedron is decomposed into the tetrahedra below: cone from $A_1$ over the opposite triangle $B$ and
over a triangulation of the opposite quadrilateral $B_2B_3A_3A_2$. For an ordered tetrahedron $(a,b,c,d)$,
write $D(a,b,c,d)=\det(b-a,c-a,d-a)$. Direct substitution of the vertices gives:

| ordered tetrahedron | chart Z | chart P |
|---|---|---|
| $(B_1,B_2,B_3,A_1)$ | $-(1-t_1)$ | $1$ |
| $(A_1,B_2,B_3,A_3)$ | $-t_1(1-t_3)$ | $(1-\rho)pq$ |
| $(A_1,B_2,A_3,A_2)$ | $-t_1t_3(1-t_2)$ | $(1-\rho)(1-\rho pq)p$ |

For genuine prism-type shapes these signs are respectively negative in Z and positive in P, since
$0<t_1\le t_2\le t_3<1$ and $0\le\rho\le1/2$, $p,q>0$. At a pyramid boundary ($t_3=1$ or $q=0$), the
middle tetrahedron is flat and contributes zero. The signs do not reverse; the polynomial volume and
first-moment formulas extend to these boundaries by continuity. Thus the signs selected at a sample
point in the code remain valid, including for genuine pyramids.

## Appendix C. Certificate files

* **Theorem B**, `lb_joint/certs_V*.jsonl.gz`. One JSON record per line:
  * graph records: edges, triangles, roots, automorphisms;
  * certificate records `cert` (Proposition 4.2: integer $\mu$ and, per triangle, integer $y^T$);
  * records `certG` (Proposition 4.3: per apex $q$ an integer symmetric $W_q$, $\mu^{(q)}$ and blocks).

  The semantics of every field is stated in the docstring of `verify_joint.py`.
* **Theorem A**, `penta/certs_{Z,P}.json.gz`. Per chart: $\rho_0$, the grid, and one bisection tree per
  cell. A node holds optional certificates (tree id in $[0,375)$, $y^0\in\mathbb R^4$,
  $Y\in\mathbb R^{4\times3}$ as binary floats) and either two children with a bisection axis, or nothing
  (a leaf), or the skip flag (chart Z only).

## Data availability

All code, certificate files and logs, together with the sources of this manuscript, are archived on
Zenodo as version v1.0 of `monostatic-polyhedra`:
[doi:10.5281/zenodo.24286668](https://doi.org/10.5281/zenodo.24286668). The code is released under the
MIT licence; the manuscript, figures, certificates and data under CC BY 4.0 (see `LICENSE.md`).

## References

* [ADD25] G. Almádi, R. J. MacG. Dawson, G. Domokos, *Building a monostable tetrahedron*,
  arXiv:2506.19244 (2025).
* [BDKR22] S. Bozóki, G. Domokos, F. Kovács, K. Regős, *Mono-unstable polyhedra with point masses have at
  least 8 vertices*, Int. J. Solids Struct. 234–235 (2022).
* [BDPR24] S. Bozóki, G. Domokos, D. Papp, K. Regős, *The smallest mono-unstable, homogeneous convex
  polyhedron has at least 7 vertices*, arXiv:2401.17906 (2024).
* [Bez11] A. Bezdek, *On stability of polyhedra*, talk, Workshop on Discrete Geometry, Fields Institute,
  Toronto, 2011.
* [CG66] J. H. Conway, R. K. Guy, *Stability of polyhedra* (Problem 66-12), SIAM Rev. 8 (1966) 381.
* [CGG69] J. H. Conway, M. Goldberg, R. K. Guy, Problem 66-12 (solution and comments), SIAM Rev. 11
  (1969) 78–82, doi:10.1137/1011014.
* [Daw85] R. J. MacG. Dawson, *Monostatic simplexes*, Amer. Math. Monthly 92 (1985) 541–546; with
  W. Finbow and P. Mak, *Monostatic simplexes II*, Geom. Dedicata 70 (1998); with W. Finbow,
  *Monostatic simplexes III*, Geom. Dedicata 84 (2001).
* [DKLRV20] G. Domokos, F. Kovács, Z. Lángi, K. Regős, P. T. Varga, *Balancing polyhedra*, Ars Math.
  Contemp. 19 (2020) 95–124.
* [DPR94] G. Domokos, J. Papadopulos, A. Ruina, *Static equilibria of planar, rigid bodies: is there
  anything new?*, J. Elasticity 36 (1994) 59–66.
* [Joh17] F. Johansson, *Arb: efficient arbitrary-precision midpoint-radius interval arithmetic*, IEEE
  Trans. Comput. 66 (2017) 1281–1292.
* [Lán21] Z. Lángi, *On some problems of Conway and Guy on monostable polyhedra*, talk at the 8th
  European Congress of Mathematics, Portorož, 21 June 2021 (slides).
* [Lán22] Z. Lángi, *A solution to some problems of Conway and Guy on monostable polyhedra*, Bull. London
  Math. Soc. 54 (2022) 501–516, arXiv:2008.02090.
* [PRDB23] D. Papp, K. Regős, G. Domokos, S. Bozóki, *The smallest mono-unstable convex polyhedron with
  point masses has 8 faces and 11 vertices*, Eur. J. Oper. Res. 310 (2023) 511–517.
* [Res14] A. Reshetov, *A unistable polyhedron with 14 faces*, Int. J. Comput. Geom. Appl. 24 (2014)
  39–59.
* [Sch26] T. Schettini Gherardini, *An explicit mono-monostatic polyhedron*, arXiv:2609.07827 (2026).
* [She68] G. C. Shephard, *Twenty problems on convex polyhedra, Part I*, Math. Gaz. 52 (1968) 136–147.
* [Wag36] K. Wagner, *Bemerkungen zum Vierfarbenproblem*, Jber. Deutsch. Math.-Verein. 46 (1936) 26–32.
