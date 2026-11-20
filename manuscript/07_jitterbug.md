# Fuller's Jitterbug and the carrier graph {#sec:jitterbug}

Fuller describes the Jitterbug as a hinged assembly of the equilateral
triangles of the vector equilibrium. His illustration follows the open
cuboctahedron through a regular icosahedral configuration to an
octahedral limit [@fuller1975synergetics, §§460.01–460.08 and fig. 460.08].
Verheyen studies the mechanisms and their motions, and Stachel
constructs this branch by dividing the edges of an octahedron and
rescaling so that the attached triangles remain rigid
[@verheyen1989jitterbug; @stachel1994heureka, pp. 447–449]. We use a
centred coordinate chart for the same branch. The locations of the
joints, the bars, and the edges of the convex hull give three
different ways to read the assembly.

## One linkage, three landmark shapes {#sec:jitterbug-rigid-path}

Let $t$ be the ratio of the smaller to the larger nonzero coordinate.
Before uniform scaling the joints are the cyclic permutations of
$(0,\pm1,\pm t)$. For a chosen bar length $\ell$ set
$$
a=\frac{\ell}{\sqrt{2(1+t^2-t)}},\qquad b=ta,\qquad 0\le t\le1.
$$
For each sign triple $s$ the three positions
$(0,s_y a,s_z b)$, $(s_x b,0,s_z a)$ and $(s_x a,s_y b,0)$
form a triangle. Their pairwise squared distances are
$2(a^2+b^2-ab)=\ell^2$. Adjacent sign triples share a joint. The
triangle planes keep normals $s$, with centroids $(a+b)s/3$.
As the assembly contracts, triangles rotate about the body-diagonal
axes and their centroids slide along those axes. A coordinate reflection
exchanges the two branches obtained by interchanging the cyclic convention.

::: {.proposition #prop:jitterbug-rigid-motion}
The fixed-length linkage has {{JITTERBUG_RIGID_TRIANGLES}} rigid
triangles and {{JITTERBUG_RIGID_BARS}} labelled bars. Its initial joints
at $t=1$ form the cuboctahedral FCC shell up to scale. At
$t=(\sqrt5-1)/2$ their convex hull is a regular icosahedron, which has
the incidence of the OPH carrier. At $t=0$ the joint pairs coalesce
to the {{JITTERBUG_OCTAHEDRAL_JOINTS}} vertices of a regular octahedron
and the bars occupy {{JITTERBUG_OCTAHEDRAL_EDGES}} distinct edge
positions.
:::

The `jitterbug-rigid-motion` probe checks each landmark in
[@prop:jitterbug-rigid-motion] in exact arithmetic. It derives the
convex hull from the vertices and compares the two nondegenerate
landmark graphs with the IVM shell and OPH carrier. The limit
describes joint and bar coincidence. The separate continuation
described by Fuller and analysed as another motion by Stachel lies
beyond this branch [@fuller1979synergetics2, §938.14; @stachel1994heureka].

## The convex hull acquires edges before it becomes regular {#sec:jitterbug-graph}

For $0<t<1$, the four joints that bounded the square towards the
positive $x$ axis are $(a,\pm b,0)$ and $(b,0,\pm a)$. They are now
skew. The planes $x\pm(1-t)z=a$ support the two triangles that share
the short diagonal joining $(a,b,0)$ and $(a,-b,0)$. Sign changes and
cyclic permutations give the occupied triangles and the new hull
triangles. Together they exhaust the triangular facets of the hull
by the spherical Euler count. There are six such diagonals and each
joint lies on exactly one.

::: {.proposition #prop:jitterbug-graph-obstruction}
At every interior point of the branch the convex hull includes the
fixed bar graph and a perfect matching of short diagonals. The short
diagonal reaches the bar length exactly at $t=(\sqrt5-1)/2$.
At that point all {{JITTERBUG_ICOSAHEDRAL_EDGES}} hull edges have the
same length, including {{JITTERBUG_MATCHING_EDGES}} absent from the
physical linkage.
:::

The diagonal-to-bar squared length ratio is $2t^2/(1+t^2-t)$.
It is greater than one before the icosahedral configuration and less
than one after it. The `jitterbug-graph-obstruction` probe computes this
ratio and the three incidences in [@prop:jitterbug-graph-obstruction].
The matching is already visible as hull edges before the ratio reaches
one. At the regular icosahedron, all hull edges are nearest pairs.
Fuller explicitly identifies the additional connections required to
brace this configuration [@fuller1979synergetics2, §938.13].
The fixed-bar construction also differs from a central-sphere
deformation. Keeping the outer spheres touching a fixed central sphere
requires a configuration-dependent rescaling
[@kusner2018spheres, §5.5].

## Which seam precision follows the motion? {#sec:jitterbug-inference}

Keep unit conductance on the bars. Assign conductance $w\ge0$ to
each edge in the designated matching. This is a separate choice of
seam model. In the intrinsic scalar Gaussian, the precision on the
centred vertex space is the resulting Laplacian $L(w)$. For weighted
fair binary repair, sample seams in proportion to those same weights.
The Gaussian and repair descriptions then share the Laplacian as in
[@prop:lag-closure-transfer].

::: {.proposition #prop:jitterbug-weighted-precision}
The matching-pair contrast variances for the bar-only network and the
unit-conductance icosahedral network are
{{JITTERBUG_MATCHING_VARIANCES}}, respectively. The centred precision
determinant increases by a factor of
{{JITTERBUG_PRECISION_RATIO}}. With bar-only precision it is constant
throughout the motion.
:::

Signed coordinate symmetries split the Laplacian into even and odd
parts. The trace-free even blocks give eigenvalues {{JITTERBUG_EVEN_DIAGONAL_EIGENVALUE}} with
multiplicity two and $4+2w$ with multiplicity three. Each odd coordinate
character acts by the same symmetric matrix, in the order of the
small and large coordinate supports,
$$
K(w)=\begin{pmatrix}4+2w&-2\\-2&4\end{pmatrix}.
$$
This gives $\det{}' L(w)=36(4+2w)^3(12+8w)^3$. Differentiating the log
determinant adds the variances of the matching contrasts. They are
equal by symmetry, giving each variance
$(4+2w)^{-1}+4/(12+8w)$. The entropy decrease of the centred
Gaussian in [@prop:jitterbug-weighted-precision] is half the log of
the determinant ratio. The `jitterbug-weighted-precision` probe also
computes the determinant and contrasts directly by grounding the
weighted graph.

::: {.proposition #prop:jitterbug-coordinate-closure}
For the branch from the cuboctahedron to the regular icosahedron,
the geometric coordinate readout has exactly Gaussian Markov lag
covariances if and only if $w=1/t-t$. This choice goes from the
bar-only to the unit-weight carrier precision at the two landmarks.
With bar length two, the factors relating squared chord length to
lowest-mode contrast variance at the two landmarks are
{{JITTERBUG_CHORD_FACTORS}}, respectively.
:::

In each odd block the geometric coordinate has coefficient vector
$(t,1)$. It is an eigenvector of $K(w)$ exactly at the weight stated
in [@prop:jitterbug-coordinate-closure], with eigenvalue $4-2t$.
At that weight the geometric coordinates span the lowest positive
eigenspace. Their equal squared column norms times the eigenvalue
give the chord-to-variance factor. At fixed $t$ every other weight
leaves a nonzero lag-two covariance defect, by
[@prop:lag-closure-transfer]. The `jitterbug-coordinate-closure` probe
checks the invariance condition and the defect using the weighted
vertex operator and rational interior witnesses. This is a
covariance match for the stated weighted protocol. Higher-order
repair statistics still require their own comparisons.

![Fuller's Jitterbug at a common bar length. The solid bars keep their incidences, and the added hull diagonals are drawn separately. The icosahedral frame has {{JITTERBUG_RIGID_BARS}} bars and {{JITTERBUG_MATCHING_EDGES}} additional geometric seams. The lower plots distinguish the geometric diagonal length from the chosen conductance that preserves Gaussian lag covariance.](../output/figures/fig_jitterbug.png){#fig:jitterbug width=100%}
