# Transfer criteria beyond the two shells {#sec:transfer-criteria}

The preceding comparisons use one OPH carrier and the first IVM shell.
Here we apply the same local seam rule, or the same Gaussian construction,
to other finite graphs and cellulations. The hypotheses below identify
which conclusions transfer. The larger close-packed examples use the
$D_3$ lattice of the IVM [@fuller1975synergetics]. The higher-dimensional
root shells and periodic square cellulation are mathematical comparison models.

## Distance regularity, resistance and repair {#sec:transfer-distance}

On a finite, connected, simple graph with uniformly sampled seams and
a fair coin, binary repair moves one of the two units to a vacant
neighbour with probability $1/(2|E|)$ per directed hop. It waits on the
other draws. Thus a variance record is strongly lumpable exactly when,
for every source resistance class, all pairs have the same number of hops
into each other resistance class. This is the row-sum criterion
of [@kemeny1976finite, theorem 6.3.2]. In particular, a resistance class
that is a graph automorphism orbit meets this requirement. If a class
merges orbits, the aggregated rows must agree.

::: {.proposition #prop:distance-regular-repair}
For a finite connected regular graph, the fair binary repair kernel is
strongly lumpable by graph distance exactly when the graph is
distance-regular. On that domain, the unit intrinsic Gaussian contrast
variance is strictly increasing with graph distance and gives the same
quotient.
:::

Write $D_j$ for the distance matrices and $A=D_1$ for adjacency.
For distinct vertices $x,y$, the number of hops from $\{x,y\}$ to a
pair at distance $j>0$ is $(AD_j+D_jA)_{xy}$. The number of draws that
wait depends only on source distance: the total number of directed
hops from a pair is twice the degree, less two if that pair is adjacent.
If distance is strongly lumpable, each $AD_j+D_jA$ belongs to the span
of the distance matrices, including its diagonal, which is constant
by regularity. Its only possible terms have indices $j-1,j,j+1$.
The coefficient at $j+1$ is positive whenever that class exists,
since a geodesic gives a hop inward. Starting with $D_0=I$ and $D_1=A$,
induction makes each $D_j$ a polynomial in $A$. The two products now
commute, and the coefficients of $AD_j$ are the intersection numbers.
Conversely, those numbers give the requisite hop totals. Biggs's
potential representation makes the resistance increase strictly
from one distance class to the next [@biggs1993potential]. The
`distance-regular-repair` probe verifies [@prop:distance-regular-repair]
on the four-dimensional cube. The grounded Gaussian solve and the
intersection-number potential both give {{TRANSFER_CUBE_RESISTANCES}}
for the successive contrast variances.

::: {.proposition #prop:resistance-collision}
For every $m\ge 2$, equal Gaussian contrast variance on $K_{m,m+1}$
can conceal different one-step repair probabilities. The same failure
occurs on a regular periodic FCC graph even when the record also gives
graph distance and common-neighbour count.
:::

In the bipartite graph, the two parts have $m$ and $m+1$ vertices.
The within-part resistances are $2/(m+1)$ and $2/m$. The cross-part
resistance is also $2/(m+1)$, from $(m+(m+1)-1)/(m(m+1))$.
These formulas follow by decomposing the Laplacian into the mean-zero
subspaces on each part and the one-dimensional contrast between their
means. A pair in the smaller part cannot reach the larger-part class in
one hop. A cross-part pair can do so with probability $1/(2(m+1))$.
The `bipartite-contrast-collision` probe independently solves the
grounded systems and counts the native draws on several members of
this family.

For the periodic FCC example, use coordinates modulo four, keep the
even-sum vertices and join nearest-neighbour FCC steps. Its
{{TRANSFER_PERIODIC_FCC_VERTICES}} vertices have common degree
{{TRANSFER_PERIODIC_FCC_DEGREE}}. The two displacements
{{TRANSFER_PERIODIC_FCC_DISPLACEMENTS}} both have graph distance two and
four common neighbours.
Their common contrast variance is {{TRANSFER_PERIODIC_FCC_COLLISION}}.
Their one-attempt probabilities of reaching variance
{{TRANSFER_PERIODIC_FCC_TARGET}} are {{TRANSFER_PERIODIC_FCC_JUMPS}}.
The `periodic-fcc-contrast-collision` probe verifies
[@prop:resistance-collision] by comparing the full grounded resistance
table against a finite Fourier sum on the periodic group, then
enumerating the two-token hops.

The loss of the Markov property can persist under stationarity. In a
reversible chain let $r$ be the next-step probability of a target class
for a hidden state in the current observed class. The probability given
only that class is $\mathbb E[r]$. Adding a previous observation of the
target reweights the hidden current state by $r$ and gives
$\mathbb E[r^2]/\mathbb E[r]$. The latter is strictly larger when the
off-class rates vary. In the periodic FCC witness, the two stationary
probabilities are {{TRANSFER_PERIODIC_FCC_HISTORY}}.

## When a projected covariance closes {#sec:lag-closure-transfer}

Consider a fixed number $k$ of binary loads, with $0<k<n$ on a connected
$n$-vertex graph. Sample seam $e$ with probability $w_e/W$ for positive
weights and $W=\sum_e w_e$, then use a fair coin. Each update is an equal
mixture of identity and the transposition of the two endpoints. This is
a discrete version of the symmetric exclusion rule [@caputo2010aldous].
Uniform binary configurations are stationary. Their centred load
covariance is $\tfrac{k(n-k)}{n(n-1)}(I-J/n)$. The mean vertex operator
is $P=I-L_w/(2W)$, with $L_w$ the weighted Laplacian.

::: {.proposition #prop:lag-closure-transfer}
Let $U$ have orthonormal columns in the centred vertex space. The
stationary lag covariance of the $U$-coordinates, divided by the common
lag-zero variance, is $M_j=U^\top P^jU$. It agrees at every nonnegative integer lag
with a stationary Gaussian Markov process on those coordinates exactly
when the column space of $U$ is invariant under $L_w$.
:::

The first assertion follows by iterating the conditional mean. A
Gaussian Markov process with unit stationary covariance and one-step
cross covariance $M_1$ must have $M_j=M_1^j$. In particular,
$$
M_2-M_1^2
= U^\top P(I-UU^\top)PU
= ((I-UU^\top)PU)^\top((I-UU^\top)PU).
$$
The matrix vanishes exactly when the observed space is invariant.
Then all lags agree, and a Gaussian recursion with transition $M_1$
and innovation covariance $I-M_1^2$ supplies a match. This criterion
concerns second moments. The three-time centroid witness in
[@prop:repair-centroid-memory] still distinguishes repair from its
Gaussian match even in the invariant case.

The `linear-lag-closure` probe checks [@prop:lag-closure-transfer] using
the geometric coordinates of the $D_4$ root shell, the regular
twenty-four-cell [@cohn2007d4]. Its unit seam weights give coordinate
retention {{TRANSFER_24CELL_UNIT_LAG}}. The `single-seam-lag-obstruction`
probe changes only one weight to {{TRANSFER_24CELL_WEIGHT}}. Along
that seam's chord direction, the observed lag-two covariance is
{{TRANSFER_24CELL_WEIGHTED_LAG_TWO}}. It exceeds the square of the
lag-one covariance by {{TRANSFER_24CELL_LAG_DEFECT}}. The defect has rank
one, as follows from the rank-one Laplacian update and its nonzero
components inside and outside the coordinate space. A native enumeration
of two occupied ports and a separate weighted vertex calculation agree
on the entire lag-two defect matrix.

![Repair transfer to the cube, bipartite and periodic FCC graphs, and the twenty-four-cell. The periodic FCC witness has two one-step rates {{TRANSFER_PERIODIC_FCC_JUMPS}} from the same variance class. The weighted twenty-four-cell has a nonzero two-lag defect {{TRANSFER_24CELL_LAG_DEFECT}}.](../output/figures/fig_repair_transfer.png){#fig:repair-transfer width=100%}

## Closed surfaces leave harmonic covariance {#sec:surface-harmonic-transfer}

Use a connected closed orientable cellulation of genus $g$ with simple
face boundaries, each edge in two distinct incident faces, and a connected
dual. Let $B$ be the vertex incidence and $C$ the consistently signed face
boundary matrix. Each positive primal weight $u_e$ and dual weight $v_e$
defines a Gaussian as in [@prop:reciprocal-dual-weights]. In common
edge coordinates let $P$ and $Q$ denote their conductance-normalised
contrast covariance projectors.

::: {.proposition #prop:surface-harmonic-transfer}
The remainder $H=I-P-Q$ is positive semidefinite exactly when $u_e v_e$
is constant over the edges. In that case $H$ is the harmonic projector
of rank $2g$. For general weights, let $k=\operatorname{rank}(PQ)$.
The trace of $H$ is $2g$ and its rank is $2g+2k$.
:::

The row spaces are orthogonal exactly when
$B\operatorname{diag}(\sqrt{u_ev_e})C^\top=0$.
Cancellation at each corner of a simple face propagates the weight
product across a face and then across the connected dual, as in
[@prop:reciprocal-dual-weights]. Conversely constant product gives
orthogonality by $BC^\top=0$. Euler's formula and connectedness give
$\dim\operatorname{im}P=|V|-1$, $\dim\operatorname{im}Q=|F|-1$, and
$|E|-\dim\operatorname{im}P-\dim\operatorname{im}Q=2g$.
For $x$ in the primal row space, $x^\top Hx=-\|Qx\|^2$. Thus positivity
forces orthogonality. The orthogonal complement of the two row spaces
is the space of harmonic edge forms [@degoes2016exterior].

In the general case, $Hx=0$ implies $x=Px+Qx$, and applying the two
projectors again gives $PQx=QPx=0$. The kernel is the direct sum of
$\operatorname{im}P\cap\ker Q$ and $\operatorname{im}Q\cap\ker P$.
Each loses $k$ dimensions relative to its row space. This gives the
rank formula in [@prop:surface-harmonic-transfer] and explains why the
trace alone cannot diagnose a covariance projector.

The `surface-harmonic-covariance` probe checks a doubly periodic square
grid with {{TRANSFER_TORUS_CELLS}} vertices, edges and faces, respectively.
For unit weights the remainder is the orthogonal projection onto constant
horizontal and vertical circulations. Its diagonal is
{{TRANSFER_TORUS_UNIT_RESIDUAL}}. The probe also solves independently for
circulations with prescribed periods when reciprocal weights vary.
The `copied-surface-weight-obstruction` probe compares the changed edge
under reciprocal and copied weights. The corresponding remainder
diagonals are {{TRANSFER_TORUS_RECIPROCAL_RESIDUAL}} and
{{TRANSFER_TORUS_COPIED_RESIDUAL}}. Copying leaves the trace unchanged,
doubles the remainder rank, and gives a negative direction.

## Quadratic variance across root shells {#sec:root-quadratic-transfer}

For each $d\ge 3$, the $D_d$ shell consists of $\pm e_i\pm e_j$ for
$i<j$. Edges have squared length two. The $D_3$ shell is the
cuboctahedral IVM comparison, and the $D_4$ shell is the
twenty-four-cell [@cohn2007d4]. Evaluate trace-free real quadratic
forms on the unit vertex directions with Frobenius coefficient norm.
Equip the vertices with the unit-seam intrinsic Gaussian precision.

::: {.proposition #prop:root-quadratic-transfer}
The quadratic observation space is invariant under the shell Laplacian.
On diagonal and off-diagonal trace-free forms, the Gaussian
coefficient precisions are $2d(d-2)$ and $8(d-2)$, respectively.
Their covariance is isotropic exactly at $d=4$.
:::

The squared coordinate sums are equal in all directions, so each
trace-free form evaluates to a centred field. For its fourth moments
the sums of $x_i^4$ and $x_i^2x_j^2$ on the unit directions are $d-1$
and one. Sign changes cancel the odd products. Thus the Gram
factors on Frobenius unit forms are $d-2$ and two.
A neighbour changes one element of the support and keeps the sign
on the shared coordinate. Summing the diagonal form over neighbouring
supports and using the trace-free condition gives Laplacian
eigenvalue $2d$. A mixed coordinate product has neighbour sum zero
by sign cancellation, so its eigenvalue is the vertex degree
$4(d-2)$. This proves invariance and [@prop:root-quadratic-transfer]
by pulling the Laplacian back through the evaluation map. Projection
and restriction have the same precision on that invariant space.
The `root-shell-quadratic-precision` probe enumerates the unit-direction
frames, adds the squared seam differences, and checks each Laplacian
eigenfield for the dimensions {{TRANSFER_ROOT_DIMENSIONS}}. At the
twenty-four-cell both variances are {{TRANSFER_24CELL_QUADRATIC_VARIANCE}}.

## What a plane screens depends on the constraint {#sec:sum-constraint-transfer}

In integer coordinates the frequency-$n$ FCC boundary has even
coordinate sum and $\max(\|x\|_\infty,\|x\|_1/2)=n$.
The central planes $x_1=0$ and $x_1+x_2+x_3=0$, together with their
images under signed coordinate permutations, separate the strict sides
of any connected induced FCC graph with vertices on both sides. A neighbour
changes one coordinate by at most one and the signed sum by at most
two. The latter has even parity, so a single move cannot join opposite
strict signs in either case. On the shell at frequency $n$, the
axis plane cuts out a square perimeter with $4n$ records, and the
diagonal plane cuts out a hexagon with $6n$ records. The square
and triangular shell facets each give a connected lattice grid, and
adjacent facets meet. Thus the `fcc-plane-screen-transfer` claim uses
the same scalar seam-screening criterion as [@prop:central-plane-blanket]
at every positive integer frequency. Its second-shell record counts
are {{TRANSFER_FCC_SECOND_PLANE_RECORDS}}.

::: {.proposition #prop:sum-constraint-transfer}
Let nonempty, clamped records separate two nonempty sides in the
grounded scalar Gaussian. A global sum constraint on the free field
reintroduces cross covariance of rank one, even though there is no
crossing seam. Conditioning on a noisy field sum gives finite
conditional mutual information; conditioning on the exact sum makes
it infinite.
:::

Let $K_A,K_B$ be the two grounded precision blocks,
$u_A=K_A^{-1}\mathbf 1$, $u_B=K_B^{-1}\mathbf 1$, and
$a=\mathbf 1^\top u_A$, $b=\mathbf 1^\top u_B$. The zero-sum intrinsic
field conditioned on the records is the grounded field subject to
the remaining sum constraint. Gaussian linear conditioning
[@bolin2021sparse] gives the cross covariance
$-u_Au_B^\top/(a+b+\tau^2)$ for sum-noise variance $\tau^2$.
For $\tau^2>0$ the mutual information is
$$
I(A;B\mid S,\text{sum record})=\frac12\log
\frac{(a+\tau^2)(b+\tau^2)}{\tau^2(a+b+\tau^2)}.
$$
For the exact sum, the individual half-field covariances remain
nonsingular, while their joint law lies on a hyperplane. Relative
entropy to the product of the half-field marginals is therefore infinite.
The two sums carry the cross dependence; residuals after Gaussian
regression on each side's sum are independent across the cut.

The noisy-sum dependence can grow with shell frequency. For either
strict side $A$, Cauchy–Schwarz gives
$a\ge |A|^2/(\mathbf 1^\top K_A\mathbf 1)$. The denominator counts seams
from the side to the recorded plane. The plane has $O(n)$ records of
bounded degree while each side has $\Theta(n^2)$ sites. Thus
$a=b=\Omega(n^3)$ by antipodal symmetry. At fixed positive sum-noise
variance, the determinant ratio in the information formula is at least
$a/(2\tau^2)$, so the conditional mutual information diverges with
frequency despite the absence of crossing seams.

The `centered-gaussian-screen-obstruction` probe checks
[@prop:sum-constraint-transfer] on the first two FCC boundaries. On
the second, unit sum noise gives determinant ratios
{{TRANSFER_FCC_SECOND_SUM_RATIOS}} for the two plane families. Grounding
records without a sum constraint leaves their half-fields independent.

![Transfer in the Gaussian geometry models. The torus remainder at the strengthened edge is {{TRANSFER_TORUS_RECIPROCAL_RESIDUAL}} under reciprocal weights and {{TRANSFER_TORUS_COPIED_RESIDUAL}} under copied weights. Root shell variances coincide at the twenty-four-cell. The FCC sum constraint reintroduces cross-plane dependence despite the grounded plane cuts.](../output/figures/fig_geometry_transfer.png){#fig:geometry-transfer width=100%}
