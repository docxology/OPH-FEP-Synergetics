# Three-way correspondences and counterexamples {#sec:three-way-limits}

The same twelve-port count allows a controlled comparison. The icosahedral
carrier and the cuboctahedral first shell of the IVM have different seams,
faces and symmetry groups [@oph; @fuller1975synergetics]. Here we hold an
inference rule or a geometric operation fixed and ask which conclusions
survive the change in graph. The FCC experiments apply OPH's local rule
to a comparison graph under the specified protocol.

## Which pair records preserve repair dynamics? {#sec:repair-contrast-record}

Select a seam uniformly, then use the fair binary repair rule of
[@sec:repair-stopping-distance]. If the seam has one loaded endpoint, one
coin outcome moves the unit across it. If it has zero or two loaded
endpoints, the binary state stays fixed. Thus in a state consisting of
two occupied ports every unoccupied neighbour of either unit contributes
one directed hop, with probability $1/(2|E|)$. This also proves the symmetry of
the transition matrix and stationarity of the uniform law over pairs.

::: {.proposition #prop:repair-contrast-record}
On both graphs, the unit intrinsic Gaussian contrast variance of the two
occupied ports is a strongly lumpable record of fair binary repair. On
the FCC first shell, graph distance alone fails that criterion. The two
distance-two classes have variances {{THREE_FCC_MIDDLE_VARIANCES}} and
one-attempt probabilities {{THREE_FCC_INWARD}} of reaching adjacency,
respectively.
:::

The carrier is distance-transitive: up to graph symmetry, a pair is
specified by its distance. The FCC shell admits permutations and sign
changes of the coordinate axes. On its twelve directions with two
nonzero components, they give one pair orbit for each of the four
squared chords in [@prop:fcc-spectral-geometry]. To distinguish the two
middle orbits, pairs with the shorter chord lie in the same coordinate
plane and have orthogonal directions. Pairs with the longer chord share
one coordinate axis with opposite signs. The full variance in each
case also distinguishes the orbits. Because seam sampling and the coin
are invariant under graph symmetries, the total transition probability
into each orbit is constant on the source orbit. This is the row-sum
criterion for strong lumpability [@kemeny1976finite, theorem 6.3.2].

For a nonadjacent pair, moving into adjacency requires a hop to a common
neighbour. Each common neighbour contributes two possible hops. The two
FCC middle classes have {{THREE_FCC_COMMON_NEIGHBOURS}} common neighbours,
respectively. Merging them into one distance-two record leaves different
one-step exit probabilities for the same observed value. The
`repair-contrast-quotient` and `repair-distance-counterexample` probes
check all binary states on both graphs. They compare native attempt counts
with the graph-hop rule and compute contrast variance by exact network solves.

The FCC quotient has {{THREE_FCC_ATTEMPTS}} equiprobable draws per row, in
increasing squared-chord order. Its ordered rows of draw counts are
{{THREE_FCC_QUOTIENT}}, with stationary masses {{THREE_FCC_STATIONARY}}.
Even under stationarity, the next-step chance of adjacency from distance
two changes from {{THREE_FCC_STATIONARY_INWARD}} to
{{THREE_FCC_HISTORY_INWARD}} when the previous pair was adjacent. To see
this, let $r$ be the inward probability of the hidden pair type. Reversibility
weights each type by $r$ after observing previous adjacency. The new
average is $\mathbb E[r^2]/\mathbb E[r]$.

## Matching Gaussian covariance can leave observable memory {#sec:repair-centroid-memory}

Write $X_t$ for the binary load vector, $R$ for the matrix of vertex positions,
and $Y_t=R^\top X_t$ for twice the centroid of the loaded pair. If $\lambda$
is the Laplacian eigenvalue of the coordinate columns and $m=|E|$, then
$\mathbb E[Y_{t+1}\mid X_t]=\alpha Y_t$ with $\alpha=1-\lambda/(2m)$.
This follows by adding one-hop coordinate changes from the native kernel.
The coordinate columns have zero sum and equal squared norm on both graphs.
The uniform pair law gives
$\operatorname{Cov}(X_t)=c(I-J/|V|)$, with $c$ determined by the one- and
two-port marginals. This makes the centroid covariance $\sigma^2 I$.
Iterating the conditional mean gives
$\operatorname{Cov}(Y_t,Y_{t+k})=\alpha^k\sigma^2 I$ for every nonnegative
integer $k$, the stationary lag covariance of a sampled isotropic
Ornstein–Uhlenbeck process with decay rate $-\log\alpha$ [@uhlenbeck1930brownian].

::: {.proposition #prop:repair-centroid-memory}
Despite this covariance match, the stationary centroid on each graph
violates the Markov property. There is a nonzero centroid $y$ for which
$\Pr(Y_{t+1}=y\mid Y_t=0)$ and
$\Pr(Y_{t+1}=y\mid Y_t=0,Y_{t-1}=y)$ differ. On the carrier the two probabilities
are {{THREE_CARRIER_CENTROID_PROBS}}. On the FCC shell they are
{{THREE_FCC_CENTROID_PROBS}}.
:::

At zero centroid the pair must be antipodal. Those hidden states are
uniformly distributed under stationarity. The probability of jumping
to the target centroid from an antipodal pair is either zero or the
probability of a single directed hop. Conditional on the previous centroid
being the target, reversibility restricts the present antipodal pair to
one that can jump back to it. The `repair-centroid-memory` probe finds an
explicit target on each graph and checks the three-time probabilities both
by summing the native transition counts and by directly counting hops.
It also checks the first conditional moment and stationary variance over
the entire binary shell.

![The repair centroid and its Gaussian match share a normalised lag covariance on each shell. Conditioning on the previous centroid changes the return probabilities from zero. The two probabilities are {{THREE_CARRIER_CENTROID_PROBS}} on the carrier and {{THREE_FCC_CENTROID_PROBS}} on the FCC shell.](../output/figures/fig_centroid_memory.png){#fig:centroid-memory width=100%}

## Reciprocal conductances preserve planar duality {#sec:reciprocal-dual-weights}

Give each matched primal and dual edge positive conductances $g_e$ and $h_e$.
For an intrinsic Gaussian with precision $B D_g B^\top$, measure the
edge difference on seam $e$ in units multiplied by $\sqrt{g_e}$.
Here $B$ is the vertex incidence matrix and $D_g$ is diagonal. The resulting
covariance is the orthogonal projector onto the row space of
$B D_g^{1/2}$. The matching expression on the dual uses the oriented face
matrix $C$ in place of $B$.

::: {.proposition #prop:reciprocal-dual-weights}
For the two spherical polyhedral pairs of [@prop:planar-dual-covariance],
the conductance-normalised primal and dual covariance matrices sum to the
identity if and only if $g_e h_e$ is constant across all matched edges.
Copying a nonuniform conductance vector onto the dual fails the identity.
:::

Connectivity and Euler's spherical formula make the projector ranks sum
to the number of edges. The two row spaces are orthogonal exactly when
$B\operatorname{diag}(\sqrt{g_e h_e})C^\top=0$. At each vertex of a simple
face boundary, the two incident terms in the corresponding column must
cancel. This forces the products to be equal along that boundary.
Neighbouring faces share an edge, and the face graph is connected, so
the same product holds everywhere. Conversely, a constant product reduces
the cross term to a multiple of $BC^\top=0$. The diagonal part is the
weighted resistance sum of [@furrer2020resistors, theorem 5].

For the numerical check, change one conductance from unit weight to
{{THREE_DUAL_CHANGED_WEIGHT}}. At that seam the matched diagonal sums are
one under the reciprocal policy and {{THREE_DUAL_COPIED_SUMS}} under the
copied policy, in carrier and FCC order. The probes
`reciprocal-dual-weights` and `copied-dual-weight-obstruction` solve both
weighted networks and compare the modified entry against the rank-one
update from unit precision. For an edge with original resistance $r$,
the normalised diagonal after changing its conductance to $w$ is
$wr/(1+(w-1)r)$.

## The linear frame does not determine quadratic uncertainty {#sec:quadratic-mode-transfer}

Write $n_i$ for the direction of vertex $i$ on the unit sphere and let
$Q$ be a symmetric trace-free matrix, with the Frobenius inner product
on such matrices. The quadratic field is $TQ=(n_i^\top Q n_i)_i$. This is
a way to express changes that preserve the sign parity of a direction.
The centred intrinsic Gaussian has graph Laplacian $L$ as precision
on the zero-mean vertex space.

::: {.proposition #prop:quadratic-mode-transfer}
On either shell, $T$ is a linear bijection onto the centred antipodally
even vertex fields. The Gaussian coefficient precision of that part of
the vertex field is $T^*LT$. It is isotropic on the icosahedral carrier.
On the FCC shell, the variance of a unit diagonal coefficient is
{{THREE_FCC_DIAGONAL_VARIANCE}} and that of a unit off-diagonal coefficient
is {{THREE_FCC_OFF_DIAGONAL_VARIANCE}}. A rotation of one direction in
coefficient space into the other changes the variance by a ratio of
{{THREE_FCC_QUADRATIC_VARIANCE_RATIO}}.
:::

Both vertex sets satisfy $\sum_i n_i n_i^\top=(|V|/3)I$,
which gives the isotropic linear frame and centres every trace-free quadratic.
The Gram matrix $T^*T$ instead uses fourth moments. The icosahedral set
integrates fourth moments as on the uniform sphere, the relevant part
of its spherical design property [@delsarte1977codes]. Direct summation
on the unit directions gives {{THREE_CARRIER_QUADRATIC_FRAME}} times the
identity. On the FCC shell the corresponding factors are
{{THREE_FCC_QUADRATIC_FRAME}} on diagonal and off-diagonal matrices.
Both frames have full rank on the trace-free space. There are six
antipodal pairs, so the centred even space has the same dimension as
the trace-free space. Both Laplacians preserve that space and the odd
space, which gives the coefficient precision in
[@prop:quadratic-mode-transfer] by changing variables in the Gaussian
density.

Each carrier vertex has one neighbour in each other antipodal pair.
The Laplacian thus acts on the centred even space by a single
eigenvalue. On the FCC shell, diagonal forms are constant within a
coordinate-plane class and neighbour twice into each of the other two.
A mixed product $n_xn_y$ vanishes outside its coordinate plane.
Its neighbour sum is zero everywhere by cancellation of signs.
Consequently the FCC diagonal and off-diagonal eigenvalues are
{{THREE_FCC_QUADRATIC_EIGENVALUES}}. The resulting pulled-back precisions
are {{THREE_CARRIER_QUADRATIC_PRECISION}} on the carrier and
{{THREE_FCC_QUADRATIC_PRECISIONS}} on the two FCC subspaces.
Rotating the normalised form $n_x^2-n_y^2$ in its coordinate plane by
$\pi/4$ gives the normalised form $2n_xn_y$, up to sign, exhibiting the
anisotropy with equal Frobenius norms. The `quadratic-even-field` and
`quadratic-isotropy-obstruction` probes compute the exact moment and
Laplacian matrices and recompute the precision from the quadratic
differences along individual seams.

## When does a plane screen the two halves? {#sec:central-plane-blanket}

For a plane through the centre of a shell, let $S$ be its on-plane ports
and $A_+,A_-$ its two strict halves. In a scalar Gaussian with seam
precision $L$ and a positive diagonal shift, conditioning on $S$ leaves
a precision block on $A_+\cup A_-$. Its cross term is zero exactly when
no seam joins the strict halves, which is equivalent to their conditional
independence [@dempster1972covariance]. For the intrinsic calculation we
ground the on-plane ports and use the Dirichlet Gaussian on the free
coordinates, with the grounded Laplacian as precision.

::: {.proposition #prop:central-plane-blanket}
None of the carrier's {{THREE_CARRIER_PLANES}} central planes spanned by
vertex directions screens its strict halves. Each leaves
{{THREE_CARRIER_CROSS_SEAMS}} crossing seams and conditional mutual
information $\tfrac12\log({{THREE_CARRIER_SCREEN_RATIO}})$ nats at unit
intrinsic seam precision. All {{THREE_FCC_PLANES}} vertex-spanned central
planes on the FCC shell screen their strict halves.
:::

The vertex-spanned planes contain two independent vertex directions.
The census takes their exact cross products, identifies proportional
normals, and classifies each remaining port by the sign of its inner
product with that normal. On the carrier, each plane contains two
antipodal pairs. On the FCC shell, the coordinate planes contain four
ports and the planes normal to body diagonals contain six. Any other
central plane contains at most one antipodal pair. The graph remains
connected on either shell after removing any two ports, so those
other planes cannot screen. To see the remaining information on the
carrier, let $K$ be the precision after clamping $S$. Block determinant
identities for the Gaussian give
$$
I(A_+;A_-\mid S)=\tfrac12\log
\frac{\det K_{++}\det K_{--}}{\det K}.
$$
The `central-plane-screening` probe checks exact signs, positive
graph separators and the clamped determinants. Thus a geometric
equator has different screening properties on the two graphs, even
though the same seam-precision rule applies.

![Four three-way tests on the carrier and FCC shell. A shows the one-step inward probabilities at distance two, including the distinct FCC rates {{THREE_FCC_INWARD}}. B compares weighted variance sums under reciprocal and copied conductances. Copied sums are {{THREE_DUAL_COPIED_SUMS}}. C shows the FCC quadratic variance split between {{THREE_FCC_DIAGONAL_VARIANCE}} and {{THREE_FCC_OFF_DIAGONAL_VARIANCE}}. D contrasts the {{THREE_CARRIER_PLANES}} crossed carrier planes and {{THREE_FCC_PLANES}} screening FCC planes.](../output/figures/fig_three_way.png){#fig:three-way width=100%}
