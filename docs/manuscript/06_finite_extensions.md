# Exact extensions of the finite models {#sec:finite-extensions}

The registry also tests what the basic correspondences imply for first
passage, coarse-graining, record placement and dual geometry. These
extensions use the same carrier and FCC shell as the preceding sections.
[@tbl:ontology-exact-extensions] identifies the objects and assumptions for
these reductions.

## First binary entry and distance reduction {#sec:repair-stopping-distance}

::: {#prop:repair-first-entry .proposition}
Choose a seam uniformly on each attempt and use a fair coin. Start with both
units at a port $i$ and let $T$ be the first attempt that reaches the binary
shell. Then $T$ is geometric on the positive integers with hazard
{{ENTRY_HAZARD}}, mean {{ENTRY_MEAN}}, and variance {{ENTRY_VARIANCE}} in
squared attempts. The two occupied ports at entry are $i$ and a uniformly
chosen neighbour, independently of $T$.
:::

Only a seam incident to $i$ splits the double. Either coin gives the same
split. All other attempts leave the double at $i$, so the transient block
is a scalar multiple of the identity and each attempt succeeds with
probability $\deg(i)/\lvert E\rvert$. At a successful attempt all its
incident seams have equal probability, for every value of $T$. The
first-entry construction applies by stopping the chain at the first binary
configuration [@kemeny1976finite, §3.2]. The `repair-first-entry`
probe compares the complete transient and entry blocks from native attempts
with the incidence formula in [@prop:repair-first-entry].

A uniform doubled start gives each seam two equally weighted ways to enter
its adjacent pair. The entry distribution $\nu$ is therefore the binary
equilibrium $u$ conditioned on adjacency. Both distributions give each port
occupancy probability {{ENTRY_PORT_MARGINAL}}. Their total variation
distance is {{ENTRY_TV}} and $D(\nu\Vert u)=\ln({{ENTRY_SUPPORT_RATIO}})$.
At every finite attempt from a doubled start, the full-shell distribution
still assigns positive probability to states of zero equilibrium mass.
Its relative entropy to the binary equilibrium on that full space is
infinite.

The first adjacent pair also gives a simple Bayesian observation of the
original doubled port. Conditional on the pair, its endpoints are the two
equally likely initial ports. The time law is common to all initial ports
and independent of the entered pair, so it adds no information about the
original port. The mutual information for a uniform doubled start is
$\ln(\lvert V\rvert/2)$ nats.

::: {#prop:repair-distance-lumping .proposition}
On the same protocol, the graph distance between the two loaded units is
itself a Markov chain for every initial law. Order the distance classes
from zero, which denotes the doubled configurations. Their sizes are
{{DISTANCE_CLASS_SIZES}}. The quotient attempt matrix is
$$
Q=\frac{1}{60}
\begin{pmatrix}
50&10&0&0\\
0&56&4&0\\
0&4&54&2\\
0&0&10&50
\end{pmatrix}.
$$
Its stationary masses are {{DISTANCE_STATIONARY_MASSES}}. On positive
distances the nonconstant eigenvalues are $(25\pm\sqrt{3})/30$.
:::

At positive distance the intersection array counts the inward and outward
seams for each unit. Each such move uses one of the two coins, and the seam
joining two occupied neighbours leaves their loads unchanged. The doubled
class uses the entry law of [@prop:repair-first-entry]. The
`repair-distance-lumping` probe compares these counts with all native
transition rows and checks the eigenvalue coefficients. This is the
row-sum criterion for strong lumpability [@kemeny1976finite, theorem 6.3.2].
The probe checks the incoming coefficients after weighting by class size.
If $K$ is the full kernel and the rows of $U$ give uniform laws in the
corresponding distance classes, this verifies $UK=QU$. The dynamics
preserves uniformity conditional on the current class. For such laws on
positive distances, relative entropy to binary equilibrium equals the
relative entropy of their distance masses to the stationary masses. The
distance gap $(5-\sqrt{3})/30$ in [@prop:repair-distance-lumping] exceeds the
full binary gap of [@sec:pr01].

## Grounded precision counts rooted forests {#sec:grounded-forest-precision}

::: {#prop:grounded-forest-precision .proposition}
On the unit-weight scalar carrier, let $S$ be a nonempty record set and $I$
its complement. Combine the Laplacian quadratic with a proper Gaussian prior on the
record coordinates. Let $F(S)$ count spanning seam forests having exactly one
record in each component. The conditional Gaussian precision is $L_{II}$
and its determinant is $F(S)$. Thus the conditional covariance has determinant
$1/F(S)$.
:::

With the record values clamped, the quadratic in the free ports has Hessian
$L_{II}$. Deleting the record rows of the vertex–seam incidence matrix gives
$B_I$ and $L_{II}=B_IB_I^{\top}$. Cauchy–Binet sums the squares of the
incidence minors. A cycle makes a minor singular, as does a component devoid
of roots. A forest with one root per component gives determinant $\pm1$ by
successive leaf elimination. This is the principal-minor case of the
all-minors matrix-tree theorem [@chaiken1982minors].

The `grounded-forest-precision` probe checks [@prop:grounded-forest-precision]
for all {{FOREST_PAIRS_CHECKED}} pairs and {{FOREST_TRIPLES_CHECKED}} triples
of records. An independent formula computes the same minors from the
contrast variances in [@prop:resistance]. Ground one selected record $a$.
For the remaining records, the contrast covariance has entries
$G_{ij}=(R_{ai}+R_{aj}-R_{ij})/2$. By the complementary-minor identity,
$F(S)=\tau\det G$, where $\tau$ is the spanning-tree count. A separate test
counts the forests by assigning each free port to a rooted component and
enumerating the possible trees within each component.

## Greedy clamping can miss the best triple {#sec:greedy-record-placement}

::: {#prop:greedy-record-placement .proposition}
A rule that starts with any single record and maximises the grounded
determinant at every addition chooses an antipodal pair. Every triple
extending that pair leaves higher conditional entropy than a triple whose
three mutual graph distances are all two.
:::

At a fixed record-set size the conditional entropy is
$\bigl(\lvert I\rvert\log(2\pi e)-\log F(S)\bigr)/2$. For a possible
new record $v$, the ratio $F(S\cup\{v\})/F(S)$ is the current conditional
variance at $v$, so the greedy rule clamps a most uncertain port. Relative to
an antipodal pair, any third port is adjacent to one endpoint and at
distance two from the other. Every resulting triple has determinant
{{FOREST_GREEDY_TRIPLE}}. Exhaustive comparison instead gives a maximum of
{{FOREST_BEST_TRIPLE}} at triples with every mutual distance two and a
minimum of {{FOREST_WORST_TRIPLE}} at facial triangles. By
[@prop:greedy-record-placement], the greedy entropy excess is
$\tfrac12\log({{FOREST_ENTROPY_RATIO}})$ nats. The `greedy-record-placement`
probe checks this obstruction for the carrier with unit weights and exact
clamping.

## Planar dual contrasts fill the edge space {#sec:planar-dual-covariance}

The oriented carrier seams also specify a graph of neighbouring faces. Two
faces are adjacent in this dual graph when they share a seam. The
cuboctahedral shell and its rhombic dodecahedral Voronoi cell supply a second
case with polygonal faces.

::: {#prop:planar-dual-covariance .proposition}
For each of these spherical graph pairs, assign unit precision to every
edge of each graph. The edge-difference covariance matrices of the two
intrinsic Gaussian vertex fields add to the identity in matched, signed
edge coordinates.
:::

Let $B$ be the vertex incidence and let the rows of $C$ be the oriented
face boundaries in the same seam coordinates. Every seam appears with
opposite signs in its two faces, so $C$ is also the incidence matrix of the
dual graph. Boundary of boundary gives $BC^{\top}=0$. Connectivity and
Euler's spherical identity make the two row-space dimensions sum to the
edge count. The two Gaussian contrast covariances are the orthogonal
projections onto those row spaces, which gives
[@prop:planar-dual-covariance].

The `planar-dual-covariance` probe compares {{DUAL_MATRIX_ENTRIES}} entries
from separately grounded exact inverses and verifies equal spanning-tree
cofactors within each pair. The primal and dual edge resistances are
{{DUAL_CARRIER_RESISTANCES}} for the carrier pair and
{{DUAL_IVM_RESISTANCES}} for the IVM pair. These diagonals give the unit
case of the resistance sum formula [@furrer2020resistors, theorem 5].

## The FCC shell splits its middle distance class {#sec:fcc-spectral-geometry}

::: {#prop:fcc-spectral-geometry .proposition}
Use the cuboctahedral FCC neighbour shell, scaled to edge length two, and
the centred intrinsic Gaussian with its unit graph Laplacian as precision.
The squared chord of every pair is {{FCC_CHORD_FACTOR}} times the
contribution to contrast variance from the first positive eigenspace.
Pairs two graph steps apart give two squared chords,
{{FCC_DISTANCE_TWO_CHORDS}}, and two full variances,
{{FCC_DISTANCE_TWO_RESISTANCES}}, in corresponding order.
:::

The shell graph is the cube's line graph. For the unsigned vertex–edge
incidence matrix $B$ of the cube, its shell Laplacian is $6I-B^{\top}B$.
The cube has $BB^{\top}=3I+A_{\mathrm{cube}}$. On cube vertices labelled
by three signs, products of the signs are adjacency eigenvectors with
eigenvalues $3-2k$, where $k$ is the number of signs in the product.
The nonzero eigenvalues of $B^{\top}B$ are those of $BB^{\top}$, with
additional zero eigenvalues supplied by the kernel of $B$. This identifies
the first positive shell eigenvalue and its multiplicity. In the
edge-two coordinates $X$, the identities $LX=2X$ and $X^{\top}X=16I$
then give the scale of [@prop:fcc-spectral-geometry]. The
`cuboctahedral-spectral-chords` probe checks this on the integer quadray
embedding by rescaling its coordinate Gram.

For the full variance the probe grounds a vertex and solves exactly,
then checks the contrast of every pair against polynomial spectral
projectors. The Kirchhoff index is {{FCC_KIRCHHOFF_INDEX}}. Contrast
variance here is effective resistance as in [@prop:resistance; @klein1993resistance].
The `cuboctahedral-distance-split` probe identifies both pair classes
at the same graph distance. The icosahedral carrier of
[@prop:carrier-spectral-chords] gives one variance per distance class,
so its distance-only formula needs that graph's structure.

![Exact extensions of the finite models. A follows the distribution of graph
distance from a doubled start through the repair quotient. B compares
grounded determinants over triple record layouts and highlights the greedy
and optimal choices. C shows the diagonal of the matched edge covariance
decomposition on the carrier and FCC dual pairs. D compares full Gaussian
contrast variance and lowest-mode variance over the FCC shell's four pair
classes. Generated from the exact model censuses.](../../output/figures/fig_extensions.png){#fig:extensions width=100%}

[@fig:extensions] connects the exact transition coefficients, Gaussian record
design, planar dual covariance and shell geometry in the same finite models.
