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
At each pose with distinct joints, $0<t\le1$, the geometric coordinate
readout has exactly Gaussian Markov lag covariances if and only if
$w=1/t-t$. The weight is zero at the cuboctahedron, is unit at the regular
icosahedron and rises without bound as the octahedral pairs coalesce.
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

## Past the regular icosahedron {#sec:jitterbug-after-icosahedron}

Write $\tau=(\sqrt{5}-1)/2$. Below $\tau$ the twelve joints are still
distinct and the hull retains the carrier incidence. The six short
diagonals, which used to compete with the rigid bars for shortest
distance, now form the entire nearest-pair graph.

::: {.proposition #prop:jitterbug-nearest-only-obstruction}
For $0<t<\tau$, the numbers of hull and nearest-pair edges are
{{JITTERBUG_CONTRACTED_EDGE_COUNTS}}, respectively. The nearest-pair
graph consists of disjoint dimers. Fair repair on those seams preserves
each pair total, and a Gaussian using those seams leaves
{{JITTERBUG_CENTRED_NULLITY}} null directions after global centering.
Collapsing the nearest pairs and rebuilding nearest pairs after collapsing
give edge counts {{JITTERBUG_METRIC_FIRST_VS_COLLAPSE_FIRST_EDGES}},
respectively.
:::

For joints with the same zero coordinate, the short sign flip has raw
squared chord $4t^2$. The other sign flips give the pair of squared
chords $\bigl(4,\,4(1+t^2)\bigr)$.
For different zero coordinates, the shared nonzero slot has signs that
either agree or disagree, giving $2(1+t^2-t)$ and $2(1+t^2+t)$.
These cases exhaust all pairs. The bars realise the smaller of the
cross-plane chords. This proves the switch to the matching at $t<\tau$.
The matching Laplacian has one zero mode per dimer. Coalescing the dimer
edges sends each of them to a loop, while recomputing distances on the
distinct sites finds the octahedral edges.

`jitterbug-nearest-only-obstruction` checks the chord classes against
pair enumeration and tests both orders of operation in
[@prop:jitterbug-nearest-only-obstruction] against the endpoint geometry.

Sphere contact imposes a different threshold. Place spheres at the
labelled joints with radii half the common centre-to-joint distance and
put a sphere of that same radius at the origin. The squared separation
needed for outer spheres to touch, before rescaling, is $1+t^2$.
Every chord other than the short flip stays above that threshold on the
open post-icosahedral branch. For the short flip the squared margin is
$3t^2-1$. Outer contacts first occur at
$t^2={{JITTERBUG_SPHERE_CONTACT_THRESHOLD_SQUARED}}$, and the same pairs
overlap below that point while the joints are distinct. This extends
the distance comparison used in the cuboctahedron-to-icosahedron
sphere protocol [@kusner2018spheres, §5.5].

::: {.proposition #prop:jitterbug-weight-policy-obstruction}
Past the regular pose, a unit-weight hull fails the coordinate lag
test. Inverse-square weighting closes at $t=1/2$ and fails at every
other strictly post-icosahedral pose. At $t=1/2$ the tuned weight is
{{JITTERBUG_HALF_POSE_WEIGHT}} and the normalised one-attempt
retention of a coordinate readout is {{JITTERBUG_HALF_POSE_RETENTION}}.
At the same pose the unit-weight hull defect is
{{JITTERBUG_UNIT_HULL_HALF_DEFECT}}. At $t=1/3$ the inverse-square
defect is {{JITTERBUG_INVERSE_SQUARE_THIRD_DEFECT}}.
:::

Keep the bars at unit conductance. From the chord ratio, weighting
diagonals by inverse squared length relative to the bars gives
$w_{\mathrm{length}}=(1+t^2-t)/(2t^2)$. Its difference from the tuned
weight factors as $(2t-1)(t^2+t-1)/(2t^2)$. A unit hull meets the tuning
condition precisely at the regular pose. In contrast, the positive
weights from $w=1/t-t$ keep the coordinate space invariant at all
distinct-joint poses. The `jitterbug-weight-policy-obstruction` probe
checks the lag defect by projecting full weighted transition matrices.

## Two reductions at the octahedral limit {#sec:jitterbug-octahedral-reductions}

Let $D$ copy one value at each octahedral site onto its pair of joints.
The matching Laplacian is zero on the range of $D$ and acts as twice
the identity on pair differences. Sending $w$ to infinity confines
the centred Gaussian to fields of the form $x=Du$.

::: {.proposition #prop:jitterbug-octahedral-gaussian-limit}
In the common site values $u$, the limiting Gaussian precision is
$2L_{\mathrm{oct}}$, where $L_{\mathrm{oct}}$ gives each octahedral
edge unit conductance. The full contrast variances for adjacent and
opposite sites are {{JITTERBUG_OCTAHEDRAL_FULL_CONTRASTS}}. In the
lowest coordinate modes the unit-bar squared chord factor is
{{JITTERBUG_OCTAHEDRAL_UNIT_BAR_CHORD_FACTOR}}.
:::

Each octahedral edge represents two bars between the corresponding
pairs, so $D^{\mathsf T}L_{\mathrm{bars}}D=2L_{\mathrm{oct}}$. In the
orthonormal basis $D/\sqrt{2}$ the same restriction is
$L_{\mathrm{oct}}$. At finite weight the pair averages have a
well-defined marginal, found by eliminating the pair differences
rather than clamping them. In orthonormal pair coordinates the marginal
precision eigenvalue on each of the three odd modes is
$4-4/(4+2w)\to4$ as $w\to\infty$. The even block keeps eigenvalue
{{JITTERBUG_EVEN_DIAGONAL_EIGENVALUE}} on the two trace-free axis modes.
In this basis, Gaussian elimination shows that the limiting
cross-covariance and the pair-difference covariance vanish.

`jitterbug-octahedral-gaussian-limit` checks the quotient and contrasts
with the exact network matrices.

::: {.proposition #prop:jitterbug-pair-load-memory}
Projecting the labelled binary repair process to the pair loads at
any finite $w$ produces memory even in stationarity. At the tuned
half-ratio pose the chance of returning from a specified split
load to a specified double load is
{{JITTERBUG_STATIONARY_RETURN_PROBABILITIES_HALF_POSE}} when only the
current split load is given and when the previous double load is
also given, respectively.
:::

Choose adjoining pairs such that the two bars between them share an
endpoint in the first pair. Given one token in each pair, only half
of the compatible label states can move to double occupancy in that
first pair. Their transition probability is $1/(2W)$, where
$W=24+6w$ is total conductance. Uniform microstate stationarity gives
$1/(4W)$ when only the current split is given. By reversibility,
knowing the previous double load selects exactly the label states
able to return, giving $1/(2W)$. The `jitterbug-pair-load-memory` probe
checks both the unequal transitions within a fibre and this stationary
conditioning witness. A Gaussian limit on field values therefore
does not supply a finite-weight Markov quotient for pair loads.

::: {.proposition #prop:jitterbug-averaged-repair-limit}
On the clock $s=k/W$, the tuned limit gives a Markov generator on
{{JITTERBUG_PAIR_LOAD_COUNT}} pair-load states obtained from
{{JITTERBUG_MICROSTATE_COUNT}} binary label states. For neighbouring
sites with occupancies $n_i$ and $n_j$, the rate of moving one token
from $i$ to $j$ is $n_i(2-n_j)/4$. The mean operator is
$-L_{\mathrm{oct}}/2$, and the normalised coordinate covariance
decays exponentially at rate {{JITTERBUG_MASS_CLOCK_DECAY_RATE}}.
:::

On the accelerated clock, within-pair attempts dominate and mix the
label arrangements uniformly at each fixed vector of pair totals.
Let $Q_B$ and $Q_M$ be the sums of fair single-seam update operators
minus the identity on the bars and matching, respectively. Then the
per-attempt operator is $P=I+(Q_B+wQ_M)/W$. The orthogonal projection
$\Pi$ onto functions of pair totals averages uniformly within each
fibre. On its complement the limiting matching-only operator is a
strict contraction because the matching updates connect each fibre
and have self-transitions. The off-diagonal blocks between these
spaces are of order $W^{-1}$. Over $k=\lfloor Ws\rfloor$ attempts the
fast block contracts geometrically, and the sum of its feedback into
the slow block is of order $W^{-1}$. The remaining slow block gives
the semigroup with generator $\Pi Q_B\Pi$.

For either bar between adjacent pairs the probability of a token at
its source and a vacancy at its target, averaged within the fibre,
is $(n_i/2)(1-n_j/2)$. There are two such bars, and each fair update
moves the token with probability one half when its endpoints differ.
This gives the stated rate. In the difference between flows in the
two directions the product $n_i n_j$ cancels, giving the mean
operator. The `jitterbug-averaged-repair-limit` probe enumerates the
microstate updates and checks every entry of the generator against
this formula. Up to the clock factor, the resulting process is symmetric
partial exclusion with two slots per site [@giardina2009duality, §3.2].
Its stationary pair-load probabilities are proportional to the number of
compatible binary states. Their covariance equals
{{JITTERBUG_CENTRED_COUNT_COVARIANCE}} times the centred projector. The
mean operator therefore determines all stationary linear count lag
covariances, with decay on each mode of $L_{\mathrm{oct}}$ at half its
eigenvalue. A Gaussian Markov process can reproduce these covariances
with the same drift and matched stationary variance. The occupancy
distribution is discrete, and its covariance differs from the static
field covariance in [@prop:jitterbug-octahedral-gaussian-limit].

The single-attempt lag factor on the tuned coordinate modes is
$\rho(t)=1-(4-2t)/(2W)$. At fixed $k$, the correlation $\rho(t)^k$ tends
to one as $t$ tends to zero. For each fixed positive $t$,
$\rho(t)^k$ decays to zero as $k$ grows, so the two limits do not
commute. Keeping $k/W$ fixed
instead gives $\exp(-2s)$, the same decay as the averaged octahedral
coordinate modes. This separates a geometric endpoint, a static
precision limit and a limit of repair dynamics.

![Fuller's Jitterbug at a common bar length. The solid bars keep their incidences, and the added hull diagonals are drawn separately. Beyond the regular icosahedron the diagonals are shorter than the bars. The lower plots compare the geometric length ratio, tuned conductance and the inverse-square length policy.](../output/figures/fig_jitterbug.png){#fig:jitterbug width=100%}
