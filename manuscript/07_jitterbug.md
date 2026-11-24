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
different ways to read the assembly. [@fig:jitterbug] compares the bars and
hull across poses at a common bar length. [@tbl:ontology-jitterbug-linkage]
states the scope of the geometric, Gaussian and repair comparisons.

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

The two alternatives in [@prop:jitterbug-weight-policy-obstruction] keep
the bars at unit conductance. From the chord ratio, weighting
diagonals by inverse squared length relative to the bars gives
$w_{\mathrm{length}}=(1+t^2-t)/(2t^2)$. Its difference from the tuned
weight factors as $(2t-1)(t^2+t-1)/(2t^2)$. A unit hull meets the tuning
condition precisely at the regular pose. In contrast, the positive
weights from $w=1/t-t$ keep the coordinate space invariant at all
distinct-joint poses. The `jitterbug-weight-policy-obstruction` probe
checks the lag defect by projecting full weighted transition matrices.

## Moving poses and the missing coordinate {#sec:jitterbug-moving}

The posewise result in [@prop:jitterbug-coordinate-closure] compares repeated
attempts at a fixed geometry. Consider instead three observation times at
poses $a,b,c$ with $0<t\leq1$ for each pose ratio $t$, allowing repeats.
Use one fair binary repair attempt at each departure pose, first $a$ and
then $b$. At each attempt the bars have unit conductance
and the closing seams use the tuned weight $w(t)=1/t-t$. Each repair kernel
preserves the uniform law on the two-token shell, so the same uniform law
holds at all three times when used initially. The protocol sets the poses
and attempt times independently of the repair outcomes.

::: {.proposition #prop:jitterbug-moving-closure}
Let $Y_t$ be the sum of the two occupied unscaled joint positions at pose $t$.
For the protocol above, write $C_{ij}$ for their cross-covariance divided by
{{JITTERBUG_MOVING_COVARIANCE_SCALE}}. The normalised Gaussian Markov covariance residual at the middle
pose is
$$
C_{02}-C_{01}C_{11}^{-1}C_{12}
=\rho(a)\sigma(b)\,\frac{(a-b)(c-b)}{1+b^2}\,I_3,
$$
where $\rho(t)=1-(4-2t)/(2W_t)$,
$\sigma(t)=1-(4+2/t)/(2W_t)$, and $W_t=24+6w(t)$.
It is strictly negative on the diagonal during strict contraction
$a>b>c>0$, and also during strict expansion $0<a<b<c$. A reversal
with unequal legs gives a positive diagonal residual. Holding either
leg fixed gives zero.
:::

To obtain [@prop:jitterbug-moving-closure], the stationary covariance of the
binary occupancy vector is a quarter of that scale times the centred
projector. A coordinate
axis uses two disjoint sign-supported vectors, each of squared norm four. Its
coefficients are $b_t=(t,1)$ and its orthogonal companion is $v_t=(1,-t)$.
In this two-dimensional block, tuned mean repair has eigenvectors $b_t$ and
$v_t$ with eigenvalues $\rho(t)$ and $\sigma(t)$. Conditioning a Gaussian
on the middle coordinate subtracts its projection onto $b_b$. The residual
therefore contains $b_a^{\mathsf T}v_b v_b^{\mathsf T}b_c/(1+b^2)$,
with the two stated retention factors. Both factors are positive. For
positions scaled to a fixed bar length, the residual is multiplied by the
positive scale factors of the first and last poses. If the departure poses
are held for $m$ and $n$ attempts, the two factors become $\rho(a)^m$ and
$\sigma(b)^n$. The sign persists for finite attempt counts, including
zero. It follows from the changing observation subspace.

The `jitterbug-moving-closure-obstruction` probe enumerates all
{{JITTERBUG_MICROSTATE_COUNT}} binary states and both weighted transition
matrices. At the three post-icosahedral poses {{JITTERBUG_MOVING_POSES}},
the normalised residual is {{JITTERBUG_MOVING_DEFECT}}. The reversal control
gives {{JITTERBUG_REVERSAL_DEFECT}}. The strict expansion control gives a
negative residual, and holding either leg gives zero.
A Gaussian process on the measured coordinates would require zero
conditional covariance between the endpoints given the middle observation.
[@fig:jitterbug-observation-limits] shows the signed test as the final pose varies.

::: {.proposition #prop:jitterbug-augmented-lag}
Retaining both sign-supported slots for each of the three spatial axes
gives a fixed {{JITTERBUG_ODD_READOUT_DIMENSION}}-dimensional readout.
For any sequence of nonnegative closing weights, this readout admits a
Gaussian Markov process with the same linear lag covariances as fair repair.
Any fixed linear readout space that contains all three geometric coordinates
at two distinct poses must contain these six dimensions.
:::

Both the bar and matching kernels leave the six slots in
[@prop:jitterbug-augmented-lag] invariant. After normalising their common
stationary covariance to the identity, the transition restricted to the
slots is a symmetric contraction $A_i$ at step $i$. A Gaussian recursion
with transition $A_i$ and independent innovations of covariance $I-A_i^2$
preserves the identity as its marginal covariance. Its cross-covariances
are the same products of the restricted transitions that compute the
repair cross-covariances. This uses a time-dependent linear state-space
representation [@kalman1960filtering]. The
`jitterbug-augmented-lag-realization` probe verifies invariance for the two
constituent kernels and the two-step covariance identity in the augmented
readout. At each pose, including the companion $v_t$ alongside $b_t$ recovers
the fixed six-slot basis. The resulting agreement concerns second moments.


The normalised formulas also hold for any uniform fixed binary load with
$1\leq k\leq11$. The centred occupancy covariance scales by
$k(12-k)/(12\cdot11)$, and the linear mean repair is unchanged. The
enumerated probes use $k=2$.

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

To obtain the return witness in [@prop:jitterbug-pair-load-memory], choose
adjoining pairs such that the two bars between them share an
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

On the clock in [@prop:jitterbug-averaged-repair-limit], within-pair attempts
dominate and mix the label arrangements uniformly at each fixed vector of
pair totals.
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
[@fig:ga-jitterbug] puts the geometric endpoint alongside the two model limits.

![Fuller's Jitterbug at a common bar length. The bars keep their incidences,
and dashed lines mark the hull diagonals. Beyond the regular icosahedron the
diagonals are shorter than the bars. The lower plots show the length ratio
and compare the conductance that preserves coordinate lag covariance with
inverse-square weighting and the dotted unit baseline.](../output/figures/fig_jitterbug.png){#fig:jitterbug width=100%}

## Finite pair means and their hidden contrasts {#sec:jitterbug-finite-projection}

The limit in [@prop:jitterbug-octahedral-gaussian-limit] gives hard equality
within each closing pair. At finite weight one can instead observe the
arithmetic means of its two field values. Let $U$ embed orthonormal pair
means and $V$ embed pair differences. Work on the centred label space and
write the precision blocks in these bases as $A=U^{\mathsf T}LU$,
$B=V^{\mathsf T}LU$, and $H=V^{\mathsf T}LV$. The finite Gaussian marginal
precision in orthonormal mean coordinates is $A-B^{\mathsf T}H^{-1}B$.
The hard equality restriction has precision $A$. For the associated fair
repair, let $\Delta=B^{\mathsf T}B/(4W^2)$ be the normalised lag-two covariance
defect of pair totals, by [@prop:lag-closure-transfer]. For the fixed-shell
fair repair and any complementary orthonormal split of the centred label
space with positive $H$, the Schur correction $A-S$, where
$S=A-B^{\mathsf T}H^{-1}B$, has the same nullspace as $\Delta$. In positive
semidefinite order,
$$
\frac{4W^2}{\lambda_{\max}(H)}\Delta
\preceq A-S \preceq
\frac{4W^2}{\lambda_{\min}(H)}\Delta.
$$
At the Jitterbug pairs $H$ is scalar, so both bounds coincide.

::: {.proposition #prop:jitterbug-pair-projection}
For finite matching weight $w\geq0$, put $h=4+2w$ and let $\Pi_-$ project onto
the three odd opposite-site modes of the octahedron. The precision blocks
satisfy $A=L_{\mathrm{oct}}$, $H=hI$, and $B^{\mathsf T}B=4\Pi_-$.
The marginal precision in literal pair means is
$2L_{\mathrm{oct}}-(8/h)\Pi_-$. It has positive opposite-site entries
$4/h$. The normalised repair lag-two defect of pair totals is
$\Pi_-/W^2$, with $W=24+6w$.
:::

For [@prop:jitterbug-pair-projection], each unit octahedral edge represents
two bars, and their action on pair sums gives $A=L_{\mathrm{oct}}$ in the
orthonormal basis. The three odd modes couple with coefficient magnitude
two to their pair differences. The two centred even modes have no such
coupling. The hidden block has eigenvalue $h$ on all six difference modes.
Literal means are $1/\sqrt{2}$ times the orthonormal coordinates, which
doubles precision. The odd projector has opposite-site entry minus one half,
so the finite marginal's zero-row-sum precision cannot be realised with
nonnegative octahedral edge conductances. Ordinary Kron reduction removes
selected vertex potentials and preserves the Laplacian sign pattern
[@dorfler2013kron]. Here the eliminated coordinates are differences in a
rotated mean-and-difference basis.

At the tuned half-ratio pose, $w={{JITTERBUG_HALF_POSE_WEIGHT}}$, the
opposite-site entry is {{JITTERBUG_PAIR_MEAN_OPPOSITE_PRECISION}}.
Adjacent and opposite marginal variances are
{{JITTERBUG_PAIR_MEAN_CONTRAST_VARIANCES}}, compared with the hard-limit
values in [@prop:jitterbug-octahedral-gaussian-limit]. The odd lag-two
covariance defect at the same pose is {{JITTERBUG_PAIR_ODD_LAG_DEFECT}}.
The `jitterbug-pair-projection-defect` probe checks the blocks and Gaussian
contrasts using the exact twelve-vertex Green matrix. It checks the lag
defect by aggregating the two-step counts of all binary microstates. The
static correction decays with inverse closing weight, while the one-attempt
clock lag defect decays with inverse weight squared. The separate finite
pair-load memory result is tested in [@prop:jitterbug-pair-load-memory].

![Two tests of Jitterbug readouts. Left, with the two departure poses fixed, contraction gives a negative residual in the Gaussian Markov covariance identity. A hold gives zero, and reversal gives a positive residual. Right, the finite pair-mean precision has a positive opposite-site entry and the pair-total lag covariance has a positive defect on odd opposite-site modes and zero on even modes. The highlighted conductance is the tuned weight at the middle pose.](../output/figures/fig_jitterbug_observation_limits.png){#fig:jitterbug-observation-limits width=100%}

The same tests extend to prescribed observations on other finite graphs.
See [@prop:moving-readout-history] for the history criterion and
[@prop:balanced-pair-transfer] for the balance condition on fixed pairs.
The tetrahedral control in [@prop:switching-pair-obstruction] shows that
changing the pairing can undo an exact fixed-partition count quotient.
