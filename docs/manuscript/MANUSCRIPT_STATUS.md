# Observers, blankets, and screens — outline

This is the prose skeleton. Every number in it is written as a bound key from
`manuscript_variables/values.py`, so the text cannot quote a value the code does
not produce. When a sentence needs a number, add or reuse a `BoundValue` and
render [docs/MANUSCRIPT_VALUES.md](../MANUSCRIPT_VALUES.md).

Each section names the registry entries it rests on. A section may not claim
more than its entries' statuses.

The cover is on page 1. The five graphical abstract panels in
`00_00_graphical_abstract.md` occupy pages 2 to 6. The Jitterbug poster on
page 6 follows the linkage to octahedral coalescence. The contents starts on
page 7. Each entry has one of three pairwise scopes or the three-way scope.
Use `oph-fep-synergetics verify --pair` to select a scope and `--kind` to
select a claim kind. The [generated concordance](../CONCORDANCE.md)
records where prose invokes each numbered figure, table, equation or formal
statement.

---

## 1. One object, three languages

**Entries:** `carrier-port-chart`, `sheaf-agreement`, `spectral-gap`,
`shell-counts`, `kissing-twelve`, `quadray-ivm`, `tree-number`,
`gmrf-independence`, `sheaf-cohomology`, `icosahedral-symmetry`,
`eigenspace-irreducibility`, `law-equivariance`, `kl-lattice-metric`,
`subgraph-census`, `angle-defect`, `effective-resistance`.

The carrier is a fixed finite object: `carrier.ports` ports, `carrier.seams`
seams, `carrier.faces` faces, every port meeting `carrier.degree` triangles.
Its port chart is not taken on trust — the adjacency is derived from OPH's own
twenty oriented faces, and the derivation is asserted against its consequence
(the triangles of the derived graph are exactly those faces).

Patch agreement — the requirement that neighbouring patches agree on their
overlap — is the global-section condition of a cellular sheaf on that complex.
The sheaf Laplacian is the carrier Laplacian, so the quadratic energy of
disagreement is the carrier's own quadratic form.

Twelve around one is where the three programs touch, and it is not one object:
the carrier is the icosahedral kissing arrangement, with slack between
neighbours, while the IVM's first shell is the cuboctahedral one,
`ivm.neighbor_moves` sites at squared distance `ivm.neighbor_squared_distance`,
with contacts between outer spheres. Both attain the kissing number. The
identities of this section are about the carrier itself, and the
twelve-around-one conflation is an obstruction in section 3.

## 2. Agreement as a minimisation

**Entries:** `variational-equivalence`, `free-energy-dictionary`,
`repair-is-conditioning`, `two-record-resistance`,
`binary-conditional-variance`, `carrier-spectral-chords`.

In the linear-Gaussian carrier model, three descriptions of "settle to agreement"
coincide: the harmonic extension of the records (an interior solve), the
constrained relative-entropy minimiser (a separate KKT solve), and the limit of
the perception flow. The free energy at the posterior, computed directly as
energy minus entropy, equals the surprisal; away from the posterior it exceeds
it by exactly the divergence.

Record surprisal uses the normalised marginal with a proper prior. Bayesian
surprise measures the change of posterior and is bracketed by the update
residual and the extreme interior eigenvalues. For two carrier records, the
Schur precision governing the contrast is set by the effective resistance.

The Gaussian claims are `exact-on-domain`. An asymmetric spring gives a
non-Gaussian counterexample to equality of posterior mean and energy minimum.

## 3. Where the three programs disagree

**Entries:** `pr01-obstruction`, `free-energy-name-collision`,
`mediator-obstruction`, `stationary-blanket`, `solenoidal-freedom`,
`d-stability`, `schedule-path`, `reduction-criterion`,
`existence-without-reachability`, `directed-gradient`,
`four-cycle-local-minimum`, `cuboctahedron-vs-icosahedron`,
`shell-vs-free-energy`, `vertex-shells`, `twelve-around-one-families`.

These are results, not caveats. Each is an executable demonstration that an
identification that looks natural in prose fails as stated.

- **Native repair permits motion at constant load energy.** The
  `shell.zero_descent` changes at mismatch one can continue indefinitely.
  On the total-load-two shell the ungated relation is confluent with no normal
  forms, while the broken-seam gate terminates with multiple reachable normal
  forms from a doubled state. On the binary shell the ensemble relative
  entropy to the uniform law is nonincreasing under the uniform schedule.
- **Three values at fixed records.** Minimum interface energy, the raw Gaussian
  record potential and normalised surprisal have different meanings. The
  interior determinant contributes to the raw potential and cancels in the
  normalised probability.
- **Eliminating an interior mediator is not free.** The reduced precision is
  `(1/3) v v^T + epsilon I`; a sum of pairwise semidefinite terms reproduces it
  if and only if `epsilon >= mediator.threshold`. Below the threshold the
  reduction is not merely inaccurate, it is unavailable.
- **A stationary blanket does not bound the dynamics.** Two states with equal
  blanket statistics can differ in internal motion. Stationarity of a
  distribution is not conditional independence of the mechanism.
- **One stationary law, many laws.** Antisymmetric drift leaves the stationary
  density unchanged while changing the spectrum and the entropy production.
- **A free-energy form needs stability.** The Lyapunov pair exists with
  `d_stability.lyapunov_determinant` for the stable case; at `k = 4` the Routh
  first column `d_stability.routh_k4` has a negative entry and the scaled
  matrix has two right-half-plane roots, so the form is gone while the rest is
  not.
- **Finite phases are path functions.** Full and half-damped iteration of the
  same affine rule share every answer, yet their centered contrasts after
  twelve steps differ by exactly `schedule.affine_contrast`. Agreement and
  dissipation are state functions; a finite schedule is not.
- **Dropping a record costs.** Two histories can share the present observation
  and need different futures; the worst-case error floor is `reduction_criterion`.
- **Existence is not reachability.** The clipped two-node map has the origin as
  its only equilibrium and a source there. The relaxation orbit stays outside
  `two_node.relaxation_min_radius` with period `two_node.relaxation_period`,
  and the synchronous map stays outside `two_node.sync_min_radius`.
- **A forward response is not a gradient.** For a directed network the direct
  nudging response and the adjoint differ; the direct route recovers
  `directed_gradient.direct` where the gradient needs
  `directed_gradient.adjoint`, and only the adjoint is the gradient.
- **A local minimum can be incompatible.** A four-cycle assignment can violate
  two edges while every single-bit flip also violates two.
- **Stable is not D-stable.** The stable control matrix admits a Lyapunov (free
  energy) form; one positive local gain, `k = 4`, removes it with
  `d_stability.right_half_plane_count` right-half-plane roots.
- **Twelve vertices is not one solid.** The neighbour graph is
  `cuboctahedron.degree`-regular with tetravolume
  `cuboctahedron.tetravolume`; the carrier is `icosahedron.degree`-regular with
  tetravolume `icosahedron.tetravolume`. Both are kissing arrangements, the
  cuboctahedral one with outer contacts and the icosahedral one with gaps
  (`kissing.icosahedron_edge_ratio_squared`). A shared vertex count is not an
  identification.
- **The shell index and the free energy determine neither each other.** The
  shell is the cuboctahedral gauge (checked through shell
  `shell.gauge_checked_through`); the first five shells carry
  `shell.free_energy_level_counts` distinct levels and the levels
  `shell.shared_free_energy_levels` are each shared by two shells. The free
  energy does bound the shell from above; what fails is the identification.

## 4. What each program supplies that the others do not

**Entries:** `equilibrium-propagation`, `replicator-potential`,
`native-repair-law`, `certificate-species`, `quadray-ivm`,
`icosahedron-tetravolume`, `packing-density`, `tetravolume-arithmetic`,
`dual-cells`, `conjugation-norm`.

- **OPH** supplies a concrete, finite, conserving law on a fixed carrier, with
  an exact energy identity and a sharp certificate — a spectral gap, not a
  sup-norm contraction.
- **Energy-based learning**, which the FEP reaches through predictive coding,
  supplies the parameter gradient the repair law does not define. The
  two-settle contrast is exactly `(theta - target) / (1 + beta)` and tends to
  the exact loss gradient as the nudge vanishes. This is a complementarity, not
  an equivalence: it extends the family rather than identifying a member.
- **Synergetics** supplies exact integer geometry: the tetravolume table
  `ivm.volume_ratios`, every entry derived from an explicit vertex set, and the
  coordination-shell sequence
  `ivm.shell_counts` shared with the Goldberg and geodesic counts.
- **Synergetics also supplies the volume of the carrier's own solid**, which its
  table does not contain: `icosahedron.tetravolume` tetravolumes, whose square is
  `icosahedron.tetravolume_squared` against the exact `175 + 75 sqrt 5`, and the
  twelve-around-one cell, whose `packing.cell_tetravolume` and packing density
  `packing.density` are derived from the same quadratic form rather than read
  from the table.
- **The three models compute different quantities.** The lattice
  displacement divergence, mean repair gap and normalised volumes use
  separate units and arithmetic conventions (`tetravolume-arithmetic`).
  The planar-dual result below gives a separate three-way identity.

The replicator potential is included as a comparison class: selection climbs a
potential, so the repair language covers it without sharing its mechanism.

## Exact extensions

**Entries:** `repair-first-entry`, `repair-distance-lumping`,
`grounded-forest-precision`, `greedy-record-placement`,
`planar-dual-covariance`, `cuboctahedral-spectral-chords`,
`cuboctahedral-distance-split`.

The extensions distinguish a first-passage distribution from equilibrium,
give an exact graph-distance reduction, compare record layouts using rooted
forest counts, and test spectral geometry and planar duality on two solids.

## 5. Scope and limits

**Entries:** none — this section states what is *not* claimed.

- The isomorphism claims are `exact-on-domain` where a domain restriction
  applies, and the domain is stated in the entry itself.
- OPH's quantum event algebra and continuum limit are cited, not implemented.
- The inference results concern the stated Gaussian models, the asymmetric
  seam control and finite conditional repair. The treatment of general
  nonlinear variational dynamics remains open.
- Fuller's apparatus appears only in the pieces the correspondence needs.
- Scale separation is not addressed: OPH's patches and an FEP blanket are not
  mapped onto one another here.
- The `35/264` port-zero occupation and the `13/15` and `1/60` transition
  probabilities are reproduced exactly from the protocols in OPH's public model
  (`reduction-criterion`, `pr01-obstruction`).

See [docs/SCOPE.md](../SCOPE.md) for the full statement.

## Three-way correspondences and counterexamples

**Entries:** `repair-contrast-quotient`, `repair-distance-counterexample`,
`reciprocal-dual-weights`, `copied-dual-weight-obstruction`,
`quadratic-even-field`, `quadratic-isotropy-obstruction`,
`central-plane-screening`, `repair-centroid-memory`.

The comparison applies the same fair binary repair law to the two
twelve-vertex graphs. Exact covariance records preserve the dynamics after
coarse-graining, while graph distance on the FCC shell and the observed
centroid on either shell lose information about the next step. Three-time
probabilities distinguish the stationary centroid from its Gaussian
lag-covariance match. Matched planar-dual covariances require reciprocal
conductances. The quadratic coefficient split and the central-plane census
give separate geometric counterexamples to transfer of Gaussian isotropy
and screening. The chapter and structured concordance specify each
protocol and cite its executable probe. The two figures show the
conditional-probability and geometry results directly.

## Transfer criteria beyond the two shells

**Entries:** `distance-regular-repair`, `bipartite-contrast-collision`,
`periodic-fcc-contrast-collision`, `linear-lag-closure`,
`single-seam-lag-obstruction`, `surface-harmonic-covariance`,
`copied-surface-weight-obstruction`, `root-shell-quadratic-precision`,
`fcc-plane-screen-transfer`, `centered-gaussian-screen-obstruction`,
`moving-readout-history`, `consecutive-readout-obstruction`, `balanced-pair-transfer`,
`pair-precision-sign-obstruction`, `pair-count-lumpability-obstruction`,
`switching-pair-obstruction`.

The transfer chapter gives conditions for a resistance quotient,
Gaussian lag covariance closure and reciprocal duality on closed
surfaces. The cube, periodic FCC and strengthened 24-cell examples
show which data determine these properties. The root-shell calculation
separates diagonal and mixed modes. The FCC plane proof extends beyond
the first shell and compares grounded, noisy-sum and exact-sum conditioning.
The moving square tests all recorded splits. Paired paths separate the
precision sign and lag tests. A square and a tetrahedron distinguish
fixed-pair linear agreement, full count lumpability and switching memory.

## Fuller's Jitterbug and the carrier graph

**Entries:** `jitterbug-rigid-motion`, `jitterbug-graph-obstruction`,
`jitterbug-weighted-precision`, `jitterbug-coordinate-closure`,
`jitterbug-nearest-only-obstruction`, `jitterbug-weight-policy-obstruction`,
`jitterbug-octahedral-gaussian-limit`, `jitterbug-pair-load-memory`,
`jitterbug-averaged-repair-limit`, `jitterbug-moving-closure-obstruction`,
`jitterbug-augmented-lag-realization`, `jitterbug-pair-projection-defect`.

The coordinate path keeps the eight triangles rigid. The hull already
has triangulated openings on the contracting branch before the short
diagonals reach the bar length. The Gaussian comparison distinguishes
bar-only precision from adding conductance to the short matching.
The latter gives a unique weight law for the linear coordinate modes
at each fixed pose. Past regularity, nearest-only seams disconnect the graph and
standard weight policies can break lag closure. At coalescence, the
hard-constraint Gaussian adds the parallel bar conductances, while
finite-weight pair counts retain memory. Averaging fast matching repair
on the accelerated clock gives an octahedral occupancy generator.

The moving-pose result tests Gaussian covariance closure across two successive
departure-pose updates. A fixed six-slot readout restores the covariance
match for arbitrary nonnegative weight schedules. Finite pair means and
pair totals expose the same hidden coupling in Gaussian precision and
repair lag covariances. The new figure and revised poster show these tests.
