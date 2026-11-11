# Observers, blankets, and screens — outline

This is the prose skeleton. Every number in it is written as a bound key from
`manuscript_variables/values.py`, so the text cannot quote a value the code does
not produce. When a sentence needs a number, add or reuse a `BoundValue` and
render [docs/MANUSCRIPT_VALUES.md](../docs/MANUSCRIPT_VALUES.md).

Each section names the registry entries it rests on. A section may not claim
more than its entries' statuses.

The rendered edition adds a cover page (page 1) and a four-panel graphical
abstract (pages 2 to 5) written in `00_00_graphical_abstract.md`, with the
contents page re-emitted after the panels. The registry carries a pairwise scope
per entry, and `oph-fep-synergetics verify --pair` executes one scope at a time;
`--kind` narrows a run to the identifications, the obstructions or the
complementarities. `fig:network` draws the carrier's network readings — the
resistance ladder, the shells, the subgraph census and the dual pair — beside
the registry and triangle figures.

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

The same twelve-around-one configuration is the kissing-number-twelve
arrangement and the IVM's first shell: `ivm.neighbor_moves` sites at squared
distance `ivm.neighbor_squared_distance`. This is a genuine identity of one
object seen three ways, and it is the reason the three programs can be put on
one page at all.

## 2. Agreement as a minimisation

**Entries:** `variational-equivalence`, `free-energy-dictionary`.

On the linear-Gaussian corner, three descriptions of "settle to agreement"
coincide: the minimiser of the interface energy over the interior, the harmonic
extension of the records, and the posterior mean. The probe reports all three
and the residual of the perception flow, and the excess free energy vanishes at
the posterior.

The dictionary that makes this exact: surprisal is what no repair removes,
Bayesian surprise is what repair removes, and the update residual pins the
Bayesian surprise between the extreme interior eigenvalues.

**Status note.** Both entries are `exact-on-domain`. The nonlinear and
non-Gaussian cases are outside this work; see [docs/SCOPE.md](../docs/SCOPE.md).

## 3. Where the three programs disagree

**Entries:** `pr01-obstruction`, `free-energy-name-collision`,
`mediator-obstruction`, `stationary-blanket`, `solenoidal-freedom`,
`d-stability`, `schedule-path`, `reduction-criterion`,
`existence-without-reachability`, `directed-gradient`,
`four-cycle-local-minimum`, `cuboctahedron-vs-icosahedron`,
`shell-vs-free-energy`, `vertex-shells`, `twelve-around-one-families`.

These are results, not caveats. Each is an executable demonstration that an
identification that looks natural in prose fails as stated.

- **The native law is not a descent method.** Of `shell.transitions` native
  transitions, `shell.strict_descent` lower the quadratic energy,
  `shell.zero_descent` move at mismatch one without lowering it, and
  `shell.noop` do nothing. A one-unit walk keeps energy `1/2` while changing
  state. OPH records this as its own obstructed premise, and this repository
  reproduces it exhaustively.
- **"Free energy" names two different objects.** The interior-minimised
  interface energy and the variational free energy at the posterior differ by a
  model-dependent constant. Writing one for the other is a gauge choice.
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
  `directed_gradient.adjoint`.
- **A local minimum can be incompatible.** A four-cycle assignment can violate
  two edges while every single-bit flip also violates two.
- **Twelve vertices is not one solid.** The neighbour graph is
  `cuboctahedron.degree`-regular with tetravolume
  `cuboctahedron.tetravolume`; the carrier is `icosahedron.degree`-regular with
  tetravolume `icosahedron.tetravolume`. A shared vertex count is not an
  identification.
- **The shell index is not a function of the free energy.** The first five
  shells carry `shell.free_energy_level_counts` distinct levels and
  `shell.shared_free_energy_levels` of them are shared between shells, so a
  bounded free energy does not bound the shell.

## 4. What each program supplies that the others do not

**Entries:** `equilibrium-propagation`, `replicator-potential`,
`native-repair-law`, `certificate-species`, `quadray-ivm`,
`icosahedron-tetravolume`, `packing-density`, `tetravolume-arithmetic`,
`dual-cells`, `conjugation-norm`.

- **OPH** supplies a concrete, finite, conserving law on a fixed carrier, with
  an exact energy identity and a sharp certificate — a spectral gap, not a
  sup-norm contraction.
- **FEP** supplies the parameter gradient the repair law does not define. The
  two-settle contrast is exactly `(theta - target) / (1 + beta)` and tends to
  the exact loss gradient as the nudge vanishes. This is a complementarity, not
  an equivalence: it extends the family rather than identifying a member.
- **Synergetics** supplies exact integer geometry: rational tetravolumes
  including `ivm.volume_ratios`, and the shell sequence
  `ivm.shell_counts` shared with the Goldberg and geodesic counts.
- **Synergetics also supplies the volume of the carrier's own solid**, which its
  table does not contain: `icosahedron.tetravolume` tetravolumes, whose square is
  `icosahedron.tetravolume_squared` against the exact `175 + 75 sqrt 5`, and the
  twelve-around-one cell, whose `packing.cell_tetravolume` and packing density
  `packing.density` are derived from the same quadratic form rather than read
  from the table.
- **No single programme supplies the exchange rate.** One lattice step has three
  exact costs in three number fields (`tetravolume-arithmetic`), which is why the
  only three-way entry is a complementarity rather than an identification.

The replicator potential is included as a comparison class: selection climbs a
potential, so the repair language covers it without sharing its mechanism.

## 5. Scope and limits

**Entries:** none — this section states what is *not* claimed.

- The isomorphism claims are `exact-on-domain` where a domain restriction
  applies, and the domain is stated in the entry itself.
- OPH's quantum event algebra and continuum limit are cited, not implemented.
- FEP outside the linear-Gaussian corner is out of scope.
- Fuller's apparatus appears only in the pieces the correspondence needs.
- Scale separation is not addressed: OPH's patches and an FEP blanket are not
  mapped onto one another here.
- `35/264` port-zero occupation and the `13/15` and `1/60` move probabilities
  are `not-reproduced`, because the source text does not contain the
  record/selector protocol they need.

See [docs/SCOPE.md](../docs/SCOPE.md) for the full statement.
