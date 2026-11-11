# Method: a registry whose every claim is executed, scoped and re-runnable {#sec:method}

## The unit of work

The paper's unit is a **correspondence**: a claim about the same object made in
two or three vocabularies, together with executable evidence and an explicit
pairwise scope.

::: {#def:correspondence .definition}
A **correspondence** is a record $c = (\kappa, \rho, \sigma, \pi)$ where $\kappa$
names the claim in the OPH, FEP and Synergetics vocabularies, $\rho$ is the
*pairwise scope*, one of `oph-fep`, `fep-synergetics`, `synergetics-oph` or
`three-way`, naming the two or three programmes the claim relates, $\sigma$ is a
status from the fixed vocabulary below, and $\pi$ is a *probe*: a total function
from no arguments to a mapping from names to values. The probe **passes** when
every value it reports satisfies the assertion the claim makes.
:::

The scope $\rho$ is not bookkeeping. A correspondence is a claim that two
vocabularies are talking about one object, and which two they are is part of the
claim: the same icosahedral symmetry is a statement about OPH and Synergetics,
while the Markov blanket of the carrier precision is a statement about OPH and
the FEP. [@sec:correspondences] reports the three pairwise ledgers separately for
that reason, and `oph-fep-synergetics verify --pair synergetics-oph` executes one
ledger on its own.

Making the probe return values rather than a boolean is the single design choice
that keeps the registry honest. A boolean can be made true by weakening a
tolerance; a mapping of named values has to be *right*, and it can be compared
against a source that states a different value.

| Status | Meaning |
|:-------|:--------|
| `exact` | the probe reproduces the claim exactly: integers, rationals, or closed-form surds |
| `exact-on-domain` | the identification holds on a stated domain and the probe asserts the boundary |
| `obstructed` | the identification fails, and the probe exhibits the failure precisely |
| `not-reproduced` | the source states a quantity its text does not determine; recorded, never invented |

: The status vocabulary. It is closed: a new status requires a new definition.
{#tbl:status}

## The registry, grouped by pairwise scope

[@tbl:registry] is the registry, reproduced from
`src/oph_fep_synergetics/analysis/registry.py`. The grouping is the paper's
argument in miniature: {{REGISTRY_ISOMORPHISMS}} claims survive as
identifications, {{REGISTRY_DISAGREEMENTS}} fail as identifications while
remaining informative, and {{REGISTRY_COMPLEMENTARITIES}} are complements in
which one program supplies what another lacks without either reducing to the
other. By scope, {{REGISTRY_OPH_FEP}} claims relate OPH to the FEP,
{{REGISTRY_SYNERGETICS_OPH}} relate Synergetics to OPH,
{{REGISTRY_FEP_SYNERGETICS}} relate the FEP to Synergetics, and
{{REGISTRY_THREE_WAY}} needs all three vocabularies at once.

| Key | Kind | Status | Pair | What it tests |
|:-----------------------------------|:----------------|:-------------------|:-----------------|:-------------------------------|
| `carrier-port-chart` | isomorphism | `exact` | `synergetics-oph` | derived adjacency reproduces the twenty oriented faces |
| `spectral-gap` | isomorphism | `exact` | `oph-fep` | $\lambda_1 = 5-\sqrt{5}$ and normalised gap $(5-\sqrt{5})/60$ |
| `shell-counts` | isomorphism | `exact` | `oph-fep` | $\binom{13}{11}$ shell size and the binary sub-shell |
| `native-repair-law` | isomorphism | `exact` | `oph-fep` | OPH's integer rule reimplemented, with its closed-form energy drop |
| `variational-equivalence` | isomorphism | `exact-on-domain` | `oph-fep` | harmonic-extension residual equals a Schur complement |
| `free-energy-dictionary` | isomorphism | `exact-on-domain` | `oph-fep` | numerical, Bayesian and surprisal free energy coincide to curvature |
| `sheaf-agreement` | isomorphism | `exact` | `oph-fep` | the cellular sheaf Laplacian is the carrier Laplacian |
| `certificate-species` | isomorphism | `exact` | `oph-fep` | the certificate species of the shell and its spectral gap |
| `tree-number` | isomorphism | `exact` | `oph-fep` | the reduced Laplacian determinant is 5 184 000, an integer tree count |
| `gmrf-independence` | isomorphism | `exact` | `oph-fep` | independence on the seams and nowhere else; blanket five, variance $7/36$ |
| `sheaf-cohomology` | isomorphism | `exact` | `oph-fep` | $\dim H^0 = 1$, $\dim H^1 = 19$, and $\dim H^0 - \dim H^1 = 12 - 30$ |
| `icosahedral-symmetry` | isomorphism | `exact` | `synergetics-oph` | order $120$, transitive on ports, seams and faces: $60$ rotations, $60$ reflections |
| `eigenspace-irreducibility` | isomorphism | `exact` | `synergetics-oph` | eigenspaces of dimensions 3, 5 and 3 irreducible under the icosahedral group |
| `law-equivariance` | isomorphism | `exact` | `synergetics-oph` | $120 \times 78$ pairs preserve the outcome multiset: 4 orbits, 3 multisets |
| `kl-lattice-metric` | isomorphism | `exact` | `fep-synergetics` | one IVM step costs four nats and squared norm $8k^2$ costs $4k^2$ |
| `pr01-obstruction` | disagreement | `obstructed` | `oph-fep` | the native law is not confluent |
| `free-energy-name-collision` | disagreement | `obstructed` | `oph-fep` | "free energy" names inequivalent functionals |
| `mediator-obstruction` | disagreement | `obstructed` | `oph-fep` | pairwise boundary reduction needs a threshold |
| `stationary-blanket` | disagreement | `obstructed` | `oph-fep` | a stationary blanket is not an attractor |
| `solenoidal-freedom` | disagreement | `obstructed` | `oph-fep` | stationary blankets leave entropy production free |
| `d-stability` | disagreement | `obstructed` | `oph-fep` | stability is not generic in the coupling |
| `schedule-path` | disagreement | `obstructed` | `oph-fep` | iteration schedules are not interchangeable |
| `reduction-criterion` | disagreement | `obstructed` | `oph-fep` | competing reduction criteria disagree |
| `existence-without-reachability` | disagreement | `obstructed` | `oph-fep` | existence of an equilibrium does not imply reachability |
| `directed-gradient` | disagreement | `obstructed` | `oph-fep` | direct and adjoint gradients differ |
| `four-cycle-local-minimum` | disagreement | `obstructed` | `oph-fep` | the four-cycle has spurious local minima |
| `cuboctahedron-vs-icosahedron` | disagreement | `obstructed` | `synergetics-oph` | twelve vertices is not one solid: tetravolume 20 against $18.5123\ldots$ |
| `shell-vs-free-energy` | disagreement | `obstructed` | `fep-synergetics` | shell two spans three free energies, shell four seven; 36 and 52 are shared |
| `quadray-ivm` | complementarity | `exact` | `synergetics-oph` | IVM coordinates supply the carrier's geometry |
| `equilibrium-propagation` | complementarity | `exact` | `oph-fep` | the energy-based learning identity holds exactly |
| `replicator-potential` | complementarity | `exact` | `oph-fep` | the replicator field is a metric gradient |
| `kissing-twelve` | complementarity | `exact` | `synergetics-oph` | twelve neighbours at squared distance eight |
| `icosahedron-tetravolume` | complementarity | `exact` | `synergetics-oph` | the carrier's solid is $18.5123\ldots$, squared exactly $175 + 75\sqrt{5}$ |
| `packing-density` | complementarity | `exact` | `synergetics-oph` | Voronoi cell $14/12/24$, tetravolume six, density $\pi/(3\sqrt{2})$ |
| `tetravolume-arithmetic` | complementarity | `exact` | `three-way` | four nats, the rational gap, and an integral table beside a non-integral carrier |

: The {{REGISTRY_ENTRIES}} correspondences, in registry order:
{{REGISTRY_ISOMORPHISMS}} isomorphisms, {{REGISTRY_DISAGREEMENTS}}
disagreements and {{REGISTRY_COMPLEMENTARITIES}} complementarities, scoped
`oph-fep` {{REGISTRY_OPH_FEP}}, `synergetics-oph`
{{REGISTRY_SYNERGETICS_OPH}}, `fep-synergetics`
{{REGISTRY_FEP_SYNERGETICS}}, `three-way` {{REGISTRY_THREE_WAY}}. Every row is
executed by `oph-fep-synergetics verify`.
{#tbl:registry}

The registry is the source of truth for the paper. The prose in
[@sec:correspondences], [@sec:disagreements] and [@sec:complementarities] says
what the probes found; it does not add claims the probes do not carry.
[@fig:scope-matrix] shows the same {{REGISTRY_ENTRIES}} rows as a matrix of
scope, kind and status, so that the shape of the registry is visible without
reading the table.

![The registry as a scope matrix. Each row is one correspondence in registry
order; the columns are its pairwise scope, its kind, its status, and whether the
entry carries an executable probe, with cells coloured by category. The figure is
drawn from `analysis/registry.py`, so an entry cannot appear in the matrix
without appearing in `verify`.](../output/figures/fig_scope_matrix.png){#fig:scope-matrix width=90%}

## Reading a source without copying it

The three source programs are read at pinned commits, and the pins are data, not
prose: `data/sources.json` records the commit, the branch, the licence terms and,
for each borrowed identifier, the file, the line range and the symbol. The
command `oph-fep-synergetics sources --vault-root <root>` re-reads the live
checkouts and reports `match`, `moved`, `absent` or `unavailable` for each pin.

Three rules follow from the licences, and they are enforced rather than
recommended.

OPH is dual-licensed: its software, including the Lean development and the model
code, is Apache-2.0, while its papers, book, documentation, figures, data and
ledgers are CC BY-NC-SA 4.0 [@oph]. This package therefore reimplements the
software behaviour with attribution and reads no data or ledger from OPH. The
carrier adjacency is the clearest case: it is *derived* from OPH's twenty
oriented faces and then checked against them, rather than transcribed from any
table, and no OPH data file is read at run time.

FEP_Formal is CC BY 4.0 and is cited by DOI and by topic identifier
[@fepformal]. The catalogue is read through its machine-readable metadata
(`config/catalogue_metadata.yaml`), which is the schema-2 roster sealed at
`fep-001` through `fep-{{FEP_FORMAL_TOPICS}}`; the generated projection
`config/topics.yaml` is never authored or edited, because FEP_Formal declares it
a projection and regenerates it.

QuadMath is Apache-2.0 [@quadmath]. The constants this package needs are
reimplemented and the source is cited by symbol: the neighbour moves, the
embedding, the shell formula and the tetravolume convention. This paper quotes no
QuadMath prose.

Two corrections in the provenance record are worth reporting, because both were
caught by checking rather than by care. An earlier draft of this package
attributed the carrier table to a Lean file that does not contain it; the
attribution was false and the table was replaced by a derivation. And an earlier
draft cited a QuadMath commit that had been superseded; the live checkout had
moved, and the pin was updated. Both are recorded in `docs/PROVENANCE.md`. A pin
is a point-in-time reading, and the discipline that matters is that the reading
is re-done rather than remembered.

## Why exactness, not agreement

It would be easier to report that three literatures "converge". The registry does
not do that, for a reason that [@sec:disagreements] makes concrete: two
programmes can agree on a number while disagreeing on the object the number is
about. Existence without reachability [@prop:existence] is the clearest example.
Both programmes can be right that an equilibrium exists, and both can be right
about its value, while one of them is wrong to infer that a trajectory reaches
it. Agreement on the value would have hidden the disagreement about the dynamics.

The registry therefore distinguishes *what is computed* from *what is claimed
about it*, and the statuses are assigned to claims. This is why the paper reports
more disagreements than a synthesis-oriented reading would, and why the
disagreements carry the same weight as the agreements.
