# Method: a registry that executes {#sec:method}

## The unit of work

The paper's unit is a **correspondence**: a claim about the same object made in
two or three vocabularies, together with executable evidence.

::: {#def:correspondence .definition}
A **correspondence** is a record $c = (\kappa, \sigma, \pi)$ where $\kappa$ names
the claim in the OPH, FEP and Synergetics vocabularies, $\sigma$ is a status from
the fixed vocabulary below, and $\pi$ is a *probe*: a total function from no
arguments to a mapping from names to values. The probe **passes** when every
value it reports satisfies the assertion the claim makes.
:::

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

## The registry

[@tbl:registry] is the registry, reproduced from
`src/oph_fep_synergetics/analysis/registry.py`. The grouping is the paper's
argument in miniature: eight claims survive as identifications, eleven fail as
identifications while remaining informative, and four are complements in which
one program supplies what another lacks without either reducing to the other.

| Key | Kind | Status | What it tests |
|:-----------------------------------|:----------------|:-------------------|:-------------------------------|
| `carrier-port-chart` | isomorphism | `exact` | derived adjacency reproduces the twenty oriented faces |
| `spectral-gap` | isomorphism | `exact` | $\lambda_1 = 5-\sqrt{5}$ and normalised gap $(5-\sqrt{5})/60$ |
| `shell-counts` | isomorphism | `exact` | $\binom{13}{11}$ shell size and the binary sub-shell |
| `native-repair-law` | isomorphism | `exact` | OPH's integer rule reimplemented, with its closed-form energy drop |
| `variational-equivalence` | isomorphism | `exact-on-domain` | harmonic-extension residual equals a Schur complement |
| `free-energy-dictionary` | isomorphism | `exact-on-domain` | numerical, Bayesian and surprisal free energy coincide to curvature |
| `sheaf-agreement` | isomorphism | `exact` | the cellular sheaf Laplacian is the carrier Laplacian |
| `certificate-species` | isomorphism | `exact` | the certificate species of the shell and its spectral gap |
| `pr01-obstruction` | disagreement | `obstructed` | the native law is not confluent |
| `free-energy-name-collision` | disagreement | `obstructed` | "free energy" names inequivalent functionals |
| `mediator-obstruction` | disagreement | `obstructed` | pairwise boundary reduction needs a threshold |
| `stationary-blanket` | disagreement | `obstructed` | a stationary blanket is not an attractor |
| `solenoidal-freedom` | disagreement | `obstructed` | stationary blankets leave entropy production free |
| `d-stability` | disagreement | `obstructed` | stability is not generic in the coupling |
| `schedule-path` | disagreement | `obstructed` | iteration schedules are not interchangeable |
| `reduction-criterion` | disagreement | `obstructed` | competing reduction criteria disagree |
| `existence-without-reachability` | disagreement | `obstructed` | existence of an equilibrium does not imply reachability |
| `directed-gradient` | disagreement | `obstructed` | direct and adjoint gradients differ |
| `four-cycle-local-minimum` | disagreement | `obstructed` | the four-cycle has spurious local minima |
| `quadray-ivm` | complementarity | `exact` | IVM coordinates supply the carrier's geometry |
| `equilibrium-propagation` | complementarity | `exact` | the energy-based learning identity holds exactly |
| `replicator-potential` | complementarity | `exact` | the replicator field is a metric gradient |
| `kissing-twelve` | complementarity | `exact` | twelve neighbours at squared distance eight |

: The twenty-three correspondences: {{REGISTRY_ISOMORPHISMS}} isomorphisms,
{{REGISTRY_DISAGREEMENTS}} disagreements, {{REGISTRY_COMPLEMENTARITIES}}
complementarities. Every row is executed by `oph-fep-synergetics verify`.
{#tbl:registry}

The registry is the source of truth for the paper. The prose in
[@sec:correspondences], [@sec:disagreements] and [@sec:complementarities] says
what the probes found; it does not add claims the probes do not carry.

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
