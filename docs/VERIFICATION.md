# Verification

Run the registry probes to reproduce the scientific results. The test suite
checks that the code, documentation and rendered PDF agree with their stated
contracts.

## The registry

```
.venv/bin/oph-fep-synergetics verify
.venv/bin/oph-fep-synergetics verify --json
.venv/bin/oph-fep-synergetics verify --key <key> [--key <key> ...]
.venv/bin/oph-fep-synergetics verify --pair <pair>
.venv/bin/oph-fep-synergetics verify --kind <kind>
```

Each of the 75 entries carries an executable probe. Its return value
`(passed, values)` records both the verdict and the quantities checked.

Exit codes:

| code | meaning |
| --- | --- |
| `0` | every selected probe passed |
| `1` | at least one selected probe failed |
| `2` | a filter was unknown or the filters selected no probes |

When combining filters, the probes must match every specified filter.

## The suite

```
.venv/bin/python -m pytest tests/ -q
.venv/bin/python -m pytest tests/ -q --cov --cov-report=term-missing
```

460 tests. The suite runs offline from a fresh process and exercises the actual
scientific code paths.

| file | what it pins down |
| --- | --- |
| `test_carrier.py` | 12/30/20, five-regular, the derived port chart, the exact spectrum and gap, the 78/66 shell, the isomorphism witness, and the rejection cases |
| `test_repair.py` | conservation and the energy identity over all 4680 transitions, the odd-unit coin, the argument guards, the gated schedule reaching rest, and an ungated move that lowers nothing |
| `test_native_dynamics.py` | ungated confluence and nontermination, gated normal forms, binary-chain stationarity, gap and ensemble relaxation, equality of the conditional-variance forms, first binary entry, and the distance and contrast-variance quotients |
| `test_sheaf.py` | laplacians, coboundaries, global sections, harmonic extension (including the fully-clamped case), record injectivity, diffusion, and every constructor guard |
| `test_free_energy.py` | the KL closed forms, the Schur complement, the update-residual form of Bayesian surprise, the variational bound, rooted-forest enumeration for the record layouts, greedy clamping and the interior-free model |
| `test_planar_duality.py` | complementary Gaussian edge covariances from centred pseudoinverses and the dual graphs against exact geometric incidence |
| `test_weighted_duality.py` | all exact entries under reciprocal conductances, a copied-weight defect, and a separate weighted projector and resistance check |
| `test_quadratic_modes.py` | exact quadratic frame and precision factors, checked independently by rotating and evaluating explicit forms |
| `test_geometric_screening.py` | exact central-plane and separator counts, and residual conditional information checked through covariance conditioning |
| `test_centroid_dynamics.py` | full lag-covariance match and three-time event probabilities for the stationary centroid |
| `test_contrast_quotients.py` | exact four-cube quotient, complete bipartite resistance collisions and a periodic FCC Fourier check |
| `test_root_shells.py` | quadratic precision formula across dimensions and exact lag-two covariance leakage from one strengthened seam |
| `test_surface_duality.py` | toroidal harmonic periods and exact ranks under reciprocal and copied conductances |
| `test_constraint_screening.py` | exact FCC plane cuts, joint Gaussian hard and noisy sum conditioning, and a general weighted-graph witness |
| `test_jitterbug.py` | rigid triangle lengths, supporting hulls, endpoint graph isomorphism and octahedral bar coincidence |
| `test_jitterbug_inference.py` | full weighted Laplacian invariants, exact coordinate projection and independent binary-kernel covariance |
| `test_controls.py` | the descent census, the Routh column in closed form, the Lyapunov pair and its determinant, the mediator threshold, the affine contrast as an exact rational, the orbits, the four-cycle trap, and the replicator potential |
| `test_geometry.py` | the twelve moves and their metric, the gauge, the exact xyz round-trip over shell sites, the shell sequence, the exact volumes, the vector equilibrium, and the subdivision identities |
| `test_geometry_guards.py` | the icosahedral adjacency matrix and every range guard |
| `test_invariants.py` | the reduced-Laplacian determinant in exact rationals, the spanning-tree number against the spectral product, the reciprocal spectrum sum, the GMRF precision and its independence pattern, Markov blankets, and the guesses the invariants reject |
| `test_network.py` | the Laplacian, distances and intersection array, exact pair resistances and their Kirchhoff index, the spectral chord comparison, the subgraph census, and isomorphism under relabelling |
| `test_lattice.py` | the cube's line graph and its isomorphism with the cuboctahedron, and the FCC site and interstice classes against a Cartesian nearest-site enumeration |
| `test_symmetry.py` | the automorphism group by exhaustive search, closure, orbit counts, the rotation/reflection split from the face orientations, eigenspace irreducibility, and the equivariance of the descent law |
| `test_solids.py` | oriented faces, the icosahedron's volume against its closed form and its exact squared ratio, and the rejection of an unoriented face list |
| `test_exact.py` | arithmetic and exact signs in `Q(sqrt 5)`, exact convex hulls and their facet counts, the six-entry tetravolume table derived from vertex sets, the icosahedron's tetravolume as `sqrt 2` times an exact element with norm 50, and the two kissing arrangements compared exactly |
| `test_independent_routes.py` | the second route to every claim that used to be checked against itself: surface cohomology beside graph cohomology, the cuboctahedral gauge beside the breadth-first shells, the exact fixed-point set of the clipped map, the KKT minimiser beside the harmonic extension, energy minus entropy beside surprisal plus divergence, transition-level equivariance, the replicator's asymmetric negative control, and the reciprocal spectrum sum from an exact inverse |
| `test_registry.py` | the registry contract, the pairwise scopes, the selection commands, the status/kind vocabulary invariant, plus an exact-value assertion per probe |
| `test_provenance.py` | the pin file, the licence split, identifier lookup, and `match`/`moved`/`absent`/`unavailable` against a real temporary repository |
| `test_manuscript_variables.py` | every bound value, cross-checked against the function that produces it |
| `test_cli.py` | every subcommand, the JSON payloads, the exit codes, and the console-script declaration |
| `test_public_api.py` | that every exported name resolves, and that no prefix of a pinned commit appears anywhere in the package |
| `test_docs.py` | current generated documents and ontology supplement, registry mentions in the main prose, and working relative links in the guides |
| `test_manuscript.py` | the authored and generated sections the renderer discovers, token and citation coverage, an explicit label on each top-level heading, registry table order, and correspondence between figure links and the accessibility registry when artefacts exist |
| `test_figures.py` | that every figure builder exists and is orchestrated, and that a sample renders non-trivial PNGs plus a well-formed registry (skipped without the `figures` extra) |
| `test_render_contract.py` | render diagnostics and the front matter, author attribution and ontology supplement in the tracked PDF |

The provenance test builds its own temporary repository by writing a commit
object directly with `hash-object`; it reads and writes no configured Git
identity.

## Coverage

`pyproject.toml` enforces `fail_under = 90` with branch coverage over
`src/oph_fep_synergetics`. The command above reports the current result and any
uncovered statements or branches. GitHub Actions runs the same gate on the
oldest supported Python and on Python 3.14. It installs the runtime, test,
and figure dependencies from `uv.lock` and checks the packaged source register
from outside the checkout.

```sh
uv lock --check
uv sync --locked --extra dev --extra figures
```

## Generated documents

```
.venv/bin/python scripts/render_docs.py
.venv/bin/python scripts/render_docs.py --check
```

`docs/CORRESPONDENCES.md`, `docs/MANUSCRIPT_VALUES.md` and `docs/CONCORDANCE.md`
are rendered from the registry, the manuscript value bindings and the text.
The ontology supplement is generated from `data/ontology_concordance.json`.
Its probe and source references are checked against the registry and source
register. `test_docs.py` checks that generated content is current.

The registry document rounds measured values for stable rendering across
numerical runtimes. The JSON verifier keeps the raw values and the probes keep
their original pass/fail tolerances.

## Live source pins

```
.venv/bin/oph-fep-synergetics sources --vault-root /path/to/vault
```

This command checks the source checkouts visible at the supplied vault root.
`test_cli.py` uses `OPH_FEP_VAULT_ROOT` when set and otherwise looks for the
author's local vault layout. If that root is unavailable, the live test skips.
A temporary-directory test independently exercises the `absent` result.

A `moved` result is reported as `moved`. It is a reading, not a failure.

## The manuscript gates

The PDF is produced by the `docxology/template` pipeline; the commands are in
[README.md](../README.md#the-manuscript), and `scripts/build.py` runs them in
order, ending in the contract check. The manuscript's own gates run in the
suite rather than in LaTeX, so a defect is reported by `pytest`:

| gate | what it rejects |
| --- | --- |
| `test_every_token_resolves_to_a_bound_value` | a `{{TOKEN}}` no `BoundValue` publishes (the hydration script also exits 1) |
| `test_the_published_value_set_is_used_except_where_documented` | a bound value the prose never quotes; the documented exception set is asserted to be exact |
| `test_every_citation_resolves_and_every_entry_is_cited` | an undefined citation (also a hard render failure) and an uncited bibliography entry |
| `test_every_top_level_heading_carries_a_section_label` | a heading the renderer would annotate, which it writes back into the source tree |
| `test_the_registry_table_lists_exactly_the_registered_keys` | a manuscript registry table that drifts from `analysis/registry.py` |
| `test_figure_links_match_the_accessibility_registry` | the wrong number of manuscript figure links or, after generating the figures, a missing registry entry or PNG |
| `test_figures.py` | a builder that is not orchestrated, or a registry entry without alt text |
| `test_render_contract.py` | PDF front matter, panel captions, registry keys, printed reference placeholders and the 0.6 in left-margin limit, plus LaTeX log diagnostics |

The [tracked PDF](../output/pdf/OPH-FEP-Synergetics_combined.pdf) places the
cover on page 1, the four graphical-abstract panels and their captions on
pages 2 to 5, and the contents on page 6. Its measured left text margin is
0.500 in. The render check reports the current page and caption counts and
looks for every live registry key in the extracted text.

The suite checks the tracked PDF on a fresh checkout. A local
`scripts/build.py` run also checks the LaTeX log for unresolved tokens,
overfull boxes and undefined citations or references. Rasterise the pages
that contain changed figures at page size and inspect text, markers and
legends for overlap. The figure registry binds each generated PNG to an
alternative description for the tagged PDF.

## Tolerances and platforms

The numerical contract behind the numbers is exactness by construction. Every counting claim
the registry calls exact is computed in Python integers and `fractions.Fraction`, and every
volume claim in `fractions.Fraction` or in the field `Q(sqrt 5)` (`geometry.exact.Surd5`), so a
pass there is not a tolerance but an equality. The float claims are numpy float64 values reported as
measured, alongside the closed form when one is available. Each probe sets its
own tolerance for a numerical comparison. The suite runs on arm64 macOS in
development and on Linux in GitHub Actions. Last-bit floating-point values can
differ across numerical libraries without changing a probe verdict.

## The three source quantities

The `13/15` and `1/60` transition probabilities and the `35/264` port-zero
occupation from *Agreement and Surprise* are reproduced exactly. The protocols
are in OPH's public model: the one-attempt state kernel in
`model/common_source/NATIVE_FIELDS.md` and the two-event observer in
`model/observer/native.py`. `pr01-obstruction` enumerates the sixty draws out of
`(1, 1, 0, ..., 0)` and checks them against a seam count; `reduction-criterion`
enumerates the 264 paths and checks them against a case analysis.
[docs/SCOPE.md](SCOPE.md) states both protocols.
