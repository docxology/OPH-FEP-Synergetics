# Verification

Everything this repository claims is checked by running it. There are two
surfaces: the registry (`oph-fep-synergetics verify`) and the test suite.

## The registry

```
.venv/bin/oph-fep-synergetics verify
.venv/bin/oph-fep-synergetics verify --json
.venv/bin/oph-fep-synergetics verify --key <key> [--key <key> ...]
.venv/bin/oph-fep-synergetics verify --pair <pair>
.venv/bin/oph-fep-synergetics verify --kind <kind>
```

Each of the 44 entries carries an executable probe. A probe returns `(passed,
values)`; `values` is what the pass is made of, so a number that changes without
changing the verdict is still visible in the report.

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

420 tests. The suite runs offline from a fresh process and exercises the actual
scientific code paths.

| file | what it pins down |
| --- | --- |
| `test_carrier.py` | 12/30/20, five-regular, the derived port chart, the exact spectrum and gap, the 78/66 shell, the isomorphism witness, and the rejection cases |
| `test_repair.py` | conservation and the energy identity over all 4680 transitions, the odd-unit coin, the argument guards, the gated schedule reaching rest, and an ungated move that lowers nothing |
| `test_native_dynamics.py` | ungated confluence and nontermination, the gated normal forms, and binary-chain stationarity, gap and ensemble relaxation |
| `test_sheaf.py` | laplacians, coboundaries, global sections, harmonic extension (including the fully-clamped case), record injectivity, diffusion, and every constructor guard |
| `test_free_energy.py` | the KL closed forms, the Schur complement, the update-residual form of Bayesian surprise, the variational bound attained at the posterior, and the interior-free model |
| `test_controls.py` | the descent census, the Routh column in closed form, the Lyapunov pair and its determinant, the mediator threshold, the affine contrast as an exact rational, the orbits, the four-cycle trap, and the replicator potential |
| `test_geometry.py` | the twelve moves and their metric, the gauge, the exact xyz round-trip over shell sites, the shell sequence, the exact volumes, the vector equilibrium, and the subdivision identities |
| `test_geometry_guards.py` | the icosahedral adjacency matrix and every range guard |
| `test_invariants.py` | the reduced-Laplacian determinant in exact rationals, the spanning-tree number against the spectral product, the reciprocal spectrum sum, the GMRF precision and its independence pattern, Markov blankets, and the guesses the invariants reject |
| `test_network.py` | the graph layer: the Laplacian's rows and degrees, distances and the intersection array, the exact effective resistance of every port pair as a rational, the ladder and its Kirchhoff index, the subgraph census, and an isomorphism test that is invariant under relabelling |
| `test_lattice.py` | the cube's corners and its line graph, and the isomorphism between that line graph and the cuboctahedron |
| `test_symmetry.py` | the automorphism group by exhaustive search, closure, orbit counts, the rotation/reflection split from the face orientations, eigenspace irreducibility, and the equivariance of the descent law |
| `test_solids.py` | oriented faces, the icosahedron's volume against its closed form and its exact squared ratio, and the rejection of an unoriented face list |
| `test_exact.py` | arithmetic and exact signs in `Q(sqrt 5)`, exact convex hulls and their facet counts, the six-entry tetravolume table derived from vertex sets, the icosahedron's tetravolume as `sqrt 2` times an exact element with norm 50, and the two kissing arrangements compared exactly |
| `test_independent_routes.py` | the second route to every claim that used to be checked against itself: surface cohomology beside graph cohomology, the cuboctahedral gauge beside the breadth-first shells, the exact fixed-point set of the clipped map, the KKT minimiser beside the harmonic extension, energy minus entropy beside surprisal plus divergence, transition-level equivariance, the replicator's asymmetric negative control, and the reciprocal spectrum sum from an exact inverse |
| `test_registry.py` | the registry contract, the pairwise scopes, the selection commands, the status/kind vocabulary invariant, plus an exact-value assertion per probe |
| `test_provenance.py` | the pin file, the licence split, identifier lookup, and `match`/`moved`/`absent`/`unavailable` against a real temporary repository |
| `test_manuscript_variables.py` | every bound value, cross-checked against the function that produces it |
| `test_cli.py` | every subcommand, the JSON payloads, the exit codes, and the console-script declaration |
| `test_public_api.py` | that every exported name resolves, and that no prefix of a pinned commit appears anywhere in the package |
| `test_docs.py` | that the three generated documents are current, that the outline and the sections name every registered claim, and that the guides' relative links resolve |
| `test_manuscript.py` | the eleven section files the renderer discovers, every token resolving to a bound value and every bound value being quoted, citations equal to the bibliography in both directions, a `{#sec:…}` label on every top-level heading, the manuscript's registry table against the registry, and the figure links against the accessibility registry |
| `test_figures.py` | that every figure builder exists and is orchestrated, and that a sample renders non-trivial PNGs plus a well-formed registry (skipped without the `figures` extra) |

The provenance test builds its own temporary repository by writing a commit
object directly with `hash-object`; it reads and writes no configured Git
identity.

## Coverage

`pyproject.toml` enforces `fail_under = 90` with branch coverage over
`src/oph_fep_synergetics`. The command above reports the current result and any
uncovered statements or branches. GitHub Actions runs the same gate on the
oldest supported Python and on Python 3.14.

## Generated documents

```
.venv/bin/python scripts/render_docs.py
.venv/bin/python scripts/render_docs.py --check
```

`docs/CORRESPONDENCES.md`, `docs/MANUSCRIPT_VALUES.md` and `docs/CONCORDANCE.md`
are rendered from the registry, from `MANUSCRIPT_VALUES` and from the manuscript
text itself. `test_docs.py` runs the same render and compares byte-for-byte, so
a claim, a value or a cross-reference that changes without a re-render fails
the suite.

The registry document rounds measured values for stable rendering across
numerical runtimes. The JSON verifier keeps the raw values and the probes keep
their original pass/fail tolerances.

## Live source pins

```
.venv/bin/oph-fep-synergetics sources --vault-root /path/to/vault
```

This is a live check against the real checkouts and is therefore
environment-dependent. It is part of the acceptance run, not part of the
portable test suite: `test_cli.py` skips the live check when the vault is not
mounted, and a temporary-directory variant asserts the `absent` path instead.
The live-check test reads its root from the `OPH_FEP_VAULT_ROOT` environment
variable when set, falling back to the author's local checkout layout; without
a mounted root the test skips.

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
| `test_figure_links_match_the_accessibility_registry` | a figure link without a registry entry, or a missing PNG |
| `test_figures.py` | a builder that is not orchestrated, or a registry entry without alt text |
| `test_render_contract.py` | PDF front matter, panel captions, registry keys, printed reference placeholders and the 0.6 in left-margin limit, plus LaTeX log diagnostics |

The [tracked PDF](../output/pdf/OPH-FEP-Synergetics_combined.pdf) places the
cover on page 1, the four graphical-abstract panels and their captions on
pages 2 to 5, and the contents on page 6. Its measured left text margin is
0.500 in. The render check reports the current page and caption counts and
looks for every live registry key in the extracted text.

The suite checks the tracked PDF on a fresh checkout. To also check the
LaTeX log for overfull boxes, undefined citations or references, and unresolved
tokens, run `scripts/build.py`. It runs the renderer and then checks the PDF
against the log from that build.

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
