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

Each of the 42 entries carries an executable probe. A probe returns `(passed,
values)`; `values` is what the pass is made of, so a number that changes without
changing the verdict is still visible in the report.

Exit codes:

| code | meaning |
| --- | --- |
| `0` | every selected probe passed |
| `1` | at least one selected probe failed |
| `2` | an unknown key, pair or kind was requested |

An unknown key is an error rather than an empty selection, so a typo cannot be
mistaken for a clean run.

## The suite

```
.venv/bin/python -m pytest tests/ -q
.venv/bin/python -m pytest tests/ -q --cov --cov-report=term-missing
```

345 tests, no mocks, no network, no fixtures that fake a result. The suite runs
from a fresh process against the real code paths.

| file | what it pins down |
| --- | --- |
| `test_carrier.py` | 12/30/20, five-regular, the derived port chart, the exact spectrum and gap, the 78/66 shell, the isomorphism witness, and the rejection cases |
| `test_repair.py` | conservation and the energy identity over all 4680 transitions, the odd-unit coin, the argument guards, the gated schedule reaching rest, and an ungated move that lowers nothing |
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
| `test_registry.py` | the registry contract, the pairwise scopes, the selection commands, the status/kind vocabulary invariant, plus an exact-value assertion per probe |
| `test_provenance.py` | the pin file, the licence split, identifier lookup, and `match`/`moved`/`absent`/`unavailable` against a real temporary repository |
| `test_manuscript_variables.py` | every bound value, cross-checked against the function that produces it |
| `test_cli.py` | every subcommand, the JSON payloads, the exit codes, and the console-script declaration |
| `test_public_api.py` | that every exported name resolves, and that no commit hash escapes the pin file |
| `test_docs.py` | that the three generated documents are current, that the outline and the sections name every registered claim, and that the guides' relative links resolve |
| `test_manuscript.py` | the eleven section files the renderer discovers, every token resolving to a bound value and every bound value being quoted, citations equal to the bibliography in both directions, a `{#sec:…}` label on every top-level heading, the manuscript's registry table against the registry, and the figure links against the accessibility registry |
| `test_figures.py` | that every figure builder exists and is orchestrated, and that a sample renders non-trivial PNGs plus a well-formed registry (skipped without the `figures` extra) |

The provenance test builds its own temporary repository by writing a commit
object directly with `hash-object`; it reads and writes no configured Git
identity.

## Coverage

`pyproject.toml` enforces `fail_under = 90` with branch coverage over
`src/oph_fep_synergetics`. The current result is **96.34%**. The remaining lines
are defensive or unreachable branches (an assertion that cannot fire on a
well-formed carrier, a guard behind a guard), not untested behaviour.

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
| `test_render_contract.py` | a render whose cover, graphical-abstract panels or contents page are not on the pages the manuscript promises, a panel page that does not carry its own figure caption, a measured left text margin wider than the 0.5 in the manuscript claims, a registry key that does not reach the printed text, or a log with an overfull box, an undefined reference or an unresolved token |

The rendered result was checked the way a reader checks it, on the artefact
itself: 43 pages, the cover on page 1, the four graphical-abstract panels on
pages 2 to 5, the contents page at 6, and a 0.500 in left text margin measured
from the PDF's own text bounding boxes. Each panel page carries its own figure
caption, so no panel floated away from its figure. The body carries nineteen
numbered figures plus the cover image, fourteen numbered equations, three
definitions, thirty-two numbered propositions and eight tables, with `0`
undefined citations and `0` unresolved `??` cross-references in
`_combined_manuscript.tex`, `0` overfull boxes in `_combined_manuscript.log`, and
all 42 registry keys present and unbroken in the extracted text of the PDF. `scripts/check_render.py` performs
that reading and exits non-zero when any of it drifts; `test_render_contract.py`
runs it against the rendered PDF when one is present.

## What is deliberately not verified

The `35/264` port-zero occupation and the `13/15` and `1/60` move probabilities
from *Agreement and Surprise* need the record/selector protocol specification
that the paper does not give. They are recorded as `not-reproduced` in
[docs/SCOPE.md](SCOPE.md) and are not asserted anywhere in this repository.
