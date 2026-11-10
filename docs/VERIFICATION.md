# Verification

Everything this repository claims is checked by running it. There are two
surfaces: the registry (`oph-fep-synergetics verify`) and the test suite.

## The registry

```
.venv/bin/oph-fep-synergetics verify
.venv/bin/oph-fep-synergetics verify --json
.venv/bin/oph-fep-synergetics verify --key <key> [--key <key> ...]
```

Each of the 23 entries carries an executable probe. A probe returns `(passed,
values)`; `values` is what the pass is made of, so a number that changes without
changing the verdict is still visible in the report.

Exit codes:

| code | meaning |
| --- | --- |
| `0` | every selected probe passed |
| `1` | at least one selected probe failed |
| `2` | an unknown key was requested |

An unknown key is an error rather than an empty selection, so a typo cannot be
mistaken for a clean run.

## The suite

```
.venv/bin/python -m pytest tests/ -q
.venv/bin/python -m pytest tests/ -q --cov --cov-report=term-missing
```

234 tests, no mocks, no network, no fixtures that fake a result. The suite runs
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
| `test_registry.py` | the registry contract plus an exact-value assertion per probe |
| `test_provenance.py` | the pin file, the licence split, identifier lookup, and `match`/`moved`/`absent`/`unavailable` against a real temporary repository |
| `test_manuscript_variables.py` | every bound value, cross-checked against the function that produces it |
| `test_cli.py` | every subcommand, the JSON payloads, the exit codes, and the console-script declaration |
| `test_public_api.py` | that every exported name resolves, and that no commit hash escapes the pin file |
| `test_docs.py` | that the two generated documents are current and that the guides' relative links resolve |

The provenance test builds its own temporary repository by writing a commit
object directly with `hash-object`; it reads and writes no configured Git
identity.

## Coverage

`pyproject.toml` enforces `fail_under = 90` with branch coverage over
`src/oph_fep_synergetics`. The current result is **96.15%**. The remaining lines
are defensive or unreachable branches (an assertion that cannot fire on a
well-formed carrier, a guard behind a guard), not untested behaviour.

## Generated documents

```
.venv/bin/python scripts/render_docs.py
.venv/bin/python scripts/render_docs.py --check
```

`docs/CORRESPONDENCES.md` and `docs/MANUSCRIPT_VALUES.md` are rendered from the
registry and from `MANUSCRIPT_VALUES`. `test_docs.py` runs the same render and
compares byte-for-byte, so a claim or a value that changes without a re-render
fails the suite.

## Live source pins

```
.venv/bin/oph-fep-synergetics sources --vault-root /path/to/vault
```

This is a live check against the real checkouts and is therefore
environment-dependent. It is part of the acceptance run, not part of the
portable test suite: `test_cli.py` skips the live check when the vault is not
mounted, and a temporary-directory variant asserts the `absent` path instead.

A `moved` result is reported as `moved`. It is a reading, not a failure.

## What is deliberately not verified

The `35/264` port-zero occupation and the `13/15` and `1/60` move probabilities
from *Agreement and Surprise* need the record/selector protocol specification
that the paper does not give. They are recorded as `not-reproduced` in
[docs/SCOPE.md](SCOPE.md) and are not asserted anywhere in this repository.
