# Provenance and reuse

Three programs are read, cited and pinned. Nothing in this repository is
relicensed, and no source prose or data is copied.

## The pins

`data/sources.json` (schema 1) records one entry per source:

- `id`, `name`, `url`
- `checkout` — the path relative to the vault root where the source was read; the path is vault-relative and meaningless outside a vault root
- `commit`, `branch` — the exact reading
- `license` — an SPDX expression, with `split` terms when a repository licenses
  different scopes differently
- `reuse_policy` — what this project may and may not take
- `cited_identifiers` — every symbol this package relies on, as
  `kind` / `path` / `lines` / `symbol` / `role`

| id | source | licence |
| --- | --- | --- |
| `oph` | FloatingPragma's Observer Patch Holography | `Apache-2.0 AND CC-BY-NC-SA-4.0` |
| `fep_formal` | FEP_Formal | `CC-BY-4.0` |
| `quadmath` | QuadMath | `Apache-2.0` |

OPH's licence is a genuine split and is recorded as one:

- **Apache-2.0** covers the software, including the Lean development and the
  model code, evidenced at `code/LICENSE`.
- **CC BY-NC-SA 4.0** covers papers, documentation, figures, data and ledgers,
  evidenced at `LICENSE`.

That means: reimplement the software with attribution, do not copy paper, book,
documentation or figure text, and do not reuse the data or the ledgers.

FEP_Formal is CC BY 4.0: attribution is required, and the DOI and the exact
topic identifiers used must be cited.

QuadMath is Apache-2.0: reuse with notice and attribution is permitted. This
package nonetheless reimplements the constants it needs rather than vendoring
the source.

## Checking the pins

```
oph-fep-synergetics sources
oph-fep-synergetics sources --vault-root /path/to/vault
oph-fep-synergetics sources --json
```

With no `--vault-root` the command reports the pins and marks every check
`unavailable`: no location is guessed. A vault root is an ordinary local
directory beneath which the pinned checkouts live; with one given, each pin is
compared against `git -C <checkout> rev-parse HEAD` and reported as:

| status | meaning |
| --- | --- |
| `match` | the checkout is at the pinned commit |
| `moved` | the checkout resolves to a different commit |
| `absent` | the checkout has no resolvable `HEAD` |
| `unavailable` | no vault root was given, or no checkout was consulted |

**A moved pin is a reading, not a failure.** `check_pins` reports it as `moved`
and `verify` does not treat it as a broken claim. A pin records what was read on
a date; the source is free to move on.

Nothing in the check mutates the sources: it is one `rev-parse` per checkout,
with the repository's `stdin` closed.

## How a cited identifier is used

Each `cited_identifiers` entry names where in the source a convention came from
and what role it plays here. Three examples, one per source:

- `oph` / `model/observer/native.py:45-52` / `native` — the native seam law this
  package reimplements on exact integers.
- `oph` / `.../self_readback_repair_closure.py:93` / `ORIENTED_BASE_FACES` — the
  twenty oriented faces from which the carrier adjacency is *derived*, with the
  derivation asserted against its own consequence (the derived graph's triangles
  must be exactly those faces).
- `quadmath` / `src/quadmath/core/quadray.py:175` / `IVM_NEIGHBOR_MOVES` — the
  twelve neighbour moves, reconstructed as the sorted distinct permutations of
  `(2, 1, 1, 0)` and asserted to have squared distance eight.

The commit actually read is deliberately **not** restated in code or prose. It
lives in `data/sources.json`, is printed by `sources`, and a test fails if any
prefix of seven or more hex digits of a pinned commit appears anywhere in the
package. That prevents the drift that a hand-copied hash, full or abbreviated,
invites. Where a source content-addresses a file itself, as OPH does for the
oriented faces in `evidence/native_geometric_source_20260925/spec.json`, the
code cites that address rather than a commit.

## Re-reading the pins

The pins are re-read against current upstream rather than remembered. Each
re-pin resolves the upstream `main` with `git ls-remote`, brings the local
checkout to it, and re-reads every cited identifier at the new commit: the file,
the line range, the symbol and the paraphrased role. A pin moves only when every
identifier still says what the package relies on, and the line ranges move with
it. The pins in `data/sources.json` were last re-read this way on 2026-11-13; at
that reading every OPH file the package cites was unchanged except
`model/PREMISES.md`, whose PR-01 row was byte-identical, and the Lean citation
was corrected from three comment lines to the theorem itself
(`Primitives.lean:1198`).

## What this package does not take

- No source prose, in any form.
- No source data or ledgers, including OPH's retained evidence archives.
- No vendored source files; the geometry conventions are reimplemented.
- No claim that a source endorses any result here.

## Unpublished and private sources

Besides the three pinned repositories, the manuscript cites three unpublished
texts: the author's own working dossier and integration proposal
(`dossier2026`, `proposal2026`), and one third-party private draft,
*Agreement and Surprise: Local Repair, Dynamics and the Limits of Equilibrium*
(Bernhard Mueller, Pragma Research, private draft, 22 October 2026). The draft
is cited the way its source states its claims: none of its prose is quoted, it
is not redistributed here, and the two quantities this package could not
reproduce from its text are recorded `not-reproduced`
([SCOPE.md](SCOPE.md)) beside the package's own measurements.

## The render pipeline

The manuscript renders through the `docxology/template` research monorepo
(Pandoc, LaTeX, the `formalism.lua` filter). This repository vendors none of
that infrastructure and refers to it only by template-relative paths; the
probes, the suite and the CLI need nothing from it.

## Licence of this project

MIT, with the third-party citations above retained. See [LICENSE](../LICENSE).
The MIT licence covers the code and the manuscript prose in `manuscript/`; the third-party
cited works keep their own terms.

## Editorial log

Three provenance corrections are recorded here because they changed what the
package claims:

1. An earlier draft of `core/carrier.py` described its adjacency table as
   transcribed from a Lean file that does not contain it. The adjacency is now
   derived from OPH's own oriented faces, and a probe asserts the derivation
   against the faces rather than against a prose claim.
2. The geometry docstrings carried a stale QuadMath commit. The commit is now
   named in one place only, and a test enforces that.
3. Two bibliography identifiers were wrong when checked against Crossref (the
   DOI of Biehl, Pollock and Kanai's critique, and the DOI and year of Da Costa
   et al.'s discrete-state synthesis); both were corrected, and every DOI in
   `manuscript/references.bib` now resolves to the metadata the entry states.
