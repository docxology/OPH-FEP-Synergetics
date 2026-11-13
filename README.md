# OPH-FEP-Synergetics

Executable correspondences, obstructions and complementarities across three
programs that all claim to describe how a bounded part of a world comes to
agree with itself:

- **Observer Patch Holography** (OPH) — an executable observer on a twelve-port
  icosahedral carrier with a native integer repair law.
- **The Free Energy Principle** (FEP) — variational free energy, Bayes-optimal
  perception, and conditional independence at a Markov blanket.
- **Fuller's Synergetics** — quadray coordinates, the isotropic vector matrix,
  the cuboctahedral vector equilibrium, and rational tetravolumes.

The project's purpose is to say precisely which of the three-way identifications
survive an executable test, which are obstructed, and which are complementary
rather than equivalent. Every claim in this repository is a probe that either
passes or fails, and the failing direction is reported as a result rather than
hidden as a gap.

**The paper:** [*Observers, Blankets, and Screens*](output/pdf/OPH-FEP-Synergetics_combined.pdf)
(PDF, rendered from `manuscript/` by `scripts/build.py`; every number in it is
bound to the function that produces it).

## What runs

```
oph-fep-synergetics verify                 # all 42 correspondence probes
oph-fep-synergetics verify --json          # the same, as a machine-readable report
oph-fep-synergetics verify --key pr01-obstruction
oph-fep-synergetics verify --pair synergetics-oph   # one pairwise ledger
oph-fep-synergetics verify --kind disagreement       # one kind of claim
oph-fep-synergetics sources                # pins, licences and reuse policy
oph-fep-synergetics sources --vault-root /path/to/vault   # check live checkouts
oph-fep-synergetics values                 # every number the manuscript may quote
oph-fep-synergetics carrier                # carrier counts, gap, shell sizes
oph-fep-synergetics geometry               # IVM moves, shell counts, volumes
```

`verify` exits non-zero when any probe fails, and exits `2` on a misspelled key
or pair so that a typo cannot pass as an empty run.

## Install

```
uv venv --python 3.15 .venv
uv pip install --python .venv/bin/python -e ".[dev,figures]"
```

Any Python >= 3.12 works; without `uv`, `python3 -m venv .venv &&
.venv/bin/pip install -e ".[dev,figures]"` is equivalent. The only runtime
dependency is `numpy`; `figures` (matplotlib) and `dev` (pytest, pytest-cov) are
optional extras. `.[dev]` alone runs the suite (the figure tests skip without
matplotlib); the `figures` extra is what `scripts/make_figures.py` and the
manuscript pipeline need.

## Layers

| layer | module | what it holds |
| --- | --- | --- |
| `core` | `core/carrier.py` | OPH's twelve-port carrier, its port chart, Laplacian, spectrum and load shell |
| | `core/repair.py` | the native seam-repair law, its conservation law and its exact energy decrease |
| | `core/sheaf.py` | cellular sheaves, agreement energy, harmonic extension, record injectivity, cohomology dimensions |
| | `core/free_energy.py` | the linear-Gaussian free-energy dictionary and the three named surprises |
| | `core/controls.py` | the obstruction and dynamics controls: descent census, D-stability, mediators, orbits, free-energy forms |
| | `core/invariants.py` | the Laplacian's counting identities: spanning trees, reciprocal spectrum, GMRF precision, Markov blankets |
| | `core/symmetry.py` | the carrier's automorphism group, rotations versus reflections, orbit counts, law equivariance |
| | `core/network.py` | the carrier as a resistive network: exact effective resistances, the Kirchhoff index, shells around a port, and the subgraph census |
| `geometry` | `geometry/quadray.py` | quadray coordinates, the A3/IVM embedding and exact tetravolumes |
| | `geometry/lattice.py` | IVM shells, the cuboctahedral vector equilibrium, Voronoi cell and packing density, Goldberg/geodesic counts |
| | `geometry/subdivision.py` | icosahedral subdivision and the dual count identities |
| | `geometry/solids.py` | oriented faces, polyhedron volume, and the icosahedron's tetravolume against its closed form |
| `analysis` | `analysis/registry.py` | the 42 registered correspondences, each with a pairwise scope and an executable probe |
| | `analysis/provenance.py` | source pins, licence terms and reuse policy |
| `manuscript_variables` | `manuscript_variables/values.py` | every number the prose may quote, derived from the package where a derivation exists |

## The registry

42 entries: **17 isomorphisms**, **15 disagreements**, **10 complementarities**,
scoped to the pairs `oph-fep` (24), `synergetics-oph` (14), `fep-synergetics` (2)
and `three-way` (2). All 42 probes pass. `verify --kind` narrows a run to one of
the three kinds and `--pair` to one of the four scopes. The full table, with the OPH / FEP /
Synergetics reading of each entry and the values its probe asserts, is generated
into [docs/CORRESPONDENCES.md](docs/CORRESPONDENCES.md).

The load-bearing results are:

| claim | status | evidence |
| --- | --- | --- |
| OPH's carrier is the kissing-twelve configuration, and its port chart is derived from OPH's own oriented faces | exact | 12 ports, 30 seams, 20 faces, five-regular; the 20 triangles of the derived graph are exactly OPH's 20 oriented faces |
| The carrier gap is the icosahedral gap | exact | `lambda_1 = 5 - sqrt(5)` in both readings (the generator of OPH's mean repair operator, and the FEP's precision), and the normalised gap `lambda_1 / (2|E|)` is exactly `(5 - sqrt(5)) / 60` |
| The native repair law conserves load and lowers `sum(x^2)/2` by `(d^2 - (|d| mod 2)) / 4` | exact | checked on all 4680 shell transitions |
| Patch agreement is the global-section condition of a cellular sheaf, and the sheaf Laplacian is the carrier Laplacian | exact | `carrier_sheaf(carrier).laplacian() == carrier.laplacian()` |
| On the linear-Gaussian domain, constrained relative-entropy minimisation and predictive-coding perception are the same object | exact on that domain | the harmonic extension, a separate KKT constrained minimiser and the perception flow coincide; the free energy computed as energy minus entropy equals the surprisal at the posterior |
| **PR-01**: the native law is not confluent, and accepted repairs need not lower a measure | **obstructed** | 120 strict descents, 600 zero-descent mismatch-one moves, 3960 no-ops out of 4680 transitions; a one-unit walk changes state at constant energy `1/2` |
| **"Free energy" names two different objects** | **obstructed** | the interior-minimised interface energy and the variational free energy at the posterior differ by a model-dependent constant |
| **Eliminating an interior mediator is not free** | **obstructed** | the reduced precision is `(1/3) v v^T + epsilon I`; a sum of pairwise semidefinite terms reproduces it iff `epsilon >= 1/3` |
| **A stationary blanket need not bound the dynamics** | **obstructed** | with a zero precision entry the internal drift is `2q eta + (q - 1) b - 2 mu`, so equal blanket statistics do not fix the internal motion |
| **Agreement can exist and be unreachable** | **obstructed** | the clipped two-node map has the origin as its only equilibrium and a source there; the relaxation orbit stays at radius `>= 0.84` with period `~7.58` |
| **D-stability is not free** | **obstructed** | Lyapunov form `P = (1/3)[[11,4,-8],[4,5,-2],[-8,-2,9]]`, `det P = 115/27`; at `k = 4` the Routh first column is `(1, 6, -2, 24)` and the scaled matrix has two right-half-plane roots |
| **Finite phases are path functions** | **obstructed** | full and half-damped iteration share every answer, yet their centered contrasts after twelve steps differ by exactly `921633999 / 40000000000` |
| Quadray geometry supplies the exact integer side of the correspondence | exact | 12 neighbour moves at squared distance 8, shells `10k^2 + 2`, rational volumes `1, 3, 4, 6, 20, 96` |
| Equilibrium propagation supplies the parameter gradient the repair law does not define | exact | the two-settle contrast is `(theta - target) / (1 + beta)` exactly and tends to the exact gradient |
| The carrier counts its own spanning trees | exact | the reduced Laplacian's determinant is `5184000` — by Kirchhoff's matrix-tree theorem the number of spanning trees in OPH's reading and, as the determinant of the same matrix read as a Gaussian precision with one port pinned, the normalisation of the FEP's linear-Gaussian interface integral — equal to the product of the eleven nonzero eigenvalues over `12` |
| The Gaussian reading makes the seam graph the independence graph | exact | ports are conditionally dependent exactly across seams; pairwise precision `-1` at every seam, each port's Markov blanket is its five neighbours, and the variance of a port's deviation from the carrier mean is exactly `7/36` |
| The agreement sheaf has one global section and nineteen obstructions, and the faces fill them all | exact | on the graph `dim H^0 = 1`, `dim H^1 = 19` (the cycle rank) with Euler characteristic `12 - 30`; with the twenty faces the face coboundary has rank 19, the surface's Betti numbers are `(1, 0, 1)` and its Euler characteristic is `2`, and OPH's oriented faces sum to a fundamental cycle |
| The symmetry group closes the descent census | exact | order 120 = 60 rotations + 60 reflections, transitive on ports/seams/faces; every one of the 1123200 transitions (group element x configuration x oriented seam x coin) commutes with the group, giving 4 configuration orbits and 3 outcome multisets |
| Information cost is half the synergetic metric | exact | at unit covariance in the integer embedding one IVM step costs 4 nats and the sites `k m` cost `4k^2`; the four is the unit of length, the ratios are the content |
| **Twelve vertices is not one solid** | **obstructed** | the neighbour graph is four-regular with tetravolume 20; the carrier is five-regular with tetravolume `18.512295868219...` |
| **The shell index and the free energy determine neither each other** | **obstructed** | the coordination shell is the cuboctahedral gauge `max(|u|_inf, |u|_1 / 2)`, a polyhedral norm, checked against every breadth-first shell through 12; the free energy on shell `k` ranges from `2k^2` (`2k^2 + 2` for odd `k`) to `4k^2`, attained only at the twelve vertices, and the levels 36 and 52 are each shared by two shells |
| Synergetics' tetravolume table is derived, not transcribed | exact | `1, 3, 4, 6, 20, 96` come out of explicit integer vertex sets by exact convex-hull volumes in the rationals |
| The carrier's own solid has an exact algebraic volume | exact | computed in exact field arithmetic as `sqrt 2 (15/2 + (5/2) sqrt 5) = 18.512295868219...`, squared exactly `175 + 75 sqrt 5`; its conjugate is `2.7009...` and the product is exactly `50` |
| Both twelve-around-one arrangements attain the kissing number, differently | exact | squared edge over squared circumradius is exactly `1` for the cuboctahedron (rigid) and `2 - 2 sqrt 5 / 5` for the icosahedron (slack): the carrier is the loose arrangement |
| The twelve-around-one cell is the densest packing's cell | exact | 14 vertices, 12 rhombic faces, 24 edges, tetravolume exactly 6, ball tetravolume exactly `pi sqrt 2`, density exactly `pi / (3 sqrt 2)` — optimal among lattice packings (Gauss) and all packings (Hales) |

## Scope

This package reimplements the *software* it needs and cites everything else. It
does not implement OPH's quantum event algebra or continuum limit, FEP's
nonlinear or non-Gaussian cases, or Fuller's full synergetic apparatus. Four
cautions from the review are load-bearing here and are recorded in
[docs/SCOPE.md](docs/SCOPE.md): disambiguate "free energy"; respect OPH's
premise discipline; do not let numerological agreement stand in for a frozen
prediction; and keep the scale separations explicit.

Two quantities in the source paper are **not independently reproducible** from
its text and are recorded as `not-reproduced` rather than verified: the
`35/264` port-zero occupation and the `13/15` and `1/60` move probabilities,
which need the record/selector protocol specification.

## Provenance and licences

Three source programs are pinned by commit in [data/sources.json](data/sources.json)
and reported by `oph-fep-synergetics sources`. Nothing here is relicensed.

| source | licence | reuse |
| --- | --- | --- |
| [Observer Patch Holography](https://github.com/FloatingPragma/observer-patch-holography) | Apache-2.0 for software, CC BY-NC-SA 4.0 for papers, docs, figures, data and ledgers | reimplement the software with attribution; do not copy prose, and do not reuse the data or ledgers |
| [FEP_Formal](https://github.com/ActiveInferenceInstitute/fep_formal) | CC BY 4.0 | cite the DOI and the exact topic identifiers used |
| [QuadMath](https://github.com/docxology/QuadMath) | Apache-2.0 | reuse with notice and attribution |

The conventions, the embedding, the neighbour moves and the tetra-volume rule
follow QuadMath. The port chart follows OPH's `ORIENTED_BASE_FACES` and is
*derived* from it rather than transcribed, so the derivation is checkable. The
comments in this repository record the commit actually read only in
`data/sources.json`; the package's own `sources` command is the single place a
commit is reported, and a test fails if any prefix of seven or more hex digits
of a pinned commit appears in the package.

### Third-party notice

Reimplemented here, with attribution: OPH's twelve-port carrier and its derived
port chart, the native seam-repair law, and the sheaf, symmetry, network and
controls machinery built on them — OPH's Apache-2.0 *software* scope; the
linear-Gaussian free-energy dictionary of the FEP — cited, CC BY 4.0, DOI
`10.5281/zenodo.19699233`; and QuadMath's quadray conventions, neighbour moves,
shell counts and tetravolume arithmetic — Apache-2.0, DOI
`10.5281/zenodo.16887791`. No source prose, documentation, figure, data file or
ledger is copied or redistributed, and no source endorses any result here. Each
source keeps its own terms, recorded with commit pins in
[data/sources.json](data/sources.json) and restated in
[docs/PROVENANCE.md](docs/PROVENANCE.md).

This project is released under the MIT licence; see [LICENSE](LICENSE).

## The manuscript

The paper, *Observers, Blankets, and Screens*, is a project of the
`docxology/template` research monorepo, which owns the Pandoc and LaTeX
pipeline. The source of truth is `manuscript/`: the render reads `config.yaml`,
`references.bib`, the eleven sections and `output/figures/`, and writes
`output/pdf/OPH-FEP-Synergetics_combined.pdf`.

```
# 1. the twenty-three figures and their accessibility registry
.venv/bin/python scripts/make_figures.py

# 2. the tokens and the resolved tree, with the template checkout's interpreter
cd /path/to/template
TEMPLATE_REPO_ROOT=$PWD .venv/bin/python /path/to/project/scripts/z_generate_manuscript_variables.py

# 3. the PDF (hydration already ran, so the stage's own hydration is skipped)
.venv/bin/python scripts/pipeline/stage_03_render.py \
    --project working/OPH-FEP-Synergetics --skip-manuscript-hydration

# 4. the render contract, read back out of the artefact
.venv/bin/python scripts/check_render.py
```
All four steps run in order with one command:

```
.venv/bin/python scripts/build.py --deliver
```

The template checkout comes from `--template-root`, `$TEMPLATE_ROOT`, or the
sibling `template/` directory beside this project; `--deliver` refreshes the
byte-identical copy in `~/Downloads/` and prints both sha256 digests.

`projects/working/OPH-FEP-Synergetics` is a symlink to this checkout, which the
pipeline resolves to its real path, so the PDF lands in this project's
`output/pdf/`. The stage's own hydration call is skipped because it expects a
project interpreter that passes the template's bytecode audit; step 2 does the
same work with the interpreter that does.

The rendered edition has a deliberate front matter. Page 1 is the cover, drawn
from `output/figures/fig_cover.png`; pages 2 to 5 are the four panels of the
graphical abstract, one panel per page, placed with `\floatplacement{figure}{H}`
so that a panel cannot float away from its caption; the contents page follows the
panels, which is why `preamble.md` keeps the renderer's automatic
`\tableofcontents` from landing between the cover and the panels and re-emits it
from the end of `00_00_graphical_abstract.md`. The margins are 0.5 in, because
this is dense reference matter with wide tables.

Three disciplines hold the manuscript together, and each one is a test:

- **Every number is a token.** The prose carries `{{TOKEN}}`, never a literal,
  and `scripts/z_generate_manuscript_variables.py` fails if a token survives
  hydration. A test asserts that the tokens the prose uses and the values the
  package publishes are the same set, so a number cannot be quoted without being
  bound and a bound value cannot silently fall out of use.
- **Every formalism is numbered by the renderer.** Definitions and propositions
  are fenced `::: {.proposition #prop:key}` blocks referenced as `[@prop:key]`,
  and the template's `formalism.lua` owns the counters, so no statement is
  hand-numbered and the numbering cannot drift between editions.
- **Every citation resolves.** Pandoc's `--natbib` resolves `[@key]` against
  `manuscript/references.bib`. A test asserts that the cited keys and the
  bibliography are the same set: no undefined citation, and no uncited entry.

## Layout

```
data/sources.json                    pinned sources, licences, cited identifiers
docs/CORRESPONDENCES.md              generated: the 42 entries and their probe values
docs/MANUSCRIPT_VALUES.md            generated: every bound value
docs/CONCORDANCE.md                  generated: where each claim is stated and each figure appears
docs/PROVENANCE.md                   the pins and the reuse rules
docs/SCOPE.md                        what is reimplemented, what is cited, what is out
docs/VERIFICATION.md                 the gates and how to run them
manuscript/                          config.yaml, preamble.md, references.bib and eleven sections
scripts/render_docs.py               regenerates the three generated documents
scripts/make_figures.py              renders the twenty-three figures and their accessibility registry
scripts/build.py                     runs the four manuscript steps, ending in the render contract
scripts/z_generate_manuscript_variables.py   resolves every token into output/manuscript/
src/oph_fep_synergetics/             the package
tests/                               the suite (388 tests)
```

## Tests

```
.venv/bin/python -m pytest tests/ -q --cov --cov-report=term-missing
```

388 tests, no mocks, no network. The suite checks exact values rather than
smoke behaviour: the 78-state shell and all 4680 transitions, the exact
rational contrast, the closed-form Routh column, the Lyapunov determinant
`115/27`, the mediator threshold at `epsilon = 1/3`, the IVM round-trip over
shell sites, the resistance ladder's exact fractions against a from-scratch solve,
and the exact volume ratios. Coverage is 96.58% with a 90% gate
enforced by `pyproject.toml`.

## Regenerating the documents

```
.venv/bin/python scripts/render_docs.py           # write
.venv/bin/python scripts/render_docs.py --check   # fail when stale
```

A test runs the `--check` path, so the committed documents cannot drift from the
registry they describe.

## Contributing and citation

Issues and pull requests are welcome on the public remote,
[docxology/OPH-FEP-Synergetics](https://github.com/docxology/OPH-FEP-Synergetics).
A pull request that changes a claim, a number or a figure has to keep the gates
green: [docs/VERIFICATION.md](docs/VERIFICATION.md) lists them, and
[AGENTS.md](AGENTS.md) records the house rules — the two most important being
that manuscript numbers are bound tokens and that the generated documents are
re-rendered, never hand-edited.

To cite this repository use the machine-readable
[CITATION.cff](CITATION.cff), which GitHub renders under "Cite this
repository"; the release history lives in [CHANGELOG.md](CHANGELOG.md).
