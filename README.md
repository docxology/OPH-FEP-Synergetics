# OPH-FEP-Synergetics

This package compares OPH's finite observer, a Gaussian model for variational
inference, and Fuller's lattice geometry. Its executable registry reports
which correspondences hold, where a proposed identification fails, and which
constructions can be used together.

- **Observer Patch Holography** (OPH) — an executable observer on a twelve-port
  icosahedral carrier with a native integer repair law.
- **The Free Energy Principle** (FEP) — variational inference and the
  conditional-independence structure of a Markov blanket.
- **Fuller's Synergetics** — quadray coordinates, the isotropic vector matrix,
  the cuboctahedral vector equilibrium, and rational tetravolumes.

![Three pairwise examples. The linear-Gaussian model gives the same posterior mean by three routes. The two kissing arrangements have different outer contacts. A coordination shell spans several Gaussian divergence levels.](docs/overview.svg)

Each edge shows one example from the [executable correspondence ledger](docs/CORRESPONDENCES.md).

[*Observers, Blankets, and Screens*](output/pdf/OPH-FEP-Synergetics_combined.pdf)
is the paper by Daniel Ari Friedman and Bernhard Mueller. Its numerical
results are bound to computations or cited sources in the manuscript build.

## What runs

```
oph-fep-synergetics verify                 # all 44 correspondence probes
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

`verify` exits `1` when a probe fails and `2` when a filter is unknown or
the combination of filters selects no probes.

## Install

```
uv sync --locked --extra dev --extra figures
uv run --locked --extra dev --extra figures oph-fep-synergetics verify
```

The package requires Python 3.12 or newer and NumPy. The `dev` extra adds
pytest and pytest-cov. The `figures` extra adds Matplotlib for the figure suite
and manuscript build. Add `--no-editable` to `uv sync` for a regular install.
The installed `sources` command reads the bundled source register from any
working directory.

For a resolver-based installation using pip, run the following with a supported
Python interpreter. Pip resolves the dependency bounds in `pyproject.toml`.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -e ".[dev,figures]"
```

## Layers

| layer | module | what it holds |
| --- | --- | --- |
| `core` | `core/carrier.py` | OPH's twelve-port carrier, its port chart, Laplacian, spectrum and load shell |
| | `core/repair.py` | the native seam-repair law, its conservation law and its exact energy decrease |
| | `core/native_dynamics.py` | ungated and gated rewrite relations, the binary attempt kernel, and ensemble relaxation |
| | `core/sheaf.py` | cellular sheaves, agreement energy, harmonic extension, record injectivity, cohomology dimensions |
| | `core/free_energy.py` | the linear-Gaussian free-energy dictionary and the three named surprises |
| | `core/controls.py` | the obstruction and dynamics controls: descent census, D-stability, mediators, orbits, free-energy forms |
| | `core/invariants.py` | the Laplacian's counting identities: spanning trees, reciprocal spectrum, GMRF precision, Markov blankets |
| | `core/symmetry.py` | the carrier's automorphism group, orientation classes, orbit counts and equivariance |
| | `core/network.py` | the carrier as a resistive network: exact effective resistances, the Kirchhoff index, shells around a port, and the subgraph census |
| `geometry` | `geometry/quadray.py` | quadray coordinates, the A3/IVM embedding and exact tetravolumes |
| | `geometry/lattice.py` | IVM shells, the cuboctahedral vector equilibrium, Voronoi cell and packing density, Goldberg/geodesic counts |
| | `geometry/subdivision.py` | icosahedral subdivision and the dual count identities |
| | `geometry/solids.py` | oriented faces, polyhedron volume, and the icosahedron's tetravolume against its closed form |
| `analysis` | `analysis/registry.py` | the 44 registered correspondences, each with a pairwise scope and an executable probe |
| | `analysis/provenance.py` | source pins, licence terms and reuse policy |
| `manuscript_variables` | `manuscript_variables/values.py` | every number the prose may quote, derived from the package where a derivation exists |

## The registry

44 entries: **19 isomorphisms**, **15 disagreements**, **10 complementarities**,
scoped to the pairs `oph-fep` (26), `synergetics-oph` (14), `fep-synergetics` (2)
and `three-way` (2). All 44 probes pass. `verify --kind` narrows a run to one of
the three kinds and `--pair` to one of the four scopes. The full table, with the OPH / FEP /
Synergetics reading of each entry and the values its probe asserts, is generated
into [docs/CORRESPONDENCES.md](docs/CORRESPONDENCES.md).

Selected results follow.

| claim | status | evidence | source or derivation |
| --- | --- | --- | --- |
| OPH's carrier port chart is derived from its oriented faces | exact | 12 ports, 30 seams and 20 faces. Each port has five neighbours, and the graph triangles recover the face vertex sets. | OPH’s oriented faces, with the adjacency derived here |
| The icosahedral Laplacian determines the mean port relaxation gap | exact | `lambda_1 = 5 - sqrt(5)`. With uniform seam selection and a fair coin, the mean repair operator is `I - L / (2\|E\|)`, with gap `(5 - sqrt(5)) / 60`. The same `L` supplies the blocks of the Gaussian model. | [Agreement and Surprise][agreement-and-surprise], §6.2, for the mean gap |
| The native repair law conserves load and changes `sum(x^2)/2` by `-(d^2 - (\|d\| mod 2)) / 4` | exact | checked on all 4680 transitions on the shell of total load two | [OPH][oph-source], native repair law |
| On the binary shell each native repair is exact Bayesian conditioning and a variational step | exact on that domain | in all 3960 cases the repair's law is the uniform reference conditioned on the datum the repair leaves visible; it minimises `F[q] = E_q[ln q - ln pi(s, o)]` over the fibre with minimum the surprisal `log 33` or `log 66`; off the shell the identity fails on the 60 seams that meet a doubled port | [Agreement and Surprise][agreement-and-surprise], §§6.2 and 7.3; [OPH synthesis][oph-synthesis] |
| The mean of the native law is one step of the carrier Laplacian | exact | averaging the integer rule over the 60 seam-and-coin draws gives `x - L x / 60` in exact rationals on all 78 configurations | [OPH][oph-source], native fields; exact enumeration here |
| Patch agreement is the global-section condition of a cellular sheaf, and the sheaf Laplacian is the carrier Laplacian | exact | `carrier_sheaf(carrier).laplacian() == carrier.laplacian()` | [Agreement and Surprise][agreement-and-surprise], §7.2 |
| On the linear-Gaussian carrier, constrained relative-entropy minimisation and predictive-coding perception give the posterior mean | exact on that domain | harmonic extension, a separate KKT solve and perception flow agree. With a proper record prior, variational free energy equals normalised surprisal at the posterior. An asymmetric seam separates its energy minimum from its posterior mean. | [Agreement and Surprise][agreement-and-surprise], Theorem 7, items 2–3; carrier, chain and asymmetric spring checks here |
| **PR-01**. Ungated repair is confluent and can cycle. The broken-seam gate terminates and admits multiple reachable normal forms. | **obstructed** | 78 states, no ungated normal forms, 66 gated normal forms and 5 reachable from each doubled state. The binary attempt kernel separately has the same computed gap as the mean port operator. Individual binary paths preserve load energy and ensemble relative entropy is nonincreasing. | [OPH][oph-source], PR-01; exhaustive reachability and kernel checks here |
| **Three distinct values at fixed records** | **obstructed** | minimum interface energy, the unnormalised record potential and normalised surprisal differ. The interior determinant cancels in the normalised probability. | [Agreement and Surprise][agreement-and-surprise], Theorem 7, item 3; equal-Schur chain control here |
| Two-record Gaussian marginals are determined by effective resistance and the record prior | exact | all 66 carrier pairs satisfy the Schur identity; seam, distance-two and antipodal resistance are 11/30, 7/15 and 1/2 | Schur elimination and an independent rational network solve |
| **Eliminating an interior mediator is not free** | **obstructed** | the reduced precision is `(1/3) v v^T + epsilon I`; a sum of pairwise semidefinite terms reproduces it iff `epsilon >= 1/3` | [Agreement and Surprise][agreement-and-surprise], Proposition 8; rank-one construction here |
| **Distinct drifts preserve the same stationary blanket** | **obstructed** | with a zero precision entry the internal drift is `2q eta + (q - 1) b - 2 mu`, which varies with `q` at fixed stationary covariance | [Agreement and Surprise][agreement-and-surprise], Proposition 11; Da Costa et al.’s stationary blanket example |
| **Agreement can exist and be unreachable** | **obstructed** | the clipped two-node map has the origin as its only equilibrium and a source there; the relaxation orbit stays at radius `>= 0.84` with period `~7.58` | [Agreement and Surprise][agreement-and-surprise], Proposition 15; exact fixed-point census and orbit measurement here |
| **D-stability is not free** | **obstructed** | Lyapunov form `P = (1/3)[[11,4,-8],[4,5,-2],[-8,-2,9]]`, `det P = 115/27`; at `k = 4` the Routh first column is `(1, 6, -2, 24)` and the scaled matrix has two right-half-plane roots | [Agreement and Surprise][agreement-and-surprise], §4.1 and Proposition 10; general Routh column here |
| **Finite phases are path functions** | **obstructed** | full and half-damped iteration share every answer, yet their centered contrasts after twelve steps differ by exactly `921633999 / 40000000000` | [Agreement and Surprise][agreement-and-surprise], Proposition 6 and §4.4; exact contrast here |
| Quadray geometry supplies the exact integer side of the correspondence | exact | 12 neighbour moves at squared distance 8, shells `10k^2 + 2`, rational volumes `1, 3, 4, 6, 20, 96` | QuadMath conventions; lattice counts and volumes derived here |
| Equilibrium propagation supplies the parameter gradient the repair law does not define | exact | the two-settle contrast is `(theta - target) / (1 + beta)` exactly and tends to the exact gradient | [Agreement and Surprise][agreement-and-surprise], Theorem 21; Scellier and Bengio |
| The carrier counts its own spanning trees | exact | the reduced Laplacian's determinant is `5184000` — by Kirchhoff's matrix-tree theorem the number of spanning trees in OPH's reading and, as the determinant of the same matrix read as a Gaussian precision with one port pinned, the normalisation of the FEP's linear-Gaussian interface integral — equal to the product of the eleven nonzero eigenvalues over `12` | Kirchhoff’s matrix-tree theorem applied here to the carrier |
| A proper Gaussian with precision `L + I/n` has the seam graph's conditional dependencies | exact | seam precision is `-1` and each port's blanket is its five neighbours. The centered intrinsic field separately gives port variance `7/36` | Carrier Gaussian precisions and intrinsic variance calculated here |
| The faces fill the cycle space of the carrier graph | exact | For the scalar sheaf, the graph has a one-dimensional space of global sections and `dim H^1 = 19`. Adding the faces gives Betti numbers `(1, 0, 1)`. The oriented faces sum to a fundamental cycle. | Cellular sheaf construction from [Agreement and Surprise][agreement-and-surprise], §7.2; carrier cohomology calculated here |
| The symmetry group preserves the native repair law | exact | The group has 120 elements, with 60 rotations and 60 orientation-reversing symmetries. Equivariance is checked over 1123200 combinations of group element, configuration, oriented seam and coin. The shell has 4 configuration orbits and 3 outcome multisets. | Carrier automorphisms and native transitions enumerated here |
| Gaussian divergence measures half the squared Euclidean displacement | exact | At unit covariance, one IVM neighbour displacement gives 4 nats. Displacement `k m` along a neighbour direction `m` gives `4k^2` nats from the origin. | Equal-covariance Gaussian divergence applied to QuadMath’s embedding |
| **Twelve vertices is not one solid** | **obstructed** | the neighbour graph is four-regular with tetravolume 20; the carrier is five-regular with tetravolume `18.512295868219...` | [OPH][oph-source] and QuadMath geometry compared here |
| **The shell index and the free energy determine neither each other** | **obstructed** | the coordination shell is the cuboctahedral gauge `max(\|u\|_inf, \|u\|_1 / 2)`, a polyhedral norm, checked against every breadth-first shell through 12; the free energy on shell `k` ranges from `2k^2` (`2k^2 + 2` for odd `k`) to `4k^2`, attained only at the twelve vertices, and the levels 36 and 52 are each shared by two shells | The Gaussian cost and cuboctahedral gauge compared here |
| Synergetics' tetravolume table is derived, not transcribed | exact | `1, 3, 4, 6, 20, 96` come out of explicit integer vertex sets by exact convex-hull volumes in the rationals | QuadMath conventions; exact convex-hull volumes derived here |
| The carrier's own solid has an exact algebraic volume | exact | computed in exact field arithmetic as `sqrt 2 (15/2 + (5/2) sqrt 5) = 18.512295868219...`, squared exactly `175 + 75 sqrt 5`; its conjugate is `2.7009...` and the product is exactly `50` | Carrier geometry; exact field computation here |
| Both twelve-around-one arrangements attain the kissing number, differently | exact | squared edge over squared circumradius is exactly `1` for the cuboctahedron, where each outer sphere touches four others, and `2 - 2 sqrt 5 / 5` for the icosahedron, where each outer sphere touches only the centre | Carrier and vector equilibrium compared here |
| The twelve-around-one cell is the densest packing's cell | exact | 14 vertices, 12 rhombic faces, 24 edges, tetravolume exactly 6, ball tetravolume exactly `pi sqrt 2`, density exactly `pi / (3 sqrt 2)` — optimal among lattice packings (Gauss) and all packings (Hales) | Exact Voronoi cell calculated here; packing optimality from Gauss and Hales |

[agreement-and-surprise]: https://philpapers.org/rec/MUEAAS-2
[oph-source]: https://github.com/FloatingPragma/observer-patch-holography
[oph-synthesis]: https://philpapers.org/rec/MUEFOC

## Scope

The executable scope is the finite OPH carrier, specified Gaussian models and
controls, and selected lattice and volume constructions. The domain of a
claim includes its reference distribution, coordinate scale, and repair
schedule where those assumptions enter the proof. Extensions to nonlinear
inference, continuum OPH, and empirical comparisons remain open. See the
[scope and limits](docs/SCOPE.md).

Three quantities that *Agreement and Surprise* states for the native law are
reproduced exactly from the protocols in OPH's public model. One attempt is a
uniform seam and a fair coin: from `(1, 1, 0, ..., 0)` the configuration stays
with probability `13/15` and reaches each of eight neighbouring configurations with probability
`1/60` (`pr01-obstruction`). The two-event observer that chooses its second
seam from its first record ends with a port-zero occupation of `35/264`,
against `1/6` with the selector disconnected (`reduction-criterion`).

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

The carrier and native seam-repair law reimplement OPH's Apache-2.0 software.
The linear-Gaussian dictionary follows [Agreement and Surprise][agreement-and-surprise].
The FEP_Formal catalogue provides topic references under CC BY 4.0 (DOI
`10.5281/zenodo.19699233`). The quadray conventions and lattice arithmetic
follow Apache-2.0 QuadMath (DOI `10.5281/zenodo.16887791`). Each source keeps
its own terms, recorded with commit pins in
[data/sources.json](data/sources.json) and restated in
[docs/PROVENANCE.md](docs/PROVENANCE.md).

This project is released under the MIT licence; see [LICENSE](LICENSE).

## The manuscript

The Pandoc and LaTeX pipeline from `docxology/template` builds the paper
from `manuscript/`. Run the build from this repository after installing the
figure and test extras.

```sh
.venv/bin/python scripts/build.py --template-root /path/to/template
```

The script generates figures and numerical bindings, renders the manuscript,
and checks the PDF and its log. It selects the template checkout from
`--template-root`, `$TEMPLATE_ROOT`, or the sibling `template/` directory.
An optional `--deliver` also copies the verified PDF to `~/Downloads/` and
prints both SHA-256 digests.

`projects/working/OPH-FEP-Synergetics` in the template checkout must point to
the checkout you want to build. The script checks that link, uses the
template interpreter to resolve the manuscript bindings, and writes the PDF
in this project's `output/pdf/`.

The rendered edition has a deliberate front matter. Page 1 is the cover, drawn
from `output/figures/fig_cover.png`; pages 2 to 5 are the four panels of the
graphical abstract, one panel per page, placed with `\floatplacement{figure}{H}`
so that a panel cannot float away from its caption; the contents page follows the
panels, which is why `preamble.md` keeps the renderer's automatic
`\tableofcontents` from landing between the cover and the panels and re-emits it
from the end of `00_00_graphical_abstract.md`. The margins are 0.5 in, because
this is dense reference matter with wide tables.

Three disciplines hold the manuscript together, and each one is a test:

- **Numerical results have bindings.** The prose quotes computed or sourced
  results as `{{TOKEN}}`. The build checks that every token resolves and the
  tests compare the manuscript's tokens with the package's bindings.
- **Every formalism is numbered by the renderer.** Definitions and propositions
  are fenced `::: {.proposition #prop:key}` blocks referenced as `[@prop:key]`,
  and the template's `formalism.lua` owns the counters. Cross-references
  follow each statement's label when its number changes.
- **Every citation resolves.** Pandoc's `--natbib` resolves `[@key]` against
  `manuscript/references.bib`. A test asserts that the cited keys and the
  bibliography are the same set: no undefined citation, and no uncited entry.

## Layout

```
data/sources.json                    pinned sources, licences, cited identifiers
docs/CORRESPONDENCES.md              generated: the 44 entries and their probe values
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
tests/                               the suite (420 tests)
```

## Tests

```
.venv/bin/python -m pytest tests/ -q --cov --cov-report=term-missing
```

420 tests. The suite checks exact values rather than
smoke behaviour: the 78-state shell and all 4680 transitions, the exact
rational contrast, the closed-form Routh column, the Lyapunov determinant
`115/27`, the mediator threshold at `epsilon = 1/3`, the IVM round-trip over
shell sites, the resistance ladder's exact fractions against a from-scratch solve,
and the exact volume ratios. The coverage gate is 90%, set in `pyproject.toml`.

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
repository". See [CHANGELOG.md](CHANGELOG.md) for changes by edition.
