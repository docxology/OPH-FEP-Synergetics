# Reproducibility {#sec:reproducibility}

## What this manuscript is made of

The manuscript uses the `docxology/template` rendering pipeline [@template].
Pandoc produces LaTeX, and BibTeX processes
[`references.bib`](https://github.com/docxology/OPH-FEP-Synergetics/blob/9cf3cd59fdcdc848fb4cfe19ba1604ff3bcb06ba/manuscript/references.bib). The `pandoc-crossref` and `formalism.lua`
filters use stable labels to generate numbering and cross-references.
Numbers can change when statements move between editions.

Quoted numerical results are brace-delimited tokens resolved at build time.
`manuscript_variables/values.py` binds each token to its computation or source.
The build fails if a token remains unresolved, and the tests verify that the
manuscript and value table use the same set of bindings.

[@tbl:surfaces] identifies the source or build step for each surface. The
figure PNGs, accessibility registry and standalone poster share the `figures`
subdirectory of `output`.

| Surface | Authored | Generated |
|:-------------------------------------------|:-------------|:------------------------------------------|
| Main manuscript sections | yes | never |
| `data/ontology_concordance.json` | yes | never |
| `docs/manuscript/90_ontology_concordance.md` | no | `scripts/render_docs.py` |
| `docs/manuscript/references.bib` | yes | never |
| `src/oph_fep_synergetics/` | yes | never |
| `uv.lock` | no | `uv lock` |
| Figure PNGs | no | `scripts/make_figures.py` |
| `figure_registry.json` | no | `scripts/make_figures.py` |
| `OPH-FEP-Synergetics_Jitterbug_Poster.pdf` | no | `scripts/make_figures.py` |
| `OPH-FEP-Synergetics_combined.pdf` in `output/pdf/` | no | `scripts/build.py` with the template checkout |
| `docs/CORRESPONDENCES.md`, `docs/MANUSCRIPT_VALUES.md`, `docs/CONCORDANCE.md` | no | `scripts/render_docs.py` |
| manuscript tokens | no | `scripts/z_generate_manuscript_variables.py` |

: Which surfaces are written by hand and which are produced. A generated surface
is never edited to make a check pass; the producer is changed instead.
{#tbl:surfaces}

The ontology concordance in [@sec:ontology-concordance] maps terms to their
model-specific meanings, applicability and relevant registry probes.
[@tbl:ontology-dynamics] distinguishes the repair, inference and dynamical
quantities that the probes report.

## The three source pins

The three programmes are read at pinned commits. The commit itself lives in
`data/sources.json` and is reported by `oph-fep-synergetics sources`, which also
reports whether the live checkout still matches it. [@tbl:pins] lists the
licences and identifier counts.

| Source | Repository | Licence | Recorded |
|:-----------------|:-------------------------------------------|:-------------------------|:-----------|
| Observer Patch Holography | `FloatingPragma/observer-patch-holography` | Apache-2.0 (software) and CC BY-NC-SA 4.0 (papers, data, ledgers) | commit, branch, thirteen cited identifiers |
| FEP_Formal | `ActiveInferenceInstitute/fep_formal` | CC BY 4.0, DOI `10.5281/zenodo.19699233` | commit, branch, four cited identifiers |
| QuadMath | `docxology/QuadMath` | Apache-2.0, DOI `10.5281/zenodo.16887791` | commit, branch, five cited identifiers |

: The three pins. The `Recorded` column is the number of identifiers the pin file
attributes to each source; each identifier names a file, a line range and a
symbol, so that a reader can check the reading rather than trust it.
{#tbl:pins}

The pin records the commit used for the reading. With `--vault-root`,
`oph-fep-synergetics sources` reports `match` for the pinned commit, `moved`
for a different commit, `absent` for a path without Git metadata, and
`unavailable` for an unreadable checkout or HEAD. Without a vault root the
command reports the recorded pins and leaves the live checks unavailable.

## The commands

Run the commands in [@tbl:commands] from the project root. File paths for
the vault and template checkout are examples.

| Step | Command |
|:------------------|:----------------------------------------------------------|
| Install the locked environment | `uv sync --locked --extra dev --extra figures` |
| Run the tests and coverage gate | `.venv/bin/python -m pytest tests/ --cov --cov-report=term-missing` |
| Verify all probes | `.venv/bin/oph-fep-synergetics verify` |
| Filter by scope | `.venv/bin/oph-fep-synergetics verify --pair synergetics-oph` |
| Filter by claim kind | `.venv/bin/oph-fep-synergetics verify --kind disagreement` |
| Export probe results | `.venv/bin/oph-fep-synergetics verify --json` |
| Compare pins against checkouts | `.venv/bin/oph-fep-synergetics sources --vault-root /path/to/vault` |
| Check the generated documents | `.venv/bin/python scripts/render_docs.py --check` |
| Build figures and the poster | `.venv/bin/python scripts/make_figures.py` |
| Build and check the manuscript PDF | `.venv/bin/python scripts/build.py --template-root /path/to/template` |

: Reproduction commands. Scope and kind filters select subsets of the probe
registry. {#tbl:commands}

The manuscript build already runs the figure step. Adding `--deliver`
copies the checked manuscript PDF to the `Downloads` folder in your home
directory and verifies that the copy matches the build output. The standalone
poster stays with the generated figure PNGs.

## Environment

NumPy is the runtime dependency. The test and figure extras add `pytest`,
`pytest-cov` and Matplotlib. `pyproject.toml` declares the supported Python
range and dependency bounds. `uv.lock` records the Python dependency
resolution, and `uv sync --locked` checks its consistency with the project
metadata. Platform and interpreter markers determine which locked packages
are installed. A `pip` install resolves the declared bounds independently.

The suite runs offline once the dependencies and tools are installed.
The exact geometry uses integer and rational arithmetic and the quadratic
field implemented in the package. Gaussian calculations and numerical sheaf
ranks use NumPy. The build prints the interpreter and NumPy versions for
the run.

Rendering also uses the separate `docxology/template` environment, Pandoc,
`pandoc-crossref`, and LuaLaTeX with PDF tagging support [@template].
The render check uses Poppler's `pdftotext` to inspect the generated PDF. The
figure registry supplies the embedded alternative descriptions. In the template checkout,
`projects/working/OPH-FEP-Synergetics` must resolve to the project being
built. The Python lock covers this project. The rendered PDF also depends on
the template revision and the external toolchain.

## What a reader should check first

Start with `oph-fep-synergetics verify` for the probe verdicts and computed
values. Run `oph-fep-synergetics sources --vault-root /path/to/vault` to
compare the recorded pins with local checkouts. Run the test suite to check
the package and manuscript bindings, and use `scripts/build.py` to validate
a fresh render.
