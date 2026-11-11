# Reproducibility: pins, commands, and the binding of every number to a function {#sec:reproducibility}

## What this manuscript is made of

The manuscript is a project of the `docxology/template` research monorepo: LaTeX
is produced by Pandoc, citations are resolved with `--natbib` against
[`references.bib`](references.bib), equations are numbered by `pandoc-crossref`,
and formal statements are numbered by the repository's `formalism.lua` filter, so
that a definition or a proposition has the same number in every rendered edition
[@template].

Every number in the prose is a brace-delimited token resolved at build time from
`manuscript_variables/values.py`, which binds each value
to the function or source that produces it. A token that the generator does not
emit fails the build, and a value that changes changes the prose. The binding is
the reason the paper can claim that its numbers are not transcribed.

| Surface | Authored | Generated |
|:-------------------------------------------|:-------------|:------------------------------------------|
| `manuscript/*.md` | yes | never |
| `references.bib` | yes | never |
| `src/`, `geometry/`, `analysis/` | yes | never |
| `output/figures/*.png` | no | `scripts/make_figures.py` |
| `output/figures/figure_registry.json` | no | `scripts/make_figures.py` |
| `docs/CORRESPONDENCES.md`, `docs/MANUSCRIPT_VALUES.md` | no | `scripts/render_docs.py` |
| manuscript tokens | no | `scripts/z_generate_manuscript_variables.py` |

: Which surfaces are written by hand and which are produced. A generated surface
is never edited to make a check pass; the producer is changed instead.
{#tbl:surfaces}

## The three source pins

The three programmes are read at pinned commits. The commit itself lives in
`data/sources.json` and is reported by `oph-fep-synergetics sources`, which also
reports whether the live checkout still matches it.

| Source | Repository | Licence | Recorded |
|:-----------------|:-------------------------------------------|:-------------------------|:-----------|
| Observer Patch Holography | `FloatingPragma/observer-patch-holography` | Apache-2.0 (software) and CC BY-NC-SA 4.0 (papers, data, ledgers) | commit, branch, seven cited identifiers |
| FEP_Formal | `ActiveInferenceInstitute/FEP_Formal` | CC BY 4.0, DOI `10.5281/zenodo.19699233` | commit, branch, four cited identifiers |
| QuadMath | `docxology/quadmath` | Apache-2.0, DOI `10.5281/zenodo.16887791` | commit, branch, five cited identifiers |

: The three pins. The `Recorded` column is the number of identifiers the pin file
attributes to each source; each identifier names a file, a line range and a
symbol, so that a reader can check the reading rather than trust it.
{#tbl:pins}

A pin is a reading, not a promise. `oph-fep-synergetics sources` reports `match`
when the live checkout is at the pinned commit, `moved` when it has advanced,
`absent` when the checkout is not present, and `unavailable` when Git cannot be
queried. A `moved` result is not a failure: it says the reading is historical.

## The commands

```bash
# the suite, and the coverage gate (90% branch minimum)
.venv/bin/python -m pytest tests/ --cov --cov-report=term-missing

# every correspondence, in one run or one pairwise scope at a time
.venv/bin/oph-fep-synergetics verify
.venv/bin/oph-fep-synergetics verify --pair synergetics-oph
.venv/bin/oph-fep-synergetics verify --kind disagreement
.venv/bin/oph-fep-synergetics verify --json

# the pins, checked against the live checkouts
.venv/bin/oph-fep-synergetics sources --vault-root /path/to/vault

# the figures and their accessibility registry
.venv/bin/python scripts/make_figures.py

# the generated reference documents
.venv/bin/python scripts/render_docs.py --check

# the manuscript tokens, then the PDF; both run from the template checkout
cd /path/to/template
TEMPLATE_REPO_ROOT=$PWD .venv/bin/python /path/to/project/scripts/z_generate_manuscript_variables.py
.venv/bin/python scripts/pipeline/stage_03_render.py --project working/OPH-FEP-Synergetics --skip-manuscript-hydration
```

: The reproduction commands, in order. The plain `verify` runs all
{{REGISTRY_ENTRIES}} probes; `--pair` runs one ledger, of
{{REGISTRY_OPH_FEP}}, {{REGISTRY_SYNERGETICS_OPH}},
{{REGISTRY_FEP_SYNERGETICS}} and {{REGISTRY_THREE_WAY}} claims respectively, and
`--json` writes the same report as machine-readable values. The last two commands
run from the `docxology/template` checkout, because the rendering infrastructure
lives there. {#tbl:commands}

## Environment

The package depends on `numpy` and on nothing else at run time; the test and
figure extras add `pytest`, `pytest-cov` and `matplotlib`. The suite runs on
Python 3.15 with `numpy` 2.5 and has no network dependency; the provenance check
reads the local checkouts with `git rev-parse` and reports `absent` when they are
not mounted.

The rendering toolchain is Pandoc 3.13 with `pandoc-crossref` and a TeX Live
installation. The template's own documentation records the full toolchain
[@template]; this section records only what the manuscript needs.

## What a reader should check first

A reader who wants to test the paper rather than read it should run three
commands. `oph-fep-synergetics verify` reproduces every number in
[@sec:correspondences], [@sec:disagreements] and [@sec:complementarities] and
exits non-zero on any failure. `oph-fep-synergetics sources` shows whether the
three readings are still current. And `pytest` shows whether the values the prose
quotes are still the values the code returns, which is the property that makes
this paper different from a summary of three literatures.
