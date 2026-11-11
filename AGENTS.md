# OPH-FEP-Synergetics — agent guidance

Mapping the executable parts of Observer Patch Holography, the Free Energy
Principle and Fuller's Synergetics onto one another, and reporting honestly
where they are isomorphic, where they disagree, and where they merely
complement each other.

Read [README.md](README.md) first, then [docs/SCOPE.md](docs/SCOPE.md). The
scope limits there are part of the deliverable, not caveats to be trimmed.

## Where things go

| you are adding | put it in |
| --- | --- |
| a numerical or structural claim | `src/oph_fep_synergetics/core/` or `geometry/`, as a function that returns its values |
| a correspondence, an obstruction or a complementarity | `analysis/registry.py`, as a `Correspondence` with a probe |
| a number the prose will quote | `manuscript_variables/values.py`, as a `BoundValue` with its source |
| a source pin or a licence term | `data/sources.json` |
| a manuscript section | `manuscript/<nn>_<name>.md`, with a `{#sec:…}` label on its top-level heading |
| a figure | a builder in `scripts/make_figures.py`, plotted only from package functions |
| a citation | `manuscript/references.bib`, cited from a section and nowhere else |
| a test | `tests/`, named for the layer it exercises |

`docs/CORRESPONDENCES.md` and `docs/MANUSCRIPT_VALUES.md` are **generated**.
Never edit them by hand: change the registry or the values and re-render.

## Invariants

1. **The registry is the source of truth.** Every isomorphisms/disagreement/
   complementarity claim is one `Correspondence` with an executable probe that
   asserts exact values. A claim without a probe does not belong in the
   registry.
2. **Statuses come from a fixed vocabulary**: `exact`, `exact-on-domain`,
   `obstructed`, `not-reproduced`. Use `exact-on-domain` when a claim holds only
   in the linear-Gaussian (or otherwise stated) corner, and `obstructed` when
   the result is a failure of the identification rather than a proof of it.
   Never upgrade a status because a neighbouring claim passed.
3. **A probe reports values, not just a boolean.** `Correspondence.run()` returns
   `CorrespondenceResult` with `values`; tests assert those values so a silent
   numeric change is caught.
4. **Do not restate a source commit in code or prose.** The commit actually read
   lives in `data/sources.json` and is reported by `oph-fep-synergetics sources`.
   A test fails when a 40-hex commit appears in the package without appearing in
   the pin file.
5. **Never copy source prose, and never reuse OPH's data or ledgers.** OPH's
   software is Apache-2.0 and may be reimplemented with attribution; its papers,
   book, documentation, figures, data and ledgers are CC BY-NC-SA 4.0. Cite the
   identifier (file, line, symbol), reimplement the behaviour, and do not paste
   text or data.
6. **Derive, do not transcribe, where a derivation is possible.** The carrier
   adjacency is derived from OPH's `ORIENTED_BASE_FACES` and the derivation is
   asserted: the 20 triangles of the derived graph must be exactly those faces.
   Prefer this pattern whenever a source table can be checked against its own
   consequences.
7. **Report a moved pin as moved, not as a failure.** A pin is a point-in-time
   reading. `check_pins` distinguishes `match`, `moved`, `absent` and
   `unavailable`.
8. **Do not claim more than the code does.** Measurements are labelled as
   measurements: the source states the `0.84` bound and the `7.58` period, and
   this package states its own measured `0.8401` and `7.575` alongside them.
   Where the source cannot be reproduced from its text, record
   `not-reproduced` and say why.
9. **The prose quotes tokens, never numbers.** Every numeric claim in
   `manuscript/` is a `{{TOKEN}}` bound to a `BoundValue`. The hydration script
   fails when a token survives, and a test asserts that the tokens the prose
   uses and the values the package publishes are the same set. Never paste a
   literal number into a section, and never widen a tolerance so that a token
   prints what the sentence expects.
10. **Formal statements are numbered by the renderer.** A definition or a
    proposition is a fenced block, `::: {.proposition #prop:key}`, referenced as
    `[@prop:key]`. The template's `formalism.lua` owns the counters, so nothing
    hand-numbers a statement and the numbering cannot drift between editions.
11. **Citations are exact in both directions.** Every `[@key]` resolves in
    `manuscript/references.bib`, and every entry there is cited somewhere. A
    test asserts the two sets are equal.
12. **Every top-level heading carries a `{#sec:…}` label.** The render
    annotates a heading that lacks one and writes the annotation back into the
    source, so keeping the labels explicit is what makes a render idempotent.
13. **The front matter is placed on purpose, and the tests check the page
    positions.** Page 1 is the cover, pages 2 to 5 are the four panels of
    `manuscript/00_00_graphical_abstract.md` (one panel per page, held in place
    by `\floatplacement{figure}{H}`), and the contents page is re-emitted by
    `\ophContentsPage` at the end of that file because `preamble.md` keeps the
    renderer's automatic `\tableofcontents` from landing between the cover and
    the panels. Changing the margin, the panel figure height fraction or a
    panel's length without re-rendering moves a panel off its page.

## Commands

```
.venv/bin/python -m pytest tests/ -q                  # the suite
.venv/bin/python -m pytest tests/ -q --cov --cov-report=term-missing
.venv/bin/oph-fep-synergetics verify                  # the registry, all 35 probes
.venv/bin/oph-fep-synergetics verify --pair oph-fep   # one pairwise ledger
.venv/bin/oph-fep-synergetics sources --vault-root /Volumes/GitVault
.venv/bin/python scripts/render_docs.py               # regenerate the two documents
.venv/bin/python scripts/render_docs.py --check       # fail when stale
.venv/bin/python scripts/make_figures.py              # the sixteen figures and their registry
```

The manuscript renders through the `docxology/template` checkout, whose
pipeline owns Pandoc and LaTeX:

```
cd /path/to/template
TEMPLATE_REPO_ROOT=$PWD .venv/bin/python /path/to/project/scripts/z_generate_manuscript_variables.py
.venv/bin/python scripts/pipeline/stage_03_render.py \
    --project working/OPH-FEP-Synergetics --skip-manuscript-hydration
```

`projects/working/OPH-FEP-Synergetics` is a symlink to this project. Hydration
runs with the template's interpreter because the stage's own hydration subprocess
requires an interpreter that passes the template's bytecode audit.

## Definition of done

- `pytest` is green and coverage stays at or above the 90% gate.
- `oph-fep-synergetics verify` exits 0 with all probes passing.
- `scripts/render_docs.py --check` exits 0.
- Any new number in prose exists as a `BoundValue` and appears in
  `docs/MANUSCRIPT_VALUES.md`.
- Any new correspondence names its OPH, FEP and Synergetics reading and its
  source.
- The manuscript renders: the PDF exists, its LaTeX log has no overfull box, and
  its `.tex` has no unresolved `??` cross-reference.
- The test count in `README.md` and `docs/VERIFICATION.md` matches the collected
  count; `test_docs.py` fails when it does not.

## Do not

- Do not add a dependency beyond `numpy` without a reason recorded here.
- Do not commit `output/`, `.venv/`, coverage artefacts or `__pycache__`.
- Do not weaken the coverage gate or the strictness of `verify` to make a run
  pass. Fix the claim or restate it at a weaker status.
- Do not hand-number a definition, a proposition or an equation; the renderer
  owns the counters, and a hand number cannot survive an edit above it.
- Do not edit anything under `output/`: regenerate the figures, the resolved
  tokens and the PDF with the scripts and the template pipeline.
- Do not push this repository anywhere; it is a local working copy.
