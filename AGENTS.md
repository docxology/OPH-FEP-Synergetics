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
| a terminology comparison for the ontology supplement | `data/ontology_concordance.json` |
| a manuscript section | `manuscript/<nn>_<name>.md`, with a `{#sec:…}` label on its top-level heading |
| a figure | a builder in `scripts/make_figures.py`, plotted only from package functions |
| a citation | `manuscript/references.bib`, cited from a section and nowhere else |
| a test | `tests/`, named for the layer it exercises |

`docs/CORRESPONDENCES.md`, `docs/MANUSCRIPT_VALUES.md`, `docs/CONCORDANCE.md`
and `manuscript/90_ontology_concordance.md` are **generated**. Edit their data
or manuscript sources and regenerate them with `scripts/render_docs.py`.

## Invariants

1. **The registry is the source of truth.** Every isomorphism, disagreement or
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
   A test fails when any prefix of seven or more hex digits of a pinned commit
   appears in the package. Re-pin by re-reading every cited identifier at the
   new upstream commit, not by editing the hash.
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
    `[@prop:key]`. The template's `formalism.lua` owns the counters.
    Cross-references follow the labels when statement numbers change.
11. **Citations are exact in both directions.** Every `[@key]` resolves in
    `manuscript/references.bib`, and every entry there is cited somewhere. A
    test asserts the two sets are equal. Authored sections also refer to every
    numbered figure, table, equation and formal statement. The generated
    concordance reports the source sections for these references.
12. **Every top-level heading carries a `{#sec:…}` label.** The render
    annotates a heading that lacks one and writes the annotation back into the
    source, so keeping the labels explicit is what makes a render idempotent.
13. **Exact where the geometry allows.** Counts, ranks, determinants and volumes
    are computed in `fractions.Fraction` or in `geometry.exact.Surd5`
    (`Q(sqrt 5)`), not in floating point; a float measurement may stand beside
    an exact value as a second route, never in place of one. A table entry is
    derived from an explicit vertex set, not typed.
14. **A probe never checks a function against itself.** An identity needs two
    routes that share no code path (an interior solve and a KKT solve; surprisal
    plus divergence and energy minus entropy; an enumeration and a closed form).
    A literal such as `"integral": False` is not evidence.
15. **Say which object a word names.** OPH's load shell fixes total load. The
    lattice coordination shell fixes the cuboctahedral gauge. The carrier
    Laplacian determines the mean repair update. The native repair energy is
    half the sum of squared port loads. The proper Gaussian model places its
    conditional dependencies on the seams. Use the relevant term in every
    claim and caption.
16. **The front matter is placed on purpose, and the tests check the page
    positions.** Page 1 is the cover, pages 2 to 6 are the five panels of
    `manuscript/00_00_graphical_abstract.md` (one panel per page, held in place
    by `\floatplacement{figure}{H}`), and the contents page is re-emitted by
    `\ophContentsPage` at the end of that file because `preamble.md` keeps the
    renderer's automatic `\tableofcontents` from landing between the cover and
    the panels. Regenerate and check the PDF after changing the margin, the
    panel figure height fraction or a panel's length.

## Commands

```
uv run --no-project --python 3.12 python scripts/check_local_ci.py  # both local CI versions
.venv/bin/python -m pytest tests/ -q                  # the suite
.venv/bin/python -m pytest tests/ -q --cov --cov-report=term-missing
.venv/bin/oph-fep-synergetics verify                  # the registry, all 80 probes
.venv/bin/oph-fep-synergetics verify --pair oph-fep   # one pairwise ledger
.venv/bin/oph-fep-synergetics verify --kind disagreement   # one kind of claim
.venv/bin/oph-fep-synergetics sources --vault-root /path/to/vault
.venv/bin/python scripts/render_docs.py               # regenerate the reference docs and supplement
.venv/bin/python scripts/render_docs.py --check       # fail when stale
.venv/bin/python scripts/make_figures.py              # thirty figures, vector Jitterbug poster and registry
.venv/bin/python scripts/build.py --template-root /path/to/template
```

The render check needs Poppler's `pdftotext`. Add `--deliver` to copy the
checked combined PDF to `~/Downloads/`. The standalone poster is generated
in `output/figures/`.

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

- The [local CI gate](docs/VERIFICATION.md#local-ci) passes on Python 3.12 and
  3.14 for the exact candidate commit. It installs the locked package, checks the
  installed source register, runs the tests at the 90% coverage gate, verifies
  all probes and checks the generated documents.
- Any new number in prose exists as a `BoundValue` and appears in
  `docs/MANUSCRIPT_VALUES.md`.
- Any new correspondence names its OPH, FEP and Synergetics reading and its
  source.
- Rebuild the manuscript when its rendering inputs change. The PDF has no
  unresolved reference placeholders, and its LaTeX log has no overfull boxes or
  undefined references or citations. The printed captions match the numbered
  manuscript inventory.
  Check the cover and changed figure pages at page size. The tagged PDF uses
  the alternative descriptions in the figure registry and cover configuration.
- The test count in `README.md` and `docs/VERIFICATION.md` matches the collected
  count; `test_docs.py` fails when it does not.

## Do not

- Do not add a dependency beyond `numpy` without a reason recorded here.
- Do not commit `output/`, `.venv/`, coverage artefacts or `__pycache__`. The
  one tracked artefact is the release PDF,
  `output/pdf/OPH-FEP-Synergetics_combined.pdf`: it is re-rendered with
  `scripts/build.py` and must pass `scripts/check_render.py` before it is
  committed, so the published PDF is always the one the sources produce.
- Do not weaken the coverage gate or the strictness of `verify` to make a run
  pass. Fix the claim or restate it at a weaker status.
- Do not hand-number a definition, a proposition or an equation; the renderer
  owns the counters, and a hand number cannot survive an edit above it.
- Do not edit anything under `output/`: regenerate the figures, the resolved
  tokens and the PDF with the scripts and the template pipeline.
- Publish only to the configured
  [docxology/OPH-FEP-Synergetics](https://github.com/docxology/OPH-FEP-Synergetics)
  remote through the maintainer's approved workflow. In GitVault, follow the
  owned-project publication procedure in the vault's operations guide. Complete
  both local CI checks before pushing, then verify that the remote points to the
  tested commit. Project CI runs locally. Managed Dependabot operates separately.
  Never force-push.
