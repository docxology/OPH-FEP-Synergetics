# References {#sec:references}

The PDF bibliography is generated from [`references.bib`](references.bib).
Pandoc converts the manuscript citations to LaTeX commands, and BibTeX resolves
their keys.

The `oph`, `fepformal` and `quadmath` entries identify the software sources.
`data/sources.json` records the commits, licences and cited identifiers.
`oph-fep-synergetics sources --vault-root /path/to/vault` compares those pins
with local checkouts.

Run the template's citation validator from the `docxology/template` checkout
to check bibliography syntax and required fields.

```bash
uv run python -m infrastructure.reference.citation.cli validate \
    projects/working/OPH-FEP-Synergetics/manuscript/references.bib --strict
```

This repository checks that every cited key has a bibliography entry and each
entry is cited. `scripts/check_render.py` also checks that citations resolve
in the rendered PDF.
