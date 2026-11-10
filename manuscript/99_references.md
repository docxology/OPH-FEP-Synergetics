# References {#sec:references}

The bibliography is [`references.bib`](references.bib), read by Pandoc during the
PDF render. The build invokes Pandoc with `--natbib`, so every `[@key]` citation
in the manuscript is rewritten to the corresponding LaTeX command and resolved
against the bib file.

Three entries are provenance records rather than literature: `oph`, `fepformal`
and `quadmath` identify the pinned checkouts, their licences and their DOIs. A
reader who wants the exact commits should read `data/sources.json` or run
`oph-fep-synergetics sources`, which reports the live status of each pin rather
than repeating a snapshot here.

The bibliography is validated for syntax and for the required fields of each
entry type by the template's citation validator:

```bash
uv run python -m infrastructure.reference.citation.cli validate \
    projects/working/OPH-FEP-Synergetics/manuscript/references.bib --strict
```
