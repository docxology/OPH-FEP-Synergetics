"""The manuscript as it is rendered: tokens, citations, labels and the registry table.

These are the gates the PDF render enforces indirectly. Each one is checked here
so that a defect is reported by the suite rather than by a failed LaTeX run.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from oph_fep_synergetics.analysis.registry import CORRESPONDENCES
from oph_fep_synergetics.manuscript_variables import MANUSCRIPT_TOKENS

PROJECT = Path(__file__).resolve().parents[1]
MANUSCRIPT = PROJECT / "manuscript"
SECTIONS = sorted(MANUSCRIPT.glob("[0-9]*.md"))

#: Tokens that are part of the published value set but are not inlined in the
#: prose. Each entry states why, and the set is asserted to be exact, so a token
#: cannot quietly stop being quoted.
UNUSED_TOKENS: dict[str, str] = {}

CODE_SPAN = re.compile(r"`[^`\n]*`")
TOKEN = re.compile(r"\{\{([A-Z][A-Z0-9_]*)\}\}")
CITATION = re.compile(r"@([A-Za-z][A-Za-z0-9_]*)")
REFERENCE_PREFIXES = ("sec", "eq", "fig", "tbl", "def", "prop", "thm", "lem", "cor", "rmk", "exa")
FIGURE_LINK = re.compile(
    r"\]\(\.\./output/figures/(fig_[a-z0-9_]+\.png)\)\{#fig:([a-z0-9_-]+)"
)
SECTION_HEADING = re.compile(r"# .+ \{#sec:[a-z0-9_-]+\}")


def _raw_text() -> str:
    """Every rendered section, exactly as the hydration step receives it."""
    return "\n".join(path.read_text(encoding="utf-8") for path in SECTIONS)


def _text_without_code_spans() -> str:
    """The same text with code spans removed.

    A citation inside backticks documents the syntax rather than resolving to an
    entry, which is how the manuscript shows `[@key]` and `[@sec:…]`.
    """
    return "\n".join(CODE_SPAN.sub("", path.read_text(encoding="utf-8")) for path in SECTIONS)


def _prose_lines(path: Path) -> list[str]:
    """The lines of a section that are not inside a fenced code block."""
    lines: list[str] = []
    in_fence = False
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("```"):
            in_fence = not in_fence
            continue
        if not in_fence:
            lines.append(line)
    return lines


def test_the_section_files_are_the_ones_the_renderer_discovers():
    expected = {
        "00_00_graphical_abstract.md",
        "00_abstract.md",
        "01_introduction.md",
        "02_framework.md",
        "03_method.md",
        "04_correspondences.md",
        "05_disagreements.md",
        "06_complementarities.md",
        "07_scope.md",
        "08_reproducibility.md",
        "99_references.md",
    }
    assert {path.name for path in SECTIONS} == expected


def test_every_token_resolves_to_a_bound_value():
    used = set(TOKEN.findall(_raw_text()))
    assert used, "the manuscript quotes no bound value"
    assert used <= set(MANUSCRIPT_TOKENS)


def test_the_published_value_set_is_used_except_where_documented():
    used = set(TOKEN.findall(_raw_text()))
    assert set(MANUSCRIPT_TOKENS) - used - set(UNUSED_TOKENS) == set()
    assert set(UNUSED_TOKENS) <= set(MANUSCRIPT_TOKENS)


def test_every_citation_resolves_and_every_entry_is_cited():
    cited = {
        key for key in CITATION.findall(_text_without_code_spans()) if key not in REFERENCE_PREFIXES
    }
    bib = set(
        re.findall(
            r"^@[a-zA-Z]+\{([^,]+),",
            (MANUSCRIPT / "references.bib").read_text(encoding="utf-8"),
            re.M,
        )
    )
    assert cited == bib


def test_every_top_level_heading_carries_a_section_label():
    for path in SECTIONS:
        for line in _prose_lines(path):
            if line.startswith("# "):
                assert SECTION_HEADING.fullmatch(line), f"{path.name}: {line}"


def test_the_manuscript_status_notes_are_not_a_rendered_section():
    # The renderer treats any top-level *.md as a section, so the outline lives
    # under a name the discovery step excludes.
    assert (MANUSCRIPT / "MANUSCRIPT_STATUS.md").is_file()
    assert not (MANUSCRIPT / "OUTLINE.md").exists()


def test_config_is_token_free_and_pdf_only():
    text = CODE_SPAN.sub("", (MANUSCRIPT / "config.yaml").read_text(encoding="utf-8"))
    assert not TOKEN.search(text)
    assert "pdf: true" in text
    for fmt in ("html", "slides", "docx", "epub"):
        assert f"{fmt}: false" in text


def test_the_registry_table_lists_exactly_the_registered_keys():
    lines = (MANUSCRIPT / "03_method.md").read_text(encoding="utf-8").splitlines()
    header = next(
        index for index, line in enumerate(lines) if line.startswith("| Key | Kind | Status |")
    )
    rows: list[str] = []
    for line in lines[header + 2 :]:
        if not line.startswith("|"):
            break
        rows.append(line.split("|")[1].strip().strip("`"))
    assert rows == [entry.key for entry in CORRESPONDENCES]


def test_figure_links_match_the_accessibility_registry():
    links: dict[str, str] = {}
    for path in SECTIONS:
        for match in FIGURE_LINK.finditer(path.read_text(encoding="utf-8")):
            links[match.group(2)] = match.group(1)
    assert len(links) == 18, sorted(links)
    registry = PROJECT / "output" / "figures" / "figure_registry.json"
    if not registry.is_file():  # the figures are generated artefacts
        return
    payload = json.loads(registry.read_text(encoding="utf-8"))
    entries = {entry["label"]: entry for entry in payload["figures"]}
    linked = {f"fig:{key}": name for key, name in links.items()}
    # The cover is the one registry entry no section links: `config.yaml` places
    # it as the page-1 cover art through `paper.cover.image`, and its alt text
    # lives in the same registry so that the tagged render and the EPUB cover
    # share one source.
    registry_only = {"fig:cover": "fig_cover.png"}
    assert {label: entry["filename"] for label, entry in entries.items()} == linked | registry_only
    for label in [*linked, *registry_only]:
        assert entries[label]["alt"], label
    for filename in [*linked.values(), *registry_only.values()]:
        assert (registry.parent / filename).stat().st_size > 10_000, filename
