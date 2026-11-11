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


#: Literal numbers that are structure, citation metadata or layout rather than
#: claim values. Each pattern says what it excuses; any other number in the
#: manuscript must arrive as a {{TOKEN}} bound to a BoundValue.
LITERAL_ALLOWLIST = tuple(
    re.compile(pattern)
    for pattern in (
        r"@[\w:-]+",                          # citations and cross-references: [@sec:pr01]
        r"width=\d+%",                        # figure layout
        r"\\setcounter\{page\}\{\d+\}",       # page plumbing
        r"\bPR-01\b",                         # OPH's own premise label
        r"\bA3\b",                            # the root lattice's name
        r"Panel \d",                          # panel headings
        r"^\s*\d+\.\s",                       # ordered-list markers
        r"Items \d+ to \d+",                  # the contributions cross-reference
        r"\|\s*\d+\s*\|",                     # the adjacency table's port-index columns
        r"port\s+0\s+meets\s+ports\s+1,2,3,4,6",  # the common-labelling example's port labels
        r"2, ?1, ?1, ?0",                     # the neighbour-move generator tuple (2,1,1,0)
        r"\{#sec:[a-z0-9_-]+\}",              # section labels in headings
        r"\\binom\{13\}\{11\}",               # the shell-size formula's arguments
        r"\bCC BY-NC-SA 4\.0\b",
        r"\bCC BY 4\.0\b",
        r"\bApache-2\.0\b",
        r"\bMIT\b",
        r"10\.5281/zenodo\.\d+",              # DOIs
        r"schema-\d",                         # FEP_Formal's roster schema version
    )
)

#: A math span that carries a letter or a command is a formula, not a value; a
#: digits-only span is a claim and must be a token instead.
def _numeric_fragments() -> list[str]:
    """Every numeric literal outside tokens, code, formulas and the allowlist."""
    fragments: list[str] = []
    for path in SECTIONS:
        fenced = False
        display_math = False
        math_carry = False
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            stripped = line.strip()
            if stripped.startswith("```"):
                fenced = not fenced
                continue
            if fenced:
                continue
            if display_math:
                if "$$" in stripped:
                    display_math = False
                continue
            if stripped.startswith("$$"):
                display_math = True
                continue
            if math_carry:
                if line.count("$") % 2 == 1:
                    math_carry = False
                continue
            text = CODE_SPAN.sub("", line)
            if text.count("$") % 2 == 1:
                # an unclosed span: its tail is a formula continuing next line
                head, sep, _ = text.rpartition("$")
                text = head + sep
                math_carry = True
            pieces: list[str] = []
            for part in re.split(r"(\$[^$]*\$)", text):
                if len(part) > 1 and part.startswith("$") and part.endswith("$"):
                    if re.search(r"[A-Za-z\\]", part[1:-1]):
                        continue
                    pieces.append(f" {part[1:-1]} ")
                else:
                    pieces.append(part)
            text = TOKEN.sub(" ", "".join(pieces))
            for pattern in LITERAL_ALLOWLIST:
                text = pattern.sub(" ", text)
            for match in re.finditer(r"\d+(?:\.\d+)?", text):
                fragments.append(
                    f"{path.name}:{number}: {match.group(0)} :: {stripped[:70]}"
                )
    return fragments


def test_no_hard_coded_numbers_in_the_manuscript():
    """Every number the prose quotes arrives as a bound token; nothing is pasted."""
    fragments = _numeric_fragments()
    assert fragments == [], "\n".join(fragments[:25])


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
    assert len(links) == 19, sorted(links)
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
