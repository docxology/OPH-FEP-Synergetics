"""The rendered PDF contract, checked on the pure parsing helpers and on the real render.

The release PDF is tracked under `output/`. The parsing helpers also use
synthetic page text; the LaTeX log is available only after a local rebuild.
"""

from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import pytest

PROJECT = Path(__file__).resolve().parents[1]
PDF = PROJECT / "output" / "pdf" / "OPH-FEP-Synergetics_combined.pdf"


def _load(script: str):
    spec = importlib.util.spec_from_file_location(script, PROJECT / "scripts" / f"{script}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


check_render = _load("check_render")


def _pages(*texts: str) -> list[str]:
    return list(texts)


def test_panel_pages_are_read_in_physical_order():
    pages = _pages(
        "Observers, Blankets, and Screens\nCover page",
        "1.1 Panel 1 — the pairwise ledger",
        "1.2 Panel 2 — the carrier",
        "1.3 Panel 3 — the dynamics",
        "1.4 Panel 4 — the accounting",
        "1.5 Panel 5 — Fuller's Jitterbug and the OPH carrier\nFigure 5: poster",
        "Contents\n1 Graphical abstract 2",
    )
    assert check_render.panel_page_map(pages) == {
        "Panel 1": 2,
        "Panel 2": 3,
        "Panel 3": 4,
        "Panel 4": 5,
        "Panel 5": 6,
    }
    assert check_render.contents_page(pages) == 7


def test_a_panel_that_moves_is_reported_by_page():
    pages = _pages("cover", "Panel 2", "Panel 1", "Panel 3", "Panel 4", "Contents")
    found = check_render.panel_page_map(pages)
    assert found["Panel 1"] == 3
    assert found["Panel 2"] == 2


def test_a_missing_panel_is_absent_rather_than_guessed():
    pages = _pages("cover", "Panel 1", "Panel 2", "Panel 3", "Contents")
    assert "Panel 4" not in check_render.panel_page_map(pages)


def test_contents_must_not_be_the_cover_and_is_absent_when_unprinted():
    assert check_render.contents_page(_pages("cover", "Panel 1")) is None
    assert check_render.contents_page(_pages("cover", "Table of Contents\n")) is None


def test_caption_pages_count_each_page_separately():
    pages = _pages(
        "cover with no caption",
        "Panel 1\nFigure 1: the ledger triangle\nmore text",
        "Panel 2\nFigure 2: the carrier",
        "Figure 3: a float that left its panel",
        "Panel 3",
        "Panel 4",
    )
    counts = check_render.caption_pages(pages)
    assert counts[1] == 0
    assert counts[2] == 1
    assert counts[4] == 1
    assert sum(counts.values()) == 3


def test_caption_inventory_requires_every_number_once():
    pages = _pages("  Figure 1: earlier panel\nTable 2: second table", "Figure 1: duplicate\n  Table 1: first table")
    assert check_render.caption_inventory(pages, "Table", 2) == (True, "table captions: 2/2")
    assert check_render.caption_inventory(pages, "Figure", 2) == (
        False, "figure captions: 2/2 (missing [2]; repeated [1])",
    )


def test_missing_keys_are_reported_in_registry_order():
    absent = check_render.missing_keys(
        "carrier-port-chart and spectral-\n  gap\n\ntree-num\fber",
        ("carrier-port-chart", "spectral-gap", "tree-number"),
    )
    assert absent == ["tree-number"]


def test_log_findings_catch_each_failure_mode():
    findings = check_render.log_findings(
        "Overfull \\hbox (12.0pt too wide) in paragraph\n"
        "Overfull \\vbox (1.0pt too high) detected\n"
        "LaTeX Warning: Citation `oph' on page 3 undefined on input line 40.\n"
        "LaTeX Warning: Reference `tbl:registry' on page 5 undefined on input line 60.\n"
        "the value {{CARRIER_PORTS}} was not resolved\n"
    )
    assert findings["overfull"] == 2
    assert findings["undefined_citations"] is True
    assert findings["undefined_references"] is True
    assert findings["unresolved_tokens"] == ["{{CARRIER_PORTS}}"]
    assert check_render.log_findings("LaTeX Warning: There were undefined references.")[
        "undefined_references"
    ]


def test_a_clean_log_reports_nothing():
    findings = check_render.log_findings("Output written on combined.pdf (40 pages).\n")
    assert findings == {
        "overfull": 0,
        "undefined_citations": False,
        "undefined_references": False,
        "unresolved_tokens": [],
    }


def test_missing_explicit_template_root_does_not_fall_back(tmp_path, monkeypatch, capsys):
    build = _load("build")
    fallback = tmp_path / "fallback"
    for name in (".venv/bin/python", "scripts/pipeline/stage_03_render.py"):
        file = fallback / name
        file.parent.mkdir(parents=True, exist_ok=True)
        file.touch()
    monkeypatch.setenv("TEMPLATE_ROOT", str(fallback))

    with pytest.raises(SystemExit) as exit_info:
        build.main(["--template-root", str(tmp_path / "missing")])
    assert exit_info.value.code == 2
    assert "no usable template checkout" in capsys.readouterr().err


def test_standalone_pdf_rejects_printed_reference_markers(tmp_path, monkeypatch):
    pdf = tmp_path / "test.pdf"
    pdf.touch()
    pages = [
        "Observers, Blankets, and Screens",
        *(f"Panel {index}\nFigure {index}: panel" for index in range(1, 6)),
        "Contents",
        "See Proposition ?? and [@def:missing]. Citation [?]. Raw label {#tbl:commands}. Is this punctuation??",
        "\n".join(f"Figure {index}: caption" for index in range(6, 17)),
    ]
    monkeypatch.setattr(check_render, "pdf_text", lambda _: pages)
    monkeypatch.setattr(check_render, "text_left_margin", lambda *_: 36.0)

    passed, lines = check_render.check(pdf, log=None, keys=())

    assert not passed
    assert "printed unresolved references or citations: 4" in lines


@pytest.mark.skipif(not PDF.is_file(), reason="no rendered PDF in this checkout")
def test_the_rendered_pdf_keeps_the_front_matter_and_ontology_supplement():
    passed, lines = check_render.check(PDF, log=None)
    assert passed, "\n".join(lines)
    pages = check_render.pdf_text(PDF)
    poster_page = pages[5]
    assert "Jitterbug and the OPH carrier" in poster_page
    assert "geometric hull diagonals" in " ".join(poster_page.split())
    assert "accelerated clock" in " ".join(poster_page.split())
    cover = " ".join(pages[0].split())
    first_author = cover.index("Daniel Ari Friedman")
    second_author = cover.index("Bernhard Mueller")
    assert first_author < second_author
    assert "Active Inference Institute" in cover[first_author:second_author]
    assert "Pragma Research, Active Inference Institute" in cover[second_author:]
    assert "0009-0002-5638-6240" in cover[second_author:]
    supplement_page = next(
        index for index, page in enumerate(pages)
        if re.search(r"(?m)^\s*\d+\s+Ontology concordance\s*$", page)
    )
    references_page = next(
        index for index, page in enumerate(pages)
        if re.search(r"(?m)^\s*\d+\s+References\s*$", page)
    )
    assert supplement_page < references_page
    ontology = json.loads((PROJECT / "data" / "ontology_concordance.json").read_text(encoding="utf-8"))
    group_pages = [
        next(index for index in range(supplement_page, references_page)
             if re.search(rf"(?m)^\s*\d+\.\d+\s+{re.escape(group['title'])}\s*$", pages[index]))
        for group in ontology["groups"]
    ]
    assert group_pages == sorted(set(group_pages))
