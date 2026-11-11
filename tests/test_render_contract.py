"""The rendered PDF contract, checked on the pure parsing helpers and on the real render.

The render itself is an artefact: it lives under `output/`, which is generated
and never committed. The parsing helpers are therefore tested against synthetic
page text, and the whole contract is exercised against the real PDF only when one
has been rendered.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

PROJECT = Path(__file__).resolve().parents[1]
PDF = PROJECT / "output" / "pdf" / "OPH-FEP-Synergetics_combined.pdf"
LOG = PROJECT / "output" / "pdf" / "_combined_manuscript.log"


def _load():
    spec = importlib.util.spec_from_file_location(
        "check_render", PROJECT / "scripts" / "check_render.py"
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


check_render = _load()


def _pages(*texts: str) -> list[str]:
    return list(texts)


def test_panel_pages_are_read_in_physical_order():
    pages = _pages(
        "Observers, Blankets, and Screens\nCover page",
        "1.1 Panel 1 — the pairwise ledger",
        "1.2 Panel 2 — the carrier",
        "1.3 Panel 3 — the dynamics",
        "1.4 Panel 4 — the accounting",
        "Contents\n1 Graphical abstract 2",
    )
    assert check_render.panel_page_map(pages) == {
        "Panel 1": 2,
        "Panel 2": 3,
        "Panel 3": 4,
        "Panel 4": 5,
    }
    assert check_render.contents_page(pages) == 6


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


def test_missing_keys_are_reported_in_registry_order():
    absent = check_render.missing_keys(
        "carrier-port-chart and spectral-gap", ("carrier-port-chart", "spectral-gap", "tree-number")
    )
    assert absent == ["tree-number"]


def test_log_findings_catch_each_failure_mode():
    findings = check_render.log_findings(
        "Overfull \\hbox (12.0pt too wide) in paragraph\n"
        "LaTeX Warning: Citation `oph' on page 3 undefined on input line 40.\n"
        "LaTeX Warning: Reference `tbl:registry' on page 5 undefined on input line 60.\n"
        "the value {{CARRIER_PORTS}} was not resolved\n"
    )
    assert findings["overfull"] == 1
    assert findings["undefined_citations"] is True
    assert findings["undefined_references"] is True
    assert findings["unresolved_tokens"] == ["{{CARRIER_PORTS}}"]


def test_a_clean_log_reports_nothing():
    findings = check_render.log_findings("Output written on combined.pdf (40 pages).\n")
    assert findings == {
        "overfull": 0,
        "undefined_citations": False,
        "undefined_references": False,
        "unresolved_tokens": [],
    }


@pytest.mark.skipif(not PDF.is_file(), reason="no rendered PDF in this checkout")
def test_the_rendered_pdf_keeps_the_front_matter_contract():
    passed, lines = check_render.check(PDF, LOG)
    assert passed, "\n".join(lines)
