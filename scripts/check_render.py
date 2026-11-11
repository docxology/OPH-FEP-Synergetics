#!/usr/bin/env python3
"""Check a rendered PDF against the front-matter and reference contract.

The manuscript makes three promises that only a rendered artefact can keep:

1. page 1 is the cover, pages 2 to 5 are the four panels of the graphical
   abstract, and the contents page follows the front matter;
2. every registry key appears in the printed text, which is what shows that the
   registry table is the registry rather than a summary of it;
3. LaTeX resolved every reference and the log carries no overfull box.

The first two are read from the PDF with ``pdftotext -layout``; the third is
read from the LaTeX log. Nothing here computes project results: it is a check on
the render, and it exits non-zero when a promise is broken.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oph_fep_synergetics.analysis.registry import CORRESPONDENCES  # noqa: E402

PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_PDF = PROJECT / "output" / "pdf" / "OPH-FEP-Synergetics_combined.pdf"
DEFAULT_LOG = PROJECT / "output" / "pdf" / "_combined_manuscript.log"
DEFAULT_TEX = PROJECT / "output" / "pdf" / "_combined_manuscript.tex"

#: The panels of the graphical abstract, in the order they must occupy pages.
PANELS = ("Panel 1", "Panel 2", "Panel 3", "Panel 4")

#: The physical pages the four panels must occupy.
PANEL_PAGES = (2, 3, 4, 5)

OVERFULL = re.compile(r"Overfull \\hbox")
UNDEFINED_CITATION = re.compile(r"Citation [`'][^']+'[^\n]*undefined|There were undefined citations")
UNRESOLVED_REFERENCE = re.compile(r"\?\?")
TOKEN = re.compile(r"\{\{[A-Z][A-Z0-9_]*\}\}")


def pdf_text(pdf: Path) -> list[str]:
    """Return the layout text of each physical page, one entry per page."""
    binary = shutil.which("pdftotext")
    if binary is None:
        raise SystemExit("pdftotext is not on PATH; install poppler to check a render")
    completed = subprocess.run(
        [binary, "-layout", str(pdf), "-"],
        capture_output=True,
        check=True,
        text=True,
    )
    return completed.stdout.split("\f")[:-1]


def panel_page_map(pages: list[str]) -> dict[str, int]:
    """Return the physical page each panel label first appears on."""
    found: dict[str, int] = {}
    for panel in PANELS:
        for number, text in enumerate(pages, start=1):
            if panel in text:
                found[panel] = number
                break
    return found


def contents_page(pages: list[str]) -> int | None:
    """Return the physical page holding the contents heading, if any."""
    for number, text in enumerate(pages, start=1):
        if re.search(r"^\s*Contents\s*$", text, re.M):
            return number
    return None


def missing_keys(text: str, keys: tuple[str, ...]) -> list[str]:
    """Return the registry keys that do not appear in the rendered text."""
    return [key for key in keys if key not in text]


def log_findings(log_text: str) -> dict[str, object]:
    """Return the LaTeX log findings a render must not produce."""
    return {
        "overfull": len(OVERFULL.findall(log_text)),
        "undefined_citations": bool(UNDEFINED_CITATION.search(log_text)),
        "undefined_references": bool(re.search(r"LaTeX Warning: Reference", log_text)),
        "unresolved_tokens": TOKEN.findall(log_text),
    }


def check(
    pdf: Path = DEFAULT_PDF,
    log: Path = DEFAULT_LOG,
    *,
    keys: tuple[str, ...] = tuple(entry.key for entry in CORRESPONDENCES),
) -> tuple[bool, list[str]]:
    """Run every check and return (passed, report lines)."""
    lines: list[str] = []
    passed = True

    if not pdf.is_file():
        return False, [f"missing render: {pdf}"]
    pages = pdf_text(pdf)
    lines.append(f"pdf: {pdf} ({len(pages)} physical pages)")

    cover = pages[0] if pages else ""
    cover_ok = "Observers, Blankets, and Screens" in cover
    lines.append(f"cover on page 1: {'ok' if cover_ok else 'FAIL'}")
    passed &= cover_ok

    found = panel_page_map(pages)
    for panel, expected in zip(PANELS, PANEL_PAGES, strict=True):
        actual = found.get(panel)
        ok = actual == expected
        passed &= ok
        lines.append(f"{panel} on page {expected}: {'ok' if ok else f'FAIL (page {actual})'}")

    contents = contents_page(pages)
    contents_ok = contents is not None and contents > PANEL_PAGES[-1]
    lines.append(
        f"contents after the panels: {'ok (page %d)' % contents if contents_ok else 'FAIL'}"
    )
    passed &= contents_ok

    text = "\n".join(pages)
    absent = missing_keys(text, keys)
    lines.append(f"registry keys printed: {len(keys) - len(absent)}/{len(keys)}")
    if absent:
        lines.append("  absent: " + ", ".join(absent))
    passed &= not absent

    if log.is_file():
        findings = log_findings(log.read_text(encoding="utf-8", errors="replace"))
        lines.append(
            "log: {overfull} overfull hbox, undefined citations {undefined_citations}, "
            "undefined references {undefined_references}, unresolved tokens {unresolved_tokens}".format(
                **findings
            )
        )
        passed &= findings["overfull"] == 0
        passed &= not findings["undefined_citations"]
        passed &= not findings["undefined_references"]
        passed &= not findings["unresolved_tokens"]
    else:
        lines.append(f"log: not found at {log}")
        passed = False

    return passed, lines


def main(argv: list[str] | None = None) -> int:
    """Run the checks and print the report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, default=DEFAULT_PDF)
    parser.add_argument("--log", type=Path, default=DEFAULT_LOG)
    arguments = parser.parse_args(argv)
    passed, lines = check(arguments.pdf, arguments.log)
    for line in lines:
        print(line)
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
