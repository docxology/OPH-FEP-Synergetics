"""The figure generator: coverage, the accessibility registry, and a render pass.

The generator is a thin orchestrator over the package's own functions, so the
figures cannot disagree with the probes. The checks here cover the properties the
suite can assert without a human looking at the images: every builder exists,
the labels and file names agree with the manuscript, and a sample of builders
really writes non-trivial PNG files.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

pytest.importorskip("matplotlib", reason="the figures extra is not installed")

import scripts.make_figures as make_figures  # noqa: E402

EXPECTED_BUILDERS = (
    "carrier_figure",
    "registry_figure",
    "descent_figure",
    "dstability_figure",
    "mediator_figure",
    "blanket_figure",
    "orbit_figure",
    "ivm_figure",
    "subdivision_figure",
    "cover_figure",
    "ga_triangle_figure",
    "ga_carrier_figure",
    "ga_descent_figure",
    "ga_synergetics_figure",
    "pair_triangles_figure",
    "scope_matrix_figure",
    "network_figure",
    "shell_gauge_figure",
)


def test_the_generator_has_every_builder_and_orchestrates_them():
    source = (Path(make_figures.__file__)).read_text(encoding="utf-8")
    for name in EXPECTED_BUILDERS:
        assert f"def {name}(" in source
    main_body = source.split("def main() -> int:")[1]
    for name in EXPECTED_BUILDERS:
        assert name in main_body


def test_a_sample_of_builders_writes_figures_and_the_registry(tmp_path, monkeypatch):
    monkeypatch.setattr(make_figures, "FIGURES", tmp_path)
    monkeypatch.setattr(make_figures, "ALT_TEXT", {})

    written = [
        make_figures.registry_figure(),
        make_figures.descent_figure(),
        make_figures.ivm_figure(),
    ]
    for path in written:
        assert path.parent == tmp_path
        assert path.stat().st_size > 10_000

    payload = json.loads(make_figures.write_registry().read_text(encoding="utf-8"))
    assert payload["schema_version"]
    assert payload["generator"] == "scripts/make_figures.py"
    assert {entry["label"] for entry in payload["figures"]} == {
        "fig:registry",
        "fig:descent",
        "fig:ivm",
    }
    for entry in payload["figures"]:
        assert entry["filename"].endswith(".png")
        assert len(entry["alt"]) > 60
        assert entry["alt"].endswith(".")


def test_the_same_builder_is_byte_deterministic(tmp_path, monkeypatch):
    """A rebuilt figure is byte-identical, so the delivered copy stays identical."""
    monkeypatch.setattr(make_figures, "FIGURES", tmp_path)
    monkeypatch.setattr(make_figures, "ALT_TEXT", {})
    first = make_figures.descent_figure().read_bytes()
    second = make_figures.descent_figure().read_bytes()
    assert second == first
