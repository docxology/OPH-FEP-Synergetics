"""The public surface of each layer stays importable and provenance-consistent."""

from __future__ import annotations

import collections
import re
from pathlib import Path

import pytest

import oph_fep_synergetics
from oph_fep_synergetics import analysis, core, geometry, manuscript_variables
from oph_fep_synergetics.analysis import DEFAULT_SOURCES

LAYERS = (core, geometry, analysis, manuscript_variables)


def test_version_is_declared():
    assert isinstance(oph_fep_synergetics.__version__, str)
    assert oph_fep_synergetics.__version__.count(".") == 2


@pytest.mark.parametrize("layer", LAYERS)
def test_every_exported_name_resolves(layer):
    assert layer.__all__
    for name in layer.__all__:
        assert getattr(layer, name, None) is not None, f"{layer.__name__}.{name} is missing"


@pytest.mark.parametrize("layer", LAYERS)
def test_exported_names_are_unique(layer):
    assert len(set(layer.__all__)) == len(layer.__all__)


def test_core_exposes_the_correspondence_surface():
    for name in (
        "IcosahedralCarrier",
        "OPH_NEIGHBOURS",
        "OPH_ORIENTED_BASE_FACES",
        "RepairStep",
        "CellularSheaf",
        "GaussianSheafModel",
        "carrier_isomorphism",
        "gated_descent_census",
        "neighbours_from_faces",
        "shell_gap",
    ):
        assert hasattr(core, name), name


def test_geometry_exposes_the_synergetic_surface():
    for name in (
        "Quadray",
        "VOLUME_RATIOS",
        "goldberg_vertex_count",
        "neighbor_moves",
        "shell_count",
        "tetra_volume",
    ):
        assert hasattr(geometry, name), name


def test_analysis_exposes_the_registry_surface():
    for name in ("CORRESPONDENCES", "by_key", "check_pins", "load_sources", "summary", "verify"):
        assert hasattr(analysis, name), name


def test_manuscript_variables_expose_the_values_surface():
    for name in ("MANUSCRIPT_VALUES", "BoundValue", "as_dict"):
        assert hasattr(manuscript_variables, name), name


def test_provenance_hashes_are_pinned_in_one_place():
    """A commit cited in code must also appear in the committed pin file."""
    if not DEFAULT_SOURCES.is_file():
        pytest.skip("the pin file is not part of this checkout")
    pinned = DEFAULT_SOURCES.read_text(encoding="utf-8")
    package_root = Path(oph_fep_synergetics.__file__).parent
    cited = collections.Counter(
        digest
        for path in sorted(package_root.rglob("*.py"))
        for digest in re.findall(r"\b[0-9a-f]{40}\b", path.read_text(encoding="utf-8"))
    )
    assert cited, "the OPH carrier faces should cite their source commit"
    assert [digest for digest in cited if digest not in pinned] == []
