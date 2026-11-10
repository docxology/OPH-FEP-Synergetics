"""The source-bound values the manuscript may quote, checked against the code."""

from __future__ import annotations

import pytest

from oph_fep_synergetics.core import (
    affine_path_contrast_float,
    all_shell_configurations,
    binary_shell_configurations,
    d_stability_free_energy_form,
    gated_descent_census,
    icosahedral_carrier,
    path_multiplicity_census,
    path_multiplicity_values,
    routh_first_column,
    shell_gap,
)
from oph_fep_synergetics.geometry import (
    VOLUME_RATIOS,
    cumulative_count,
    neighbor_moves,
    shell_count,
)
from oph_fep_synergetics.manuscript_variables import MANUSCRIPT_VALUES, as_dict
from oph_fep_synergetics.manuscript_variables.values import _value

EXPECTED_KEYS = {
    "carrier.ports",
    "carrier.seams",
    "carrier.faces",
    "carrier.degree",
    "carrier.lambda1",
    "carrier.shell_gap",
    "shell.total_two",
    "shell.binary",
    "shell.transitions",
    "shell.strict_descent",
    "shell.zero_descent",
    "shell.noop",
    "shell.multiplicities",
    "schedule.affine_contrast",
    "d_stability.routh_k4",
    "d_stability.lyapunov_determinant",
    "d_stability.routh_onset",
    "two_node.relaxation_min_radius",
    "two_node.relaxation_period",
    "two_node.relaxation_min_radius_measured",
    "two_node.relaxation_max_radius",
    "two_node.relaxation_period_measured",
    "two_node.sync_min_radius",
    "ivm.neighbor_moves",
    "ivm.neighbor_squared_distance",
    "ivm.shell_counts",
    "ivm.cumulative_counts",
    "ivm.volume_ratios",
    "mediator.threshold",
    "directed_gradient.direct",
    "directed_gradient.adjoint",
    "fep_formal.topics",
    "paper.five_surprises",
    "registry.entries",
    "registry.isomorphisms",
    "registry.disagreements",
    "registry.complementarities",
}


def test_the_bound_value_set_is_exactly_the_documented_one():
    keys = {entry.key for entry in MANUSCRIPT_VALUES}
    assert keys == EXPECTED_KEYS
    assert len(MANUSCRIPT_VALUES) == len(keys)


def test_every_value_names_its_source():
    for entry in MANUSCRIPT_VALUES:
        assert entry.source
        assert entry.value is not None
        assert entry.exact is None or isinstance(entry.exact, str)
        assert entry.note is None or isinstance(entry.note, str)


def test_as_dict_projects_the_same_values():
    assert as_dict() == {entry.key: entry.value for entry in MANUSCRIPT_VALUES}


def test_value_lookup_by_key():
    assert _value("carrier.ports").value == 12
    with pytest.raises(KeyError):
        _value("carrier.not-a-value")


def test_carrier_values_match_the_executed_code():
    carrier = icosahedral_carrier()
    assert _value("carrier.ports").value == carrier.port_count
    assert _value("carrier.seams").value == carrier.seam_count
    assert _value("carrier.faces").value == carrier.face_count
    assert _value("carrier.degree").value == carrier.degree(0)
    assert _value("carrier.lambda1").value == pytest.approx(carrier.spectral_gap(), abs=1e-12)
    assert _value("carrier.shell_gap").value == pytest.approx(shell_gap(), abs=1e-18)


def test_shell_values_match_the_exhaustive_census():
    census = gated_descent_census()
    assert _value("shell.total_two").value == len(all_shell_configurations()) == 78
    assert _value("shell.binary").value == len(binary_shell_configurations()) == 66
    assert _value("shell.transitions").value == census["states"] * 30 * 2
    assert _value("shell.strict_descent").value == census["strict"]
    assert _value("shell.zero_descent").value == census["zero"]
    assert _value("shell.noop").value == census["noop"]
    assert tuple(_value("shell.multiplicities").value) == path_multiplicity_values()
    assert tuple(_value("shell.multiplicities").value) == tuple(sorted(path_multiplicity_census()))


def test_control_values_match_the_executed_controls():
    assert _value("schedule.affine_contrast").value == pytest.approx(
        affine_path_contrast_float(steps=12), abs=1e-18
    )
    assert tuple(_value("d_stability.routh_k4").value) == routh_first_column(4.0)
    assert _value("d_stability.lyapunov_determinant").value == pytest.approx(
        d_stability_free_energy_form()["determinant"], abs=1e-12
    )
    assert _value("mediator.threshold").value == 1.0 / 3.0
    assert tuple(_value("directed_gradient.direct").value) == (1.0, 0.0)
    assert tuple(_value("directed_gradient.adjoint").value) == (1.0, 2.0)


def test_orbit_values_separate_the_source_bound_from_the_measured():
    radius = _value("two_node.relaxation_min_radius")
    period = _value("two_node.relaxation_period")
    assert radius.value == 0.84
    assert period.value == 7.58
    assert radius.note and "measured" in radius.note
    assert period.note and "measured" in period.note
    assert _value("two_node.sync_min_radius").value == 1.0
    assert "different object" in _value("two_node.sync_min_radius").note


def test_ivm_values_match_the_geometry_module():
    assert _value("ivm.neighbor_moves").value == len(neighbor_moves())
    assert _value("ivm.neighbor_squared_distance").value == 8
    assert tuple(_value("ivm.shell_counts").value) == tuple(shell_count(k) for k in range(5))
    assert tuple(_value("ivm.cumulative_counts").value) == tuple(
        cumulative_count(k) for k in range(5)
    )
    assert _value("ivm.volume_ratios").value == {
        str(name): int(value) for name, value in VOLUME_RATIOS.items()
    }


def test_paper_values_are_the_cited_ones():
    assert _value("fep_formal.topics").value == 168
    surprises = _value("paper.five_surprises")
    assert tuple(surprises.value) == (
        "numerical",
        "bayesian",
        "surprisal",
        "critic",
        "meaningful",
    )
    assert surprises.source.endswith("Table 7")
