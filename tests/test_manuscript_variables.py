"""The source-bound values the manuscript may quote, checked against the code."""

from __future__ import annotations

import math

import pytest

from oph_fep_synergetics.analysis.registry import CORRESPONDENCES, PAIRS, by_kind, by_pair
from oph_fep_synergetics.core import (
    affine_path_contrast_float,
    all_shell_configurations,
    binary_shell_configurations,
    cohomology_dimensions,
    d_stability_free_energy_form,
    gated_descent_census,
    icosahedral_carrier,
    invariants as invariant_module,
    path_multiplicity_census,
    path_multiplicity_values,
    routh_first_column,
    shell_gap,
    symmetry as symmetry_module,
)
from oph_fep_synergetics.core.sheaf import carrier_sheaf
from oph_fep_synergetics.geometry import (
    VOLUME_RATIOS,
    cuboctahedron_adjacency,
    cumulative_count,
    generate_shell,
    icosahedron_tetravolume,
    neighbor_moves,
    packing_density,
    shell_count,
    squared_norm,
    voronoi_cell,
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
    "registry.oph_fep",
    "registry.fep_synergetics",
    "registry.synergetics_oph",
    "registry.three_way",
    "carrier.tree_number",
    "gmrf.marginal_variance",
    "gmrf.blanket_size",
    "gmrf.pairwise_precision",
    "sheaf.coboundary_rank",
    "sheaf.h0",
    "sheaf.h1",
    "sheaf.euler_characteristic",
    "symmetry.group_order",
    "symmetry.rotations",
    "symmetry.reflections",
    "symmetry.configurations",
    "symmetry.orbits",
    "symmetry.outcome_multisets",
    "icosahedron.tetravolume",
    "icosahedron.tetravolume_squared",
    "icosahedron.faces_reoriented",
    "icosahedron.degree",
    "cuboctahedron.degree",
    "cuboctahedron.tetravolume",
    "tetravolume.integral_solids",
    "voronoi.vertices",
    "voronoi.faces",
    "voronoi.edges",
    "packing.density",
    "packing.cell_tetravolume",
    "packing.ball_tetravolume",
    "shell.free_energy_level_counts",
    "shell.shared_free_energy_levels",
    "kl.one_step",
    "kl.step_squared_distance",
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


def test_registry_counts_are_derived_from_the_registry():
    assert _value("registry.entries").value == len(CORRESPONDENCES)
    kinds = by_kind()
    pairs = by_pair()
    for kind, key in (
        ("isomorphism", "registry.isomorphisms"),
        ("disagreement", "registry.disagreements"),
        ("complementarity", "registry.complementarities"),
    ):
        assert _value(key).value == len(kinds[kind])
    for pair in PAIRS:
        assert _value("registry." + pair.replace("-", "_")).value == len(pairs[pair])
    assert _value("registry.entries").value == sum(len(pairs[pair]) for pair in PAIRS)
    assert sum(
        _value(key).value
        for key in (
            "registry.isomorphisms",
            "registry.disagreements",
            "registry.complementarities",
        )
    ) == len(CORRESPONDENCES)


def test_invariant_values_match_the_executed_invariants():
    assert _value("carrier.tree_number").value == invariant_module.spanning_tree_number()
    assert _value("carrier.tree_number").value == 5_184_000
    assert invariant_module.tree_number_from_spectrum() == pytest.approx(5_184_000, abs=1e-4)
    assert float(invariant_module.marginal_variance()) == pytest.approx(7 / 36, abs=1e-15)
    assert _value("gmrf.marginal_variance").value == pytest.approx(7 / 36, abs=1e-15)
    assert set(invariant_module.blanket_sizes()) == {_value("gmrf.blanket_size").value} == {5}
    assert set(invariant_module.pairwise_precision().values()) == {
        _value("gmrf.pairwise_precision").value
    } == {-1}


def test_sheaf_values_match_the_executed_chain_complex():
    dimensions = cohomology_dimensions(carrier_sheaf(icosahedral_carrier()))
    assert _value("sheaf.coboundary_rank").value == dimensions["coboundary_rank"] == 11
    assert _value("sheaf.h0").value == dimensions["h0"] == 1
    assert _value("sheaf.h1").value == dimensions["h1"] == 19
    assert _value("sheaf.euler_characteristic").value == dimensions["euler_characteristic"] == -18
    assert dimensions["h0"] - dimensions["h1"] == dimensions["euler_characteristic"]


def test_symmetry_values_match_the_executed_group():
    order = len(symmetry_module.carrier_automorphisms())
    rotations = len(symmetry_module.rotation_subgroup())
    reflections = len(symmetry_module.reflection_subgroup())
    assert _value("symmetry.group_order").value == order == 120
    assert _value("symmetry.rotations").value == rotations == 60
    assert _value("symmetry.reflections").value == reflections == 60
    assert rotations + reflections == order
    equivariance = symmetry_module.law_equivariance()
    assert _value("symmetry.orbits").value == equivariance["orbits"] == 4
    assert _value("symmetry.outcome_multisets").value == equivariance["distinct_outcome_multisets"] == 3
    assert _value("symmetry.configurations").value == len(
        icosahedral_carrier().shell_configurations(total=2)
    ) == 78


def test_geometry_values_match_the_executed_solids_and_cell():
    ico = icosahedron_tetravolume()
    assert _value("icosahedron.tetravolume").value == ico["tetravolume"]
    assert ico["tetravolume"] == pytest.approx(2.5 * (3 * math.sqrt(2) + math.sqrt(10)), abs=1e-12)
    assert _value("icosahedron.tetravolume_squared").value == pytest.approx(175 + 75 * math.sqrt(5), abs=1e-9)
    assert _value("icosahedron.faces_reoriented").value == ico["faces_reoriented"] == 10
    assert _value("icosahedron.degree").value == icosahedral_carrier().degree(0) == 5
    assert _value("cuboctahedron.degree").value == max(sum(row) for row in cuboctahedron_adjacency()) == 4
    assert _value("cuboctahedron.tetravolume").value == float(VOLUME_RATIOS["cuboctahedron"]) == 20.0
    assert all(float(value).is_integer() for value in VOLUME_RATIOS.values())
    assert not float(ico["tetravolume"]).is_integer()
    assert _value("tetravolume.integral_solids").value == tuple(
        int(value) for value in VOLUME_RATIOS.values()
    )


def test_packing_values_match_the_executed_cell_and_density():
    cell = voronoi_cell()
    assert (
        _value("voronoi.vertices").value,
        _value("voronoi.faces").value,
        _value("voronoi.edges").value,
    ) == (cell["vertex_count"], cell["face_count"], cell["edge_count"]) == (14, 12, 24)
    assert cell["vertex_count"] - cell["edge_count"] + cell["face_count"] == 2
    density = packing_density()
    assert _value("packing.density").value == pytest.approx(density["density"], abs=1e-15)
    assert _value("packing.density").value == pytest.approx(math.pi / (3 * math.sqrt(2)), abs=1e-12)
    assert _value("packing.cell_tetravolume").value == pytest.approx(6.0, abs=1e-9)
    assert _value("packing.ball_tetravolume").value == pytest.approx(math.pi * math.sqrt(2), abs=1e-9)


def test_free_energy_and_kl_values_match_the_lattice():
    counts = tuple(
        len({0.5 * float(squared_norm(site)) for site in generate_shell(k)}) for k in range(1, 6)
    )
    assert _value("shell.free_energy_level_counts").value == counts == (1, 3, 4, 7, 7)
    carriage: dict[float, set[int]] = {}
    for k in range(1, 6):
        for site in generate_shell(k):
            carriage.setdefault(0.5 * float(squared_norm(site)), set()).add(k)
    shared = tuple(sorted(level for level, shells in carriage.items() if len(shells) > 1))
    assert _value("shell.shared_free_energy_levels").value == shared == (36.0, 52.0)
    assert _value("kl.step_squared_distance").value == float(squared_norm(neighbor_moves()[0])) == 8.0
    assert _value("kl.one_step").value == 0.5 * _value("kl.step_squared_distance").value == 4.0
