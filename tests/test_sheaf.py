"""Cellular sheaves: agreement energy, harmonic extension and record injectivity."""

from __future__ import annotations

import numpy as np
import pytest

from oph_fep_synergetics.core import (
    CellularSheaf,
    EdgeRestriction,
    carrier_sheaf,
    coboundary_matrix,
    coboundary_rank,
    cohomology_dimensions,
    icosahedral_carrier,
    linear_chain_sheaf,
)


def identity_edge(u: int = 0, v: int = 1, weight: float = 1.0) -> EdgeRestriction:
    return EdgeRestriction(
        u=u, v=v, map_u=np.array([[weight]]), map_v=np.array([[weight]])
    )


def test_chain_laplacian_is_the_path_laplacian():
    sheaf = linear_chain_sheaf(3)
    assert sheaf.node_count == 3
    assert sheaf.total_dim == 3
    assert sheaf.offsets() == (0, 1, 2, 3)
    assert np.allclose(
        sheaf.laplacian(),
        np.array([[1.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 1.0]]),
    )


def test_edge_weight_scales_the_laplacian():
    sheaf = linear_chain_sheaf(2, edge_weight=4.0)
    assert np.allclose(sheaf.laplacian(), np.array([[4.0, -4.0], [-4.0, 4.0]]))


def test_a_chain_needs_at_least_one_node():
    with pytest.raises(ValueError):
        linear_chain_sheaf(0)
    assert linear_chain_sheaf(1).edges == ()


def test_coboundary_is_the_stacked_disagreement():
    sheaf = linear_chain_sheaf(3)
    assert np.allclose(sheaf.coboundary(np.array([2.0, 2.0, 2.0])), [0.0, 0.0])
    assert np.allclose(sheaf.coboundary(np.array([0.0, 1.0, 3.0])), [1.0, 2.0])
    assert np.allclose(linear_chain_sheaf(1).coboundary(np.array([7.0])), np.zeros(0))


def test_energy_is_half_the_squared_disagreement():
    sheaf = linear_chain_sheaf(3)
    vector = np.array([0.0, 1.0, 3.0])
    disagreement = sheaf.coboundary(vector)
    assert sheaf.energy(vector) == 0.5 * float(disagreement @ disagreement)
    assert sheaf.energy(np.array([5.0, 5.0, 5.0])) == 0.0


def test_global_sections_of_a_chain_are_the_constants():
    sheaf = linear_chain_sheaf(4)
    basis = sheaf.global_section_basis()
    assert basis.shape == (4, 1)
    constant = basis[:, 0]
    assert np.allclose(np.abs(constant), np.full(4, 0.5))
    assert sheaf.is_global_section(np.array([3.0, 3.0, 3.0, 3.0]))
    assert not sheaf.is_global_section(np.array([0.0, 1.0, 3.0, 4.0]))


def test_pack_and_unpack_round_trip():
    sheaf = linear_chain_sheaf(3)
    nodes = (np.array([1.0]), np.array([2.0]), np.array([3.0]))
    packed = sheaf.pack(nodes)
    assert np.allclose(packed, [1.0, 2.0, 3.0])
    assert all(np.allclose(a, b) for a, b in zip(sheaf.unpack(packed), nodes))
    with pytest.raises(ValueError):
        sheaf.pack(nodes[:2])
    with pytest.raises(ValueError):
        sheaf.unpack(np.zeros(4))
    with pytest.raises(ValueError):
        sheaf.is_global_section(np.zeros(4))


def test_harmonic_extension_interpolates_between_records():
    sheaf = linear_chain_sheaf(3)
    nodes = sheaf.harmonic_extension({0: [0.0], 2: [2.0]})
    assert np.allclose([float(node[0]) for node in nodes], [0.0, 1.0, 2.0])
    # Clamped data with disagreeing records is never a global section; what the
    # free stalk satisfies is a vanishing gradient, not edge agreement.
    assert not sheaf.is_global_section(sheaf.pack(nodes))
    gradient = 2.0 * float(nodes[1][0]) - 0.0 - 2.0
    assert abs(gradient) < 1e-12


def test_harmonic_extension_matches_the_sheaf_laplacian_solve():
    sheaf = linear_chain_sheaf(5)
    nodes = sheaf.harmonic_extension({0: [1.0], 4: [5.0]})
    free = [1, 2, 3]
    laplacian = sheaf.laplacian()[1:4, 1:4]
    rhs = -(sheaf.laplacian()[1:4, [0, 4]] @ np.array([1.0, 5.0]))
    assert np.allclose([float(node[0]) for node in nodes][1:4], np.linalg.solve(laplacian, rhs))


def test_harmonic_extension_rejects_bad_observations():
    sheaf = linear_chain_sheaf(3)
    with pytest.raises(IndexError):
        sheaf.harmonic_extension({7: [0.0]})
    with pytest.raises(ValueError):
        sheaf.harmonic_extension({0: [0.0, 1.0]})


def test_fully_clamped_assignment_must_be_a_global_section():
    sheaf = linear_chain_sheaf(3)
    consistent = sheaf.harmonic_extension({0: [1.0], 1: [1.0], 2: [1.0]})
    assert np.allclose([float(node[0]) for node in consistent], [1.0, 1.0, 1.0])
    with pytest.raises(ValueError):
        sheaf.harmonic_extension({0: [0.0], 1: [5.0], 2: [0.0]})


def test_record_injectivity_follows_the_interior_block():
    sheaf = linear_chain_sheaf(3)
    assert sheaf.record_injectivity([0, 2]) is True
    assert sheaf.record_injectivity([0, 1, 2]) is True
    # with no records the interior block is the whole singular Laplacian, so
    # the constants survive as sections the records cannot separate
    assert sheaf.record_injectivity([]) is False
    with pytest.raises(IndexError):
        sheaf.record_injectivity([0, 9])


def test_a_degenerate_restriction_loses_injectivity():
    """A zero restriction map leaves an interior direction the records cannot see."""
    sheaf = CellularSheaf(
        stalk_dims=(1, 1),
        edges=[EdgeRestriction(0, 1, np.array([[0.0]]), np.array([[0.0]]))],
    )
    assert sheaf.record_injectivity([0]) is False
    assert sheaf.energy(np.array([1.0, 2.0])) == 0.0
    assert sheaf.is_global_section(np.array([1.0, 2.0]))


def test_diffusion_converges_to_the_harmonic_extension():
    sheaf = linear_chain_sheaf(3)
    observations = {0: [0.0], 2: [2.0]}
    final = sheaf.diffuse(np.zeros(3), observations, steps=400, step_size=0.1)
    assert np.allclose(final, [0.0, 1.0, 2.0], atol=1e-6)


def test_diffusion_rejects_bad_arguments():
    sheaf = linear_chain_sheaf(3)
    with pytest.raises(ValueError):
        sheaf.diffuse(np.zeros(3), {}, steps=-1)
    with pytest.raises(ValueError):
        sheaf.diffuse(np.zeros(2), {})


def test_sheaf_rejects_malformed_definitions():
    with pytest.raises(ValueError):
        CellularSheaf((), ())
    with pytest.raises(ValueError):
        CellularSheaf((0, 1), ())
    with pytest.raises(ValueError):
        CellularSheaf((1, 1), (), node_names=("only-one",))
    with pytest.raises(ValueError):
        CellularSheaf((1, 1), [EdgeRestriction(0, 0, np.array([[1.0]]), np.array([[1.0]]))])
    with pytest.raises(IndexError):
        CellularSheaf((1, 1), [EdgeRestriction(0, 5, np.array([[1.0]]), np.array([[1.0]]))])
    with pytest.raises(ValueError):
        CellularSheaf((1, 1), [EdgeRestriction(0, 1, np.zeros((2, 1)), np.zeros((1, 1)))])
    with pytest.raises(ValueError):
        CellularSheaf((1, 1), [EdgeRestriction(0, 1, np.zeros((1, 2)), np.zeros((1, 1)))])


def test_carrier_sheaf_reproduces_the_carrier_laplacian():
    carrier = icosahedral_carrier()
    sheaf = carrier_sheaf(carrier)
    assert sheaf.node_count == 12
    assert len(sheaf.edges) == 30
    assert np.allclose(sheaf.laplacian(), carrier.laplacian())
    assert sheaf.record_injectivity([0, 5])


def test_coboundary_matrix_of_the_carrier_has_thirty_rows_and_eleven_rank():
    sheaf = carrier_sheaf(icosahedral_carrier())
    assert coboundary_matrix(sheaf).shape == (30, 12)
    assert coboundary_rank(sheaf) == 11


def test_carrier_cohomology_dimensions():
    assert cohomology_dimensions(carrier_sheaf(icosahedral_carrier())) == {
        "nodes": 12,
        "edges": 30,
        "total_dim": 12,
        "edge_total_dim": 30,
        "coboundary_rank": 11,
        "h0": 1,
        "h1": 19,
        "euler_characteristic": -18,
        "h0_minus_h1": -18,
        "laplacian_equals_coboundary_square": True,
    }


def test_carrier_global_section_basis_is_the_constants():
    sheaf = carrier_sheaf(icosahedral_carrier())
    basis = sheaf.global_section_basis()
    assert basis.shape == (12, 1)
    assert sheaf.is_global_section(np.ones(12))


def test_chain_cohomology_dimensions():
    report = cohomology_dimensions(linear_chain_sheaf(4))
    assert report["nodes"] == 4
    assert report["edges"] == 3
    assert report["coboundary_rank"] == 3
    assert report["h0"] == 1
    assert report["h1"] == 0
    assert report["euler_characteristic"] == 1
    assert report["h0_minus_h1"] == 1
    assert report["laplacian_equals_coboundary_square"] is True
