"""SC2: exact carrier counts, spectrum, shell size and port chart."""

from __future__ import annotations

import math

import numpy as np
import pytest

from oph_fep_synergetics.core import (
    OPH_NEIGHBOURS,
    OPH_ORIENTED_BASE_FACES,
    IcosahedralCarrier,
    all_shell_configurations,
    binary_shell_configurations,
    carrier_isomorphism,
    carrier_spectral_gap,
    geometric_carrier,
    graph_signature,
    neighbours_from_faces,
    relaxation_gap_from_laplacian,
    shell_gap,
)

#: The exact first positive Laplacian eigenvalue of the icosahedron.
LAMBDA_1 = 5.0 - math.sqrt(5.0)


def test_counts_are_exact(carrier):
    assert carrier.port_count == 12
    assert carrier.seam_count == 30
    assert carrier.face_count == 20
    assert len(set(carrier.vertices)) == 12


def test_euler_characteristic_of_the_sphere(carrier):
    assert carrier.port_count - carrier.seam_count + carrier.face_count == 2


def test_every_port_meets_five_seams_and_five_faces(carrier):
    assert [carrier.degree(port) for port in range(12)] == [5] * 12
    face_incidence = [sum(1 for face in carrier.faces if port in face) for port in range(12)]
    assert face_incidence == [5] * 12


def test_adjacency_is_derived_from_ophs_oriented_faces(carrier):
    declared = neighbours_from_faces(OPH_ORIENTED_BASE_FACES)
    assert declared == OPH_NEIGHBOURS
    assert [carrier.neighbours(port) for port in range(12)] == list(OPH_NEIGHBOURS)
    # OPH's chart, not the ``0 ~ {1, 2, 3, 4, 6}`` labelling of an earlier draft.
    assert carrier.neighbours(0) == (1, 5, 7, 10, 11)


def test_faces_are_exactly_the_oph_oriented_faces(carrier):
    expected = tuple(sorted(tuple(sorted(face)) for face in OPH_ORIENTED_BASE_FACES))
    assert len(set(expected)) == 20
    assert carrier.faces == expected


def test_edges_are_sorted_unique_and_unoriented(carrier):
    assert list(carrier.edges) == sorted(carrier.edges)
    assert len(set(carrier.edges)) == 30
    assert all(a < b for a, b in carrier.edges)


def test_laplacian_is_five_identity_minus_adjacency(carrier):
    adjacency = carrier.adjacency_matrix()
    assert np.allclose(adjacency, adjacency.T)
    assert np.allclose(np.diag(adjacency), 0.0)
    assert np.allclose(carrier.laplacian(), 5.0 * np.eye(12) - adjacency)


def test_spectrum_is_the_icosahedral_spectrum(carrier):
    spectrum = np.sort(carrier.laplacian_spectrum())
    assert abs(float(spectrum[0])) < 1e-9
    expected = [LAMBDA_1] * 3 + [6.0] * 5 + [5.0 + math.sqrt(5.0)] * 3
    assert np.allclose(spectrum[1:], expected, atol=1e-9)


def test_spectral_gap_is_five_minus_sqrt_five(carrier):
    assert abs(carrier.spectral_gap() - LAMBDA_1) < 1e-12
    assert abs(carrier_spectral_gap() - LAMBDA_1) < 1e-12
    assert abs(carrier_spectral_gap(carrier) - LAMBDA_1) < 1e-12
    assert abs(carrier_spectral_gap(geometric_carrier()) - LAMBDA_1) < 1e-9


def test_shell_gap_is_the_spectral_gap_over_sixty():
    assert shell_gap() == (5.0 - math.sqrt(5.0)) / 60.0
    assert abs(shell_gap() - LAMBDA_1 / 60.0) < 1e-18
    assert abs(relaxation_gap_from_laplacian() - shell_gap()) < 1e-15
    assert abs(relaxation_gap_from_laplacian(geometric_carrier()) - shell_gap()) < 1e-9


def test_shell_of_total_two_has_seventy_eight_states(carrier):
    states = carrier.shell_configurations(total=2)
    assert len(states) == math.comb(13, 11) == 78
    assert len(set(states)) == 78
    assert all(sum(state) == 2 for state in states)
    assert all(min(state) >= 0 for state in states)


def test_binary_shell_has_sixty_six_states(carrier):
    binary = carrier.binary_shell_configurations()
    assert len(binary) == math.comb(12, 2) == 66
    assert all(sorted(state) == [0] * 10 + [1, 1] for state in binary)
    assert set(binary) < set(carrier.shell_configurations(total=2))


def test_module_wrappers_match_the_carrier(carrier):
    assert len(all_shell_configurations()) == 78
    assert len(binary_shell_configurations()) == 66
    assert all_shell_configurations() == all_shell_configurations(carrier)
    assert binary_shell_configurations() == binary_shell_configurations(carrier)


def test_shell_of_total_three_has_three_hundred_sixty_four_states(carrier):
    states = carrier.shell_configurations(total=3)
    assert len(states) == 364
    assert all(sum(state) == 3 for state in states)


def test_negative_total_is_rejected(carrier):
    with pytest.raises(ValueError):
        carrier.shell_configurations(total=-1)


def test_seam_index_round_trips(carrier):
    for index, (a, b) in enumerate(carrier.edges):
        assert carrier.seam_index(a, b) == index
        assert carrier.seam_index(b, a) == index


def test_seam_index_rejects_non_seams(carrier):
    with pytest.raises(ValueError):
        carrier.seam_index(0, 0)
    with pytest.raises(ValueError):
        carrier.seam_index(0, 2)


def test_neighbours_from_faces_rejects_malformed_faces():
    with pytest.raises(ValueError):
        neighbours_from_faces(((0, 0, 1),))
    with pytest.raises(ValueError):
        neighbours_from_faces(((0, 1, 12),))


def test_geometric_realisation_is_the_same_graph(carrier):
    geometric = geometric_carrier()
    assert graph_signature(geometric)[0] == graph_signature(carrier)[0]
    assert np.allclose(geometric.laplacian_spectrum(), carrier.laplacian_spectrum(), atol=1e-9)
    witness = carrier_isomorphism(carrier, geometric)
    assert witness is not None
    assert sorted(witness) == list(range(12))
    image = {(min(witness[a], witness[b]), max(witness[a], witness[b])) for a, b in carrier.edges}
    assert image == set(geometric.edges)


def test_isomorphism_of_the_carrier_with_itself_is_the_identity(carrier):
    witness = carrier_isomorphism(carrier, carrier)
    assert witness is not None
    for port in range(12):
        assert set(carrier.neighbours(port)) == {witness[other] for other in carrier.neighbours(port)}


def test_isomorphism_returns_none_for_a_triangle_free_carrier(carrier):
    # K(6, 6) minus a perfect matching is 5-regular with 30 seams and no triangle,
    # so it passes the degree checks and still cannot be the icosahedron.
    left, right = range(6), range(6, 12)
    bipartite = IcosahedralCarrier(
        vertices=carrier.vertices,
        edges=tuple(sorted((a, b) for a in left for b in right if b - 6 != a)),
        faces=(),
    )
    assert bipartite.port_count == 12
    assert bipartite.seam_count == 30
    assert [bipartite.degree(port) for port in range(12)] == [5] * 12
    assert carrier_isomorphism(carrier, bipartite) is None


def test_isomorphism_rejects_a_different_seam_count(carrier):
    fewer = IcosahedralCarrier(vertices=carrier.vertices, edges=carrier.edges[:29], faces=())
    assert carrier_isomorphism(carrier, fewer) is None


def test_isomorphism_rejects_a_different_port_count(carrier):
    smaller = IcosahedralCarrier(vertices=carrier.vertices[:11], edges=(), faces=())
    assert carrier_isomorphism(carrier, smaller) is None
