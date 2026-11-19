"""SC5: the carrier's resistance, distance shells, subgraph census and isomorphism."""

from __future__ import annotations

from fractions import Fraction

import numpy as np

from oph_fep_synergetics.core import icosahedral_carrier, reciprocal_spectrum_sum
from oph_fep_synergetics.core.carrier import geometric_carrier
from oph_fep_synergetics.core.network import (
    carrier_resistance_geometry,
    distance_distribution,
    effective_resistance,
    graph_distances,
    graphs_are_isomorphic,
    intersection_array,
    kirchhoff_index,
    laplacian_matrix,
    resistance_ladder,
    subgraph_census,
)
from oph_fep_synergetics.geometry import Surd5, cuboctahedron_adjacency


def _adjacency() -> list[list[int]]:
    """Return the 12x12 carrier adjacency as plain integers."""
    return [[int(value) for value in row] for row in icosahedral_carrier().adjacency_matrix()]


def _relabel(adjacency: list[list[int]], permutation: list[int]) -> list[list[int]]:
    """Return the adjacency matrix under a vertex permutation."""
    return [[adjacency[permutation[i]][permutation[j]] for j in range(len(adjacency))] for i in range(len(adjacency))]


def test_carrier_adjacency_is_twelve_by_twelve_five_regular_with_thirty_edges():
    adjacency = _adjacency()
    assert len(adjacency) == 12
    assert all(len(row) == 12 for row in adjacency)
    assert {sum(row) for row in adjacency} == {5}
    assert sum(sum(row) for row in adjacency) // 2 == 30


def test_laplacian_matrix_is_the_degree_matrix_minus_the_adjacency():
    adjacency = _adjacency()
    laplacian = laplacian_matrix(adjacency)
    assert len(laplacian) == 12
    for i in range(12):
        for j in range(12):
            expected = sum(adjacency[i]) if i == j else -adjacency[i][j]
            assert laplacian[i][j] == expected
        assert sum(laplacian[i]) == 0


def test_graph_distances_are_symmetric_with_zero_diagonal_and_max_three():
    distances = graph_distances(_adjacency())
    assert all(distances[i][j] == distances[j][i] for i in range(12) for j in range(12))
    assert all(distances[i][i] == 0 for i in range(12))
    assert max(max(row) for row in distances) == 3


def test_distance_distribution_is_one_five_five_one():
    assert distance_distribution(_adjacency()) == (1, 5, 5, 1)


def test_intersection_array_is_five_two_one_by_one_two_five():
    assert intersection_array(_adjacency()) == ((5, 2, 1), (1, 2, 5))


def test_effective_resistance_is_symmetric_and_zero_on_the_diagonal():
    adjacency = _adjacency()
    for i in range(12):
        assert effective_resistance(adjacency, i, i) == 0
        for j in range(12):
            value = effective_resistance(adjacency, i, j)
            assert isinstance(value, Fraction)
            assert value == effective_resistance(adjacency, j, i)


def test_effective_resistance_of_an_adjacent_pair_is_eleven_thirtieths():
    adjacency = _adjacency()
    assert adjacency[0][1] == 1
    assert effective_resistance(adjacency, 0, 1) == Fraction(11, 30)


def test_effective_resistance_of_a_distance_three_pair_is_one_half():
    adjacency = _adjacency()
    assert graph_distances(adjacency)[0][3] == 3
    assert effective_resistance(adjacency, 0, 3) == Fraction(1, 2)


def test_resistance_ladder_is_one_shell_at_a_time():
    ladder = resistance_ladder(_adjacency())
    assert ladder == (
        (1, Fraction(11, 30), 30),
        (2, Fraction(7, 15), 30),
        (3, Fraction(1, 2), 6),
    )
    # The three shells partition all unordered pairs exactly once.
    assert sum(count for _, _, count in ladder) == 12 * 11 // 2


def test_carrier_chords_recover_only_the_lowest_mode_contrast_variance():
    comparison = carrier_resistance_geometry()
    assert comparison["shells"] == (
        (1, Fraction(11, 30), Surd5(4), 30),
        (2, Fraction(7, 15), Surd5(6, 2), 30),
        (3, Fraction(1, 2), Surd5(10, 2), 6),
    )
    assert comparison["lowest_eigenvalue"] == Surd5(5, -1)
    assert comparison["lowest_multiplicity"] == 3
    assert comparison["coordinate_norm_squared"] == Surd5(10, 2)
    assert comparison["chord_to_lowest_variance_factor"] == 40
    assert comparison["lowest_mode_variances"] == (
        Fraction(1, 10),
        Surd5(Fraction(3, 20), Fraction(1, 20)),
        Surd5(Fraction(1, 4), Fraction(1, 20)),
    )
    assert comparison["full_variance_proportional"] is False

    # Independently build the covariance projector by numerical eigensolve
    # and the chords from the regular carrier on the unit sphere.
    carrier = geometric_carrier()
    laplacian = carrier.laplacian()
    eigenvalues, eigenvectors = np.linalg.eigh(laplacian)
    lowest = np.isclose(eigenvalues, eigenvalues[1])
    basis = eigenvectors[:, lowest]
    lowest_covariance = basis @ basis.T / eigenvalues[1]
    full_covariance = np.linalg.pinv(laplacian)

    points = np.asarray(carrier.vertices)
    a, b = carrier.edges[0]
    points *= 2 / np.linalg.norm(points[a] - points[b])
    differences = points[:, None, :] - points[None, :, :]
    chord_squared = (differences * differences).sum(axis=2)
    lowest_diagonal = lowest_covariance.diagonal()
    lowest_variances = lowest_diagonal[:, None] + lowest_diagonal[None, :] - 2 * lowest_covariance
    full_diagonal = full_covariance.diagonal()
    full_variances = full_diagonal[:, None] + full_diagonal[None, :] - 2 * full_covariance
    np.testing.assert_allclose(chord_squared, 40 * lowest_variances, atol=1e-12)
    assert not np.allclose(chord_squared, full_variances * (4 / full_variances[a, b]))

    center = np.eye(carrier.port_count) - 1 / carrier.port_count
    chord_gram = -0.5 * center @ chord_squared @ center
    inferred_full_gram = -0.5 * center @ full_variances @ center
    np.testing.assert_allclose(inferred_full_gram, full_covariance, atol=1e-12)
    assert comparison["ranks"] == (
        np.linalg.matrix_rank(inferred_full_gram),
        np.linalg.matrix_rank(chord_gram),
    ) == (11, 3)


def test_kirchhoff_index_is_twenty_eight_and_matches_the_spectrum():
    adjacency = _adjacency()
    assert kirchhoff_index(adjacency) == Fraction(28, 1)
    assert kirchhoff_index(adjacency) == 12 * reciprocal_spectrum_sum(icosahedral_carrier())


def test_subgraph_census_of_the_icosahedron():
    assert subgraph_census(_adjacency()) == {
        "triangles": 20,
        "pentagons": 12,
        "perfect_matchings": 125,
        "hamiltonian_cycles": 1280,
    }


def test_subgraph_census_of_the_cuboctahedron():
    assert subgraph_census(cuboctahedron_adjacency()) == {
        "triangles": 8,
        "pentagons": 0,
        "perfect_matchings": 32,
        "hamiltonian_cycles": 100,
    }


def test_graphs_are_isomorphic_is_reflexive_and_permutation_invariant():
    adjacency = _adjacency()
    assert graphs_are_isomorphic(adjacency, adjacency) is True
    permutation = list(range(11, -1, -1))
    assert graphs_are_isomorphic(_relabel(adjacency, permutation), adjacency) is True


def test_graphs_are_isomorphic_rejects_a_different_order():
    adjacency = _adjacency()
    assert graphs_are_isomorphic(adjacency, cuboctahedron_adjacency()) is False
