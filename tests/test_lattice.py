"""SC5: the cube lattice, its line graph and the cuboctahedron they generate."""

from __future__ import annotations

from collections import Counter
from itertools import combinations, product

from oph_fep_synergetics.core.network import graphs_are_isomorphic
from oph_fep_synergetics.geometry import (
    cube_adjacency,
    cube_line_graph_adjacency,
    cuboctahedron_adjacency,
)
from oph_fep_synergetics.geometry.lattice import ivm_quadray_cosets


def _edges(adjacency: tuple[tuple[int, ...], ...]) -> int:
    """Return the number of undirected edges in an adjacency matrix."""
    return sum(sum(row) for row in adjacency) // 2


def test_cube_adjacency_is_eight_by_eight_three_regular_with_twelve_edges():
    adjacency = cube_adjacency()
    assert len(adjacency) == 8
    assert all(len(row) == 8 for row in adjacency)
    assert {sum(row) for row in adjacency} == {3}
    assert _edges(adjacency) == 12
    assert all(adjacency[i][j] == adjacency[j][i] for i in range(8) for j in range(8))
    assert all(adjacency[i][i] == 0 for i in range(8))


def test_cube_line_graph_adjacency_is_twelve_by_twelve_four_regular_with_twenty_four_edges():
    adjacency = cube_line_graph_adjacency()
    assert len(adjacency) == 12
    assert all(len(row) == 12 for row in adjacency)
    assert {sum(row) for row in adjacency} == {4}
    assert _edges(adjacency) == 24
    assert all(adjacency[i][j] == adjacency[j][i] for i in range(12) for j in range(12))
    assert all(adjacency[i][i] == 0 for i in range(12))


def test_the_line_graph_of_the_cube_is_the_cuboctahedron():
    assert graphs_are_isomorphic(cube_line_graph_adjacency(), cuboctahedron_adjacency()) is True


def test_the_cube_is_not_isomorphic_to_its_line_graph():
    assert graphs_are_isomorphic(cube_adjacency(), cube_line_graph_adjacency()) is False


def test_quadray_residues_classify_tetrahedral_and_octahedral_holes():
    geometry = ivm_quadray_cosets()
    assert geometry == {
        "minimum_squared_distances": (0, 3, 4, 3),
        "nearest_counts": (1, 4, 6, 4),
        "nearest_pair_squared_distance_counts": ((), ((8, 6),), ((8, 12), (16, 3)), ((8, 6),)),
    }

    def distance_squared(a: tuple[int, ...], b: tuple[int, ...]) -> int:
        return sum((x - y) ** 2 for x, y in zip(a, b))

    # Independently enumerate scaled Cartesian D3 around r times the first
    # quadray direction, including a translated representative of each class.
    for residue in range(4):
        for translation in ((0, 0, 0), (4, 2, 2)):
            centre = tuple(residue + coordinate for coordinate in translation)
            candidates = [
                site for site in product(*(range(value - 2, value + 3) for value in centre))
                if all(coordinate % 2 == 0 for coordinate in site) and sum(site) % 4 == 0
            ]
            minimum = min(distance_squared(centre, site) for site in candidates)
            nearest = [site for site in candidates if distance_squared(centre, site) == minimum]
            pair_distances = Counter(distance_squared(a, b) for a, b in combinations(nearest, 2))
            assert minimum == geometry["minimum_squared_distances"][residue]
            assert len(nearest) == geometry["nearest_counts"][residue]
            assert tuple(sorted(pair_distances.items())) == (
                geometry["nearest_pair_squared_distance_counts"][residue]
            )
