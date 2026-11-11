"""SC5: the cube lattice, its line graph and the cuboctahedron they generate."""

from __future__ import annotations

from oph_fep_synergetics.core.network import graphs_are_isomorphic
from oph_fep_synergetics.geometry import (
    cube_adjacency,
    cube_line_graph_adjacency,
    cuboctahedron_adjacency,
)


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
