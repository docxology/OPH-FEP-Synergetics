"""Guards and the remaining public helpers of the geometry layer."""

from __future__ import annotations

import pytest

from oph_fep_synergetics.geometry import (
    MAX_SHELL,
    cumulative_count,
    generate_shell,
    geodesic_edge_count,
    geodesic_face_count,
    geodesic_subdivision_counts,
    geodesic_vertex_count,
    goldberg_dual_counts,
    goldberg_face_count,
    goldberg_vertex_count,
    icosahedron_adjacency,
    sites_through_shell,
    shell_count,
)


def test_icosahedron_adjacency_is_the_twelve_by_twelve_matrix():
    adjacency = icosahedron_adjacency()
    assert len(adjacency) == 12
    assert all(len(row) == 12 for row in adjacency)
    assert all(sum(row) == 5 for row in adjacency)
    assert all(row[vertex] == 0 for vertex, row in enumerate(adjacency))
    assert all(
        row[other] == adjacency[other][vertex]
        for vertex, row in enumerate(adjacency)
        for other in range(12)
    )
    assert sum(sum(row) for row in adjacency) // 2 == 30


def test_shell_helpers_reject_a_negative_index():
    for helper in (shell_count, cumulative_count, generate_shell, sites_through_shell):
        with pytest.raises(ValueError):
            helper(-1)
    assert shell_count(0) == cumulative_count(0) == 1


def test_shell_helpers_reject_an_index_above_the_bound():
    with pytest.raises(ValueError):
        generate_shell(MAX_SHELL + 1)
    with pytest.raises(ValueError):
        sites_through_shell(MAX_SHELL + 1)
    assert shell_count(MAX_SHELL + 1) == 10 * (MAX_SHELL + 1) ** 2 + 2


def test_subdivision_helpers_reject_a_negative_frequency():
    for helper in (
        goldberg_vertex_count,
        goldberg_face_count,
        geodesic_vertex_count,
        geodesic_face_count,
        geodesic_edge_count,
    ):
        with pytest.raises(ValueError):
            helper(-1)
    with pytest.raises(ValueError):
        geodesic_subdivision_counts(0)
    with pytest.raises(ValueError):
        goldberg_dual_counts(0)
    assert goldberg_dual_counts(1)["vertices"] == 20
    assert goldberg_dual_counts(1)["faces"] == 12
