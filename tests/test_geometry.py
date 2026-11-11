"""SC4: quadray/IVM geometry, exact volumes, and subdivision counts."""

from __future__ import annotations

import itertools
import math
from fractions import Fraction

import pytest

from oph_fep_synergetics.geometry import (
    EMBEDDING,
    MAX_SHELL,
    NEIGHBOR_SQUARED_DISTANCE,
    PHI,
    VOLUME_RATIOS,
    Quadray,
    cuboctahedron_adjacency,
    cuboctahedron_edges,
    cuboctahedron_faces,
    cuboctahedron_neighbors,
    cuboctahedron_signature,
    cuboctahedron_volume,
    cumulative_count,
    dot,
    from_xyz,
    generate_shell,
    geodesic_edge_count,
    geodesic_face_count,
    geodesic_subdivision_counts,
    geodesic_vertex_count,
    goldberg_dual_counts,
    goldberg_face_count,
    goldberg_vertex_count,
    icosahedron_edges,
    icosahedron_faces,
    icosahedron_vertices,
    is_ivm_site,
    is_octahedron_site,
    ivm_shell_matches_geodesic,
    ivm_shell_matches_goldberg,
    neighbor_moves,
    neighbor_sites,
    octahedron_vertices,
    octahedron_volume,
    scaled_canonical,
    shell_count,
    shell_index,
    sites_through_shell,
    squared_distance,
    squared_norm,
    tetrahedron_vertices,
    tetrahedron_volume,
    tetra_volume,
    to_xyz,
)


# ---------------------------------------------------------------------------
# Quadray coordinates
# ---------------------------------------------------------------------------


def test_neighbour_moves_are_the_twelve_signature_permutations():
    moves = neighbor_moves()
    assert len(moves) == 12
    assert len(set(moves)) == 12
    assert moves == tuple(sorted(set(itertools.permutations((2, 1, 1, 0)))))
    assert all(sum(move) == 4 for move in moves)


def test_neighbour_moves_are_equidistant_and_never_closer():
    moves = neighbor_moves()
    assert squared_norm(moves[0]) == NEIGHBOR_SQUARED_DISTANCE == 8
    assert all(squared_norm(move) == 8 for move in moves)
    pairwise = {
        squared_distance(a, b) for a, b in itertools.combinations(moves, 2)
    }
    assert min(pairwise) == 8
    assert len(neighbor_sites()) == 12


def test_quadray_gauge_and_normalisation():
    assert Quadray(1, 2, 3, 4).normalize() == Quadray(0, 1, 2, 3)
    assert Quadray(5, 5, 5, 5).normalize() == Quadray(0, 0, 0, 0)
    assert Quadray(2, 3, 4, 5).scale(2) == Quadray(4, 6, 8, 10)
    assert (Quadray(1, 2, 3, 4) - Quadray(1, 1, 1, 1)) == Quadray(0, 1, 2, 3)
    assert (Quadray(1, 2, 3, 4) + Quadray(1, 1, 1, 1)) == Quadray(2, 3, 4, 5)
    assert Quadray(1, 2, 3, 4).as_tuple() == (1, 2, 3, 4)
    assert tuple(Quadray(1, 2, 3, 4)) == (1, 2, 3, 4)


def test_quadray_rejects_non_integers():
    with pytest.raises(TypeError):
        Quadray(1.5, 0, 0, 0)
    with pytest.raises(TypeError):
        Quadray(True, 0, 0, 0)
    assert Quadray(2.0, 0, 0, 0).as_tuple() == (2, 0, 0, 0)


def test_quadray_metric_is_the_diamond_form():
    """``<p, q> = 4 p.q - (sum p)(sum q)``, the form QuadMath uses."""
    for left, right in itertools.combinations(neighbor_moves()[:5], 2):
        expected = 4 * sum(a * b for a, b in zip(left, right)) - sum(left) * sum(right)
        assert abs(dot(left, right) - expected) < 1e-12
        assert abs(squared_distance(left, right) - (8 + 8 - 2 * expected)) < 1e-12


def test_embedding_has_four_tetrahedral_columns():
    assert len(EMBEDDING) == 3
    assert all(len(row) == 4 for row in EMBEDDING)
    assert EMBEDDING == ((1, -1, -1, 1), (1, 1, -1, -1), (1, -1, 1, -1))
    assert all(sum(row) == 0 for row in EMBEDDING)
    assert scaled_canonical(Quadray(2, 1, 1, 0), 1).as_tuple() == (2, 1, 1, 0)


def test_xyz_round_trip_is_exact_on_shell_sites():
    for (a, b, c, d) in sites_through_shell(2):
        point = Quadray(a, b, c, d)
        assert from_xyz(*to_xyz(point)) == point.normalize()


def test_xyz_projection_and_singular_embeddings():
    assert to_xyz((1, 2, 0, 3)) == (2.0, 0.0, -4.0)
    custom = ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0))
    assert from_xyz(1.0, 2.0, 3.0, custom) == Quadray(1, 2, 3, 0)
    with pytest.raises(ValueError):
        to_xyz((1, 1, 0, 0), ((1, 2), (3, 4), (5, 6)))
    with pytest.raises(ValueError):
        to_xyz((1, 1, 0, 0), ((1, 1), (1, 1), (1, 1), (1, 1)))
    with pytest.raises(ValueError):
        from_xyz(1.0, 1.0, 1.0, ((0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)))
    with pytest.raises(ValueError):
        from_xyz(1.0, 1.0, 1.0, ((1, -1, 0, 0), (1, -1, 0, 0), (1, -1, 0, 0)))


def test_quadray_rejects_wrong_length_points():
    with pytest.raises(ValueError):
        squared_norm((1, 1, 0))
    with pytest.raises(ValueError):
        to_xyz((1, 1, 0))


def test_site_predicates():
    assert is_ivm_site((2, 1, 1, 0))
    assert is_ivm_site((1, 1, 1, 1))
    assert not is_ivm_site((1, 1, 0, 0))
    assert is_octahedron_site((1, 1, 0, 0))
    assert not is_octahedron_site((1, 1, 1, 1))
    assert not is_octahedron_site((2, 1, 1, 0))


# ---------------------------------------------------------------------------
# The IVM shell sequence
# ---------------------------------------------------------------------------


def test_shell_counts_follow_ten_k_squared_plus_two():
    # QuadMath's own convention: the centre shell holds one site, and every
    # frequency k >= 1 holds 10k^2 + 2.
    assert shell_count(0) == 1
    for k in range(1, 9):
        assert shell_count(k) == 10 * k * k + 2
    assert [shell_count(k) for k in range(9)] == [1, 12, 42, 92, 162, 252, 362, 492, 642]
    assert MAX_SHELL == 32


def test_cumulative_counts_are_the_running_sums():
    cumulative = [cumulative_count(k) for k in range(6)]
    assert cumulative == [1, 13, 55, 147, 309, 561]
    for k in range(6):
        assert cumulative[k] == sum(shell_count(j) for j in range(k + 1))
    assert cumulative_count(0) == 1


def test_generated_shells_match_their_counts():
    for k in range(4):
        shell = generate_shell(k)
        assert len(shell) == shell_count(k)
        assert len(set(shell)) == shell_count(k)
        assert all(shell_index(point) == k for point in shell)
    assert len(sites_through_shell(3)) == cumulative_count(3) == 147


def test_first_shell_is_the_twelve_around_one():
    first = generate_shell(1)
    assert len(first) == 12
    assert {point.as_tuple() for point in first} == set(neighbor_moves())
    assert all(isinstance(point, Quadray) for point in first)
    assert len(cuboctahedron_neighbors()) == 12
    assert squared_norm(first[0]) == NEIGHBOR_SQUARED_DISTANCE


def test_shell_index_rejects_points_outside_the_enumeration():
    with pytest.raises(ValueError):
        shell_index((1, 1, 0, 0))


# ---------------------------------------------------------------------------
# Volumes and the vector-equilibrium shell
# ---------------------------------------------------------------------------


def test_exact_ivm_volume_ratios():
    assert tetrahedron_volume() == Fraction(1)
    assert octahedron_volume() == Fraction(4)
    assert cuboctahedron_volume() == Fraction(20)
    assert VOLUME_RATIOS["cube"] == Fraction(3)
    assert VOLUME_RATIOS["rhombic_dodecahedron"] == Fraction(6)
    assert VOLUME_RATIOS["truncated_octahedron"] == Fraction(96)
    assert VOLUME_RATIOS["tetrahedron"] == Fraction(1)
    assert VOLUME_RATIOS["cuboctahedron"] == Fraction(20)


def test_tetra_volume_uses_the_projected_determinant():
    axes = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    # the four axis directions form the unit IVM tetrahedron (determinant four)
    assert tetra_volume(*axes) == Fraction(1)
    # the origin with three axes is the primitive simplex (determinant one)
    assert tetra_volume((0, 0, 0, 0), *axes[:3]) == Fraction(1, 4)
    assert tetra_volume(*axes[:3], (0, 0, 0, 0)) == Fraction(1, 4)
    assert tetra_volume(*tetrahedron_vertices()) == Fraction(1)


def test_cuboctahedron_is_the_vector_equilibrium():
    signature = cuboctahedron_signature()
    assert signature["vertices"] == 12
    assert signature["edges"] == 24
    assert signature["degree"] == [4]
    assert signature["triangles"] == 8
    assert signature["squares"] == 6
    assert signature["squared_edge_length"] == 8
    assert signature["volume"] == Fraction(20)
    assert len(cuboctahedron_edges()) == 24
    adjacency = cuboctahedron_adjacency()
    assert len(adjacency) == 12
    assert all(len(row) == 12 for row in adjacency)
    assert all(sum(row) == 4 for row in adjacency)
    assert all(row[vertex] == 0 for vertex, row in enumerate(adjacency))
    assert all(
        row[other] == adjacency[other][vertex]
        for vertex, row in enumerate(adjacency)
        for other in range(12)
    )


def test_cuboctahedron_faces_are_eight_triangles_and_six_squares():
    faces = cuboctahedron_faces()
    assert set(faces) == {"triangles", "squares"}
    assert len(faces["triangles"]) == 8
    assert len(faces["squares"]) == 6
    assert all(len(face) == 3 for face in faces["triangles"])
    assert all(len(face) == 4 for face in faces["squares"])


def test_octahedron_has_six_vertices_of_the_primitive_shell():
    vertices = octahedron_vertices()
    assert len(vertices) == 6
    assert len(set(vertices)) == 6
    assert all(is_octahedron_site(vertex.as_tuple()) for vertex in vertices)
    assert all(sorted(vertex.as_tuple()) == [0, 0, 1, 1] for vertex in vertices)


# ---------------------------------------------------------------------------
# Icosahedron, geodesic and Goldberg subdivision
# ---------------------------------------------------------------------------


def test_icosahedron_helpers_are_the_twelve_port_carrier():
    assert len(icosahedron_vertices()) == 12
    assert len(set(icosahedron_vertices())) == 12
    assert len(icosahedron_edges()) == 30
    assert len(icosahedron_faces()) == 20
    assert abs(PHI - (1.0 + math.sqrt(5.0)) / 2.0) < 1e-18


def test_geodesic_class_two_counts():
    counts = geodesic_subdivision_counts(2)
    assert counts["vertices"] == 42
    assert counts["edges"] == 120
    assert counts["faces"] == 80
    assert counts["frequency"] == 2
    assert counts["original_vertices"] == 12
    assert counts["added_vertices"] == 30


def test_subdivision_count_formulas():
    for frequency in range(1, 6):
        assert geodesic_vertex_count(frequency) == 10 * frequency * frequency + 2
        assert goldberg_face_count(frequency) == 10 * frequency * frequency + 2
        assert geodesic_face_count(frequency) == 20 * frequency * frequency
        assert goldberg_vertex_count(frequency) == 20 * frequency * frequency
        assert geodesic_edge_count(frequency) == 30 * frequency * frequency
        assert geodesic_subdivision_counts(frequency)["vertices"] == 10 * frequency * frequency + 2
        assert goldberg_dual_counts(frequency)["vertices"] == 20 * frequency * frequency


def test_subdivision_complexes_obey_euler_and_are_dual():
    for frequency in range(1, 5):
        counts = geodesic_subdivision_counts(frequency)
        assert counts["vertices"] - counts["edges"] + counts["faces"] == 2
        dual = goldberg_dual_counts(frequency)
        assert dual["vertices"] - dual["edges"] + dual["faces"] == 2
        assert dual["vertices"] == counts["faces"]
        assert dual["faces"] == counts["vertices"]
        assert dual["edges"] == counts["edges"]


def test_ivm_shell_sizes_are_the_subdivision_sequence():
    for frequency in range(1, 5):
        assert ivm_shell_matches_geodesic(frequency) is True
        assert ivm_shell_matches_goldberg(frequency) is True
        assert shell_count(frequency) == geodesic_vertex_count(frequency)
        assert shell_count(frequency) == goldberg_face_count(frequency)
