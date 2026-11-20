"""Keep the Jitterbug's bar, convex-hull and metric graphs distinct."""

from fractions import Fraction
from math import dist

import pytest

from oph_fep_synergetics.core.carrier import icosahedral_carrier
from oph_fep_synergetics.core.network import graphs_are_isomorphic
from oph_fep_synergetics.geometry.jitterbug import (
    ICOSAHEDRON_PARAMETER,
    jitterbug_census,
    jitterbug_matching_edges,
    jitterbug_rigid_edges,
    jitterbug_rigid_triangles,
    jitterbug_snapshot,
    jitterbug_vertices,
    scaled_jitterbug_vertices,
)
from oph_fep_synergetics.geometry.lattice import cuboctahedron_adjacency


def test_jitterbug_keeps_its_rigid_equilateral_panels_and_six_short_diagonals():
    bars = jitterbug_rigid_edges()
    matching = jitterbug_matching_edges()
    assert len(jitterbug_rigid_triangles()) == 8
    assert len(bars) == 24
    assert len(matching) == 6
    assert set(bars).isdisjoint(matching)
    assert sorted(port for pair in matching for port in pair) == list(range(12))
    for parameter in (0, Fraction(3, 5), ICOSAHEDRON_PARAMETER, 0.75, 1):
        points = scaled_jitterbug_vertices(parameter)
        assert len(points) == 12
        assert all(dist(points[left], points[right]) == pytest.approx(2.0) for left, right in bars)
    with pytest.raises(ValueError, match="between zero and one"):
        jitterbug_vertices(Fraction(5, 4))
    with pytest.raises(ValueError, match="finite and positive"):
        scaled_jitterbug_vertices(Fraction(3, 4), edge_length=0)


def test_jitterbug_hull_triangulates_before_all_thirty_seams_become_equal():
    census = jitterbug_census()
    for snapshot in census.values():
        assert snapshot["bar_lengths_match"]
        assert snapshot["matching_lengths_match"]
        assert snapshot["rigid_orientations_match"]
        assert snapshot["hull_matches_designated_graph"]
    cube = census["cuboctahedron"]
    middle = census["intermediate"]
    regular = census["icosahedron"]
    contracted = census["contracted"]
    assert cube["facet_counts"] == ((3, 8), (4, 6))
    assert cube["hull_edges"] == cube["bar_edges"] == cube["shortest_pairs"]
    for interior in (middle, regular, contracted):
        assert interior["facet_counts"] == ((3, 20),)
        assert interior["bar_edge_count"] == 24
        assert interior["matching_edge_count"] == 6
        assert interior["hull_edge_count"] == 30
    assert middle["shortest_pairs"] == middle["bar_edges"]
    assert middle["pair_squared_length_counts"] == (
        (Fraction(13, 8), 24), (Fraction(9, 4), 6), (Fraction(4), 6),
        (Fraction(37, 8), 24), (Fraction(25, 4), 6),
    )
    assert regular["shortest_pairs"] == regular["hull_edges"]
    assert regular["matching_squared_length"] == regular["bar_squared_length"]
    assert tuple(count for _, count in regular["pair_squared_length_counts"]) == (30, 30, 6)
    assert contracted["shortest_pairs"] == contracted["matching_edges"]

    # The supporting-plane result recovers the two separately constructed
    # endpoint graphs, up to relabelling.
    references = (cuboctahedron_adjacency(), icosahedral_carrier().adjacency_matrix().tolist())
    for snapshot, reference in zip((cube, regular), references, strict=True):
        edges = set(snapshot["hull_edges"])
        adjacency = tuple(
            tuple(int(tuple(sorted((left, right))) in edges) for right in range(12))
            for left in range(12)
        )
        assert graphs_are_isomorphic(adjacency, reference)


def test_jitterbug_shortest_chords_differ_from_outer_sphere_contacts():
    census = jitterbug_census()
    assert census["cuboctahedron"]["contact_pairs"] == jitterbug_rigid_edges()
    assert census["icosahedron"]["shortest_pair_count"] == 30
    for name in ("intermediate", "icosahedron", "contracted"):
        assert census[name]["contact_pair_count"] == census[name]["overlap_pair_count"] == 0
    assert census["overlap_control"]["overlap_pairs"] == jitterbug_matching_edges()


def test_jitterbug_octahedral_limit_retains_labels_and_doubles_each_bar():
    collapsed = jitterbug_snapshot(0)
    assert len(jitterbug_vertices(0)) == 12
    assert collapsed["vertex_count"] == 6
    assert collapsed["coincident_joint_pairs"] == 6
    assert sorted(collapsed["joint_to_vertex"]) == [vertex for vertex in range(6) for _ in range(2)]
    assert collapsed["facet_counts"] == ((3, 8),)
    assert collapsed["matching_edges"] == ()
    assert collapsed["bar_edge_count"] == collapsed["hull_edge_count"] == 12
    assert collapsed["hull_edges"] == collapsed["bar_edges"] == collapsed["shortest_pairs"]
    assert collapsed["bar_multiplicities"] == (2,) * 12
    collapsed["vertex_count"] = 0
    assert jitterbug_snapshot(0)["vertex_count"] == 6
