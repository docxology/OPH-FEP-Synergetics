"""Exact arithmetic: Q(sqrt 5), exact convex hulls, and the derived tetravolume table."""

from __future__ import annotations

import math
from fractions import Fraction

import pytest

from oph_fep_synergetics.geometry import (
    IVM_TETRA_UNIT,
    PHI_EXACT,
    VOLUME_RATIOS,
    Surd5,
    convex_facets,
    convex_volume,
    derived_volume_ratios,
    exact_icosahedron_vertices,
    icosahedron_tetravolume,
    icosahedron_tetravolume_exact,
    icosahedron_vertices,
    ivm_solid_vertices,
    kissing_arrangements,
    tetravolume_galois_norm,
    tetravolume_ratio_squared_exact,
)


def test_surd5_field_operations_are_exact():
    a = Surd5(Fraction(3), Fraction(1))
    b = Surd5(Fraction(1, 2), Fraction(-2))
    assert a + b == Surd5(Fraction(7, 2), Fraction(-1))
    assert a - b == Surd5(Fraction(5, 2), Fraction(3))
    assert a * b == Surd5(Fraction(3, 2) - 10, Fraction(-6) + Fraction(1, 2))
    assert (a / b) * b == a
    assert a.conjugate() == Surd5(3, -1)
    assert a.norm() == Fraction(4)
    assert -a == Surd5(-3, -1)
    assert 2 - a == Surd5(-1, -1)
    assert float(a) == pytest.approx(3 + math.sqrt(5))


def test_surd5_sign_is_decided_exactly():
    assert Surd5(0, 0).sign() == 0
    assert Surd5(3, -1).sign() == 1  # 3 > sqrt 5
    assert Surd5(2, -1).sign() == -1  # 2 < sqrt 5
    assert Surd5(-3, 1).sign() == -1
    assert Surd5(0, -1).sign() == -1
    assert PHI_EXACT * PHI_EXACT == PHI_EXACT + 1


def test_surd5_division_by_zero_raises():
    with pytest.raises(ZeroDivisionError):
        Surd5(1, 1) / Surd5(0, 0)


def test_surd5_prints_its_exact_spelling():
    assert str(Surd5(Fraction(15, 2), Fraction(5, 2))) == "15/2 + 5/2 sqrt 5"
    assert str(Surd5(Fraction(15, 2), Fraction(-5, 2))) == "15/2 - 5/2 sqrt 5"
    assert str(Surd5(4)) == "4"


def test_the_unit_tetrahedron_has_volume_eight_thirds():
    assert convex_volume(ivm_solid_vertices()["tetrahedron"]) == IVM_TETRA_UNIT


def test_every_table_entry_is_derived_from_its_vertex_set():
    assert derived_volume_ratios() == VOLUME_RATIOS
    assert all(value.denominator == 1 for value in derived_volume_ratios().values())


@pytest.mark.parametrize(
    ("name", "vertices", "facets"),
    [
        ("tetrahedron", 4, 4),
        ("cube", 8, 6),
        ("octahedron", 6, 8),
        ("rhombic_dodecahedron", 14, 12),
        ("cuboctahedron", 12, 14),
        ("truncated_octahedron", 24, 14),
    ],
)
def test_each_solid_has_its_vertex_and_facet_count(name, vertices, facets):
    solid = ivm_solid_vertices()[name]
    assert len(solid) == vertices
    assert len(convex_facets(solid)) == facets


def test_facets_are_ordered_outward():
    # A mis-ordered facet would raise inside convex_volume; the cube is the
    # smallest case with non-triangular facets.
    cube = ivm_solid_vertices()["cube"]
    assert convex_volume(cube) == Fraction(8)
    assert all(len(facet) == 4 for facet in convex_facets(cube))


def test_degenerate_inputs_raise():
    with pytest.raises(ValueError):
        convex_facets([(0, 0, 0), (1, 0, 0), (0, 1, 0)])
    with pytest.raises(ValueError):
        convex_facets([(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)])
    with pytest.raises(ValueError):
        # an interior point is not a vertex of the hull
        convex_facets([(0, 0, 0), (4, 0, 0), (0, 4, 0), (0, 0, 4), (1, 1, 1)])
    with pytest.raises(ValueError):
        convex_facets([(0.5, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])


def test_exact_icosahedron_matches_the_float_vertices():
    for exact, measured in zip(exact_icosahedron_vertices(), icosahedron_vertices()):
        assert tuple(float(value) for value in exact) == pytest.approx(measured)


def test_icosahedron_tetravolume_is_sqrt_two_times_an_exact_element():
    report = icosahedron_tetravolume_exact()
    assert report["edge_squared"] == 4
    assert report["edges"] == 30
    assert report["facets"] == 20
    assert report["facet_sizes"] == (3,)
    assert report["volume"] == Surd5(10, Fraction(10, 3))
    assert report["sqrt2_coefficient"] == Surd5(Fraction(15, 2), Fraction(5, 2))
    assert report["square"] == Surd5(175, 75)
    assert report["norm"] == 50
    assert report["square_norm"] == 2500
    assert report["tetravolume"] == pytest.approx(
        float(icosahedron_tetravolume()["tetravolume"]), abs=1e-12
    )
    assert report["conjugate_tetravolume"] == pytest.approx(2.700907567377264, abs=1e-12)


def test_the_solids_constants_are_computed_not_typed():
    assert tetravolume_ratio_squared_exact() == (Fraction(175), Fraction(75))
    assert tetravolume_galois_norm() == Fraction(50)
    assert icosahedron_tetravolume()["integral"] is False


def test_the_two_kissing_arrangements_are_distinguished_exactly():
    report = kissing_arrangements()
    assert report["cuboctahedron_edge_over_radius_squared"] == Fraction(1)
    assert report["icosahedron_edge_over_radius_squared"] == Surd5(2, Fraction(-2, 5))
    assert report["icosahedron_ratio_float"] == pytest.approx(2 - 2 * math.sqrt(5) / 5)
    assert report["icosahedron_has_slack"] is True
    assert report["cuboctahedron_is_rigid"] is True
