"""SC4: polyhedron volumes, the icosahedron's tetravolume, and the IVM Voronoi cell."""

from __future__ import annotations

import math
from fractions import Fraction

import pytest

from oph_fep_synergetics.core import icosahedral_carrier
from oph_fep_synergetics.geometry import (
    angle_defects,
    icosahedron_faces,
    icosahedron_tetravolume,
    icosahedron_vertices,
    packing_density,
    quadratic_norm,
    tetravolume_galois_norm,
    voronoi_cell,
    voronoi_cell_tetravolume,
    voronoi_cell_volume,
)
from oph_fep_synergetics.geometry import solids


def test_shortest_edge_of_the_icosahedron_is_two():
    assert solids.shortest_edge(icosahedron_vertices()) == pytest.approx(2.0, abs=1e-12)


def test_regular_tetrahedron_volume():
    assert solids.regular_tetrahedron_volume(2.0) == pytest.approx(math.sqrt(2) / 12 * 8, abs=1e-12)


def test_icosahedron_volume():
    assert solids.icosahedron_volume(2.0) == pytest.approx(5 / 12 * (3 + math.sqrt(5)) * 8, abs=1e-9)


def test_orient_faces_outward_flips_ten_of_twenty():
    vertices = icosahedron_vertices()
    oriented, flips = solids.orient_faces_outward(vertices, icosahedron_faces())
    assert flips == 10
    assert len(oriented) == 20


def test_raw_icosahedron_faces_are_rejected_as_unoriented():
    with pytest.raises(ValueError):
        solids.polyhedron_volume(icosahedron_vertices(), icosahedron_faces())


def test_oriented_volume_matches_both_the_closed_form_and_the_measurement():
    vertices = icosahedron_vertices()
    oriented, _ = solids.orient_faces_outward(vertices, icosahedron_faces())
    volume = solids.polyhedron_volume(vertices, oriented)
    assert volume == pytest.approx(solids.icosahedron_volume(2.0), abs=1e-9)
    assert volume == pytest.approx(float(icosahedron_tetravolume()["volume"]), abs=1e-9)


def test_tetravolume_ratio_squared_is_exact():
    assert solids.tetravolume_ratio_squared_exact() == (Fraction(175), Fraction(75))


def test_icosahedron_tetravolume_signature():
    report = icosahedron_tetravolume()
    assert report["vertices"] == 12
    assert report["edges"] == 30
    assert report["faces"] == 20
    assert report["faces_reoriented"] == 10
    assert report["edge_spread"] == 0.0
    assert report["edge"] == pytest.approx(2.0, abs=1e-12)
    assert report["tetravolume"] == pytest.approx(2.5 * (3 * math.sqrt(2) + math.sqrt(10)), abs=1e-9)
    assert report["tetravolume_squared"] == pytest.approx(175 + 75 * math.sqrt(5), abs=1e-9)
    assert report["exact_squared_parts"] == (175, 75)
    assert report["exact_squared_residual"] < 1e-9
    assert report["integral"] is False


def test_voronoi_cell_signature():
    cell = voronoi_cell()
    assert cell["vertex_count"] == 14
    assert cell["face_count"] == 12
    assert cell["edge_count"] == 24
    assert cell["euler_characteristic"] == 2
    assert cell["all_faces_quadrilateral"] is True
    assert cell["all_faces_rhombic"] is True
    assert cell["max_side_spread"] == 0.0
    assert cell["feasible"] is True
    assert cell["constraint_half"] == 4.0
    assert cell["max_constraint_value"] == pytest.approx(4.0, abs=1e-12)
    assert cell["orientation_flips"] == 0
    assert cell["inradius"] == pytest.approx(math.sqrt(2), abs=1e-12)
    assert cell["volume"] == pytest.approx(16.0, abs=1e-9)
    assert cell["pyramid_volume"] == pytest.approx(16.0, abs=1e-9)
    assert cell["volume_residual"] < 1e-9


def test_voronoi_cell_vertices_have_two_radii():
    radii = sorted(math.dist((0.0, 0.0, 0.0), vertex) for vertex in voronoi_cell()["vertices"])
    assert len(radii) == 14
    assert all(abs(radius - math.sqrt(3)) < 1e-9 for radius in radii[:8])
    assert all(abs(radius - 2.0) < 1e-12 for radius in radii[8:])


def test_voronoi_cell_volumes():
    assert voronoi_cell_volume() == pytest.approx(16.0, abs=1e-9)
    assert voronoi_cell_tetravolume() == pytest.approx(6.0, abs=1e-9)


def test_packing_density_is_the_ivm_density():
    report = packing_density()
    closed = math.pi / (3 * math.sqrt(2))
    assert report["density"] == pytest.approx(closed, abs=1e-15)
    assert report["closed_form"] == closed
    assert report["density_residual"] < 1e-15
    assert report["cell_tetravolume"] == pytest.approx(6.0, abs=1e-9)
    assert report["ball_tetravolume"] == pytest.approx(math.pi * math.sqrt(2), abs=1e-9)
    assert report["inradius_is_half_the_prime_vector"] is True
    assert report["cell_tetravolume_matches_table"] is True
    assert report["vertices"] == 14
    assert report["faces"] == 12
    assert report["edges"] == 24
    assert report["volume_residual"] < 1e-9


def test_angle_defects_of_the_carrier_are_a_third_at_every_vertex():
    carrier = icosahedral_carrier()
    report = angle_defects(carrier.faces, 12)
    assert report["faces"] == 20
    assert report["per_vertex"] == (Fraction(1, 3),) * 12
    assert report["total"] == Fraction(4, 1)
    assert report["faces_per_vertex"] == (5,) * 12
    assert report["euler_characteristic_times_two"] == 4
    euler_characteristic = Fraction(report["euler_characteristic_times_two"], 2)
    assert report["total"] == 2 * euler_characteristic


def test_quadratic_norm_of_the_tetravolume_squared_parts():
    parts = solids.tetravolume_ratio_squared_exact()
    assert parts == (Fraction(175), Fraction(75))
    assert quadratic_norm(parts) == Fraction(2500)
    assert quadratic_norm((Fraction(175), Fraction(75))) == Fraction(2500)


def test_tetravolume_galois_norm_is_fifty():
    assert tetravolume_galois_norm() == Fraction(50)
    assert tetravolume_galois_norm() * tetravolume_galois_norm() == quadratic_norm(
        solids.tetravolume_ratio_squared_exact()
    )
    assert quadratic_norm(solids.tetravolume_ratio_squared_exact()) == 50**2
