"""Polyhedron volumes in tetrahedral units, and the icosahedron's irrational one.

Synergetics measures volumes in tetrahedra rather than in cubes: the prime
vector ``d`` (the IVM's nearest-neighbour distance, and the edge of the
tetrahedron whose four vertices are the four quadray axes) is the edge of the
unit tetravolume.  A solid's **tetravolume** is therefore its volume divided by
the volume of the regular tetrahedron with the same edge, and in that accounting
the IVM's solids are integral -- tetrahedron 1, cube 3, octahedron 4, rhombic
dodecahedron 6, cuboctahedron 20, truncated octahedron 96.

The icosahedron is not.  It is the carrier's own solid, and its tetravolume is

    T = (5 / 2) (3 sqrt 2 + sqrt 10) = 18.512295868219...

whose square is exactly ``175 + 75 sqrt 5``.  The irrationality is the price of
five-fold symmetry, and it is why the icosahedron cannot appear in a table that
is otherwise integral.

Volumes here are measured, not quoted: :func:`orient_faces_outward` fixes the
orientation of a face list that was enumerated without one -- ten of the
icosahedron's twenty faces come out inward-facing, and an unoriented signed sum
over the twenty faces cancels instead of measuring anything -- and
:func:`polyhedron_volume` applies the divergence theorem to the oriented
surface.  The closed form is compared against that measurement rather than
substituted for it.

This module is deliberately free of ``numpy``: the geometry of the quadrays,
the IVM and the subdivision is plain arithmetic over floats and integers.
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Iterable, Sequence

from .subdivision import icosahedron_edges, icosahedron_faces, icosahedron_vertices

from .vectors import (
    TIGHT_TOLERANCE,
    Vector,
    add,
    cross,
    dot,
    norm,
    scale,
    subtract,
)

__all__ = [
    "angle_defects",
    "icosahedron_tetravolume",
    "icosahedron_volume",
    "orient_faces_outward",
    "polyhedron_volume",
    "quadratic_norm",
    "regular_tetrahedron_volume",
    "shortest_edge",
    "tetravolume_galois_norm",
    "tetravolume_ratio_squared_exact",
]


def _centroid(vertices: Sequence[Sequence[float]]) -> Vector:
    count = len(vertices)
    return (
        sum(vertex[0] for vertex in vertices) / count,
        sum(vertex[1] for vertex in vertices) / count,
        sum(vertex[2] for vertex in vertices) / count,
    )


def shortest_edge(vertices: Sequence[Sequence[float]]) -> float:
    """Return the least distance between two vertices.

    For a regular solid the least pairwise distance is the edge length, so this
    returns the edge without being told which pairs are edges.
    """
    least = math.inf
    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            distance = norm(subtract(vertices[i], vertices[j]))
            if TIGHT_TOLERANCE < distance < least:
                least = distance
    if not math.isfinite(least):
        raise ValueError("a vertex set needs at least two distinct vertices")
    return least


def orient_faces_outward(
    vertices: Sequence[Sequence[float]],
    faces: Iterable[Sequence[int]],
) -> tuple[tuple[tuple[int, int, int], ...], int]:
    """Return the faces oriented so every normal points away from the centroid.

    A face list produced by enumerating triangles carries no orientation, and an
    unoriented signed sum over such a list cancels.  Each face is rotated so that
    its normal points away from the solid's centroid, and the number of faces
    that had to be flipped is returned beside the result.
    """
    points = [tuple(float(value) for value in vertex) for vertex in vertices]
    centre = _centroid(points)
    oriented: list[tuple[int, int, int]] = []
    flips = 0
    for face in faces:
        a, b, c = (int(index) for index in face)
        normal = cross(subtract(points[b], points[a]), subtract(points[c], points[a]))
        outward = subtract(
            scale(add(add(points[a], points[b]), points[c]), 1.0 / 3.0), centre
        )
        if dot(normal, outward) < 0.0:
            b, c = c, b
            flips += 1
        oriented.append((a, b, c))
    return tuple(oriented), flips


def polyhedron_volume(
    vertices: Sequence[Sequence[float]],
    faces: Iterable[Sequence[int]],
) -> float:
    """Return the volume enclosed by an oriented closed triangulated surface.

    The divergence theorem turns the volume into a sum over faces; every face
    must be counter-clockwise seen from outside, which
    :func:`orient_faces_outward` establishes.  A face whose normal points inward
    raises instead of contributing a negative term, so a mixed orientation cannot
    cancel silently into a wrong number.
    """
    points = [tuple(float(value) for value in vertex) for vertex in vertices]
    centre = _centroid(points)
    total = 0.0
    for face in faces:
        a, b, c = (int(index) for index in face)
        left = subtract(points[a], centre)
        middle = subtract(points[b], centre)
        right = subtract(points[c], centre)
        if dot(cross(left, middle), right) < 0.0:
            raise ValueError("face orientation points inward; orient the faces outward first")
        total += dot(left, cross(middle, right))
    return abs(total) / 6.0


def regular_tetrahedron_volume(edge: float) -> float:
    """Return the regular tetrahedron's volume, ``sqrt(2) / 12 * edge^3``."""
    return math.sqrt(2.0) / 12.0 * edge**3


def icosahedron_volume(edge: float) -> float:
    """Return the regular icosahedron's volume, ``5/12 * (3 + sqrt 5) * edge^3``."""
    return 5.0 / 12.0 * (3.0 + math.sqrt(5.0)) * edge**3


def tetravolume_ratio_squared_exact() -> tuple[Fraction, Fraction]:
    """Return ``(175, 75)``: the squared tetravolume is exactly ``175 + 75 sqrt 5``.

    Computed, not quoted: :func:`~oph_fep_synergetics.geometry.exact.icosahedron_tetravolume_exact`
    measures the icosahedron's volume in the field ``Q(sqrt 5)`` from its exact
    vertices, the tetravolume is ``sqrt 2`` times an element of that field, and
    its square is returned here as the pair of rational coefficients.
    """
    from .exact import icosahedron_tetravolume_exact

    square = icosahedron_tetravolume_exact()["square"]
    return square.parts()  # type: ignore[union-attr]


def icosahedron_tetravolume() -> dict[str, object]:
    """Measure the icosahedron's tetravolume and compare it with the exact value.

    The measurement orients the twenty faces outward (ten of them need it),
    applies the divergence theorem to the oriented surface and divides by the
    equal-edge tetrahedron's volume.  The comparison is therefore between a
    measurement and an algebraic identity, not between two spellings of one
    formula.
    """
    vertices = icosahedron_vertices()
    faces = icosahedron_faces()
    edge_list = icosahedron_edges()
    lengths = [norm(subtract(vertices[a], vertices[b])) for a, b in edge_list]
    oriented, flips = orient_faces_outward(vertices, faces)
    volume = polyhedron_volume(vertices, oriented)
    edge = sum(lengths) / len(lengths)
    tetraunit = regular_tetrahedron_volume(edge)
    ratio = volume / tetraunit
    closed = 5.0 / 2.0 * (3.0 * math.sqrt(2.0) + math.sqrt(10.0))
    exact_parts = tetravolume_ratio_squared_exact()
    exact_squared = float(exact_parts[0]) + float(exact_parts[1]) * math.sqrt(5.0)
    return {
        "vertices": len(vertices),
        "edges": len(edge_list),
        "faces": len(faces),
        "faces_reoriented": flips,
        "edge": edge,
        "edge_spread": max(lengths) - min(lengths),
        "volume": volume,
        "closed_form_volume": icosahedron_volume(edge),
        "tetraunit": tetraunit,
        "tetravolume": ratio,
        "closed_form_tetravolume": closed,
        "tetravolume_residual": abs(ratio - closed),
        "tetravolume_squared": ratio * ratio,
        "exact_squared": exact_squared,
        "exact_squared_parts": (int(exact_parts[0]), int(exact_parts[1])),
        "exact_squared_residual": abs(ratio * ratio - exact_squared),
        # Decided by exact arithmetic, not asserted: the tetravolume is sqrt 2 times
        # an element of Q(sqrt 5) whose sqrt 5 part is nonzero.
        "integral": _exact_tetravolume_is_integral(),
    }


def _exact_tetravolume_is_integral() -> bool:
    from .exact import icosahedron_tetravolume_exact

    coefficient = icosahedron_tetravolume_exact()["sqrt2_coefficient"]
    # sqrt 2 * (a + b sqrt 5) is rational only if a = b = 0; it is an integer
    # only then.  A nonzero coefficient therefore means "not integral".
    return coefficient.a == 0 and coefficient.b == 0  # type: ignore[union-attr]


def angle_defects(
    faces: Iterable[Sequence[int]], vertex_count: int
) -> dict[str, object]:
    """Return the combinatorial angle defect of a triangulated closed surface.

    At a vertex where ``k`` equilateral triangles meet, the surface angle sum is
    ``k * pi / 3``, so the defect is ``2 pi - k pi / 3``. Its total is ``2 pi``
    times the Euler characteristic -- Gauss-Bonnet, read off a face list rather
    than quoted. Defects are returned as exact multiples of ``pi``, which is the
    form in which the carrier's is a rational: ``1/3`` at every port and ``4`` in
    total, against ``2 * chi = 4`` for a triangulated sphere.
    """
    triangles = [tuple(face) for face in faces]
    for face in triangles:
        if len(face) != 3 or len(set(face)) != 3:
            raise AssertionError("the defect formula here assumes a triangle-only surface")
    incident = [0] * vertex_count
    for face in triangles:
        for vertex in face:
            if not 0 <= vertex < vertex_count:
                raise AssertionError("a face names a vertex outside the surface")
            incident[vertex] += 1
    per_vertex = tuple(Fraction(2, 1) - Fraction(count, 3) for count in incident)
    return {
        "faces": len(triangles),
        "faces_per_vertex": tuple(incident),
        "per_vertex": per_vertex,
        "total": sum(per_vertex, Fraction(0)),
        "euler_characteristic_times_two": 2 * vertex_count - len(triangles),
    }


def quadratic_norm(parts: Sequence[Fraction]) -> Fraction:
    """Return the norm ``a^2 - 5 b^2`` of ``a + b sqrt 5``.

    The conjugation ``sqrt 5 -> -sqrt 5`` is the Galois symmetry of the field the
    icosahedron's volume lives in, and the norm is the integer it lands on.
    """
    return parts[0] * parts[0] - 5 * parts[1] * parts[1]


def tetravolume_galois_norm() -> Fraction:
    """Return ``50``: the icosahedron's tetravolume times its conjugate's.

    The tetravolume is ``T = sqrt 2 * c`` with ``c = 15/2 + (5/2) sqrt 5`` exactly,
    and the conjugation ``sqrt 5 -> -sqrt 5`` fixing ``sqrt 2`` sends it to
    ``sqrt 2 * c'``, so ``T T' = 2 N(c) = 2 (225/4 - 125/4) = 50``.  The value is
    computed in exact field arithmetic and cross-checked against the square: the
    norm of ``T^2 = 175 + 75 sqrt 5`` must be ``50^2``.  A counting invariant over
    the same carrier -- the spanning-tree number -- is a rational integer and so
    is fixed by the same conjugation; the volume is not, which is the arithmetic
    difference between the two kinds of claim this package makes.
    """
    from .exact import icosahedron_tetravolume_exact

    report = icosahedron_tetravolume_exact()
    norm = Fraction(report["norm"])  # type: ignore[arg-type]
    if norm * norm != report["square_norm"]:
        raise AssertionError("the norm of the squared tetravolume must be the square of the norm")
    return norm
