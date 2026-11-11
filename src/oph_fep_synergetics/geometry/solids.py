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

__all__ = [
    "icosahedron_tetravolume",
    "icosahedron_volume",
    "orient_faces_outward",
    "polyhedron_volume",
    "regular_tetrahedron_volume",
    "shortest_edge",
    "tetravolume_ratio_squared_exact",
]

Vector = tuple[float, float, float]


def _subtract(left: Sequence[float], right: Sequence[float]) -> Vector:
    return (left[0] - right[0], left[1] - right[1], left[2] - right[2])


def _add(left: Sequence[float], right: Sequence[float]) -> Vector:
    return (left[0] + right[0], left[1] + right[1], left[2] + right[2])


def _scale(vector: Sequence[float], factor: float) -> Vector:
    return (vector[0] * factor, vector[1] * factor, vector[2] * factor)


def _dot(left: Sequence[float], right: Sequence[float]) -> float:
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2]


def _cross(left: Sequence[float], right: Sequence[float]) -> Vector:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _norm(vector: Sequence[float]) -> float:
    return math.sqrt(_dot(vector, vector))


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
            distance = _norm(_subtract(vertices[i], vertices[j]))
            if 1e-12 < distance < least:
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
        normal = _cross(_subtract(points[b], points[a]), _subtract(points[c], points[a]))
        outward = _subtract(
            _scale(_add(_add(points[a], points[b]), points[c]), 1.0 / 3.0), centre
        )
        if _dot(normal, outward) < 0.0:
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
        left = _subtract(points[a], centre)
        middle = _subtract(points[b], centre)
        right = _subtract(points[c], centre)
        if _dot(_cross(left, middle), right) < 0.0:
            raise ValueError("face orientation points inward; orient the faces outward first")
        total += _dot(left, _cross(middle, right))
    return abs(total) / 6.0


def regular_tetrahedron_volume(edge: float) -> float:
    """Return the regular tetrahedron's volume, ``sqrt(2) / 12 * edge^3``."""
    return math.sqrt(2.0) / 12.0 * edge**3


def icosahedron_volume(edge: float) -> float:
    """Return the regular icosahedron's volume, ``5/12 * (3 + sqrt 5) * edge^3``."""
    return 5.0 / 12.0 * (3.0 + math.sqrt(5.0)) * edge**3


def tetravolume_ratio_squared_exact() -> tuple[Fraction, Fraction]:
    """Return ``(175, 75)``: the squared tetravolume is exactly ``175 + 75 sqrt 5``.

    Exact in the rationals and in the quadratic field: the tetravolume ratio is
    ``5 (3 + sqrt 5) / sqrt 2``, whose square is ``(25 / 2) (14 + 6 sqrt 5)``,
    using ``(3 + sqrt 5)^2 = 14 + 6 sqrt 5`` rather than a decimal expansion.
    """
    base = (Fraction(14, 1), Fraction(6, 1))
    scale = Fraction(25, 2)
    return (scale * base[0], scale * base[1])


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
    lengths = [_norm(_subtract(vertices[a], vertices[b])) for a, b in edge_list]
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
        "integral": False,
    }
