"""Exact volumes: the rationals for the IVM's solids, ``Q(sqrt 5)`` for the carrier's.

Synergetics' volume table is a claim about exact numbers, so this module
computes it without floating point.  Every solid is given by an explicit vertex
set; its facets are found by an exact supporting-plane test, ordered
counter-clockwise as seen from outside by an exact angular comparison, and its
volume is a sum of signed tetrahedra whose signs are checked rather than
assumed.  Nothing is transcribed: the table entries ``1, 3, 4, 6, 20, 96`` come
out of the vertex sets, and the icosahedron's tetravolume comes out as an
element of ``sqrt(2) * Q(sqrt 5)``.

Two number systems are enough.

* The IVM's solids have integer vertices in the default quadray embedding,
  where the prime vector (the nearest-neighbour distance) has squared length
  ``8``.  The regular tetrahedron of that edge has volume ``8/3`` exactly, so
  every tetravolume in the table is a rational number: three eighths of a
  rational volume.
* The regular icosahedron with edge ``2`` has vertices ``(0, +-1, +-phi)`` and
  their cyclic permutations, with ``phi = (1 + sqrt 5) / 2``, so its volume is
  an element of ``Q(sqrt 5)``.  The regular tetrahedron of edge ``2`` has volume
  ``2 sqrt 2 / 3``, so the tetravolume is ``sqrt 2`` times an element of
  ``Q(sqrt 5)``, and the Galois conjugation ``sqrt 5 -> -sqrt 5`` (fixing
  ``sqrt 2``) acts on it exactly.

Provenance: the table being reproduced is QuadMath's polyhedra table
(``docxology/QuadMath``, Apache-2.0); the commit actually read is pinned once in
``data/sources.json``.  The constructions here are this package's own.
"""

from __future__ import annotations

import functools
import itertools
from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Sequence, Union

__all__ = [
    "IVM_TETRA_UNIT",
    "PHI_EXACT",
    "Surd5",
    "convex_facets",
    "convex_volume",
    "derived_volume_ratios",
    "exact_icosahedron_vertices",
    "icosahedron_tetravolume_exact",
    "ivm_solid_vertices",
    "kissing_arrangements",
]

Rational = Union[int, Fraction]


@dataclass(frozen=True, eq=False)
class Surd5:
    """An element ``a + b sqrt 5`` of the quadratic field ``Q(sqrt 5)``.

    The coefficients are exact rationals.  Signs are decided exactly, by
    comparing ``a^2`` with ``5 b^2``, which can only be equal when both
    coefficients vanish because ``sqrt 5`` is irrational.
    """

    a: Fraction
    b: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        object.__setattr__(self, "a", Fraction(self.a))
        object.__setattr__(self, "b", Fraction(self.b))

    @staticmethod
    def of(value: "Surd5 | Rational") -> "Surd5":
        return value if isinstance(value, Surd5) else Surd5(Fraction(value), Fraction(0))

    def __add__(self, other: "Surd5 | Rational") -> "Surd5":
        other = Surd5.of(other)
        return Surd5(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self) -> "Surd5":
        return Surd5(-self.a, -self.b)

    def __sub__(self, other: "Surd5 | Rational") -> "Surd5":
        return self + (-Surd5.of(other))

    def __rsub__(self, other: "Surd5 | Rational") -> "Surd5":
        return Surd5.of(other) - self

    def __mul__(self, other: "Surd5 | Rational") -> "Surd5":
        other = Surd5.of(other)
        return Surd5(self.a * other.a + 5 * self.b * other.b, self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def conjugate(self) -> "Surd5":
        """The Galois conjugate ``a - b sqrt 5``."""
        return Surd5(self.a, -self.b)

    def norm(self) -> Fraction:
        """The field norm ``a^2 - 5 b^2``, the product with the conjugate."""
        return self.a * self.a - 5 * self.b * self.b

    def __truediv__(self, other: "Surd5 | Rational") -> "Surd5":
        other = Surd5.of(other)
        denominator = other.norm()
        if denominator == 0:
            raise ZeroDivisionError("division by zero in Q(sqrt 5)")
        numerator = self * other.conjugate()
        return Surd5(numerator.a / denominator, numerator.b / denominator)

    def sign(self) -> int:
        """Return ``-1``, ``0`` or ``1`` exactly."""
        a_sign = (self.a > 0) - (self.a < 0)
        b_sign = (self.b > 0) - (self.b < 0)
        if b_sign == 0:
            return a_sign
        if a_sign == 0 or a_sign == b_sign:
            return b_sign
        return a_sign if self.a * self.a > 5 * self.b * self.b else b_sign

    def is_rational(self) -> bool:
        return self.b == 0

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Surd5):
            return self.a == other.a and self.b == other.b
        if isinstance(other, (int, Fraction)):
            return self.b == 0 and self.a == other
        return NotImplemented

    def __hash__(self) -> int:
        # A rational element hashes like the rational, so equal values hash alike.
        return hash(self.a) if self.b == 0 else hash((self.a, self.b))

    def __float__(self) -> float:
        return float(self.a) + float(self.b) * 5.0**0.5

    def parts(self) -> tuple[Fraction, Fraction]:
        return (self.a, self.b)

    def __str__(self) -> str:
        if self.b == 0:
            return str(self.a)
        sign = "+" if self.b > 0 else "-"
        return f"{self.a} {sign} {abs(self.b)} sqrt 5"


#: The golden ratio as an exact element of ``Q(sqrt 5)``.
PHI_EXACT = Surd5(Fraction(1, 2), Fraction(1, 2))

#: Volume of the regular tetrahedron whose edge is the IVM's prime vector
#: (squared length ``8``): ``(2 sqrt 2)^3 sqrt 2 / 12 = 8/3``.
IVM_TETRA_UNIT = Fraction(8, 3)

Scalar = Union[Fraction, Surd5]
Point = tuple[Scalar, Scalar, Scalar]


def _sign(value: Scalar | int) -> int:
    if isinstance(value, Surd5):
        return value.sign()
    return (value > 0) - (value < 0)


def _sub(left: Sequence[Scalar], right: Sequence[Scalar]) -> Point:
    return (left[0] - right[0], left[1] - right[1], left[2] - right[2])


def _dot(left: Sequence[Scalar], right: Sequence[Scalar]) -> Scalar:
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2]


def _cross(left: Sequence[Scalar], right: Sequence[Scalar]) -> Point:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _centroid(points: Sequence[Sequence[Scalar]]) -> Point:
    count = len(points)
    total: list[Scalar] = [Fraction(0), Fraction(0), Fraction(0)]
    for point in points:
        total = [total[i] + point[i] for i in range(3)]
    return (total[0] * Fraction(1, count), total[1] * Fraction(1, count), total[2] * Fraction(1, count))


def _as_exact(point: Iterable[object]) -> Point:
    values = []
    for value in point:
        if isinstance(value, Surd5):
            values.append(value)
        elif isinstance(value, float):
            if value != int(value):
                raise ValueError("an exact vertex needs integer, rational or Q(sqrt 5) coordinates")
            values.append(Fraction(int(value)))
        else:
            values.append(Fraction(value))  # type: ignore[arg-type]
    if len(values) != 3:
        raise ValueError("a vertex has three coordinates")
    return (values[0], values[1], values[2])


def convex_facets(vertices: Sequence[Sequence[object]]) -> tuple[tuple[int, ...], ...]:
    """Return the facets of the convex hull of ``vertices``, each ordered outward.

    A plane through three vertices is a facet plane exactly when every vertex
    lies on one closed side of it; the facet is the set of vertices on the
    plane.  Each facet is then ordered counter-clockwise as seen from outside, by
    an exact comparison of angles in the facet's own plane.  Every input vertex
    must be a vertex of the hull, so that a degenerate input raises instead of
    returning a smaller solid.
    """
    points = [_as_exact(vertex) for vertex in vertices]
    if len(points) < 4:
        raise ValueError("a solid needs at least four vertices")
    found: dict[frozenset[int], Point] = {}
    for i, j, k in itertools.combinations(range(len(points)), 3):
        normal = _cross(_sub(points[j], points[i]), _sub(points[k], points[i]))
        if all(_sign(component) == 0 for component in normal):
            continue
        sides = [_sign(_dot(normal, _sub(point, points[i]))) for point in points]
        if any(side > 0 for side in sides) and any(side < 0 for side in sides):
            continue
        if all(side == 0 for side in sides):
            raise ValueError("the vertices are coplanar")
        on_plane = frozenset(index for index, side in enumerate(sides) if side == 0)
        if on_plane in found:
            continue
        if any(side > 0 for side in sides):
            normal = (-normal[0], -normal[1], -normal[2])
        found[on_plane] = normal
    covered = set().union(*found)
    if covered != set(range(len(points))):
        raise ValueError("every vertex must be a vertex of the convex hull")
    ordered = [_order_facet(points, sorted(facet), normal) for facet, normal in found.items()]
    return tuple(sorted(ordered))


def _order_facet(points: Sequence[Point], facet: Sequence[int], normal: Point) -> tuple[int, ...]:
    centre = _centroid([points[index] for index in facet])
    first = _sub(points[facet[0]], centre)
    second = _cross(normal, first)

    def coordinates(index: int) -> tuple[Scalar, Scalar]:
        offset = _sub(points[index], centre)
        return (_dot(offset, first), _dot(offset, second))

    def half(index: int) -> int:
        x, y = coordinates(index)
        return 0 if _sign(y) > 0 or (_sign(y) == 0 and _sign(x) > 0) else 1

    def compare(left: int, right: int) -> int:
        if half(left) != half(right):
            return half(left) - half(right)
        x1, y1 = coordinates(left)
        x2, y2 = coordinates(right)
        return -_sign(x1 * y2 - y1 * x2)

    ordering = sorted(facet, key=functools.cmp_to_key(compare))
    start = ordering.index(min(ordering))
    return tuple(ordering[start:] + ordering[:start])


def convex_volume(vertices: Sequence[Sequence[object]]) -> Scalar:
    """Return the exact volume of the convex hull of ``vertices``.

    The hull is cut into tetrahedra joining its centroid to a fan triangulation
    of each outward-ordered facet.  Every tetrahedron's signed volume must be
    positive; a non-positive one raises, so a mis-ordered facet cannot cancel
    silently into a wrong number.
    """
    points = [_as_exact(vertex) for vertex in vertices]
    centre = _centroid(points)
    total: Scalar = Fraction(0)
    for facet in convex_facets(points):
        anchor = _sub(points[facet[0]], centre)
        for index in range(1, len(facet) - 1):
            middle = _sub(points[facet[index]], centre)
            right = _sub(points[facet[index + 1]], centre)
            signed = _dot(anchor, _cross(middle, right))
            if _sign(signed) <= 0:
                raise AssertionError("a facet tetrahedron has non-positive signed volume")
            total = total + signed
    return total * Fraction(1, 6)


def _xyz(point: Sequence[int]) -> Point:
    """The default quadray embedding, in exact integers."""
    from .quadray import EMBEDDING

    shift = min(point)
    values = [value - shift for value in point]
    return tuple(  # type: ignore[return-value]
        Fraction(sum(weight * value for weight, value in zip(row, values))) for row in EMBEDDING
    )


def ivm_solid_vertices() -> dict[str, tuple[Point, ...]]:
    """Return the six solids of Synergetics' table as exact vertex sets.

    Every set is built from the lattice rather than typed in: the tetrahedron
    from the four quadray axes, the cube from the axes and their negatives (the
    tetrahedron with its dual), the octahedron from the pairwise axis sums, the
    cuboctahedron from the twelve neighbour moves, the rhombic dodecahedron from
    the neighbour sums over the cuboctahedron's faces divided by four (the
    Voronoi cell), and the truncated octahedron from the points one third of the
    way along every edge of the octahedron scaled by three.
    """
    from .lattice import cuboctahedron_faces, cuboctahedron_neighbors

    axes = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    tetrahedron = tuple(_xyz(axis) for axis in axes)
    negatives = tuple(_xyz(tuple(-value for value in axis)) for axis in axes)
    octahedron = tuple(
        _xyz(tuple(a + b for a, b in zip(axes[i], axes[j])))
        for i, j in itertools.combinations(range(4), 2)
    )
    neighbours = tuple(_xyz(site.as_tuple()) for site in cuboctahedron_neighbors())
    faces = cuboctahedron_faces()
    cell = tuple(
        tuple(
            sum((neighbours[index][axis] for index in face), Fraction(0)) * Fraction(1, 4)
            for axis in range(3)
        )
        for key in ("triangles", "squares")
        for face in faces[key]
    )
    large = [tuple(3 * value for value in vertex) for vertex in octahedron]
    truncated: set[Point] = set()
    for first, second in itertools.permutations(large, 2):
        squared = sum((x - y) ** 2 for x, y in zip(first, second))
        if squared == 9 * 8:
            truncated.add(
                tuple(  # type: ignore[arg-type]
                    x + (y - x) * Fraction(1, 3) for x, y in zip(first, second)
                )
            )
    return {
        "tetrahedron": tetrahedron,
        "cube": tetrahedron + negatives,
        "octahedron": octahedron,
        "rhombic_dodecahedron": cell,  # type: ignore[dict-item]
        "cuboctahedron": neighbours,
        "truncated_octahedron": tuple(sorted(truncated)),
    }


def derived_volume_ratios() -> dict[str, Fraction]:
    """Return each table solid's tetravolume, computed exactly from its vertices.

    The volume is exact and the unit tetrahedron of the prime vector has volume
    ``8/3``, so each entry is a rational number with no rounding anywhere.
    """
    return {
        name: Fraction(convex_volume(vertices)) / IVM_TETRA_UNIT  # type: ignore[arg-type]
        for name, vertices in ivm_solid_vertices().items()
    }


def exact_icosahedron_vertices() -> tuple[Point, ...]:
    """The regular icosahedron of edge ``2`` in ``Q(sqrt 5)``, in the float module's order."""
    one = Surd5(1)
    zero = Surd5(0)
    vertices: list[Point] = []
    for first in (-1, 1):
        for second in (-1, 1):
            vertices.append((zero, one * first, PHI_EXACT * second))
            vertices.append((one * first, PHI_EXACT * second, zero))
            vertices.append((PHI_EXACT * first, zero, one * second))
    return tuple(vertices)


def icosahedron_tetravolume_exact() -> dict[str, object]:
    """Return the icosahedron's tetravolume as ``sqrt 2`` times an exact element.

    With edge ``2`` the volume ``V`` lies in ``Q(sqrt 5)`` and the unit
    tetrahedron has volume ``2 sqrt 2 / 3``, so the tetravolume is
    ``T = sqrt 2 * (3 V / 4)``.  The square ``T^2 = 2 (3 V / 4)^2`` lies in
    ``Q(sqrt 5)``, and the conjugation ``sqrt 5 -> -sqrt 5`` (fixing ``sqrt 2``)
    sends ``T`` to ``sqrt 2 * (3 V / 4)'``, so ``T T' = 2 N(3 V / 4)``.
    """
    vertices = exact_icosahedron_vertices()
    squared_lengths = [
        _dot(_sub(first, second), _sub(first, second))
        for first, second in itertools.combinations(vertices, 2)
    ]
    edge_squared = min(squared_lengths, key=float)
    edges = sum(1 for value in squared_lengths if value == edge_squared)
    facets = convex_facets(vertices)
    volume = Surd5.of(convex_volume(vertices))
    coefficient = volume * Fraction(3, 4)
    square = coefficient * coefficient * 2
    conjugate = coefficient.conjugate()
    return {
        "edge_squared": edge_squared,
        "edges": edges,
        "facets": len(facets),
        "facet_sizes": tuple(sorted({len(facet) for facet in facets})),
        "volume": volume,
        "sqrt2_coefficient": coefficient,
        "tetravolume": float(coefficient) * 2.0**0.5,
        "square": square,
        "conjugate_sqrt2_coefficient": conjugate,
        "conjugate_tetravolume": float(conjugate) * 2.0**0.5,
        "norm": 2 * coefficient.norm(),
        "square_norm": square.norm(),
    }


def kissing_arrangements() -> dict[str, object]:
    """Compare the two twelve-sphere kissing arrangements exactly.

    Put the twelve sphere centres on the sphere of radius ``2r`` around a
    central sphere of radius ``r``.  The neighbours touch each other exactly
    when the arrangement's edge equals its circumradius.  The ratio
    ``edge^2 / circumradius^2`` is computed for both arrangements: the IVM's
    cuboctahedron gives exactly ``1`` (rigid: every neighbour touches four
    others), and the icosahedron gives exactly ``2 - 2 sqrt 5 / 5``, which is
    greater than one, so its twelve spheres touch the centre and none of each
    other.  That slack is what lets Fuller's jitterbug move between the two.
    """
    icosahedron = exact_icosahedron_vertices()
    ico_radius = Surd5.of(_dot(icosahedron[0], icosahedron[0]))
    ico_edge = min(
        (
            Surd5.of(_dot(_sub(first, second), _sub(first, second)))
            for first, second in itertools.combinations(icosahedron, 2)
        ),
        key=float,
    )
    ico_ratio = ico_edge / ico_radius
    neighbours = ivm_solid_vertices()["cuboctahedron"]
    cubo_radius = Fraction(_dot(neighbours[0], neighbours[0]))  # type: ignore[arg-type]
    cubo_edge = min(
        Fraction(_dot(_sub(first, second), _sub(first, second)))  # type: ignore[arg-type]
        for first, second in itertools.combinations(neighbours, 2)
    )
    return {
        "cuboctahedron_edge_over_radius_squared": cubo_edge / cubo_radius,
        "icosahedron_edge_over_radius_squared": ico_ratio,
        "icosahedron_ratio_float": float(ico_ratio),
        "icosahedron_has_slack": ico_ratio.sign() == 1 and (ico_ratio - 1).sign() == 1,
        "cuboctahedron_is_rigid": cubo_edge == cubo_radius,
    }
