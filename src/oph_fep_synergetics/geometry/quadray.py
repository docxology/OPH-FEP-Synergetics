"""Quadray coordinates on the A3 / IVM lattice.

The quadray representation stores four integers ``(a, b, c, d)`` modulo the
common offset ``q ~ q + t * (1, 1, 1, 1)``.  The canonical representative has
minimum component zero.  The default embedding maps quadray coordinates to
:math:`\\mathbb{R}^3` with rows

``(1, -1, -1, 1)``, ``(1, 1, -1, -1)``, ``(1, -1, 1, -1)``,

whose columns are the four tetrahedron directions.  The Gram matrix of the
columns is ``4I - J``, so the Euclidean inner product of two canonical
quadrays is ``4 * sum(q_i r_i) - sum(q) * sum(r)``.

Provenance: the conventions, the embedding, the neighbour moves and the
native tetra-volume rule follow ``docxology/QuadMath`` (Apache-2.0),
``src/quadmath/core/quadray.py`` and ``src/quadmath/lattice/omni_numbering.py``.
The commit actually read is pinned once in ``data/sources.json`` and reported
by ``oph-fep-synergetics sources``; it is deliberately not restated here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations
from typing import Iterable, Sequence

#: Default 3x4 embedding: rows are X, Y, Z; columns are the tetrahedron directions.
EMBEDDING: tuple[tuple[int, int, int, int], ...] = (
    (1, -1, -1, 1),
    (1, 1, -1, -1),
    (1, -1, 1, -1),
)


@dataclass(frozen=True)
class Quadray:
    """A four-integer quadray, canonically normalised on request."""

    a: int
    b: int
    c: int
    d: int

    def __post_init__(self) -> None:
        for name in ("a", "b", "c", "d"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise TypeError("Quadray components must be integers")
            if isinstance(value, float):
                if not value.is_integer():
                    raise TypeError("Quadray components must be integers")
                object.__setattr__(self, name, int(value))

    def as_tuple(self) -> tuple[int, int, int, int]:
        """Return ``(a, b, c, d)``."""
        return (self.a, self.b, self.c, self.d)

    def normalize(self) -> "Quadray":
        """Return the canonical representative with minimum component zero."""
        shift = min(self.a, self.b, self.c, self.d)
        return Quadray(self.a - shift, self.b - shift, self.c - shift, self.d - shift)

    def __add__(self, other: "Quadray") -> "Quadray":
        return Quadray(
            self.a + other.a, self.b + other.b, self.c + other.c, self.d + other.d
        )

    def __sub__(self, other: "Quadray") -> "Quadray":
        return Quadray(
            self.a - other.a, self.b - other.b, self.c - other.c, self.d - other.d
        )

    def scale(self, factor: int) -> "Quadray":
        """Return ``factor * q`` component-wise."""
        return Quadray(
            factor * self.a, factor * self.b, factor * self.c, factor * self.d
        )

    def __iter__(self):
        return iter(self.as_tuple())


def _components(point: Quadray | Sequence[int]) -> tuple[int, int, int, int]:
    values = point.as_tuple() if isinstance(point, Quadray) else tuple(point)
    if len(values) != 4:
        raise ValueError("a quadray point must have four components")
    return (int(values[0]), int(values[1]), int(values[2]), int(values[3]))


def _rows(embedding: Iterable[Iterable[int]]) -> tuple[tuple[int, ...], ...]:
    rows = tuple(tuple(int(value) for value in row) for row in embedding)
    if len(rows) != 3 or any(len(row) != 4 for row in rows):
        raise ValueError("embedding must have shape (3, 4)")
    return rows


def to_xyz(
    point: Quadray | Sequence[int],
    embedding: Iterable[Iterable[int]] = EMBEDDING,
) -> tuple[float, float, float]:
    """Map a quadray to :math:`\\mathbb{R}^3` by the 3x4 embedding.

    The default embedding removes the common gauge first.  Custom embeddings
    keep the ordinary linear map, matching QuadMath's documented boundary.
    """
    rows = _rows(embedding)
    values = _components(point)
    if rows == EMBEDDING:
        shift = min(values)
        values = tuple(value - shift for value in values)
    return tuple(  # type: ignore[return-value]
        float(sum(weight * value for weight, value in zip(row, values)))
        for row in rows
    )


def _round_half_up(value: float) -> int:
    """Round to nearest integer, ties toward +infinity."""
    lower = math.floor(value)
    return lower + int(value - lower >= 0.5)


def from_xyz(
    x: float,
    y: float,
    z: float,
    embedding: Iterable[Iterable[int]] = EMBEDDING,
) -> Quadray:
    """Return the canonical lattice quadray whose embedding is nearest ``(x, y, z)``.

    For the default gauge-compatible embedding this is the minimum-norm
    preimage (component sum zero) rounded half-up and normalised.  It is a
    lattice projection, not an exact inverse for arbitrary off-lattice points.
    """
    rows = _rows(embedding)
    scale = max(abs(value) for row in rows for value in row)
    if scale == 0:
        raise ValueError("embedding is singular")
    scaled = [[value / scale for value in row] for row in rows]
    xyz = (x / scale, y / scale, z / scale)
    gauge_compatible = all(sum(row) == 0 for row in scaled)
    if gauge_compatible:
        # ker(M) is the all-ones direction, so sum(q) = 0 is the minimum-norm solve.
        system = [list(row) for row in scaled] + [[1.0, 1.0, 1.0, 1.0]]
        rhs = [xyz[0], xyz[1], xyz[2], 0.0]
        _solve_in_place(system, rhs)
        solution = rhs
        shift = min(solution)
        return Quadray(*(_round_half_up(value - shift) for value in solution)).normalize()
    # General full-rank embedding: least-squares minimum-norm preimage.
    solution = _least_squares(scaled, xyz)
    return Quadray(*(_round_half_up(value) for value in solution)).normalize()


def _solve_in_place(matrix: list[list[float]], rhs: list[float]) -> None:
    """Solve a small dense linear system by Gaussian elimination with pivoting."""
    size = len(matrix)
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(matrix[row][column]))
        if abs(matrix[pivot][column]) <= 1e-15:
            raise ValueError("embedding is singular")
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            rhs[column], rhs[pivot] = rhs[pivot], rhs[column]
        diagonal = matrix[column][column]
        for row in range(column + 1, size):
            factor = matrix[row][column] / diagonal
            if factor == 0.0:
                continue
            for k in range(column, size):
                matrix[row][k] -= factor * matrix[column][k]
            rhs[row] -= factor * rhs[column]
    for row in range(size - 1, -1, -1):
        total = rhs[row] - sum(
            matrix[row][column] * rhs[column] for column in range(row + 1, size)
        )
        rhs[row] = total / matrix[row][row]


def _least_squares(matrix: list[list[float]], xyz: tuple[float, float, float]) -> list[float]:
    """Minimum-norm least-squares solution of ``matrix @ q = xyz``."""
    # Normal equations are acceptable here because the embedding is tiny and
    # the fallback is only used for non-gauge custom maps.
    normal = [[sum(matrix[r][i] * matrix[r][j] for r in range(3)) for j in range(4)] for i in range(4)]
    rhs = [sum(matrix[r][i] * xyz[r] for r in range(3)) for i in range(4)]
    # Solve the 4x4 normal equations; the system is singular along the gauge
    # kernel, so add the minimum-norm regulariser.
    for i in range(4):
        normal[i][i] += 1e-12
    _solve_in_place(normal, rhs)
    return rhs


def tetra_volume(
    p0: Quadray | Sequence[int],
    p1: Quadray | Sequence[int],
    p2: Quadray | Sequence[int],
    p3: Quadray | Sequence[int],
) -> Fraction:
    """Exact IVM tetra-volume ``|det(projected edges)| / 4``.

    ``project(v) = (v.a - v.d, v.b - v.d, v.c - v.d)``.  The four axis
    directions are the unit tetrahedron (projected determinant four, volume
    one); the origin together with three axes is the primitive simplex
    (projected determinant one, volume ``1/4``).
    """
    origin = _components(p0)
    edges = []
    for point in (p1, p2, p3):
        values = _components(point)
        edge = tuple(value - base for value, base in zip(values, origin))
        edges.append((edge[0] - edge[3], edge[1] - edge[3], edge[2] - edge[3]))
    a, b, c = edges
    determinant = (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )
    return Fraction(abs(determinant), 4)


def dot(
    left: Quadray | Sequence[int],
    right: Quadray | Sequence[int],
    embedding: Iterable[Iterable[int]] = EMBEDDING,
) -> int:
    """Exact Euclidean inner product under the embedding.

    For the default embedding the closed form is
    ``4 * sum(l_i r_i) - sum(l) * sum(r)``.
    """
    a = _components(left)
    b = _components(right)
    rows = _rows(embedding)
    if rows == EMBEDDING:
        return 4 * sum(x * y for x, y in zip(a, b)) - sum(a) * sum(b)
    return sum(
        sum(weight * value for weight, value in zip(row, a))
        * sum(weight * value for weight, value in zip(row, b))
        for row in rows
    )


def squared_norm(point: Quadray | Sequence[int]) -> int:
    """Exact squared Euclidean norm under the default embedding."""
    return dot(point, point)


def squared_distance(
    left: Quadray | Sequence[int],
    right: Quadray | Sequence[int],
) -> int:
    """Exact squared Euclidean distance under the default embedding."""
    a = _components(left)
    b = _components(right)
    difference = Quadray(*(x - y for x, y in zip(a, b)))
    return squared_norm(difference)


def is_ivm_site(point: Quadray | Sequence[int]) -> bool:
    """Whether the canonical quadray is an IVM lattice site (sum divisible by 4)."""
    return sum(_components(point)) % 4 == 0


def is_octahedron_site(point: Quadray | Sequence[int]) -> bool:
    """Whether the canonical quadray lies on the octahedral-hole sublattice.

    Those sites satisfy ``sum(q) = 2 (mod 4)``, so they interleave with the IVM
    sites (``sum(q) = 0 (mod 4)``) and are never IVM sites themselves.  The six
    vertices ``e_i + e_j`` are the octahedron's corners and have tetravolume
    four at the unit tetrahedron's edge length.
    """
    return sum(_components(point)) % 4 == 2


def neighbor_moves() -> tuple[tuple[int, int, int, int], ...]:
    """The 12 distinct permutations of ``(2, 1, 1, 0)``, sorted."""
    return tuple(sorted(set(permutations((2, 1, 1, 0)))))


def neighbor_sites() -> tuple[Quadray, ...]:
    """The 12 nearest IVM neighbours as canonical quadrays."""
    return tuple(Quadray(*move) for move in neighbor_moves())


def scale(point: Quadray | Sequence[int], factor: int) -> Quadray:
    """Return ``factor * point`` in the canonical gauge."""
    return Quadray(*_components(point)).scale(factor).normalize()
