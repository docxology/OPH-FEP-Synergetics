"""The A3 / IVM close-packing lattice and its cuboctahedral shell.

The isotropic vector matrix (IVM) is the FCC close packing in quadray
coordinates.  Shell zero is the single centre site; every frequency
``k >= 1`` carries ``10 k^2 + 2`` sites, and the first shell is the
twelve-around-one cuboctahedron (vector equilibrium).  The
same quadratic form counts the vertices of a Goldberg polyhedron of frequency
``k``, so the IVM shell sequence and the Goldberg vertex sequence are one
sequence seen through two constructions.

Provenance: ``docxology/QuadMath`` (Apache-2.0),
``src/quadmath/lattice/omni_numbering.py`` and the polyhedra table in
``quadmath/markdown/03_quadray_methods.md``.  The commit actually read is
pinned once in ``data/sources.json`` and reported by
``oph-fep-synergetics sources``; it is deliberately not restated here.
"""

from __future__ import annotations

import itertools
from fractions import Fraction
from typing import Iterator

from .quadray import (
    EMBEDDING,
    Quadray,
    neighbor_moves,
    neighbor_sites,
    squared_distance,
    tetra_volume,
    to_xyz,
)

#: Deepest frequency shell this module enumerates.  QuadMath supports 32.
MAX_SHELL: int = 32

#: Neighbour squared distance in the default embedding.
NEIGHBOR_SQUARED_DISTANCE: int = 8

#: Declared tetrahedral volume ratios from the QuadMath polyhedra table.
VOLUME_RATIOS: dict[str, Fraction] = {
    "tetrahedron": Fraction(1, 1),
    "cube": Fraction(3, 1),
    "octahedron": Fraction(4, 1),
    "rhombic_dodecahedron": Fraction(6, 1),
    "cuboctahedron": Fraction(20, 1),
    "truncated_octahedron": Fraction(96, 1),
}

_SHELL_CACHE: dict[int, tuple[tuple[Quadray, ...], ...]] = {}


def shell_count(k: int) -> int:
    """Return the number of IVM sites on frequency shell ``k``.

    Shell zero is the lone centre; ``k >= 1`` carries ``10 k^2 + 2`` sites.
    """
    if k < 0:
        raise ValueError("shell index must be nonnegative")
    return 1 if k == 0 else 10 * k * k + 2


def cumulative_count(k: int) -> int:
    """Return the total number of IVM sites through shell ``k``."""
    if k < 0:
        raise ValueError("shell index must be nonnegative")
    return 1 + sum(10 * j * j + 2 for j in range(1, k + 1))


def _build_shells(max_shell: int) -> tuple[tuple[Quadray, ...], ...]:
    """Enumerate shells ``0..max_shell`` by breadth-first expansion."""
    if max_shell < 0:
        raise ValueError("max_shell must be nonnegative")
    origin = Quadray(0, 0, 0, 0)
    seen: set[Quadray] = {origin}
    shells: list[tuple[Quadray, ...]] = [(origin,)]
    frontier: set[Quadray] = {origin}
    moves = neighbor_sites()
    for _ in range(max_shell):
        next_frontier: set[Quadray] = set()
        for site in frontier:
            for move in moves:
                candidate = (site + move).normalize()
                if candidate not in seen:
                    seen.add(candidate)
                    next_frontier.add(candidate)
        shells.append(tuple(sorted(next_frontier, key=lambda q: q.as_tuple())))
        frontier = next_frontier
    return tuple(shells)


def _shells(max_shell: int) -> tuple[tuple[Quadray, ...], ...]:
    cached = _SHELL_CACHE.get(max_shell)
    if cached is not None:
        return cached
    built = _build_shells(max_shell)
    _SHELL_CACHE[max_shell] = built
    return built


def generate_shell(k: int) -> tuple[Quadray, ...]:
    """Return the canonical sites on frequency shell ``k`` in lexicographic order."""
    if k < 0:
        raise ValueError("shell index must be nonnegative")
    if k > MAX_SHELL:
        raise ValueError(f"shell index above the supported bound {MAX_SHELL}")
    return _shells(k)[k]


def sites_through_shell(k: int) -> tuple[Quadray, ...]:
    """Return all sites from the centre through shell ``k`` in shell-major order."""
    if k < 0:
        raise ValueError("shell index must be nonnegative")
    if k > MAX_SHELL:
        raise ValueError(f"shell index above the supported bound {MAX_SHELL}")
    return tuple(site for shell in _shells(k) for site in shell)


def shell_index(point: Quadray | tuple[int, int, int, int]) -> int:
    """Return the graph-distance shell of a canonical IVM site.

    The shell is found by expanding shells until the point appears.  It is a
    graph-distance definition, not a fixed Euclidean radius.
    """
    target = Quadray(*point).normalize() if not isinstance(point, Quadray) else point.normalize()
    for index in range(MAX_SHELL + 1):
        if target in _shells(index)[index]:
            return index
    raise ValueError(f"point {target.as_tuple()} is outside the enumerated shells")


# ---------------------------------------------------------------------------
# Twelve-around-one cuboctahedron
# ---------------------------------------------------------------------------


def cuboctahedron_neighbors() -> tuple[Quadray, ...]:
    """The twelve nearest IVM neighbours, the vector-equilibrium shell."""
    return neighbor_sites()


def cuboctahedron_adjacency() -> tuple[tuple[int, ...], ...]:
    """Return the 12x12 adjacency matrix of the cuboctahedron graph.

    Two neighbours are adjacent exactly when their squared distance is eight.
    """
    sites = cuboctahedron_neighbors()
    adjacency = []
    for i, left in enumerate(sites):
        row = []
        for j, right in enumerate(sites):
            if i == j:
                row.append(0)
                continue
            row.append(1 if squared_distance(left, right) == NEIGHBOR_SQUARED_DISTANCE else 0)
        adjacency.append(tuple(row))
    return tuple(adjacency)


def cuboctahedron_edges() -> tuple[tuple[int, int], ...]:
    """Return the 24 edges of the cuboctahedron graph."""
    adjacency = cuboctahedron_adjacency()
    return tuple(
        (i, j)
        for i in range(len(adjacency))
        for j in range(i + 1, len(adjacency))
        if adjacency[i][j]
    )


def _cliques(adjacency: tuple[tuple[int, ...], ...], size: int) -> Iterator[tuple[int, ...]]:
    count = len(adjacency)
    for combination in itertools.combinations(range(count), size):
        if all(
            adjacency[a][b]
            for a, b in itertools.combinations(combination, 2)
        ):
            yield combination


def cuboctahedron_faces() -> dict[str, tuple[tuple[int, ...], ...]]:
    """Return the cuboctahedron's 8 triangular and 6 square faces as index tuples."""
    adjacency = cuboctahedron_adjacency()
    triangles = tuple(_cliques(adjacency, 3))
    squares = []
    for combination in itertools.combinations(range(len(adjacency)), 4):
        if all(adjacency[a][b] for a, b in itertools.combinations(combination, 2)):
            # A 4-clique would be a tetrahedral face pattern; the cuboctahedron
            # has none.  Genuine squares are 4-cycles without diagonals.
            continue
        ordered = _square_cycle(combination, adjacency)
        if ordered is not None:
            squares.append(ordered)
    return {
        "triangles": tuple(sorted(triangles)),
        "squares": tuple(sorted(squares)),
    }


def _square_cycle(
    vertices: tuple[int, ...], adjacency: tuple[tuple[int, ...], ...]
) -> tuple[int, ...] | None:
    """Return a 4-cycle ordering of four vertices, or ``None``."""
    for ordering in itertools.permutations(vertices):
        if all(
            adjacency[ordering[i]][ordering[(i + 1) % 4]]
            for i in range(4)
        ) and not adjacency[ordering[0]][ordering[2]] and not adjacency[ordering[1]][ordering[3]]:
            return ordering
    return None


def cuboctahedron_volume() -> Fraction:
    """Return the cuboctahedron's tetrahedral volume as a sum of 20 tetrahedra.

    Each face triangle with the origin is a unit tetrahedron.  The eight
    triangular faces contribute eight; each of the six square faces splits
    into two triangles and contributes two more.
    """
    origin = Quadray(0, 0, 0, 0)
    sites = cuboctahedron_neighbors()
    faces = cuboctahedron_faces()
    total = Fraction(0, 1)
    for triangle in faces["triangles"]:
        total += tetra_volume(origin, sites[triangle[0]], sites[triangle[1]], sites[triangle[2]])
    for square in faces["squares"]:
        total += tetra_volume(origin, sites[square[0]], sites[square[1]], sites[square[2]])
        total += tetra_volume(origin, sites[square[0]], sites[square[2]], sites[square[3]])
    return total


def octahedron_vertices() -> tuple[Quadray, ...]:
    """The six octahedral sites ``e_i + e_j`` in the canonical gauge."""
    axes = [Quadray(1, 0, 0, 0), Quadray(0, 1, 0, 0), Quadray(0, 0, 1, 0), Quadray(0, 0, 0, 1)]
    return tuple(
        (axes[i] + axes[j]).normalize()
        for i, j in itertools.combinations(range(4), 2)
    )


def octahedron_volume() -> Fraction:
    """Return the octahedron's tetrahedral volume as a sum of eight tetrahedra.

    The octahedron's six vertices form three opposite pairs; choosing one
    vertex from each pair gives eight orthant tetrahedra of volume ``1/2``.
    """
    origin = Quadray(0, 0, 0, 0)
    vertices = octahedron_vertices()
    pairs = _opposite_pairs(vertices)
    total = Fraction(0, 1)
    for choice in itertools.product(*pairs):
        total += tetra_volume(origin, choice[0], choice[1], choice[2])
    return total


def _opposite_pairs(
    vertices: tuple[Quadray, ...],
) -> tuple[tuple[Quadray, Quadray], ...]:
    """Pair each octahedron vertex with its antipode in the canonical gauge."""
    remaining = set(vertices)
    pairs = []
    while remaining:
        vertex = next(iter(remaining))
        remaining.remove(vertex)
        antipode = Quadray(*vertex.as_tuple()).scale(-1).normalize()
        if antipode not in remaining:
            raise AssertionError("octahedron vertices are not centrally symmetric")
        remaining.remove(antipode)
        pairs.append((vertex, antipode))
    return tuple(pairs)


def tetrahedron_vertices() -> tuple[Quadray, ...]:
    """The unit IVM tetrahedron: the origin plus three neighbour moves."""
    return (
        Quadray(0, 0, 0, 0),
        Quadray(2, 1, 0, 1),
        Quadray(2, 1, 1, 0),
        Quadray(2, 0, 1, 1),
    )


def tetrahedron_volume() -> Fraction:
    """Return the unit IVM tetrahedron's volume, exactly one."""
    p0, p1, p2, p3 = tetrahedron_vertices()
    return tetra_volume(p0, p1, p2, p3)


# ---------------------------------------------------------------------------
# Goldberg / geodesic counts
# ---------------------------------------------------------------------------


def goldberg_vertex_count(frequency: int) -> int:
    """Return the vertex count ``20 f^2`` of a frequency-``f`` Goldberg polyhedron.

    The Goldberg polyhedron is the dual of the geodesic icosahedron, so its
    vertices correspond to geodesic faces.  Frequency one gives the regular
    dodecahedron: 20 vertices and 12 pentagonal faces.
    """
    if frequency < 0:
        raise ValueError("frequency must be nonnegative")
    return 20 * frequency * frequency


def goldberg_face_count(frequency: int) -> int:
    """Return the face count ``10 f^2 + 2`` of a frequency-``f`` Goldberg polyhedron.

    The twelve pentagons are frequency-independent and ``10 (f^2 - 1)``
    hexagons complete the shell.  Frequency one gives the dodecahedron's 12
    pentagons.  This is the quadratic form shared with the IVM shell count.
    """
    if frequency < 0:
        raise ValueError("frequency must be nonnegative")
    return 10 * frequency * frequency + 2


def geodesic_vertex_count(frequency: int) -> int:
    """Return the vertex count ``10 f^2 + 2`` of a geodesic icosahedron."""
    if frequency < 0:
        raise ValueError("frequency must be nonnegative")
    return 10 * frequency * frequency + 2


def geodesic_face_count(frequency: int) -> int:
    """Return the ``20 f^2`` triangular face count of a geodesic icosahedron."""
    if frequency < 0:
        raise ValueError("frequency must be nonnegative")
    return 20 * frequency * frequency


def geodesic_edge_count(frequency: int) -> int:
    """Return the ``30 f^2`` edge count of a geodesic icosahedron."""
    if frequency < 0:
        raise ValueError("frequency must be nonnegative")
    return 30 * frequency * frequency


def cuboctahedron_signature() -> dict[str, object]:
    """A compact, checkable signature of the twelve-around-one structure."""
    adjacency = cuboctahedron_adjacency()
    faces = cuboctahedron_faces()
    return {
        "vertices": len(adjacency),
        "edges": len(cuboctahedron_edges()),
        "degree": sorted({sum(row) for row in adjacency}),
        "triangles": len(faces["triangles"]),
        "squares": len(faces["squares"]),
        "squared_edge_length": NEIGHBOR_SQUARED_DISTANCE,
        "volume": cuboctahedron_volume(),
    }
