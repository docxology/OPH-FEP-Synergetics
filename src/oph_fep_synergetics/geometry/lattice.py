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
import math
from fractions import Fraction
from typing import Iterator

from .quadray import (
    EMBEDDING,
    Quadray,
    neighbor_moves,
    neighbor_sites,
    squared_distance,
    squared_norm,
    tetra_volume,
    to_xyz,
)
from .solids import (
    _add,
    _cross,
    _dot,
    _norm,
    _scale,
    _subtract,
    orient_faces_outward,
    polyhedron_volume,
    regular_tetrahedron_volume,
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


# ---------------------------------------------------------------------------
# The Voronoi cell of the IVM, and the packing density it implies
# ---------------------------------------------------------------------------


def voronoi_cell() -> dict[str, object]:
    """Return the Voronoi cell of the IVM, built as the dual of the twelve neighbours.

    The cell is the intersection of the half-spaces ``x . n_i <= |n_i|^2 / 2``
    over the twelve nearest neighbours.  Its vertices are the duals of the
    cuboctahedron's faces -- the neighbour sum over a face, divided by four --
    which places eight vertices on the triangular faces and six on the square
    ones.  Nothing about the cell is transcribed: its vertex, edge and face
    counts, its face shapes and Euler's characteristic are read off the
    construction, and every vertex is checked against all twelve half-spaces.

    Each face is ordered counter-clockwise as seen from outside, and
    :func:`~oph_fep_synergetics.geometry.solids.polyhedron_volume` raises if any
    of them points inward, so the measured volume cannot be a cancellation
    artefact.  A pyramid decomposition over the same faces is returned beside the
    divergence-theorem volume as an independent second measurement.
    """
    sites = cuboctahedron_neighbors()
    neighbours = [to_xyz(site) for site in sites]
    half = 0.5 * float(squared_norm(sites[0]))
    faces = cuboctahedron_faces()

    vertices: list[tuple[float, float, float]] = []
    dualized: list[tuple[int, ...]] = []
    for key in ("triangles", "squares"):
        for face in faces[key]:
            total = (0.0, 0.0, 0.0)
            for index in face:
                total = _add(total, neighbours[index])
            vertices.append(_scale(total, 0.25))
            dualized.append(tuple(face))

    cell_faces: list[tuple[int, ...]] = []
    for index, neighbour in enumerate(neighbours):
        incident = [dual for dual, face in enumerate(dualized) if index in face]
        axis = _scale(neighbour, 1.0 / _norm(neighbour))
        reference = _cross(axis, (1.0, 0.0, 0.0))
        if _norm(reference) < 1e-9:
            reference = _cross(axis, (0.0, 1.0, 0.0))
        reference = _scale(reference, 1.0 / _norm(reference))
        second = _cross(axis, reference)
        centre = _scale(
            sum((vertices[dual] for dual in incident), (0.0, 0.0, 0.0)),
            1.0 / len(incident),
        )
        cell_faces.append(
            tuple(
                sorted(
                    incident,
                    key=lambda dual: math.atan2(
                        _dot(_subtract(vertices[dual], centre), second),
                        _dot(_subtract(vertices[dual], centre), reference),
                    ),
                )
            )
        )

    edges: set[tuple[int, int]] = set()
    for face in cell_faces:
        for a, b in zip(face, face[1:] + face[:1]):
            edges.add((a, b) if a < b else (b, a))

    side_spreads = []
    for face in cell_faces:
        lengths = [
            _norm(_subtract(vertices[face[k]], vertices[face[(k + 1) % len(face)]]))
            for k in range(len(face))
        ]
        side_spreads.append(max(lengths) - min(lengths))

    triangles: list[tuple[int, int, int]] = []
    for face in cell_faces:
        for k in range(1, len(face) - 1):
            triangles.append((face[0], face[k], face[k + 1]))
    volume = polyhedron_volume(vertices, triangles)
    _, orientation_flips = orient_faces_outward(vertices, triangles)

    inradius = half / _norm(neighbours[0])
    pyramid = 0.0
    for face in cell_faces:
        area = 0.0
        for k in range(1, len(face) - 1):
            area += 0.5 * _norm(
                _cross(
                    _subtract(vertices[face[k]], vertices[face[0]]),
                    _subtract(vertices[face[k + 1]], vertices[face[0]]),
                )
            )
        pyramid += area * inradius / 3.0

    constraint_values = [
        _dot(vertex, neighbour) for vertex in vertices for neighbour in neighbours
    ]
    return {
        "vertices": tuple(vertices),
        "faces": tuple(cell_faces),
        "vertex_count": len(vertices),
        "face_count": len(cell_faces),
        "edge_count": len(edges),
        "euler_characteristic": len(vertices) - len(edges) + len(cell_faces),
        "all_faces_quadrilateral": all(len(face) == 4 for face in cell_faces),
        "max_side_spread": max(side_spreads),
        "all_faces_rhombic": max(side_spreads) < 1e-9,
        "constraint_half": half,
        "max_constraint_value": max(constraint_values),
        "feasible": max(constraint_values) <= half + 1e-9,
        "inradius": inradius,
        "orientation_flips": orientation_flips,
        "volume": volume,
        "pyramid_volume": pyramid,
        "volume_residual": abs(volume - pyramid),
    }


def voronoi_cell_volume() -> float:
    """Return the Voronoi cell's volume in the embedding's cubic units."""
    return float(voronoi_cell()["volume"])


def voronoi_cell_tetravolume() -> float:
    """Return the cell's volume in tetrahedra of the prime-vector edge.

    The prime vector is the IVM's nearest-neighbour distance, and it is also the
    edge of the tetrahedron spanned by the four quadray axes, so this ratio is
    the cell's tetravolume in Synergetics' own accounting.  It comes out exactly
    six, which is the table's rhombic dodecahedron entry -- derived here from the
    cell's own shape instead of transcribed from the table.
    """
    return float(voronoi_cell()["volume"]) / regular_tetrahedron_volume(
        math.sqrt(NEIGHBOR_SQUARED_DISTANCE)
    )


def packing_density() -> dict[str, object]:
    """Return the IVM's packing density, measured from the cell and the ball.

    The cell's inradius is half the prime vector, so the largest ball it can
    contain has radius ``d / 2``, and the density is that ball's volume divided
    by the cell's.  Both are computed here in the same units and their ratio is
    compared with the closed form ``pi / (3 sqrt 2)``.
    """
    cell = voronoi_cell()
    edge = math.sqrt(NEIGHBOR_SQUARED_DISTANCE)
    tetraunit = regular_tetrahedron_volume(edge)
    inradius = float(cell["inradius"])
    ball = 4.0 / 3.0 * math.pi * inradius**3
    density = ball / float(cell["volume"])
    closed = math.pi / (3.0 * math.sqrt(2.0))
    cell_tetravolume = float(cell["volume"]) / tetraunit
    return {
        "neighbours": len(cuboctahedron_neighbors()),
        "prime_vector": edge,
        "neighbour_squared_distance": NEIGHBOR_SQUARED_DISTANCE,
        "inradius": inradius,
        "inradius_is_half_the_prime_vector": abs(inradius - edge / 2.0) < 1e-12,
        "ball_volume": ball,
        "cell_volume": float(cell["volume"]),
        "cell_tetravolume": cell_tetravolume,
        "ball_tetravolume": ball / tetraunit,
        "density": density,
        "closed_form": closed,
        "density_residual": abs(density - closed),
        "cell_tetravolume_matches_table": abs(
            cell_tetravolume - float(VOLUME_RATIOS["rhombic_dodecahedron"])
        )
        < 1e-9,
        "vertices": cell["vertex_count"],
        "faces": cell["face_count"],
        "edges": cell["edge_count"],
        "euler_characteristic": cell["euler_characteristic"],
        "all_faces_rhombic": cell["all_faces_rhombic"],
        "volume_residual": cell["volume_residual"],
    }
