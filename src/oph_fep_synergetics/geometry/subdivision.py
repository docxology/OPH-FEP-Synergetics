"""Icosahedral subdivision and its Goldberg dual counts.

The regular icosahedron has 12 vertices, 30 edges and 20 triangular faces.  A
geodesic subdivision of frequency ``f`` refines each face into ``f^2``
triangles, so the triangulation carries ``20 f^2`` faces, ``30 f^2`` edges and
``10 f^2 + 2`` vertices.  Its Goldberg dual swaps vertices and faces, carrying
``20 f^2`` vertices and ``10 f^2 + 2`` faces.

Two counts coincide with the IVM frequency-shell count ``10 f^2 + 2``: the
geodesic vertex count and the Goldberg face count.  Frequency one gives twelve
in each case -- the icosahedron's 12 vertices and the dodecahedron's 12
pentagonal faces -- which is the same dozen as the twelve IVM neighbours around
one site.  The count coincidence is real; the frequency-one polyhedra are not
the same solid as the cuboctahedral shell.

Provenance: the icosahedron counts and the twelve-around-one correspondence
follow ``docxology/QuadMath`` (Apache-2.0) and OPH's carrier at
``Lean/Screen/PortFrameGram.lean``; the dual counts follow the standard
geodesic/Goldberg duality and are recorded as declared, executable identities.
"""

from __future__ import annotations

import itertools
import math
from typing import Iterable

from .vectors import EDGE_TOLERANCE, TIGHT_TOLERANCE

PHI = (1.0 + math.sqrt(5.0)) / 2.0


def icosahedron_vertices() -> tuple[tuple[float, float, float], ...]:
    """Return the 12 vertices of the regular icosahedron.

    Coordinates are the cyclic permutations of ``(0, +-1, +-PHI)``, on the
    sphere of radius ``sqrt(1 + PHI^2)``.
    """
    vertices: list[tuple[float, float, float]] = []
    for first in (-1.0, 1.0):
        for second in (-1.0, 1.0):
            vertices.append((0.0, first * 1.0, second * PHI))
            vertices.append((first * 1.0, second * PHI, 0.0))
            vertices.append((first * PHI, 0.0, second * 1.0))
    return tuple(vertices)


def _squared_distance(a: Iterable[float], b: Iterable[float]) -> float:
    return sum((x - y) ** 2 for x, y in zip(a, b))


def icosahedron_edges() -> tuple[tuple[int, int], ...]:
    """Return the 30 edges of the regular icosahedron."""
    vertices = icosahedron_vertices()
    pairwise = [
        (a, b, _squared_distance(vertices[a], vertices[b]))
        for a, b in itertools.combinations(range(len(vertices)), 2)
    ]
    edge_length = min(distance for _, _, distance in pairwise if distance > TIGHT_TOLERANCE)
    return tuple(
        sorted(
            (a, b)
            for a, b, distance in pairwise
            if abs(distance - edge_length) <= EDGE_TOLERANCE
        )
    )


def icosahedron_faces() -> tuple[tuple[int, int, int], ...]:
    """Return the 20 triangular faces of the regular icosahedron."""
    edges = set(icosahedron_edges())
    adjacency = {vertex: set() for vertex in range(len(icosahedron_vertices()))}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    faces = []
    for a, b, c in itertools.combinations(range(len(icosahedron_vertices())), 3):
        if (a, b) in edges and (b, c) in edges and (a, c) in edges:
            faces.append((a, b, c))
    return tuple(sorted(faces))


def icosahedron_adjacency() -> tuple[tuple[int, ...], ...]:
    """Return the 12x12 adjacency matrix of the icosahedron graph."""
    count = len(icosahedron_vertices())
    adjacency = [[0] * count for _ in range(count)]
    for a, b in icosahedron_edges():
        adjacency[a][b] = 1
        adjacency[b][a] = 1
    return tuple(tuple(row) for row in adjacency)


def geodesic_subdivision_counts(frequency: int) -> dict[str, int]:
    """Return the vertex/edge/face counts of a frequency-``f`` geodesic icosahedron.

    The subdivision keeps the 12 original vertices, adds ``10 (f^2 - 1)`` edge
    midpoints, and splits each of the 20 faces into ``f^2`` triangles.
    """
    if frequency < 1:
        raise ValueError("frequency must be a positive integer")
    faces = 20 * frequency * frequency
    edges = 30 * frequency * frequency
    vertices = 10 * frequency * frequency + 2
    return {
        "vertices": vertices,
        "edges": edges,
        "faces": faces,
        "frequency": frequency,
        "original_vertices": 12,
        "added_vertices": vertices - 12,
    }


def goldberg_dual_counts(frequency: int) -> dict[str, int]:
    """Return the counts of the Goldberg dual of a geodesic icosahedron.

    The dual has one vertex per geodesic face and one face per geodesic vertex,
    so its vertex count is ``20 f^2`` and its face count is ``10 f^2 + 2``.
    Frequency one is the regular dodecahedron: 20 vertices, 12 faces.
    """
    geodesic = geodesic_subdivision_counts(frequency)
    return {
        "vertices": geodesic["faces"],
        "faces": geodesic["vertices"],
        "edges": geodesic["edges"],
        "frequency": frequency,
    }


def ivm_shell_matches_geodesic(frequency: int) -> bool:
    """Whether the IVM shell count equals the geodesic *vertex* count."""
    from .lattice import shell_count

    return shell_count(frequency) == geodesic_subdivision_counts(frequency)["vertices"]


def ivm_shell_matches_goldberg(frequency: int) -> bool:
    """Whether the IVM shell count equals the Goldberg dual *face* count."""
    from .lattice import shell_count

    return shell_count(frequency) == goldberg_dual_counts(frequency)["faces"]
