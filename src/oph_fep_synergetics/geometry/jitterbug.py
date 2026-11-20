"""Exact joints, rigid bars and convex hulls for Fuller's Jitterbug path.

The twelve labelled joints are cyclic coordinate shifts of
``(0, ±1, ±t)``, with ``0 <= t <= 1``. Rescaling the joints by the
reciprocal bar length keeps the eight equilateral panels rigid. The
analytic hull uses the short diagonal of each former square opening
for every interior parameter. Only at the golden-ratio parameter
are these six diagonals equilateral with the twenty-four rigid bars.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from functools import cache, cmp_to_key, lru_cache
from itertools import combinations, product
from math import isfinite, prod, sqrt

from .exact import PHI_EXACT, Surd5, convex_facets

Parameter = int | Fraction | Surd5
Point = tuple[Parameter, Parameter, Parameter]
Edge = tuple[int, int]
Triangle = tuple[int, int, int]

ICOSAHEDRON_PARAMETER = PHI_EXACT - 1

_SIGNS = (-1, 1)
_JOINT_LABELS = tuple(product(range(3), _SIGNS, _SIGNS))
_JOINT_INDEX = {label: index for index, label in enumerate(_JOINT_LABELS)}


def _parameter(t: Parameter) -> Surd5:
    exact = Surd5.of(t)
    if exact.sign() < 0 or (exact - 1).sign() > 0:
        raise ValueError("Jitterbug parameter must lie between zero and one")
    return exact


def jitterbug_vertices(t: Parameter) -> tuple[Point, ...]:
    """Return all twelve exact labelled joints, retaining coincidences at zero."""
    exact = _parameter(t)
    ratio = exact.a if exact.is_rational() else exact
    vertices: list[Point] = []
    for axis, large_sign, small_sign in _JOINT_LABELS:
        point: list[Parameter] = [0, 0, 0]
        point[(axis + 1) % 3] = large_sign
        point[(axis + 2) % 3] = small_sign * ratio
        vertices.append((point[0], point[1], point[2]))
    return tuple(vertices)


def scaled_jitterbug_vertices(
    t: Parameter | float, *, edge_length: float = 2.0
) -> tuple[tuple[float, float, float], ...]:
    """Return plotting coordinates with fixed bar length along the path."""
    if not isfinite(edge_length) or edge_length <= 0:
        raise ValueError("bar length must be finite and positive")
    parameter = Fraction(t) if isinstance(t, float) else t
    exact = _parameter(parameter)
    factor = edge_length / sqrt(float(2 * (1 + exact * exact - exact)))
    return tuple(
        (factor * float(x), factor * float(y), factor * float(z))
        for x, y, z in jitterbug_vertices(exact)
    )


@cache
def jitterbug_rigid_triangles() -> tuple[Triangle, ...]:
    """Return the eight rigid triangles in outward winding order."""
    triangles = []
    for signs in product(_SIGNS, repeat=3):
        first, second, third = (
            _JOINT_INDEX[axis, signs[(axis + 1) % 3], signs[(axis + 2) % 3]]
            for axis in range(3)
        )
        # A coordinate reflection flips winding at fixed outward normal.
        triangles.append((first, second, third) if prod(signs) > 0 else (first, third, second))
    return tuple(triangles)


@cache
def jitterbug_rigid_edges() -> tuple[Edge, ...]:
    """Return the twenty-four physical bars of the triangular panels."""
    return tuple(sorted({
        edge for triangle in jitterbug_rigid_triangles() for edge in combinations(sorted(triangle), 2)
    }))


@cache
def jitterbug_matching_edges() -> tuple[Edge, ...]:
    """Return the six short-diagonal pairs, each flipping the small coordinate."""
    return tuple(
        (_JOINT_INDEX[axis, sign, -1], _JOINT_INDEX[axis, sign, 1])
        for axis, sign in product(range(3), _SIGNS)
    )


def _squared_chord(first: Point, second: Point) -> Surd5:
    delta = (left - right for left, right in zip(first, second, strict=True))
    return Surd5.of(sum(value * value for value in delta))


def _rotate_face(face: tuple[int, ...]) -> tuple[int, ...]:
    start = face.index(min(face))
    return face[start:] + face[:start]


@lru_cache(maxsize=16)
def _jitterbug_snapshot(t: Surd5) -> dict[str, object]:
    joints = jitterbug_vertices(t)
    vertices = tuple(dict.fromkeys(joints))
    joint_to_vertex = tuple(vertices.index(joint) for joint in joints)
    facets = convex_facets(vertices)
    hull_edges = tuple(sorted({
        tuple(sorted((a, b))) for face in facets for a, b in zip(face, face[1:] + face[:1], strict=True)
    }))
    bar_counts = Counter(tuple(sorted((joint_to_vertex[a], joint_to_vertex[b])))
                         for a, b in jitterbug_rigid_edges())
    bars = tuple(sorted(bar_counts))
    matching = tuple(sorted(
        tuple(sorted((joint_to_vertex[a], joint_to_vertex[b])))
        for a, b in jitterbug_matching_edges() if joint_to_vertex[a] != joint_to_vertex[b]
    ))
    expected_hull = set(bars) if t in (0, 1) else set(bars) | set(matching)
    expected_rigid_faces = {
        _rotate_face(tuple(joint_to_vertex[joint] for joint in face)) for face in jitterbug_rigid_triangles()
    }
    lengths = {
        edge: _squared_chord(vertices[edge[0]], vertices[edge[1]])
        for edge in combinations(range(len(vertices)), 2)
    }
    length_order = cmp_to_key(lambda left, right: (left - right).sign())
    minimum = min(lengths.values(), key=length_order)
    radius_squared = 1 + t * t
    bar_squared = 2 * (radius_squared - t)
    matching_squared = 4 * t * t
    shortest = tuple(edge for edge, length in lengths.items() if length == minimum)
    contacts = tuple(edge for edge, length in lengths.items() if length == radius_squared)
    overlaps = tuple(edge for edge, length in lengths.items() if (length - radius_squared).sign() < 0)
    return {
        "parameter": t,
        "vertices": vertices,
        "joint_to_vertex": joint_to_vertex,
        "vertex_count": len(vertices),
        "coincident_joint_pairs": len(joints) - len(vertices),
        "facets": facets,
        "facet_counts": tuple(sorted(Counter(map(len, facets)).items())),
        "hull_edges": hull_edges,
        "hull_edge_count": len(hull_edges),
        "bar_edges": bars,
        "bar_edge_count": len(bars),
        "bar_multiplicities": tuple(bar_counts[edge] for edge in bars),
        "matching_edges": matching,
        "matching_edge_count": len(matching),
        "shortest_pairs": shortest,
        "shortest_pair_count": len(shortest),
        "contact_pairs": contacts,
        "contact_pair_count": len(contacts),
        "overlap_pairs": overlaps,
        "overlap_pair_count": len(overlaps),
        "bar_squared_length": bar_squared,
        "matching_squared_length": matching_squared,
        "minimum_squared_length": minimum,
        "radius_squared": radius_squared,
        "pair_squared_length_counts": tuple(sorted(
            Counter(lengths.values()).items(), key=lambda row: length_order(row[0])
        )),
        "bar_lengths_match": all(lengths[edge] == bar_squared for edge in bars),
        "matching_lengths_match": all(
            _squared_chord(joints[a], joints[b]) == matching_squared for a, b in jitterbug_matching_edges()
        ),
        "rigid_orientations_match": expected_rigid_faces <= set(facets),
        "hull_matches_designated_graph": set(hull_edges) == expected_hull,
    }


def jitterbug_snapshot(t: Parameter) -> dict[str, object]:
    """Census an exact supporting hull and its chord, bar and contact graphs.

    Edges in the result use distinct-vertex indices; ``joint_to_vertex``
    relates them to the twelve physical labels. The contact and overlap
    census is for distinct centres at radius R from the origin, with
    congruent spheres of radius R/2. Coincident labels are recorded
    separately at the octahedral limit.
    """
    return _jitterbug_snapshot(_parameter(t)).copy()


def jitterbug_census() -> dict[str, dict[str, object]]:
    """Check cuboctahedral, icosahedral and octahedral phases and controls."""
    return {
        "cuboctahedron": jitterbug_snapshot(1),
        "intermediate": jitterbug_snapshot(Fraction(3, 4)),
        "icosahedron": jitterbug_snapshot(ICOSAHEDRON_PARAMETER),
        "contracted": jitterbug_snapshot(Fraction(3, 5)),
        "overlap_control": jitterbug_snapshot(Fraction(1, 2)),
        "octahedron": jitterbug_snapshot(0),
    }
