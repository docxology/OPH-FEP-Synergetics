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
SPHERE_CONTACT_PARAMETER_SQUARED = Fraction(1, 3)

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


def _component_sizes(vertex_count: int, edges: tuple[Edge, ...]) -> tuple[int, ...]:
    adjacent: list[set[int]] = [set() for _ in range(vertex_count)]
    for left, right in edges:
        adjacent[left].add(right)
        adjacent[right].add(left)
    unreached = set(range(vertex_count))
    sizes = []
    while unreached:
        frontier = [unreached.pop()]
        size = 0
        while frontier:
            size += 1
            new = adjacent[frontier.pop()] & unreached
            unreached.difference_update(new)
            frontier.extend(new)
        sizes.append(size)
    return tuple(sorted(sizes))


def _labelled_chord_coefficients() -> dict[Edge, tuple[int, int, int]]:
    """Derive each chord's coefficients in the basis 1, t, t squared."""
    initial, final = jitterbug_vertices(0), jitterbug_vertices(1)
    coefficients = {}
    for left, right in combinations(range(len(initial)), 2):
        base = tuple(int(x - y) for x, y in zip(initial[left], initial[right], strict=True))
        at_one = tuple(int(x - y) for x, y in zip(final[left], final[right], strict=True))
        slope = tuple(x - y for x, y in zip(at_one, base, strict=True))
        coefficients[left, right] = (
            sum(x * x for x in base),
            2 * sum(x * y for x, y in zip(base, slope, strict=True)),
            sum(x * x for x in slope),
        )
    return coefficients


def _sphere_threshold_gap_sign(coefficients: tuple[int, int, int]) -> int:
    """Sign of chord squared minus radius squared at t = sqrt(1/3)."""
    constant, linear, quadratic = coefficients
    # After substitution the gap is rational + linear * sqrt(1/3).
    rational = constant - 1 + (quadratic - 1) * SPHERE_CONTACT_PARAMETER_SQUARED
    rational_sign = (rational > 0) - (rational < 0)
    linear_sign = (linear > 0) - (linear < 0)
    if not rational_sign or not linear_sign or rational_sign == linear_sign:
        return rational_sign or linear_sign
    difference = rational * rational - linear * linear * SPHERE_CONTACT_PARAMETER_SQUARED
    return rational_sign * ((difference > 0) - (difference < 0))


def jitterbug_contraction_census() -> dict[str, object]:
    """Compare nearest-pair connectivity, sphere contacts and joint coalescence."""
    coefficients = _labelled_chord_coefficients()
    threshold_gaps = {edge: _sphere_threshold_gap_sign(row) for edge, row in coefficients.items()}
    checkpoints = {}
    landmarks = jitterbug_census()
    for name, snapshot in landmarks.items():
        parameter = snapshot["parameter"]
        shortest = snapshot["shortest_pairs"]
        vertex_count = snapshot["vertex_count"]
        joint_to_vertex = snapshot["joint_to_vertex"]
        measured_classes = {
            tuple(sorted((joint_to_vertex[a], joint_to_vertex[b]))):
                constant + linear * parameter + quadratic * parameter * parameter
            for (a, b), (constant, linear, quadratic) in coefficients.items()
            if joint_to_vertex[a] != joint_to_vertex[b]
        }
        checkpoints[name] = {
            "parameter": str(parameter),
            "shortest_edge_count": len(shortest),
            "shortest_component_sizes": _component_sizes(vertex_count, shortest),
            "outer_contacts": snapshot["contact_pair_count"],
            "outer_overlaps": snapshot["overlap_pair_count"],
            "hull_edge_count": snapshot["hull_edge_count"],
            "class_formula_matches": (
                dict(Counter(measured_classes.values())) == dict(snapshot["pair_squared_length_counts"])
            ),
        }

    # Quotient the labelled nearest and hull graphs from the open branch.
    # The matching becomes loops. The rigid bars instead give double edges
    # between octahedral sites and remain after suppressing loops.
    octahedron = landmarks["octahedron"]
    contracted = landmarks["contracted"]
    coalesced = octahedron["joint_to_vertex"]

    def project_edges(edges: tuple[Edge, ...]) -> tuple[Edge, ...]:
        return tuple(sorted({
            tuple(sorted((coalesced[a], coalesced[b])))
            for a, b in edges if coalesced[a] != coalesced[b]
        }))

    limit_of_nearest = project_edges(contracted["shortest_pairs"])
    vertices_at_limit = octahedron["vertex_count"]
    nearest_after_collapse = octahedron["shortest_pairs"]
    return {
        "checkpoints": checkpoints,
        "squared_chord_class_counts": tuple(sorted(Counter(coefficients.values()).items())),
        "bar_polynomial_matches": tuple(
            edge for edge, row in coefficients.items() if row == (2, -2, 2)
        ) == jitterbug_rigid_edges(),
        "matching_polynomial_matches": tuple(
            edge for edge, row in coefficients.items() if row == (0, 0, 4)
        ) == jitterbug_matching_edges(),
        "sphere_contact_threshold_squared": str(SPHERE_CONTACT_PARAMETER_SQUARED),
        "threshold_contact_edges": tuple(edge for edge, sign in threshold_gaps.items() if sign == 0),
        "threshold_overlap_edges": tuple(edge for edge, sign in threshold_gaps.items() if sign < 0),
        "coalescence": {
            "nearest_quotient_edge_count": len(limit_of_nearest),
            "nearest_quotient_component_sizes": _component_sizes(vertices_at_limit, limit_of_nearest),
            "octahedral_nearest_edge_count": len(nearest_after_collapse),
            "octahedral_nearest_component_sizes": _component_sizes(vertices_at_limit, nearest_after_collapse),
            "bar_quotient_matches_octahedron": project_edges(contracted["bar_edges"]) == nearest_after_collapse,
            "hull_quotient_matches_octahedron": project_edges(contracted["hull_edges"]) == nearest_after_collapse,
        },
    }
