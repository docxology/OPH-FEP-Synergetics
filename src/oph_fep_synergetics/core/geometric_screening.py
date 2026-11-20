"""Gaussian screening by central planes in the two twelve-site geometries."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations

from ..geometry.exact import (
    Point,
    _cross,
    _dot,
    _sign,
    _sub,
    exact_icosahedron_vertices,
    ivm_solid_vertices,
)
from ..geometry.lattice import cuboctahedron_adjacency
from .carrier import icosahedral_carrier
from .invariants import _exact_inverse, laplacian_determinant
from .network import Adjacency, graphs_are_isomorphic, laplacian_matrix


def _component_sizes(neighbours: tuple[set[int], ...], records: tuple[int, ...]) -> tuple[int, ...]:
    unvisited = set(range(len(neighbours))) - set(records)
    sizes = []
    while unvisited:
        frontier = [unvisited.pop()]
        size = 0
        while frontier:
            port = frontier.pop()
            reached = neighbours[port] & unvisited
            unvisited.difference_update(reached)
            frontier.extend(reached)
            size += 1
        sizes.append(size)
    return tuple(sorted(sizes))


def _minimum_separators(adjacency: Adjacency) -> tuple[int, tuple[tuple[tuple[int, ...], int], ...]]:
    """Enumerate record sets through the first budget that disconnects the graph."""
    ports = range(len(adjacency))
    neighbours = tuple({j for j, edge in enumerate(row) if edge} for row in adjacency)
    for budget in range(len(ports) - 1):
        sizes: Counter[tuple[int, ...]] = Counter()
        for records in combinations(ports, budget):
            components = _component_sizes(neighbours, records)
            if len(components) > 1:
                sizes[components] += 1
        if sizes:
            return budget, tuple(sorted(sizes.items()))
    raise AssertionError("the two finite geometry graphs have vertex separators")


def _hemispheres(signs: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    positive = tuple(port for port, sign in enumerate(signs) if sign > 0)
    negative = tuple(port for port, sign in enumerate(signs) if sign < 0)
    return positive, negative


def _crossing_seams(
    adjacency: Adjacency, positive: tuple[int, ...], negative: tuple[int, ...],
) -> tuple[tuple[int, int], ...]:
    return tuple((port, other) for port in positive for other in negative if adjacency[port][other])


def _plane_census(vertices: tuple[Point, ...], adjacency: Adjacency) -> dict[str, object]:
    # A plane containing independent vertex directions is uniquely determined
    # by the set of vertices in it. Cross products enumerate all such planes.
    signs_by_records = {}
    for first, second in combinations(vertices, 2):
        normal = _cross(first, second)
        if any(value != 0 for value in normal):
            signs = tuple(_sign(_dot(normal, point)) for point in vertices)
            records = tuple(port for port, sign in enumerate(signs) if sign == 0)
            signs_by_records[records] = signs

    laplacian = laplacian_matrix(adjacency)
    classes: Counter[tuple[int, int, int, int, Fraction]] = Counter()
    inverse_matches_seams = True
    for records, signs in signs_by_records.items():
        positive, negative = _hemispheres(signs)
        seams = _crossing_seams(adjacency, positive, negative)
        free = positive + negative
        covariance = _exact_inverse([[Fraction(laplacian[i][j]) for j in free] for i in free])
        cross_zero = all(
            covariance[i][j] == 0
            for i in range(len(positive)) for j in range(len(positive), len(free))
        )
        inverse_matches_seams &= cross_zero == (not seams)
        # Jacobi's identity writes the Gaussian conditional information as
        # 1/2 log(det Q_++ det Q_-- / det Q_free).
        ratio = Fraction(
            laplacian_determinant(laplacian, grounded=records + negative)
            * laplacian_determinant(laplacian, grounded=records + positive),
            laplacian_determinant(laplacian, grounded=records),
        )
        classes[len(records), len(positive), len(negative), len(seams), ratio] += 1

    # The remaining central planes can contain at most one antipodal
    # direction. Each open hemisphere of a full-dimensional centrally
    # symmetric vertex set is nonempty, so a separator of size at most two
    # would be necessary for any of those planes to screen.
    min_records, min_record_component_sizes = _minimum_separators(adjacency)
    axis_signs = tuple(_sign(point[0]) for point in vertices)
    axis_records = tuple(port for port, sign in enumerate(axis_signs) if sign == 0)
    axis_positive, axis_negative = _hemispheres(axis_signs)
    return {
        "plane_classes": tuple(
            (size, positive, negative, seams, str(ratio), count)
            for (size, positive, negative, seams, ratio), count in sorted(classes.items())
        ),
        "vertex_spanned_planes": len(signs_by_records),
        "screening_planes": sum(count for (_, _, _, seams, _), count in classes.items() if seams == 0),
        "minimum_separator_size": min_records,
        "minimum_separator_component_counts": min_record_component_sizes,
        "inverse_matches_seams": inverse_matches_seams,
        "axis_plane_partition": (axis_records, axis_positive, axis_negative),
        "axis_plane_crossing_seams": _crossing_seams(adjacency, axis_positive, axis_negative),
        "axis_plane_precision_determinants": (
            laplacian_determinant(laplacian, grounded=axis_records + axis_negative),
            laplacian_determinant(laplacian, grounded=axis_records + axis_positive),
            laplacian_determinant(laplacian, grounded=axis_records),
        ),
    }


@cache
def _central_plane_summary() -> tuple[bool, dict[str, object], dict[str, object]]:
    ico_vertices = exact_icosahedron_vertices()
    ico_adjacency = []
    for first in ico_vertices:
        row = []
        for second in ico_vertices:
            delta = _sub(first, second)
            row.append(int(_dot(delta, delta) == 4))
        ico_adjacency.append(row)
    oph = icosahedral_carrier().adjacency_matrix().astype(int).tolist()
    return (
        graphs_are_isomorphic(ico_adjacency, oph),
        _plane_census(ico_vertices, ico_adjacency),
        _plane_census(ivm_solid_vertices()["cuboctahedron"], cuboctahedron_adjacency()),
    )


def central_plane_blanket_census() -> dict[str, object]:
    """Classify central planes by Gaussian screening for both geometries.

    Each plane-class row gives records, positive and negative side sizes,
    crossing seams, the exact conditional-information determinant ratio,
    and the number of unoriented planes. Ratios use the unit intrinsic
    Gaussian grounded at the plane's recorded sites.
    """
    oph_isomorphism, icosahedron, cuboctahedron = _central_plane_summary()
    return {
        "oph_isomorphism": oph_isomorphism,
        "icosahedron": icosahedron.copy(),
        "cuboctahedron": cuboctahedron.copy(),
    }
