"""Exact centroid moments and stationary memory in the binary repair chain."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations

import numpy as np

from ..geometry.exact import Surd5, exact_icosahedron_vertices
from ..geometry.lattice import cuboctahedron_adjacency, cuboctahedron_neighbors
from ..geometry.quadray import EMBEDDING
from .carrier import icosahedral_carrier
from .native_dynamics import _attempt_counts
from .network import Adjacency, graphs_are_isomorphic, laplacian_matrix

Coordinate = tuple[Surd5 | Fraction | int, ...]


def _centroid_case(
    vertices: tuple[Coordinate, ...],
    edge_squared: int,
    target: tuple[int, int, int],
    reference_adjacency: Adjacency,
) -> dict[str, object]:
    def squared_chord(left: int, right: int) -> Surd5 | Fraction | int:
        delta = (a - b for a, b in zip(vertices[left], vertices[right], strict=True))
        return sum(value * value for value in delta)

    ports = len(vertices)
    pairs = tuple(combinations(range(ports), 2))
    adjacency = tuple(
        tuple(int(squared_chord(left, right) == edge_squared) for right in range(ports))
        for left in range(ports)
    )
    seams = tuple((left, right) for left, right in pairs if adjacency[left][right])
    states = tuple(tuple(int(port in pair) for port in range(ports)) for pair in pairs)
    counts = _attempt_counts(states, seams)
    draws = 2 * len(seams)

    coordinates = np.asarray(vertices, dtype=object)
    gram = coordinates.T @ coordinates
    laplacian = np.asarray(laplacian_matrix(adjacency), dtype=object)
    laplacian_coordinates = laplacian @ coordinates
    eigenvalue = Surd5.of((coordinates.T @ laplacian_coordinates)[0, 0]) / gram[0, 0]
    retention = 1 - eigenvalue / draws

    # Under a uniform two-subset, the occupancy covariance is the projector
    # onto the centred subspace times 2(n-2)/(n(n-1)). Compare the induced
    # coordinate covariance with the direct second moment of the full shell.
    centroids = np.asarray(states, dtype=object) @ coordinates
    covariance = Fraction(2 * (ports - 2), ports * (ports - 1)) * gram
    axis_variance = covariance[0, 0]
    isotropic = np.eye(len(target), dtype=object) * axis_variance

    zero_fibre = [i for i, point in enumerate(centroids) if not any(value != 0 for value in point)]
    target_fibre = [i for i, point in enumerate(centroids) if tuple(point) == target]
    arrive = counts[np.ix_(target_fibre, zero_fibre)].sum(axis=0)
    leave = counts[np.ix_(zero_fibre, target_fibre)].sum(axis=1)
    next_given_zero = Fraction(int(leave.sum()), draws * len(zero_fibre))
    next_given_previous_target = Fraction(int(arrive @ leave), draws * int(arrive.sum()))

    # At a zero centroid, swapping one token shifts the centroid by the
    # corresponding oriented edge vector. Count that event directly from
    # the exact vertex coordinates and incidences, without the repair matrix.
    geometric_jump_counts = []
    for index in zero_fibre:
        left, right = pairs[index]
        geometric_jump_counts.append(sum(
            tuple(x - y for x, y in zip(vertices[destination], vertices[mover], strict=True)) == target
            for mover, fixed in ((left, right), (right, left))
            for destination, edge in enumerate(adjacency[mover])
            if edge and destination != fixed
        ))

    return {
        "attempts_per_state": draws,
        "graph_matches_reference": graphs_are_isomorphic(adjacency, reference_adjacency),
        "uniform_stationary": bool(
            np.array_equal(counts, counts.T) and np.all(counts.sum(axis=1) == draws)
        ),
        "coordinate_eigenvalue": str(eigenvalue),
        "coordinate_mode_verified": bool(np.array_equal(laplacian_coordinates, coordinates * eigenvalue)),
        "retention": str(retention),
        "retention_parts": tuple(map(str, retention.parts())),
        "conditional_mean_verified": bool(np.array_equal(
            counts @ centroids, centroids * (retention * draws)
        )),
        "one_axis_variance": str(axis_variance),
        "stationary_variance_verified": bool(
            np.all(centroids.sum(axis=0) == 0)
            and np.array_equal(covariance, isotropic)
            and np.array_equal(centroids.T @ centroids, covariance * len(pairs))
        ),
        "target_centroid": target,
        "zero_fibre_size": len(zero_fibre),
        "zero_fibre_jump_counts": tuple(sorted(leave.tolist())),
        "witness_matches_edge_moves": bool(np.array_equal(leave, geometric_jump_counts)),
        "witness_matches_reversible_formula": next_given_previous_target == Fraction(
            int(leave @ leave), draws * int(leave.sum())
        ),
        "next_target_given_zero": str(next_given_zero),
        "next_target_given_previous_target_and_zero": str(next_given_previous_target),
        "stationary_centroid_non_markov": next_given_zero != next_given_previous_target,
    }


def centroid_gaussian_census() -> dict[str, dict[str, object]]:
    """Compare all-lag OU covariance and a three-time repair witness.

    The observable is the sum of the occupied vertex coordinates. The
    icosahedron has edge length two. The cuboctahedron uses the integer D3
    chart of the FCC first shell, with edge length sqrt(2). Conditional mean
    closure and isotropic stationarity give lag covariance C * alpha**k,
    which also holds for a stationary OU process sampled at integer times
    with retention alpha and invariant covariance C.
    """
    carrier = icosahedral_carrier()
    fcc_vertices = tuple(
        tuple(Fraction(sum(weight * component for weight, component in zip(
            axis, site.as_tuple(), strict=True
        )), 2) for axis in EMBEDDING)
        for site in cuboctahedron_neighbors()
    )
    return {
        "carrier": _centroid_case(
            exact_icosahedron_vertices(), 4, (0, 2, 0),
            carrier.adjacency_matrix().astype(int).tolist(),
        ),
        "cuboctahedron": _centroid_case(
            fcc_vertices, 2, (0, 1, -1), cuboctahedron_adjacency(),
        ),
    }
