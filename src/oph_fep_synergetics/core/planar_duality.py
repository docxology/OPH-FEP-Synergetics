"""Exact edge covariances on paired spherical primal and dual graphs.

For unit seam precision, the covariance of oriented Gaussian edge differences
is the orthogonal projection onto the graph's cut space. Oriented face
boundaries span the dual cut space in the same edge coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import cache

import numpy as np

from ..geometry.exact import convex_facets, ivm_solid_vertices
from ..geometry.lattice import cuboctahedron_edges
from .carrier import OPH_ORIENTED_BASE_FACES, icosahedral_carrier
from .invariants import spanning_tree_number_of_adjacency
from .network import _grounded_inverse

Edge = tuple[int, int]
Matrix = tuple[tuple[Fraction, ...], ...]


@dataclass(frozen=True)
class PlanarDuality:
    """One common edge order, with pairs written from tail to head."""

    primal_edges: tuple[Edge, ...]
    dual_edges: tuple[Edge, ...]
    primal_covariance: Matrix
    dual_covariance: Matrix
    cut_ranks: tuple[int, int]
    tree_counts: tuple[int, int]


def _graph_covariance(vertex_count: int, edges: tuple[Edge, ...]) -> tuple[Matrix, int]:
    """Compute the edge covariance and tree cofactor from this graph alone."""
    adjacency = [[0] * vertex_count for _ in range(vertex_count)]
    for tail, head in edges:
        adjacency[tail][head] = adjacency[head][tail] = 1
    inverse = _grounded_inverse(adjacency)
    zero = Fraction(0)
    green = [[zero] * vertex_count] + [[zero, *row] for row in inverse]
    covariance = tuple(
        tuple(
            green[a][c] - green[a][d] - green[b][c] + green[b][d]
            for c, d in edges
        )
        for a, b in edges
    )
    return covariance, spanning_tree_number_of_adjacency(adjacency)


def _paired_covariances(
    primal_vertices: int, primal_edges: tuple[Edge, ...],
    oriented_faces: tuple[tuple[int, ...], ...],
) -> PlanarDuality:
    """Use the oriented boundary of each face to align the dual edge signs."""
    edge_count = len(primal_edges)
    face_count = len(oriented_faces)
    rank_sum = primal_vertices + face_count - 2
    if rank_sum != edge_count:
        raise ValueError("Euler's spherical count must hold for these polyhedra")

    edge_index = {edge: index for index, edge in enumerate(primal_edges)}
    incidence = np.zeros((primal_vertices, edge_count), dtype=np.int64)
    boundary = np.zeros((face_count, edge_count), dtype=np.int64)
    for index, (tail, head) in enumerate(primal_edges):
        incidence[tail, index] = -1
        incidence[head, index] = 1
    for face_index, face in enumerate(oriented_faces):
        for a, b in zip(face, face[1:] + face[:1]):
            column = edge_index[min(a, b), max(a, b)]
            boundary[face_index, column] = 1 if a < b else -1

    dual_edges = []
    for column in boundary.T:
        tails = np.flatnonzero(column == -1)
        heads = np.flatnonzero(column == 1)
        if len(tails) != 1 or len(heads) != 1:
            raise ValueError("every edge must occur oppositely in two face boundaries")
        dual_edges.append((int(tails[0]), int(heads[0])))
    if np.any(incidence @ boundary.T):
        raise AssertionError("the boundary of a face boundary must vanish")

    paired_edges = tuple(dual_edges)
    primal, primal_trees = _graph_covariance(primal_vertices, primal_edges)
    dual, dual_trees = _graph_covariance(face_count, paired_edges)
    # Both grounded inverses exist, so their cut ranks are the node counts
    # minus one. Orthogonality and Euler's count exhaust the edge space.
    return PlanarDuality(
        primal_edges, paired_edges, primal, dual, (primal_vertices - 1, face_count - 1),
        (primal_trees, dual_trees),
    )


@cache
def carrier_planar_duality() -> PlanarDuality:
    """Pair the carrier seams with the regular dodecahedral face graph."""
    carrier = icosahedral_carrier()
    return _paired_covariances(carrier.port_count, carrier.edges, OPH_ORIENTED_BASE_FACES)


@cache
def ivm_planar_duality() -> PlanarDuality:
    """Pair cuboctahedral edges with the rhombic dodecahedral Voronoi graph."""
    vertices = ivm_solid_vertices()["cuboctahedron"]
    return _paired_covariances(len(vertices), cuboctahedron_edges(), convex_facets(vertices))


def planar_duality_census() -> dict[str, dict[str, object]]:
    """Compare every covariance entry and both exact tree counts in each pair."""
    results = {}
    for name, pair in (("carrier", carrier_planar_duality()), ("ivm", ivm_planar_duality())):
        primal, dual = pair.primal_covariance, pair.dual_covariance
        count = len(pair.primal_edges)
        matching_entries = sum(
            primal[i][j] + dual[i][j] == int(i == j)
            for i in range(count)
            for j in range(count)
        )
        results[name] = {
            "edge_count": count,
            "cut_ranks": pair.cut_ranks,
            "primal_resistances": tuple(sorted({primal[i][i] for i in range(count)})),
            "dual_resistances": tuple(sorted({dual[i][i] for i in range(count)})),
            "tree_counts": pair.tree_counts,
            "matrix_entries_checked": count * count,
            "covariances_sum_to_identity": matching_entries == count * count,
        }
    return results
