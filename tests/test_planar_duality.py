"""Planar duality pairs signed Gaussian seam contrasts."""

from __future__ import annotations

from fractions import Fraction
from itertools import product

import numpy as np

from oph_fep_synergetics.core.network import graphs_are_isomorphic
from oph_fep_synergetics.core.planar_duality import (
    carrier_planar_duality,
    ivm_planar_duality,
    planar_duality_census,
)
from oph_fep_synergetics.geometry.exact import PHI_EXACT, convex_facets, ivm_solid_vertices


def test_primal_and_dual_edge_covariances_sum_to_identity():
    assert planar_duality_census() == {
        "carrier": {
            "edge_count": 30,
            "cut_ranks": (11, 19),
            "primal_resistances": (Fraction(11, 30),),
            "dual_resistances": (Fraction(19, 30),),
            "tree_counts": (5184000, 5184000),
            "matrix_entries_checked": 900,
            "covariances_sum_to_identity": True,
        },
        "ivm": {
            "edge_count": 24,
            "cut_ranks": (11, 13),
            "primal_resistances": (Fraction(11, 24),),
            "dual_resistances": (Fraction(13, 24),),
            "tree_counts": (331776, 331776),
            "matrix_entries_checked": 576,
            "covariances_sum_to_identity": True,
        },
    }
    # Cross-check both covariance matrices via the centered spectral solve,
    # which does not reuse the grounded rational inverse.
    for pair in (carrier_planar_duality(), ivm_planar_duality()):
        for edges, rank, covariance in (
            (pair.primal_edges, pair.cut_ranks[0], pair.primal_covariance),
            (pair.dual_edges, pair.cut_ranks[1], pair.dual_covariance),
        ):
            incidence = np.zeros((rank + 1, len(edges)))
            for i, (a, b) in enumerate(edges):
                incidence[a, i], incidence[b, i] = -1, 1
            precision = incidence @ incidence.T
            centered_covariance = np.linalg.pinv(precision)
            np.testing.assert_allclose(
                np.asarray(covariance, dtype=float),
                incidence.T @ centered_covariance @ incidence,
                atol=1e-12,
            )
            assert np.linalg.matrix_rank(precision) == rank


def _adjacency(size, edges):
    adjacency = [[0] * size for _ in range(size)]
    for u, v in edges:
        adjacency[u][v] = adjacency[v][u] = 1
    return adjacency


def test_dual_face_graphs_match_geometric_dodecahedral_graphs():
    inverse_phi = PHI_EXACT - 1
    vertices = list(product((-1, 1), repeat=3))
    for short_sign, long_sign in product((-1, 1), repeat=2):
        a, b = short_sign * inverse_phi, long_sign * PHI_EXACT
        vertices.extend(((0, a, b), (a, b, 0), (b, 0, a)))
    edge = 2 * inverse_phi
    dodecahedron = [
        [
            int(sum((a - b) * (a - b) for a, b in zip(left, right)) == edge * edge)
            for right in vertices
        ]
        for left in vertices
    ]
    carrier = carrier_planar_duality()
    assert graphs_are_isomorphic(
        _adjacency(carrier.cut_ranks[1] + 1, carrier.dual_edges), dodecahedron,
    )

    cell_vertices = ivm_solid_vertices()["rhombic_dodecahedron"]
    cell_edges = [
        (a, b)
        for face in convex_facets(cell_vertices)
        for a, b in zip(face, face[1:] + face[:1])
    ]
    ivm = ivm_planar_duality()
    assert graphs_are_isomorphic(
        _adjacency(ivm.cut_ranks[1] + 1, ivm.dual_edges),
        _adjacency(len(cell_vertices), cell_edges),
    )
