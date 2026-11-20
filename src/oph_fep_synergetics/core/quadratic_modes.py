"""Exact quadratic field geometry on the carrier and the first FCC shell."""

from __future__ import annotations

from fractions import Fraction

import numpy as np

from ..geometry.exact import Surd5, exact_icosahedron_vertices, ivm_solid_vertices
from ..geometry.lattice import cuboctahedron_adjacency
from .carrier import icosahedral_carrier
from .network import graphs_are_isomorphic, laplacian_matrix


def _rational(value: Surd5 | Fraction) -> Fraction:
    exact = Surd5.of(value)
    if exact.b:
        raise AssertionError("the normalized matrix coefficient must be rational")
    return exact.a


def quadratic_mode_census() -> dict[str, dict[str, object]]:
    """Compare quadratic evaluation frames and pullback Gaussian precisions.

    Evaluate each trace-free quadratic at the unit directions of the vertices.
    The Frobenius-orthogonal basis uses two diagonal forms and the three
    off-diagonal forms, with squared norms 2, 6, 2, 2, 2. Keeping the basis
    unnormalized avoids adjoining square roots to the exact coordinate field.
    """
    carrier = exact_icosahedron_vertices()
    carrier_adj = tuple(
        tuple(
            int(sum(delta * delta for delta in (a - b for a, b in zip(left, right, strict=True))) == 4)
            for right in carrier
        )
        for left in carrier
    )
    oph_adj = icosahedral_carrier().adjacency_matrix().astype(int).tolist()
    if not graphs_are_isomorphic(carrier_adj, oph_adj):
        raise AssertionError("the exact icosahedron must realize the OPH graph")
    fcc = ivm_solid_vertices()["cuboctahedron"]
    norms = (2, 6, 2, 2, 2)
    result = {}
    for name, points, adjacency in (
        ("carrier", carrier, carrier_adj),
        ("fcc", fcc, cuboctahedron_adjacency()),
    ):
        radius_squared = sum(value * value for value in points[0])
        modes = np.asarray([
            (x * x - y * y, x * x + y * y - 2 * z * z, 2 * x * y, 2 * x * z, 2 * y * z)
            for x, y, z in points
        ], dtype=object) / radius_squared
        frame = modes.T @ modes
        laplacian = np.asarray(laplacian_matrix(adjacency))
        stiffness = laplacian @ modes
        prec = modes.T @ stiffness
        differences = np.asarray([
            modes[i] - modes[j]
            for i, row in enumerate(adjacency) for j, edge in enumerate(row)
            if i < j and edge
        ], dtype=object)
        frame_values = tuple(_rational(frame[i, i] / norm) for i, norm in enumerate(norms))
        precisions = tuple(_rational(prec[i, i] / norm) for i, norm in enumerate(norms))
        eigenvalues = tuple(p / f for p, f in zip(precisions, frame_values, strict=True))
        xyz = np.asarray(points, dtype=object)
        result[name] = {
            "linear_frame_isotropic": bool(np.array_equal(xyz.T @ xyz, np.eye(3, dtype=int) * (4 * radius_squared))),
            "quadratic_modes_centered": all(value == 0 for value in modes.sum(axis=0)),
            "frame_orthogonal": np.array_equal(frame, np.diag(frame.diagonal())),
            "precision_orthogonal": np.array_equal(prec, np.diag(prec.diagonal())),
            "edge_and_laplacian_precision_agree": np.array_equal(prec, differences.T @ differences),
            "laplacian_preserves_modes": np.array_equal(stiffness, modes * eigenvalues),
            "frame": frame_values,
            "laplacian_eigenvalues": eigenvalues,
            "precision": precisions,
            "variance": tuple(Fraction(1, 1) / value for value in precisions),
        }
    return result
