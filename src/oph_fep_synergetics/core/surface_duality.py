"""Exact Gaussian cut and harmonic covariances on a square torus."""

from __future__ import annotations

from fractions import Fraction
from functools import cache

import numpy as np

from .invariants import _exact_inverse
from .planar_duality import Edge
from .weighted_duality import _normalised_covariance


def _square_torus() -> tuple[tuple[Edge, ...], tuple[Edge, ...], np.ndarray, np.ndarray]:
    """Build the 3×3 quotient with horizontal edges before vertical edges."""
    side = 3
    sites = tuple((x, y) for y in range(side) for x in range(side))

    def vertex(x: int, y: int) -> int:
        return x % side + side * (y % side)

    edges = tuple(
        (vertex(x, y), vertex(x + dx, y + dy))
        for dx, dy in ((1, 0), (0, 1))
        for x, y in sites
    )
    incidence = np.zeros((len(sites), len(edges)), dtype=int)
    positions = {}
    for index, (tail, head) in enumerate(edges):
        incidence[tail, index] = -1
        incidence[head, index] = 1
        positions[tail, head] = (index, 1)
        positions[head, tail] = (index, -1)

    boundary = np.zeros_like(incidence)
    for index, (x, y) in enumerate(sites):
        face = (vertex(x, y), vertex(x + 1, y), vertex(x + 1, y + 1), vertex(x, y + 1))
        for tail, head in zip(face, face[1:] + face[:1], strict=True):
            edge, sign = positions[tail, head]
            boundary[index, edge] = sign
    dual = []
    for column in boundary.T:
        (tail,) = np.flatnonzero(column == -1)
        (head,) = np.flatnonzero(column == 1)
        dual.append((int(tail), int(head)))
    return edges, tuple(dual), incidence, boundary


def _rank_fraction(matrix: np.ndarray) -> int:
    """Count pivots by exact row elimination."""
    reduced = matrix.copy()
    rank = 0
    for column in range(reduced.shape[1]):
        pivot = next((row for row in range(rank, len(reduced)) if reduced[row, column]), None)
        if pivot is None:
            continue
        reduced[[rank, pivot]] = reduced[[pivot, rank]]
        for row in range(rank + 1, len(reduced)):
            if reduced[row, column]:
                factor = reduced[row, column] / reduced[rank, column]
                reduced[row, column:] -= factor * reduced[rank, column:]
        rank += 1
    return rank


def _harmonic_from_periods(
    weighted_boundary: np.ndarray, weighted_incidence: np.ndarray,
    roots: tuple[Fraction, ...],
) -> tuple[np.ndarray, np.ndarray]:
    """Find closed and coclosed representatives with unit winding periods."""
    periods = np.full((2, len(roots)), Fraction(0), dtype=object)
    # The bottom horizontal and left vertical cycles each wind once.
    for direction, cycle in enumerate(((0, 1, 2), (9, 12, 15))):
        for edge in cycle:
            periods[direction, edge] = 1 / roots[edge]

    # Remove one redundant equation per incidence. The remaining equations
    # and the two periods specify unique representatives independently of
    # the two Gaussian covariance matrices.
    equations = np.vstack((weighted_incidence[1:], weighted_boundary[1:], periods))
    inverse = np.asarray(_exact_inverse(equations.tolist()), dtype=object)
    representatives = inverse[:, -2:]
    gram = representatives.T @ representatives
    inverse_gram = np.asarray(_exact_inverse(gram.tolist()), dtype=object)
    return representatives @ inverse_gram @ representatives.T, gram


@cache
def _toroidal_census() -> tuple[dict[str, object], dict[str, dict[str, object]]]:
    primal_edges, dual_edges, incidence, boundary = _square_torus()
    vertices, edges = incidence.shape
    faces = len(boundary)
    unit = (Fraction(1),) * edges
    strengthened = (Fraction(2), *unit[1:])
    reciprocal = tuple(1 / value for value in strengthened)
    rescaled = tuple(3 * value for value in reciprocal)
    metadata = {
        "cells": (vertices, edges, faces),
        "genus": (2 - vertices + edges - faces) // 2,
        "changed_edge": primal_edges[0],
        "changed_weight": strengthened[0] ** 2,
        "entries_checked": edges**2,
        "boundary_of_boundary_zero": not np.any(incidence @ boundary.T),
    }
    identity = np.eye(edges, dtype=int)
    results = {}
    for name, primal_roots, dual_roots in (
        ("unit", unit, unit),
        ("reciprocal", strengthened, reciprocal),
        ("rescaled_reciprocal", strengthened, rescaled),
        ("copied", strengthened, strengthened),
    ):
        primal, _ = _normalised_covariance(primal_edges, vertices - 1, primal_roots)
        dual, _ = _normalised_covariance(dual_edges, faces - 1, dual_roots)
        remainder = identity - np.asarray(primal, dtype=object) - np.asarray(dual, dtype=object)
        weighted_incidence = incidence * np.asarray(primal_roots, dtype=object)
        weighted_boundary = boundary * np.asarray(dual_roots, dtype=object)
        cross_rank = _rank_fraction(weighted_incidence @ weighted_boundary.T)
        row = {
            "primal_diagonal": primal[0][0],
            "dual_diagonal": dual[0][0],
            "remainder_diagonal": remainder[0, 0],
            "remainder_trace": np.trace(remainder),
            "remainder_rank": _rank_fraction(remainder),
            "cross_pairing_rank": cross_rank,
            "is_projector": bool(np.array_equal(remainder @ remainder, remainder)),
        }
        if cross_rank == 0:
            harmonic, gram = _harmonic_from_periods(weighted_boundary, weighted_incidence, primal_roots)
            row["harmonic_periods_match"] = bool(np.array_equal(remainder, harmonic))
            row["period_gram"] = tuple(tuple(line) for line in gram.tolist())
        if name == "unit":
            direction = np.repeat(np.arange(2), vertices)
            directional = (direction[:, None] == direction[None, :]) * Fraction(1, vertices)
            row["uniform_direction_blocks_match"] = bool(np.array_equal(remainder, directional))
        results[name] = row
    return metadata, results


def toroidal_duality_census() -> dict[str, object]:
    """Return exact torus tests for reciprocal and copied seam conductances.

    The first horizontal edge has conductance four in the weighted trials.
    The unit and reciprocal remainders are compared entry by entry with the
    Gram projector of independent harmonic representatives with unit periods.
    """
    metadata, results = _toroidal_census()
    return {**metadata, **{name: row.copy() for name, row in results.items()}}
