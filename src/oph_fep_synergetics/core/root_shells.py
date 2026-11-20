"""Exact Gaussian and repair covariance checks on nearest-neighbour root shells."""

from __future__ import annotations

from fractions import Fraction
from functools import cache
from itertools import combinations, product

import numpy as np

from .native_dynamics import _attempt_counts


def _root_graph(dimension: int) -> tuple[np.ndarray, tuple[tuple[int, int], ...]]:
    """Give the integer D_d roots and edges of squared length two."""
    points = []
    for i, j in combinations(range(dimension), 2):
        for a, b in product((-1, 1), repeat=2):
            point = [0] * dimension
            point[i], point[j] = a, b
            points.append(point)
    coordinates = np.array(points, dtype=np.int64)
    first, second = np.where(np.triu(coordinates @ coordinates.T == 1, k=1))
    return coordinates, tuple(zip(first.tolist(), second.tolist(), strict=True))


def _quadratic_shell(dimension: int) -> tuple[object, ...]:
    points, edges = _root_graph(dimension)
    squares = points**2
    # Each column is twice the evaluation on the corresponding unit directions.
    diagonal = [squares[:, :k].sum(axis=1) - k * squares[:, k] for k in range(1, dimension)]
    mixed = [2 * points[:, i] * points[:, j] for i, j in combinations(range(dimension), 2)]
    fields = np.column_stack(diagonal + mixed)
    norms_sq = [k * (k + 1) for k in range(1, dimension)] + [2] * len(mixed)
    first, second = np.asarray(edges).T
    gradients = fields[first] - fields[second]
    laplacian_fields = np.zeros_like(fields)
    np.add.at(laplacian_fields, first, gradients)
    np.add.at(laplacian_fields, second, -gradients)
    gram = fields.T @ fields
    precision = gradients.T @ gradients
    expected_frame = (dimension - 2, 2)
    expected_eigenvalues = (2 * dimension, 4 * (dimension - 2))
    expected_precisions = tuple(f * e for f, e in zip(expected_frame, expected_eigenvalues, strict=True))
    frame_factors = tuple(Fraction(int(gram[i, i]), 4 * norms_sq[i]) for i in (0, dimension - 1))
    precision_factors = tuple(Fraction(int(precision[i, i]), 4 * norms_sq[i]) for i in (0, dimension - 1))

    def repeated(factors: tuple[int, int]) -> np.ndarray:
        return np.repeat(factors, [len(diagonal), len(mixed)])

    matches = (
        not np.any(fields.sum(axis=0))
        and np.array_equal(points.T @ points, np.eye(dimension, dtype=int) * (4 * (dimension - 1)))
        and np.array_equal(gram, np.diag(4 * np.array(norms_sq) * repeated(expected_frame)))
        and np.array_equal(precision, np.diag(4 * np.array(norms_sq) * repeated(expected_precisions)))
        and np.array_equal(laplacian_fields, fields * repeated(expected_eigenvalues))
    )
    return (
        dimension, len(points), len(edges), frame_factors, precision_factors,
        tuple(1 / value for value in precision_factors), matches,
    )


@cache
def root_shell_quadratics() -> tuple[tuple[object, ...], ...]:
    """Check the frame and precision formula in dimensions three through seven."""
    return tuple(_quadratic_shell(dimension) for dimension in range(3, 8))


@cache
def _lag_shell() -> tuple[tuple[str, object], ...]:
    points, edges = _root_graph(4)
    n, dimension = points.shape
    pairs = tuple(combinations(range(n), 2))
    states = tuple(tuple(int(i in pair) for i in range(n)) for pair in pairs)
    values = np.asarray(states) @ points
    counts = _attempt_counts(states, edges)
    tail, head = edges[0]
    chord = points[tail] - points[head]
    chord_norm = int(chord @ chord)
    column_norm = int((points.T @ points)[0, 0])
    scale = int((values.T @ values)[0, 0])
    laplacian = np.zeros((n, n), dtype=np.int64)
    first, second = np.asarray(edges).T
    laplacian[first, second] = laplacian[second, first] = -1
    np.fill_diagonal(laplacian, -laplacian.sum(axis=1))
    extra = len(edges)
    weighted_counts = counts + extra * _attempt_counts(states, ((tail, head),))
    edge = np.eye(n, dtype=np.int64)[tail] - np.eye(n, dtype=np.int64)[head]
    weighted_laplacian = laplacian + extra * np.outer(edge, edge)

    def verify_lags(kernel, node_laplacian, attempts):
        moved = kernel @ values
        native_lag_one = values.T @ moved
        native_lag_two = moved.T @ moved
        native_defect = scale * native_lag_two - native_lag_one @ native_lag_one
        transformed = node_laplacian @ points
        first_moment = column_norm * attempts * np.eye(dimension, dtype=np.int64) - points.T @ transformed
        node_defect = column_norm * transformed.T @ transformed - (points.T @ transformed) @ (points.T @ transformed)
        denominator = (scale * attempts)**2
        # The two derivations have different denominators. Cancel their exact ratio
        # before comparing, to keep these finite integer products in range.
        factor = scale // column_norm
        assert scale == factor * column_norm
        native_mean_matches = np.array_equal(native_lag_one, factor * first_moment)
        defects_match = np.array_equal(native_defect, factor**2 * node_defect)
        stationary = (
            np.array_equal(kernel, kernel.T) and np.all(kernel.sum(axis=1) == attempts)
            and np.array_equal(values.T @ values, scale * np.eye(dimension, dtype=np.int64))
            and not np.any(values.sum(axis=0))
        )
        return {
            "routes_match": bool(native_mean_matches and defects_match and stationary),
            "lag_one": Fraction(int(chord @ native_lag_one @ chord), scale * attempts * chord_norm),
            "lag_two": Fraction(int(chord @ native_lag_two @ chord), scale * attempts**2 * chord_norm),
            "defect_trace": Fraction(int(np.trace(native_defect)), denominator),
            "rank_one_chord": (
                bool(np.any(native_defect))
                and np.array_equal(native_defect * chord_norm, np.trace(native_defect) * np.outer(chord, chord))
            ),
        }

    unit = verify_lags(counts, laplacian, 2 * len(edges))
    weighted = verify_lags(weighted_counts, weighted_laplacian, 2 * (len(edges) + extra))
    return tuple({
        "vertices": n,
        "seams": len(edges),
        "binary_states": len(states),
        "changed_weight": extra + 1,
        "unit_coordinate_modes": np.array_equal(laplacian @ points, 4 * points),
        "unit": tuple(unit.items()),
        "weighted": tuple(weighted.items()),
    }.items())


def root_shell_lag_census() -> dict[str, object]:
    """Verify the linear covariance closure criterion on the 24-cell."""
    values = dict(_lag_shell())
    return {**values, "unit": dict(values["unit"]), "weighted": dict(values["weighted"])}
