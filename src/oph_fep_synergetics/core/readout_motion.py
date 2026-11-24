"""Exact moving-readout tests on a square with alternating conductances."""

from fractions import Fraction
from functools import cache
from itertools import combinations

import numpy as np

from .invariants import _exact_inverse
from .native_dynamics import _attempt_counts
from .surface_duality import _rank_fraction


_SQUARE = ((0, 1), (2, 3), (0, 2), (1, 3))


def _vertex_step(weights: tuple[int, ...]) -> np.ndarray:
    identity = np.eye(4, dtype=object)
    laplacian = np.zeros_like(identity)
    for (i, j), weight in zip(_SQUARE, weights, strict=True):
        direction = identity[:, i] - identity[:, j]
        laplacian += weight * np.outer(direction, direction)
    return identity - laplacian / Fraction(2 * sum(weights))


def _native_step(edge_counts: tuple[np.ndarray, ...], weights: tuple[int, ...]) -> np.ndarray:
    counts = sum(weight * block for weight, block in zip(weights, edge_counts, strict=True))
    return counts / Fraction(2 * sum(weights))


def _transport(steps: tuple[np.ndarray, ...], first: int, last: int) -> np.ndarray:
    product = np.eye(len(steps[0]), dtype=object)
    for step in steps[first:last]:
        product = product @ step
    return product


def _cross_covariances(
    steps: tuple[np.ndarray, ...], records: tuple[np.ndarray, ...], scale: Fraction
) -> dict[tuple[int, int], np.ndarray]:
    return {
        (i, j): left.T @ _transport(steps, i, j) @ records[j] / scale
        for i, left in enumerate(records)
        for j in range(i, len(records))
    }


def _inverse(matrix: np.ndarray) -> np.ndarray:
    rational = [[Fraction(cell) for cell in row] for row in matrix.tolist()]
    return np.asarray(_exact_inverse(rational), dtype=object)


def _residual(
    covariance: dict[tuple[int, int], np.ndarray], first: int, middle: int, last: int
) -> np.ndarray:
    return covariance[first, last] - (
        covariance[first, middle] @ _inverse(covariance[middle, middle]) @ covariance[middle, last]
    )


def _hidden_components(
    steps: tuple[np.ndarray, ...], frames: tuple[np.ndarray, ...], first: int, middle: int, last: int
) -> tuple[np.ndarray, np.ndarray]:
    visible = frames[middle]
    hidden = np.eye(len(visible), dtype=object) - visible @ _inverse(visible.T @ visible) @ visible.T
    left = hidden @ _transport(steps, first, middle).T @ frames[first]
    right = hidden @ _transport(steps, middle, last) @ frames[last]
    return left, right


def _minimum_common_dimension(steps: tuple[np.ndarray, ...], frames: tuple[np.ndarray, ...]) -> int:
    span = np.column_stack(frames) / Fraction(1)
    rank = _rank_fraction(span)
    while rank < len(span) - 1:
        extended = np.column_stack((span, *(step @ span for step in steps)))
        grown_rank = _rank_fraction(extended)
        if grown_rank == rank:
            break
        span, rank = extended, grown_rank
    return rank


@cache
def _moving_readout_census() -> dict[str, object]:
    # Independent routes use joint incidence and concrete native repair attempts.
    states = tuple(tuple(int(vertex in loads) for vertex in range(4)) for loads in combinations(range(4), 2))
    state_matrix = np.asarray(states, dtype=object)
    edge_counts = tuple(_attempt_counts(states, (edge,)).astype(object) for edge in _SQUARE)
    weights = ((1, 1, 2, 2), (2, 2, 1, 1), (1, 1, 2, 2))
    vertex_steps = tuple(_vertex_step(row) for row in weights)
    native_steps = tuple(_native_step(edge_counts, row) for row in weights)
    # For k=2, each sign mode has variance (4 terms) times k(4-k)/(4*3).
    mode_variance = Fraction(4, 3)
    x, y, z = (np.array(mode, dtype=object).reshape(4, 1) for mode in (
        (1, 1, -1, -1), (1, -1, 1, -1), (1, -1, -1, 1),
    ))
    three_time = (
        (x + y, x, x + y), (x + y, x, x - y), (x + y, x, x + z),
    )
    delayed = (x + y, x, x, x + y)
    changing_basis = tuple(np.column_stack(columns) for columns in (
        (x, y), (x + y, x - y), (2 * x, y - x), (x - y, 2 * y),
    ))

    transported = (x + y, 4 * x + 5 * y, x + y, 4 * x + 5 * y)
    all_frames = (*three_time, delayed, changing_basis, transported)
    covariances = tuple(_cross_covariances(vertex_steps, frames, Fraction(4)) for frames in all_frames)
    native_covariances = tuple(
        _cross_covariances(native_steps, tuple(state_matrix @ frame for frame in frames),
                           len(states) * mode_variance)
        for frames in all_frames
    )
    native_matches_transport = all(
        np.array_equal(value, native[key])
        for covariance, native in zip(covariances, native_covariances, strict=True)
        for key, value in covariance.items()
    )
    split_factorisation_matches = all(
        np.array_equal(_residual(covariance, i, j, k), (left.T @ right) / Fraction(4))
        for frames, covariance in zip(all_frames, covariances, strict=True)
        for i, j, k in combinations(range(len(frames)), 3)
        for left, right in (_hidden_components(vertex_steps, frames, i, j, k),)
    )
    hidden_norms = tuple(
        ((left.T @ left)[0, 0] / Fraction(4), (right.T @ right)[0, 0] / Fraction(4))
        for frames in three_time
        for left, right in (_hidden_components(vertex_steps, frames, 0, 1, 2),)
    )
    four = covariances[3]
    dimensions = (
        _minimum_common_dimension(vertex_steps, delayed),
        _minimum_common_dimension((*vertex_steps, _vertex_step((2, 1, 2, 2))), delayed),
    )
    return {
        "state_count": len(states),
        "covariance_scale": mode_variance,
        "defects": tuple(_residual(cov, 0, 1, 2)[0, 0] for cov in covariances[:3]),
        "hidden_component_norms": hidden_norms,
        "adjacent_defects": (_residual(four, 0, 1, 2)[0, 0], _residual(four, 1, 2, 3)[0, 0]),
        "endpoint_defects": (_residual(four, 0, 1, 3)[0, 0], _residual(four, 0, 2, 3)[0, 0]),
        "endpoint_covariance": four[0, 3][0, 0],
        "adjacent_covariances": tuple(four[i, i + 1][0, 0] for i in range(3)),
        "minimum_common_dimensions": dimensions,
        "changing_basis_closes": all(
            not np.any(_residual(covariances[4], i, j, k)) for i, j, k in combinations(range(4), 3)
        ),
        "transported_moving_closes": (
            _rank_fraction(np.column_stack(transported[:2]) / Fraction(1)) == 2
            and all(np.array_equal(vertex_steps[i] @ transported[i], factor * transported[i + 1])
                    for i, factor in enumerate((Fraction(1, 6), Fraction(10, 3), Fraction(1, 6))))
            and all(not np.any(_residual(covariances[5], i, j, k)) for i, j, k in combinations(range(4), 3))
        ),
        "native_matches_transport": native_matches_transport,
        "split_factorisation_matches": split_factorisation_matches,
        "uniform_stationary": all(np.array_equal(step, step.T) and bool(np.all(step.sum(axis=1) == 1))
                                  for step in native_steps),
    }


def moving_readout_census() -> dict[str, object]:
    """Compare both covariance routes for moving readouts on a weighted square."""
    return _moving_readout_census().copy()
