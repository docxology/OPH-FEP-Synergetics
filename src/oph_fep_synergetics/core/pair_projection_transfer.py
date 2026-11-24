"""Exact comparisons of pair-mean precision and binary pair-total repair."""

from __future__ import annotations

from collections.abc import Iterable
from fractions import Fraction
from functools import cache
from itertools import combinations, pairwise

import numpy as np

from .invariants import _exact_inverse
from .native_dynamics import _attempt_counts

Pairing = tuple[tuple[int, int], ...]


def _adjacency(size: int, edges: Iterable[tuple[int, int]]) -> np.ndarray:
    graph = np.zeros((size, size), dtype=object)
    for i, j in edges:
        graph[i, j] = graph[j, i] = Fraction(1)
    return graph


def _pair_basis(size: int, pairs: Pairing) -> tuple[np.ndarray, np.ndarray]:
    sums = np.zeros((size, len(pairs)), dtype=object)
    differences = np.zeros_like(sums)
    for column, (first, second) in enumerate(pairs):
        sums[first, column] = sums[second, column] = 1
        differences[first, column], differences[second, column] = 1, -1
    return sums, differences


def _binary_transition(graph: np.ndarray) -> tuple[np.ndarray, np.ndarray, Fraction]:
    size = len(graph)
    states = tuple(
        tuple(int(label in occupied) for label in range(size))
        for occupied in combinations(range(size), 2)
    )
    counts = np.zeros((len(states), len(states)), dtype=object)
    mass = Fraction(0)
    for i, j in combinations(range(size), 2):
        if (weight := graph[i, j]) > 0:
            counts += weight * _attempt_counts(states, ((i, j),)).astype(object)
            mass += weight
    return np.asarray(states, dtype=object), counts / (2 * mass), mass


def _count_lumpability(loads: np.ndarray, transition: np.ndarray) -> bool:
    tuples = tuple(map(tuple, loads))
    categories = tuple(sorted(set(tuples)))
    membership = np.asarray(
        [[record == category for category in categories] for record in tuples], dtype=bool
    )
    probability_rows = transition @ membership.astype(object)
    for column in membership.T:
        rows = probability_rows[column]
        if not np.all(rows == rows[0]):
            return False
    return True


def _return_probabilities(
    transition: np.ndarray, endpoint_loads: np.ndarray, middle_loads: np.ndarray
) -> tuple[Fraction, Fraction]:
    """Return the probability of a full first pair with and without its past."""
    size = endpoint_loads.shape[1]
    endpoint = np.all(endpoint_loads == (2,) + (0,) * (size - 1), axis=1)
    intermediate = np.all(middle_loads == (1, 1) + (0,) * (size - 2), axis=1)
    from_endpoint = transition[endpoint][:, intermediate].sum(axis=0)
    to_endpoint = transition[intermediate][:, endpoint].sum(axis=1)
    return (
        (from_endpoint @ to_endpoint) / sum(from_endpoint),
        sum(to_endpoint) / Fraction(len(to_endpoint)),
    )


def _uniform_cross_blocks(graph: np.ndarray, pairs: Pairing) -> bool:
    return all(len(set(graph[np.ix_(left, right)].flat)) == 1 for left, right in combinations(pairs, 2))


def _inverse(matrix: np.ndarray) -> np.ndarray:
    return np.asarray(_exact_inverse(matrix.tolist()), dtype=object)


def _pair_projection_case(graph: np.ndarray, pairs: Pairing) -> dict[str, object]:
    size, site_count = len(graph), len(pairs)
    replication, differences = _pair_basis(size, pairs)
    laplacian = np.diag(graph.sum(axis=1)) - graph
    half = Fraction(1, 2)
    visible = half * (replication.T @ laplacian @ replication)
    hidden = half * (differences.T @ laplacian @ differences)
    coupling = half * (differences.T @ laplacian @ replication)
    imbalance = differences.T @ graph @ replication
    correction = coupling.T @ _inverse(hidden) @ coupling
    site_precision = 2 * (visible - correction)
    literal_laplacian = all(site_precision[i, j] <= 0 for i, j in combinations(range(site_count), 2))

    # Invert the full centered Gaussian first, then project the covariance
    # and invert again on centered literal site means.
    joint_constant, site_constant = Fraction(1, size), Fraction(1, site_count)
    joint_covariance = _inverse(laplacian + joint_constant) - joint_constant
    projected_covariance = replication.T @ joint_covariance @ replication * Fraction(1, 4)
    direct_precision = _inverse(projected_covariance + site_constant) - site_constant

    states, transition, mass = _binary_transition(graph)
    totals = states @ replication
    centered_totals = totals - Fraction(2, site_count)
    covariance_scale = Fraction(4 * (size - 2), size * (size - 1))
    mean_factor = 1 / (len(states) * covariance_scale)
    native_covariance = centered_totals.T @ centered_totals * mean_factor
    one_step = transition @ centered_totals
    native_lag_one = centered_totals.T @ one_step * mean_factor
    native_lag_two = one_step.T @ one_step * mean_factor
    site_centering = np.eye(site_count, dtype=object) - site_constant
    lag_defect = coupling.T @ coupling / (4 * mass**2)
    return {
        "edges": tuple((i, j) for i, j in combinations(range(size), 2) if graph[i, j] > 0),
        "pairs": pairs,
        "pair_balanced": not bool(np.any(imbalance)),
        "balance_formula_matches": np.array_equal(
            coupling, half * (np.diag(imbalance.sum(axis=1)) - imbalance)
        ),
        "uniform_cross_blocks": _uniform_cross_blocks(graph, pairs),
        "static_correction_vanishes": not bool(np.any(correction)),
        "endpoint_precision": site_precision[0, -1],
        "literal_site_precision": tuple(map(tuple, site_precision)),
        "site_precision_is_laplacian": literal_laplacian,
        "lag_two_defect": tuple(map(tuple, lag_defect)),
        "lag_two_defect_trace": np.trace(lag_defect),
        "strongly_lumpable": _count_lumpability(totals, transition),
        "return_probabilities": _return_probabilities(transition, totals, totals),
        "gaussian_marginal_matches": np.array_equal(site_precision, direct_precision),
        "native_lags_match": (
            np.array_equal(transition, transition.T)
            and bool(np.all(transition.sum(axis=1) == 1))
            and np.array_equal(native_covariance, site_centering)
            and np.array_equal(native_lag_one, site_centering - visible / (2 * mass))
            and np.array_equal(native_lag_two - native_lag_one @ native_lag_one, lag_defect)
        ),
    }


@cache
def _pair_projection_transfer_census() -> dict[str, dict[str, object]]:
    path4_edges = tuple(pairwise(range(4)))
    graphs = {
        "path4": _adjacency(4, path4_edges),
        "path6": _adjacency(6, pairwise(range(6))),
        "cycle4": _adjacency(4, (*path4_edges, (0, 3))),
        "complete4": _adjacency(4, combinations(range(4), 2)),
    }
    return {
        name: _pair_projection_case(graph, tuple((i, i + 1) for i in range(0, len(graph), 2)))
        for name, graph in graphs.items()
    }


def pair_projection_transfer_census() -> dict[str, dict[str, object]]:
    """Compare the centered reductions on two paths, a cycle and a complete graph."""
    return {name: row.copy() for name, row in _pair_projection_transfer_census().items()}


@cache
def _switching_pair_census() -> dict[str, object]:
    graph = _adjacency(4, combinations(range(4), 2))
    endpoint, middle = ((0, 1), (2, 3)), ((0, 2), (1, 3))
    endpoint_basis, endpoint_differences = _pair_basis(4, endpoint)
    middle_basis, middle_differences = _pair_basis(4, middle)
    states, transition, _ = _binary_transition(graph)
    endpoint_loads, middle_loads = states @ endpoint_basis, states @ middle_basis
    x, y = endpoint_loads[:, 0] - 1, middle_loads[:, 0] - 1
    variance_mass = Fraction(x @ x)
    adjacent = (x @ transition @ y / variance_mass, y @ transition @ x / variance_mass)
    endpoint_covariance = x @ transition @ transition @ x / variance_mass
    return {
        "adjacent_covariances": adjacent,
        "endpoint_covariance": endpoint_covariance,
        "lag_two_residual": endpoint_covariance - adjacent[0] * adjacent[1],
        "return_probabilities": _return_probabilities(transition, endpoint_loads, middle_loads),
        "all_fixed_partitions_lumpable": (
            _count_lumpability(endpoint_loads, transition)
            and _count_lumpability(middle_loads, transition)
        ),
        "all_fixed_partitions_balanced": not (
            np.any(endpoint_differences.T @ graph @ endpoint_basis)
            or np.any(middle_differences.T @ graph @ middle_basis)
        ),
    }


def switching_pair_census() -> dict[str, object]:
    """Check Gaussian covariance and conditional count memory for K4's A-B-A pairing."""
    return _switching_pair_census().copy()
