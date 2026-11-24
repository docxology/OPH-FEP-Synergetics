"""Static and temporal effects of projecting onto Jitterbug joint pairs."""

from __future__ import annotations

from fractions import Fraction
from functools import cache
from itertools import combinations

import numpy as np

from ..geometry.jitterbug import jitterbug_matching_edges, jitterbug_rigid_edges, jitterbug_snapshot
from .invariants import _exact_inverse
from .jitterbug_inference import _graph_matrices
from .native_dynamics import _attempt_counts
from .network import laplacian_matrix


def jitterbug_pair_projection(weight: int | Fraction) -> dict[str, int | Fraction]:
    """Return exact finite-weight Gaussian and repair readouts on octahedral pairs.

    Gaussian precision uses literal arithmetic site means of centred joint
    fields. Repair is fair seam updating on the binary two-load shell. Lag
    covariances of the pair totals are divided by their common stationary
    variance on the centred site space.
    """
    weight = Fraction(weight)
    if weight < 0:
        raise ValueError("closing conductance must be nonnegative")

    mass = len(jitterbug_rigid_edges()) + len(jitterbug_matching_edges()) * weight
    difference_precision = 4 + 2 * weight
    softening = 8 / difference_precision
    odd_precision = 8 - softening
    even_precision = 12
    return {
        "closing_weight": weight,
        "attempt_mass": mass,
        "difference_precision": difference_precision,
        "pair_mean_odd_precision": odd_precision,
        "pair_mean_even_precision": even_precision,
        "pair_mean_opposite_precision": softening / 2,
        "pair_mean_adjacent_variance": 1 / odd_precision + Fraction(1, even_precision),
        "pair_mean_opposite_variance": 2 / odd_precision,
        "odd_precision_softening": softening,
        "odd_lag_one_retention": 1 - 2 / mass,
        "even_lag_one_retention": 1 - 3 / mass,
        "odd_lag_two_defect": 1 / mass**2,
        "precision_softening_to_lag_defect": 8 * mass**2 / difference_precision,
    }


def _native_pair_lags(
    replication: np.ndarray, weight: Fraction, mass: Fraction
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute exact covariances by enumerating updates on binary states."""
    joints, sites = replication.shape
    states = tuple(
        tuple(int(joint in pair) for joint in range(joints))
        for pair in combinations(range(joints), 2)
    )
    counts = (
        weight.denominator * _attempt_counts(states, jitterbug_rigid_edges())
        + weight.numerator * _attempt_counts(states, jitterbug_matching_edges())
    )
    pair_loads = np.asarray(states, dtype=object) @ replication - Fraction(2, sites)
    normalizer = 2 * mass * weight.denominator
    # Use a transpose for the past step so the two-lag route follows the
    # entire transition matrix even without assuming detailed balance.
    future = counts @ pair_loads / normalizer
    past = counts.T @ pair_loads / normalizer
    covariance_scale = len(states) * Fraction(10, 33)
    return (
        pair_loads.T @ pair_loads / covariance_scale,
        pair_loads.T @ future / covariance_scale,
        past.T @ future / covariance_scale,
    )


def _projection_checks(weight: Fraction) -> dict[str, int | Fraction]:
    row = jitterbug_pair_projection(weight)
    geometry = jitterbug_snapshot(0)
    vertices = np.asarray(geometry["vertices"], dtype=object)
    sites = len(vertices)
    identity = np.eye(sites, dtype=object)
    constant = Fraction(1, sites)
    centered = identity - constant
    vertex_products = vertices @ vertices.T
    odd = vertex_products / Fraction(2)
    even = centered - odd
    octahedral = np.asarray(laplacian_matrix((vertex_products == 0).astype(int).tolist()), dtype=object)

    joint_to_site = geometry["joint_to_vertex"]
    replication = identity[list(joint_to_site)]
    differences = np.zeros(replication.shape, dtype=object)
    for first, second in jitterbug_matching_edges():
        differences[first, joint_to_site[first]] = 1
        differences[second, joint_to_site[second]] = -1

    precision, covariance, _ = _graph_matrices(weight)
    pair_precision = replication.T @ precision @ replication / Fraction(2)
    coupling = differences.T @ precision @ replication / Fraction(2)
    difference_precision = differences.T @ precision @ differences / Fraction(2)

    pair_covariance = replication.T @ covariance @ replication / Fraction(4)
    pair_marginal = np.asarray(_exact_inverse((pair_covariance + constant).tolist()), dtype=object) - constant
    expected_marginal = row["pair_mean_odd_precision"] * odd + row["pair_mean_even_precision"] * even

    lag_zero, lag_one, lag_two = _native_pair_lags(replication, weight, Fraction(row["attempt_mass"]))
    lag_defect = lag_two - lag_one @ lag_one
    scaled_defect = lag_defect * row["attempt_mass"]**2
    pair_softening = 2 * pair_precision - pair_marginal
    return {
        **row,
        "matrix_pair_blocks_match": (
            np.array_equal(pair_precision, octahedral)
            and np.array_equal(difference_precision, row["difference_precision"] * identity)
            and np.array_equal(coupling.T @ coupling, 4 * odd)
        ),
        "matrix_gaussian_marginal_matches": (
            np.array_equal(pair_marginal, expected_marginal)
            and bool(np.all(pair_marginal[vertex_products == -1] == row["pair_mean_opposite_precision"]))
        ),
        "native_pair_lag_matches": (
            np.array_equal(lag_zero, centered)
            and np.array_equal(lag_one, row["odd_lag_one_retention"] * odd + row["even_lag_one_retention"] * even)
            and np.array_equal(lag_defect, row["odd_lag_two_defect"] * odd)
            and np.array_equal(scaled_defect @ scaled_defect, scaled_defect)
        ),
        "matrix_common_coupling_matches": np.array_equal(
            pair_softening, row["precision_softening_to_lag_defect"] * lag_defect
        ),
        "lag_defect_rank": int(np.trace(scaled_defect)),
    }


@cache
def _projection_samples() -> tuple[dict[str, int | Fraction], ...]:
    return tuple(_projection_checks(weight) for weight in (Fraction(0), Fraction(3, 2), Fraction(15, 4)))


def jitterbug_pair_projection_census() -> tuple[dict[str, int | Fraction], ...]:
    """Check the formulas against joint matrices and native repair enumeration."""
    return tuple(row.copy() for row in _projection_samples())
