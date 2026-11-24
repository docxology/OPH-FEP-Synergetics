"""Coordinate covariance under changes of the Jitterbug pose."""

from __future__ import annotations

from fractions import Fraction
from functools import cache
from itertools import combinations

import numpy as np

from ..geometry.exact import Surd5
from ..geometry.jitterbug import jitterbug_matching_edges, jitterbug_rigid_edges, jitterbug_vertices
from .jitterbug_inference import ExactNumber, jitterbug_coordinate_inference
from .native_dynamics import _attempt_counts


def jitterbug_motion_defect(
    first_pose: ExactNumber, middle_pose: ExactNumber, last_pose: ExactNumber
) -> Fraction | Surd5:
    """Return the coordinate's three-time Gaussian Markov covariance defect.

    Use one tuned fair repair attempt at each departure pose. The per-axis
    covariance of the unscaled sum of the two occupied joint positions is
    divided by the common odd-mode scale, 20/33, at all times. This keeps
    the formula rational at rational poses without normalising each pose.
    """
    first, middle, last = map(Surd5.of, (first_pose, middle_pose, last_pose))
    if any(t.sign() <= 0 or (t - 1).sign() > 0 for t in (first, middle, last)):
        raise ValueError("poses must lie in the noncoincident interval (0, 1]")

    rho = jitterbug_coordinate_inference(first)["lag_one_retention"]
    weight = (1 - middle * middle) / middle
    mass = len(jitterbug_rigid_edges()) + len(jitterbug_matching_edges()) * weight
    sigma = 1 - (4 + Surd5.of(2) / middle) / (2 * mass)
    defect = rho * sigma * (first - middle) * (last - middle) / (1 + middle * middle)
    return defect.a if defect.is_rational() else defect


def _covariance(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return left.T @ right / Fraction(len(left))


def _in_span(basis: np.ndarray, target: np.ndarray) -> bool:
    # All basis columns have the same norm and are mutually orthogonal.
    norm_squared = basis[:, 0] @ basis[:, 0]
    return np.array_equal(target * norm_squared, basis @ (basis.T @ target))


def _motion_sample(
    poses: tuple[Fraction, Fraction, Fraction],
    states: np.ndarray,
    large: np.ndarray,
    small: np.ndarray,
    bar_counts: np.ndarray,
    matching_counts: np.ndarray,
    scale: Fraction,
) -> dict[str, object]:
    observations = [
        states @ np.concatenate(
            (np.asarray(jitterbug_vertices(t), dtype=object), small - t * large), axis=1
        )
        for t in poses
    ]
    first, middle, last = observations
    advance = []
    stationary = True
    for t in poses[:2]:
        weight = (1 - t * t) / t
        counts = weight.denominator * bar_counts + weight.numerator * matching_counts
        attempts = 2 * (
            len(jitterbug_rigid_edges()) * weight.denominator
            + len(jitterbug_matching_edges()) * weight.numerator
        )
        stationary &= np.array_equal(counts, counts.T) and bool(np.all(counts.sum(axis=1) == attempts))
        advance.append(counts / Fraction(attempts))

    first_step, second_step = advance
    next_middle = first_step @ middle
    next_last = second_step @ last
    covariances = (
        _covariance(first, first_step @ next_last),
        _covariance(first, next_middle),
        _covariance(middle, next_last),
        _covariance(middle, middle),
    )
    cov02, cov01, cov12, cov11 = covariances
    variance = cov11[0, 0]
    coordinate_covariances = tuple(cov[:3, :3] for cov in covariances)
    obs02, obs01, obs12, _ = coordinate_covariances
    defect = obs02 - obs01 @ obs12 / variance
    normalised = defect[0, 0] / scale
    expected = jitterbug_motion_defect(*poses)
    return {
        "poses": poses,
        "coordinate_lag_defect": normalised,
        "raw_coordinate_lag_defect": defect[0, 0],
        "raw_coordinate_covariances": tuple(cov[0, 0] for cov in coordinate_covariances),
        "coordinate_isotropic": all(
            np.array_equal(cov, cov[0, 0] * np.eye(3, dtype=object))
            for cov in coordinate_covariances
        ),
        "formula_matches": np.array_equal(defect, expected * scale * np.eye(3, dtype=object)),
        "augmented_lag_closes": (
            np.array_equal(cov11, variance * np.eye(6, dtype=object))
            and np.array_equal(cov02, cov01 @ cov12 / variance)
        ),
        "augmented_mean_invariant": _in_span(first, next_middle) and _in_span(middle, next_last),
        "uniform_stationary": stationary,
    }


@cache
def _motion_census() -> dict:
    large = np.asarray(jitterbug_vertices(0), dtype=object)
    small = np.asarray(jitterbug_vertices(1), dtype=object) - large
    ports = len(large)
    states = tuple(
        tuple(int(port in occupied) for port in range(ports))
        for occupied in combinations(range(ports), 2)
    )
    bar_counts = _attempt_counts(states, jitterbug_rigid_edges()).astype(object)
    matching_counts = _attempt_counts(states, jitterbug_matching_edges()).astype(object)
    state_matrix = np.asarray(states, dtype=object)
    fixed_basis = state_matrix @ np.concatenate((large, small), axis=1)
    fixed_dimension = fixed_basis.shape[1]
    # Each signed slot has squared norm four. A uniform two-subset has
    # occupancy covariance 2(n-2)/(n(n-1)) on the centred vertex space.
    scale = Fraction(8 * (ports - 2), ports * (ports - 1))
    high, middle, low = Fraction(3, 5), Fraction(1, 2), Fraction(1, 3)
    return {
        "microstate_count": len(states),
        "odd_readout_dimension": fixed_dimension,
        "readout_covariance_scale": scale,
        "centred_fixed_odd_covariance": (
            bool(np.all(fixed_basis.sum(axis=0) == 0))
            and np.array_equal(_covariance(fixed_basis, fixed_basis),
                               scale * np.eye(fixed_dimension, dtype=object))
        ),
        "bar_and_matching_preserve_odd_basis": tuple(
            _in_span(fixed_basis, counts @ fixed_basis)
            for counts in (bar_counts, matching_counts)
        ),
        "samples": {
            name: _motion_sample(poses, state_matrix, large, small, bar_counts, matching_counts, scale)
            for name, poses in (
                ("contracting", (high, middle, low)),
                ("frozen", (middle, middle, middle)),
                ("unchanged_first", (middle, middle, low)),
                ("unchanged_second", (high, middle, middle)),
                ("expanding", (low, middle, high)),
                ("reversing", (high, middle, high)),
            )
        },
    }


def jitterbug_motion_census() -> dict:
    """Check the analytic defect against the exact two-token transition matrices."""
    census = _motion_census()
    return {**census, "samples": {name: row.copy() for name, row in census["samples"].items()}}
