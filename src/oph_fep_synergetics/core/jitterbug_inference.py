"""Gaussian and repair inference on a Jitterbug graph with closing seams."""

from __future__ import annotations

from fractions import Fraction
from functools import cache
from itertools import combinations

import numpy as np

from ..geometry.exact import Surd5
from ..geometry.jitterbug import (
    ICOSAHEDRON_PARAMETER,
    jitterbug_matching_edges,
    jitterbug_rigid_edges,
    jitterbug_vertices,
)
from .invariants import _exact_inverse, laplacian_determinant

ExactNumber = int | Fraction | Surd5


def _normalize(value: ExactNumber) -> Fraction | Surd5:
    exact = Surd5.of(value)
    return exact.a if exact.is_rational() else exact


def weighted_matching_inference(weight: ExactNumber) -> dict[str, object]:
    """Return exact scalar-Gaussian invariants for unit rods and weight ``w``.

    The six added seams have equal nonnegative conductance. Each nonzero
    spectrum row gives ``(center, sign, radical_squared, multiplicity)`` for
    an eigenvalue ``center + sign * sqrt(radical_squared)``. The entropy
    precision ratio is the centred precision determinant divided by the
    rod-only value, so half its natural log is the entropy decrease.
    """
    w = Surd5.of(weight)
    if w.sign() < 0:
        raise ValueError("closing conductance must be nonnegative")

    # The centred even field has two diagonal modes and three mixed modes.
    # In each odd coordinate-sign sector the small and large slots have
    # precision [[4 + 2w, -2], [-2, 4]].
    mixed = 4 + 2 * w
    odd_determinant = 12 + 8 * w
    factor = mixed * odd_determinant
    pseudodeterminant = 36 * factor * factor * factor
    # The derivative of log det' is the sum of six equal closing variances.
    closing_variance = (28 + 16 * w) / factor
    rod_variance = (11 - 6 * w * closing_variance) / 24  # Weighted Foster.
    return {
        "closing_weight": _normalize(w),
        "nonzero_precision_eigenvalues": (
            (6, 0, 0, 2),
            (_normalize(mixed), 0, 0, 3),
            (_normalize(4 + w), -1, _normalize(w * w + 4), 3),
            (_normalize(4 + w), 1, _normalize(w * w + 4), 3),
        ),
        "precision_pseudodeterminant": _normalize(pseudodeterminant),
        "tree_weight": _normalize(pseudodeterminant / 12),
        "closing_pair_variance": _normalize(closing_variance),
        "rod_pair_variance": _normalize(rod_variance),
        "entropy_precision_ratio": _normalize(pseudodeterminant / (36 * 48**3)),
    }


def jitterbug_coordinate_inference(t: ExactNumber, weight: ExactNumber | None = None) -> dict[str, object]:
    """Evaluate stationary coordinate-covariance closure at a fixed pose.

    On the cuboctahedron-to-icosahedron segment, sample seams proportionally
    to conductance and apply fair binary repair. The lag-two defect is the
    common diagonal of the normalised lag-two covariance minus the square
    of lag one. If the closure holds, the geometric coordinate space is
    the lowest positive precision eigenspace, and the unit-rod chord factor
    converts its Gaussian contrast contribution to squared chord length.
    With no weight supplied, choose the closing conductance giving closure.
    """
    t = Surd5.of(t)
    if (t - ICOSAHEDRON_PARAMETER).sign() < 0 or (t - 1).sign() > 0:
        raise ValueError("t must lie on the cuboctahedron-to-icosahedron segment")
    tuned = (1 - t * t) / t
    w = tuned if weight is None else Surd5.of(weight)
    if w.sign() < 0:
        raise ValueError("closing conductance must be nonnegative")
    # Leakage out of the coordinate direction (t, 1) in an odd block.
    leak = (w * t + t * t - 1) / ((1 + t * t) * (24 + 6 * w))
    defect = leak * leak
    closes = defect == 0
    eigenvalue = 4 - 2 * t
    chord_factor = 2 * (1 + t * t) * eigenvalue / (1 + t * t - t)
    return {
        "t": _normalize(t),
        "closing_weight": _normalize(w),
        "tuned_weight": _normalize(tuned),
        "coordinate_covariance_closes": closes,
        "lag_two_defect": _normalize(defect),
        "coordinate_eigenvalue": _normalize(eigenvalue) if closes else None,
        "unit_rod_chord_factor": _normalize(chord_factor) if closes else None,
    }


@cache
def _graph_matrices(weight: Fraction) -> tuple[np.ndarray, np.ndarray, list[list[int]]]:
    """Assemble precision from the geometric edge lists, independent of the spectrum."""
    rods = jitterbug_rigid_edges()
    matching = jitterbug_matching_edges()
    size = len(jitterbug_vertices(1))
    scaled = [[0] * size for _ in range(size)]
    for edges, conductance in ((rods, weight.denominator), (matching, weight.numerator)):
        for i, j in edges:
            scaled[i][i] += conductance
            scaled[j][j] += conductance
            scaled[i][j] -= conductance
            scaled[j][i] -= conductance

    precision = np.asarray(scaled, dtype=object) * (1 / Fraction(weight.denominator))
    projector = Fraction(1, size)
    inverse = _exact_inverse((precision + projector).tolist())
    green = np.asarray(inverse, dtype=object) - projector
    return precision, green, scaled


def _weighted_matrix_checks(weight: Fraction, expected: dict[str, object]) -> dict[str, bool]:
    _, green, scaled = _graph_matrices(weight)
    numerator, denominator = weight.numerator, weight.denominator
    matrix = np.asarray(scaled, dtype=object)
    power = np.eye(len(scaled), dtype=object)
    odd_sum = 8 * denominator + 2 * numerator
    odd_product = (12 * denominator + 8 * numerator) * denominator
    odd_previous, odd_current = 2, odd_sum
    traces_match = True
    # Newton's identities determine the full characteristic polynomial
    # from these exact traces of successive integer matrix powers.
    for exponent in range(1, len(scaled) + 1):
        power = power @ matrix
        target = (2 * (6 * denominator)**exponent
                  + 3 * (4 * denominator + 2 * numerator)**exponent + 3 * odd_current)
        traces_match &= power.trace() == target
        odd_previous, odd_current = odd_current, odd_sum * odd_current - odd_product * odd_previous

    cofactor = laplacian_determinant(scaled)
    unit_tree_weight = Fraction(cofactor, denominator**(len(scaled) - 1))
    return {
        "matrix_spectrum_matches": bool(traces_match),
        "matrix_tree_weight_matches": unit_tree_weight == expected["tree_weight"],
        "matrix_rod_variances_match": all(
            green[i, i] + green[j, j] - 2 * green[i, j] == expected["rod_pair_variance"]
            for i, j in jitterbug_rigid_edges()
        ),
        "matrix_closing_variances_match": all(
            green[i, i] + green[j, j] - 2 * green[i, j] == expected["closing_pair_variance"]
            for i, j in jitterbug_matching_edges()
        ),
    }


def _coordinate_matrix_checks(t: ExactNumber, weight: Fraction, expected: dict[str, object]) -> dict[str, bool]:
    precision, green, _ = _graph_matrices(weight)
    coordinates = np.asarray(jitterbug_vertices(t), dtype=object)
    gram = coordinates.T @ coordinates
    variance = gram[0, 0]
    identity = np.eye(3, dtype=object)
    assert np.array_equal(gram, identity * variance)
    repair = np.eye(len(precision), dtype=object) - precision / (2 * (24 + 6 * weight))
    step = repair @ coordinates
    lag_one = (coordinates.T @ step) / variance
    lag_two = (step.T @ step) / variance
    exact_defect = lag_two - lag_one @ lag_one
    result = {
        "matrix_lag_defect_matches": np.array_equal(exact_defect, identity * expected["lag_two_defect"]),
    }
    if expected["coordinate_covariance_closes"]:
        result["matrix_coordinate_mode_matches"] = np.array_equal(
            precision @ coordinates, coordinates * expected["coordinate_eigenvalue"],
        )
        projection = (coordinates @ coordinates.T) / variance
        low_mode = projection @ green @ projection
        i, j = jitterbug_rigid_edges()[0]
        unit_rod_squared = sum(value * value for value in coordinates[i] - coordinates[j])
        result["matrix_chord_factor_matches"] = all(
            (low_mode[i, i] + low_mode[j, j] - 2 * low_mode[i, j]) * expected["unit_rod_chord_factor"]
            == sum(value * value for value in coordinates[i] - coordinates[j]) / unit_rod_squared
            for i, j in combinations(range(len(precision)), 2)
        )
    return result


@cache
def _inference_samples() -> tuple[tuple[dict[str, object], ...], tuple[dict[str, object], ...]]:
    weighted = []
    for weight in (Fraction(0), Fraction(1, 2), Fraction(1), Fraction(2)):
        result = weighted_matching_inference(weight)
        weighted.append({**result, **_weighted_matrix_checks(weight, result)})
    coordinates = []
    for case, t, weight in (
        ("cuboctahedron", 1, Fraction(0)),
        ("interior-rods", Fraction(3, 4), Fraction(0)),
        ("interior-tuned", Fraction(3, 4), Fraction(7, 12)),
        ("interior-hull", Fraction(3, 4), Fraction(1)),
        ("icosahedral-rods", ICOSAHEDRON_PARAMETER, Fraction(0)),
        ("icosahedron", ICOSAHEDRON_PARAMETER, Fraction(1)),
    ):
        result = jitterbug_coordinate_inference(t, weight)
        coordinates.append({"case": case, **result, **_coordinate_matrix_checks(t, weight, result)})
    return tuple(weighted), tuple(coordinates)


def jitterbug_inference_census() -> dict[str, tuple[dict[str, object], ...]]:
    """Verify formulas against exact full-matrix rational and golden-pose samples."""
    weighted, coordinates = _inference_samples()
    return {
        "weighted_samples": tuple(row.copy() for row in weighted),
        "coordinate_samples": tuple(row.copy() for row in coordinates),
    }
