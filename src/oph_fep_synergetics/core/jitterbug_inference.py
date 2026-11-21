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
    jitterbug_snapshot,
    jitterbug_vertices,
)
from .invariants import _exact_inverse, laplacian_determinant
from .network import laplacian_matrix

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


def _pose_scales(t: Surd5) -> tuple[Surd5, Surd5, Surd5]:
    coordinate_eigenvalue = 4 - 2 * t
    radius_squared = 1 + t * t
    radius_to_rod_squared = radius_squared / (2 * (radius_squared - t))
    return (
        coordinate_eigenvalue,
        4 * radius_to_rod_squared * coordinate_eigenvalue,
        Fraction(20, 33) * radius_to_rod_squared,
    )


def jitterbug_coordinate_inference(t: ExactNumber, weight: ExactNumber | None = None) -> dict[str, object]:
    """Evaluate stationary coordinate-covariance closure at a fixed pose.

    Between the cuboctahedron and the open octahedral endpoint, sample seams
    proportionally to conductance and apply fair binary repair on the two-load
    shell. The lag-one retention uses coordinates normalised by stationary
    variance. The unit-rod variance uses the sum of the two occupied positions.
    The lag-two defect is the
    common diagonal of the normalised lag-two covariance minus the square
    of lag one. If the closure holds, the geometric coordinate space is
    the lowest positive precision eigenspace, and the unit-rod chord factor
    converts its Gaussian contrast contribution to squared chord length.
    With no weight supplied, choose the closing conductance giving closure.
    """
    t = Surd5.of(t)
    if t.sign() <= 0 or (t - 1).sign() > 0:
        raise ValueError("t must lie in the noncoincident interval (0, 1]")
    tuned = (1 - t * t) / t
    w = tuned if weight is None else Surd5.of(weight)
    if w.sign() < 0:
        raise ValueError("closing conductance must be nonnegative")
    # Leakage out of the coordinate direction (t, 1) in an odd block.
    mass = 24 + 6 * w
    leak = (w * t + t * t - 1) / ((1 + t * t) * mass)
    defect = leak * leak
    closes = defect == 0
    eigenvalue, chord_factor, centroid_variance = _pose_scales(t)
    projected_precision = 4 + (2 * w * t * t - 4 * t) / (1 + t * t)
    return {
        "t": _normalize(t),
        "closing_weight": _normalize(w),
        "tuned_weight": _normalize(tuned),
        "coordinate_covariance_closes": closes,
        "lag_one_retention": _normalize(1 - projected_precision / (2 * mass)),
        "unit_rod_twice_centroid_variance": _normalize(centroid_variance),
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
    pairs = tuple(combinations(range(len(precision)), 2))
    centroids = np.asarray([coordinates[i] + coordinates[j] for i, j in pairs], dtype=object)
    i, j = jitterbug_rigid_edges()[0]
    rod_squared = sum(value * value for value in coordinates[i] - coordinates[j])
    direct_stationary_variance = (centroids.T @ centroids) / (len(pairs) * rod_squared)
    result = {
        "matrix_lag_defect_matches": np.array_equal(exact_defect, identity * expected["lag_two_defect"]),
        "matrix_lag_one_matches": np.array_equal(lag_one, identity * expected["lag_one_retention"]),
        "matrix_centroid_variance_matches": np.array_equal(
            direct_stationary_variance, identity * expected["unit_rod_twice_centroid_variance"],
        ),
    }
    if expected["coordinate_covariance_closes"]:
        result["matrix_coordinate_mode_matches"] = np.array_equal(
            precision @ coordinates, coordinates * expected["coordinate_eigenvalue"],
        )
        projection = (coordinates @ coordinates.T) / variance
        low_mode = projection @ green @ projection
        result["matrix_chord_factor_matches"] = all(
            (low_mode[i, i] + low_mode[j, j] - 2 * low_mode[i, j]) * expected["unit_rod_chord_factor"]
            == sum(value * value for value in coordinates[i] - coordinates[j]) / rod_squared
            for i, j in pairs
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


def _octahedral_quotient() -> tuple[dict[str, object], np.ndarray]:
    geometry = jitterbug_snapshot(0)
    vertices = np.asarray(geometry["vertices"], dtype=object)
    size = len(vertices)
    identity = np.eye(size, dtype=object)
    replication = identity[list(geometry["joint_to_vertex"])]
    adjacency = [[0] * size for _ in range(size)]
    for (i, j), count in zip(geometry["bar_edges"], geometry["bar_multiplicities"], strict=True):
        adjacency[i][j] = adjacency[j][i] = count
    precision = np.asarray(laplacian_matrix(adjacency), dtype=object)
    constant = Fraction(1, size)
    centered_covariance = np.asarray(_exact_inverse((precision + constant).tolist()), dtype=object) - constant
    axial = (vertices @ vertices.T) * Fraction(1, 2)
    even = identity - axial - constant
    low_covariance = axial @ centered_covariance @ axial

    rod_precision, _, _ = _graph_matrices(Fraction(0))
    unit_precision, _, _ = _graph_matrices(Fraction(1))
    matching = unit_precision - rod_precision
    pair_means = replication / Fraction(2)
    within_pair = np.eye(len(replication), dtype=object) - replication @ pair_means.T
    contrasts = {}
    for i, j in combinations(range(size), 2):
        chord = sum(value * value for value in vertices[i] - vertices[j]) / Fraction(2)
        contrasts.setdefault(chord, set()).add(tuple(
            matrix[i, i] + matrix[j, j] - 2 * matrix[i, j]
            for matrix in (centered_covariance, low_covariance)
        ))
    (adjacent, axial_adjacent), = contrasts[1]
    (opposite, axial_opposite), = contrasts[2]
    normalized = [[value // 2 for value in row] for row in precision.tolist()]
    return {
        "sites": size,
        "bar_multiplicities": geometry["bar_multiplicities"],
        "pair_mean_spectrum": ((0, 1), (8, 3), (12, 2)),
        "pair_mean_pseudodeterminant": size * laplacian_determinant(precision.tolist()),
        "normalized_sum_pseudodeterminant": size * laplacian_determinant(normalized),
        "adjacent_variance": adjacent,
        "opposite_variance": opposite,
        "adjacent_coordinate_variance": axial_adjacent,
        "opposite_coordinate_variance": axial_opposite,
        "unit_rod_chord_factor": 1 / axial_adjacent,
        "matrix_spectrum_matches": bool(
            np.array_equal(precision @ axial, 8 * axial) and np.array_equal(precision @ even, 12 * even)
            and np.array_equal(precision.sum(axis=1), np.zeros(size, dtype=object))
        ),
        "matrix_hard_quotient_matches": bool(
            np.array_equal(replication.T @ rod_precision @ replication, precision)
            and np.array_equal(matching, 2 * within_pair)
        ),
        "matrix_chord_factor_matches": axial_opposite == 2 * axial_adjacent,
    }, pair_means


def _nearest_only(t: Fraction) -> dict[str, object]:
    edges = jitterbug_snapshot(t)["shortest_pairs"]
    sites = len(jitterbug_vertices(t))
    matching = edges == jitterbug_matching_edges() and len({site for edge in edges for site in edge}) == sites
    # A Laplacian on disjoint two-site edges has one zero mode per edge.
    nullity = len(edges) if matching else None
    return {
        "nearest_only_seams": len(edges),
        "nearest_only_matches_closing_pairs": matching,
        "nearest_only_nullity": nullity,
        "nearest_only_centered_nullity": nullity - 1 if nullity is not None else None,
    }


def _finite_pair_mean_checks(weight: Fraction, pair_means: np.ndarray) -> dict[str, object]:
    _, green, _ = _graph_matrices(weight)
    covariance = pair_means.T @ green @ pair_means
    odd_eigenvalue = 4 - 4 / (4 + 2 * weight)
    adjacent = 1 / (2 * odd_eigenvalue) + Fraction(1, 12)
    opposite = 1 / odd_eigenvalue
    vertices = jitterbug_snapshot(0)["vertices"]
    return {
        "pair_mean_adjacent_variance": adjacent,
        "pair_mean_opposite_variance": opposite,
        "matrix_pair_mean_variances_match": all(
            covariance[i, i] + covariance[j, j] - 2 * covariance[i, j]
            == (adjacent if sum((a - b) * (a - b) for a, b in zip(vertices[i], vertices[j], strict=True)) == 2
                else opposite)
            for i, j in combinations(range(len(vertices)), 2)
        ),
    }


def _tuned_contraction_scales(t: Fraction) -> dict[str, Fraction]:
    # Cancel explicit 1/t poles before evaluating the limits at t=0.
    scaled_weight = 1 - t * t
    scaled_mass = 24 * t + 6 * scaled_weight
    eigenvalue, chord_factor, centroid_variance = _pose_scales(Surd5.of(t))
    mixed = 4 * t + 2 * scaled_weight
    odd = 12 * t + 8 * scaled_weight
    return {
        "retention_loss_over_t": eigenvalue.a / (2 * scaled_mass),
        "attempt_mass_times_t": scaled_mass,
        "mass_clock_covariance_decay_rate": eigenvalue.a / 2,
        "unit_rod_twice_centroid_variance": centroid_variance.a,
        "unit_rod_chord_factor": chord_factor.a,
        "closing_variance_over_t": 1 / mixed + 4 / odd,
        "precision_entropy_ratio_times_t_six": (mixed * odd / 48)**3,
    }


@cache
def _contraction_samples() -> tuple[tuple[dict[str, object], ...], dict[str, object], dict[str, object]]:
    quotient, pair_means = _octahedral_quotient()
    rows = []
    half, third, quarter = Fraction(1, 2), Fraction(1, 3), Fraction(1, 4)
    for case, t, policy in (
        ("half-tuned", half, "tuned"),
        ("half-unit", half, "unit"),
        ("quarter-tuned", quarter, "tuned"),
        ("quarter-unit", quarter, "unit"),
        ("half-inverse-square", half, "inverse-square"),
        ("third-inverse-square", third, "inverse-square"),
    ):
        weights = {"tuned": 1 / t - t, "unit": Fraction(1), "inverse-square": (1 + t * t - t) / (2 * t * t)}
        weight = weights[policy]
        result = jitterbug_coordinate_inference(t, weight)
        rows.append({"case": case, **result, **_coordinate_matrix_checks(t, weight, result),
                     **_nearest_only(t), **_finite_pair_mean_checks(weight, pair_means)})
    limits = {
        **_tuned_contraction_scales(Fraction(0)),
        # At fixed unit conductance the leakage numerator tends to -1 and
        # its denominator to the thirty-seam sampling mass.
        "fixed_unit_lag_two_defect": Fraction(1, 30**2),
    }
    return tuple(rows), quotient, limits


def jitterbug_contraction_census() -> dict[str, object]:
    """Check the noncoincident continuation and the hard octahedral quotient."""
    rows, quotient, limits = _contraction_samples()
    return {
        "continuation_samples": tuple(row.copy() for row in rows),
        "octahedral_quotient": quotient.copy(),
        "tuned_limit": limits.copy(),
    }
