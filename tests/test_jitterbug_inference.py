"""Exact Gaussian and weighted-repair checks on the Jitterbug segment."""

import math
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest

from oph_fep_synergetics.core.jitterbug_inference import (
    jitterbug_contraction_census,
    jitterbug_coordinate_inference,
    jitterbug_inference_census,
    weighted_matching_inference,
)
from oph_fep_synergetics.geometry.jitterbug import ICOSAHEDRON_PARAMETER, jitterbug_vertices


def test_weighted_matching_invariants_from_the_full_laplacian():
    rows = jitterbug_inference_census()["weighted_samples"]
    observed = tuple((row["closing_weight"], row["tree_weight"], row["closing_pair_variance"],
                      row["rod_pair_variance"]) for row in rows)
    assert observed == (
        (0, 331776, F(7, 12), F(11, 24)),
        (F(1, 2), 1536000, F(9, 20), F(193, 480)),
        (1, 5184000, F(11, 30), F(11, 30)),
        (2, 33718272, F(15, 56), F(109, 336)),
    )
    assert rows[2]["entropy_precision_ratio"] == F(125, 8)
    assert rows[1]["nonzero_precision_eigenvalues"] == (
        (6, 0, 0, 2), (5, 0, 0, 3), (F(9, 2), -1, F(17, 4), 3), (F(9, 2), 1, F(17, 4), 3),
    )
    assert all(row[key] for row in rows for key in (
        "matrix_spectrum_matches", "matrix_tree_weight_matches", "matrix_rod_variances_match",
        "matrix_closing_variances_match",
    ))


def test_coordinate_closure_and_chord_factors_by_exact_matrix_projection():
    rows = {row["case"]: row for row in jitterbug_inference_census()["coordinate_samples"]}
    assert tuple(row["lag_two_defect"] for row in rows.values()) == (
        0, F(49, 360000), 0, F(1, 22500), F(1, 2880), 0,
    )
    assert tuple(row["unit_rod_chord_factor"] for row in rows.values()) == (
        8, None, F(125, 13), None, None, 10,
    )
    assert rows["interior-tuned"]["closing_weight"] == F(7, 12)
    assert rows["icosahedron"]["coordinate_eigenvalue"] == 4 - 2 * ICOSAHEDRON_PARAMETER
    assert all(row["matrix_lag_defect_matches"] for row in rows.values())
    for row in rows.values():
        if row["coordinate_covariance_closes"]:
            assert row["matrix_coordinate_mode_matches"] and row["matrix_chord_factor_matches"]
    rows["icosahedron"]["lag_two_defect"] = 1
    assert jitterbug_inference_census()["coordinate_samples"][-1]["lag_two_defect"] == 0


def _pair_transition_counts(pairs, seam_weights):
    """Sum the weighted wait and hop counts of the fair binary kernel."""
    pair_index = {pair: i for i, pair in enumerate(pairs)}
    mass = 2 * sum(seam_weights.values())
    counts = np.eye(len(pairs), dtype=object) * mass
    for state, sites in enumerate(pairs):
        occupied = set(sites)
        for seam, conductance in seam_weights.items():
            if len(occupied.intersection(seam)) == 1:
                target = tuple(sorted(occupied.symmetric_difference(seam)))
                counts[state, pair_index[target]] += conductance
                counts[state, state] -= conductance
    return counts, mass


def test_lag_defect_from_binary_transition_probabilities():
    # Construct the seams from exact endpoint distances. Rebuild the
    # 66-state kernel directly at poses on both sides of the icosahedron.
    endpoint_positions = (jitterbug_vertices(1), jitterbug_vertices(ICOSAHEDRON_PARAMETER))
    pairs = tuple(combinations(range(len(endpoint_positions[0])), 2))
    endpoint_edges = tuple({
        (i, j) for i, j in pairs
        if sum((x - y) * (x - y) for x, y in zip(points[i], points[j], strict=True)) == length
    } for points, length in zip(endpoint_positions, (2, 4 * ICOSAHEDRON_PARAMETER * ICOSAHEDRON_PARAMETER), strict=True))
    rods = endpoint_edges[0]
    closing = endpoint_edges[1] - rods
    assert len(rods) == 24 and len(closing) == 6
    for t, weights in ((F(3, 4), (F(0), F(7, 12), F(1))), (F(1, 2), (F(3, 2), F(1)))):
        positions = np.asarray(jitterbug_vertices(t), dtype=object)
        centroids = np.array([positions[i] + positions[j] for i, j in pairs], dtype=object)
        variance = (centroids.T @ centroids)[0, 0]
        for weight in weights:
            seams = {edge: F(1) for edge in rods} | {edge: weight for edge in closing}
            transition, mass = _pair_transition_counts(pairs, seams)
            assert np.array_equal(transition, transition.T)
            one_step = (transition @ centroids) / mass
            lag_one = (centroids.T @ one_step) / variance
            lag_two = (one_step.T @ one_step) / variance
            expected = jitterbug_coordinate_inference(t, weight)
            identity = np.eye(3, dtype=object)
            assert np.array_equal(lag_one, identity * expected["lag_one_retention"])
            assert np.array_equal(lag_two - lag_one @ lag_one, identity * expected["lag_two_defect"])


def test_exact_parameters_and_weight_domain():
    with pytest.raises(ValueError, match="must be nonnegative"):
        weighted_matching_inference(-1)
    with pytest.raises(ValueError, match="must be nonnegative"):
        jitterbug_coordinate_inference(1, weight=-1)
    for outside in (0, -F(1, 2), F(3, 2)):
        with pytest.raises(ValueError, match="must lie"):
            jitterbug_coordinate_inference(outside)
    other = F(4, 5)
    result = jitterbug_coordinate_inference(other)
    assert result["tuned_weight"] == F(9, 20)
    assert result["unit_rod_chord_factor"] == F(328, 35)
    assert weighted_matching_inference(ICOSAHEDRON_PARAMETER)["closing_weight"] == ICOSAHEDRON_PARAMETER


def test_post_icosahedral_inference_policies_against_exact_matrices():
    rows = {row["case"]: row for row in jitterbug_contraction_census()["continuation_samples"]}
    assert tuple((row["t"], row["closing_weight"], row["lag_two_defect"]) for row in rows.values()) == (
        (F(1, 2), F(3, 2), 0),
        (F(1, 2), 1, F(1, 22500)),
        (F(1, 4), F(15, 4), 0),
        (F(1, 4), 1, F(121, 260100)),
        (F(1, 2), F(3, 2), 0),
        (F(1, 3), F(7, 2), F(1, 32400)),
    )
    assert tuple(rows["half-tuned"][key] for key in (
        "unit_rod_chord_factor", "lag_one_retention", "unit_rod_twice_centroid_variance",
        "pair_mean_adjacent_variance", "pair_mean_opposite_variance",
    )) == (10, F(21, 22), F(50, 99), F(11, 48), F(7, 24))
    assert tuple(rows["quarter-tuned"][key] for key in (
        "unit_rod_chord_factor", "lag_one_retention", "unit_rod_twice_centroid_variance",
    )) == (F(119, 13), F(179, 186), F(170, 429))
    for row in rows.values():
        assert all(row[key] for key in (
            "matrix_lag_defect_matches", "matrix_lag_one_matches", "matrix_centroid_variance_matches",
            "matrix_pair_mean_variances_match", "nearest_only_matches_closing_pairs",
        ))
        if row["coordinate_covariance_closes"]:
            assert row["matrix_chord_factor_matches"] and row["matrix_coordinate_mode_matches"]
        assert (row["nearest_only_seams"], row["nearest_only_nullity"], row["nearest_only_centered_nullity"]) == (6, 6, 5)


def test_hard_gaussian_quotient_against_an_independent_octahedral_graph():
    quotient = jitterbug_contraction_census()["octahedral_quotient"]
    assert quotient["sites"] == 6
    assert quotient["bar_multiplicities"] == (2,) * 12
    assert quotient["pair_mean_spectrum"] == ((0, 1), (8, 3), (12, 2))
    assert (quotient["pair_mean_pseudodeterminant"], quotient["normalized_sum_pseudodeterminant"]) == (73728, 2304)
    assert (quotient["adjacent_variance"], quotient["opposite_variance"]) == (F(5, 24), F(1, 4))
    assert (quotient["adjacent_coordinate_variance"], quotient["opposite_coordinate_variance"]) == (F(1, 8), F(1, 4))
    assert quotient["unit_rod_chord_factor"] == 8
    assert all(quotient[key] for key in (
        "matrix_spectrum_matches", "matrix_hard_quotient_matches", "matrix_chord_factor_matches",
    ))

    axis_vertices = np.concatenate((np.eye(3), -np.eye(3)))
    # Opposite pairs are the only missing edges between different sites.
    adjacency = (axis_vertices @ axis_vertices.T == 0).astype(int)
    precision = 2 * (np.diag(adjacency.sum(axis=1)) - adjacency)
    covariance = np.linalg.pinv(precision, hermitian=True)
    projector = axis_vertices @ axis_vertices.T / 2
    axial = projector @ covariance @ projector
    assert np.linalg.eigvalsh(precision) == pytest.approx([0, 8, 8, 8, 12, 12], abs=1e-12)
    for i, j in combinations(range(len(axis_vertices)), 2):
        edge = bool(adjacency[i, j])
        names = ("adjacent_variance", "adjacent_coordinate_variance") if edge else (
            "opposite_variance", "opposite_coordinate_variance",
        )
        for matrix, name in zip((covariance, axial), names, strict=True):
            assert matrix[i, i] + matrix[j, j] - 2 * matrix[i, j] == pytest.approx(float(quotient[name]))


def test_tuned_limits_preserve_variance_and_require_a_rescaled_clock():
    limits = jitterbug_contraction_census()["tuned_limit"]
    assert limits == {
        "retention_loss_over_t": F(1, 3),
        "attempt_mass_times_t": 6,
        "mass_clock_covariance_decay_rate": 2,
        "unit_rod_twice_centroid_variance": F(10, 33),
        "unit_rod_chord_factor": 8,
        "closing_variance_over_t": 1,
        "precision_entropy_ratio_times_t_six": F(1, 27),
        "fixed_unit_lag_two_defect": F(1, 900),
    }
    t = F(1, 10_000)
    coordinate = jitterbug_coordinate_inference(t)
    weighted = weighted_matching_inference(coordinate["closing_weight"])
    assert float(coordinate["unit_rod_twice_centroid_variance"]) == pytest.approx(10 / 33, rel=2e-4)
    assert float(coordinate["lag_one_retention"]**int(1 / t)) == pytest.approx(math.exp(-1 / 3), rel=2e-4)
    assert float(weighted["closing_pair_variance"] / t) == pytest.approx(1, rel=2e-4)
    assert float(t**6 * weighted["entropy_precision_ratio"]) == pytest.approx(1 / 27, rel=2e-3)
