"""Exact Gaussian and weighted-repair checks on the Jitterbug segment."""

from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest

from oph_fep_synergetics.core.jitterbug_inference import (
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


def test_lag_defect_from_binary_transition_probabilities():
    # Construct the seams from exact distances at the two endpoints, and
    # move one particle across a selected occupied-to-empty seam with the
    # fair coin. This builds the 66-state transition kernel directly.
    t = F(3, 4)
    positions = np.asarray(jitterbug_vertices(t), dtype=object)
    endpoint_positions = (jitterbug_vertices(1), jitterbug_vertices(ICOSAHEDRON_PARAMETER))
    pairs = tuple(combinations(range(len(positions)), 2))
    pair_index = {pair: i for i, pair in enumerate(pairs)}
    endpoint_edges = tuple({
        (i, j) for i, j in pairs
        if sum((x - y) * (x - y) for x, y in zip(points[i], points[j], strict=True)) == length
    } for points, length in zip(endpoint_positions, (2, 4 * ICOSAHEDRON_PARAMETER * ICOSAHEDRON_PARAMETER), strict=True))
    rods = endpoint_edges[0]
    closing = endpoint_edges[1] - rods
    assert len(rods) == 24 and len(closing) == 6
    centroids = np.array([positions[i] + positions[j] for i, j in pairs], dtype=object)
    variance = (centroids.T @ centroids)[0, 0]
    for weight in (F(0), F(7, 12), F(1)):
        seams = {edge: F(1) for edge in rods} | {edge: weight for edge in closing}
        mass = 2 * sum(seams.values())
        transition = np.eye(len(pairs), dtype=object) * mass
        for state, sites in enumerate(pairs):
            occupied = set(sites)
            for seam, conductance in seams.items():
                if len(occupied.intersection(seam)) == 1:
                    target = tuple(sorted(occupied.symmetric_difference(seam)))
                    transition[state, pair_index[target]] += conductance
                    transition[state, state] -= conductance
        assert np.array_equal(transition, transition.T)
        one_step = (transition @ centroids) / mass
        lag_one = (centroids.T @ one_step) / variance
        lag_two = (one_step.T @ one_step) / variance
        expected = jitterbug_coordinate_inference(t, weight)["lag_two_defect"]
        assert np.array_equal(lag_two - lag_one @ lag_one, np.eye(3, dtype=object) * expected)


def test_exact_parameters_and_weight_domain():
    with pytest.raises(ValueError, match="must be nonnegative"):
        weighted_matching_inference(-1)
    with pytest.raises(ValueError, match="must be nonnegative"):
        jitterbug_coordinate_inference(1, weight=-1)
    for outside in (F(1, 2), F(3, 2)):
        with pytest.raises(ValueError, match="must lie"):
            jitterbug_coordinate_inference(outside)
    other = F(4, 5)
    result = jitterbug_coordinate_inference(other)
    assert result["tuned_weight"] == F(9, 20)
    assert result["unit_rod_chord_factor"] == F(328, 35)
    assert weighted_matching_inference(ICOSAHEDRON_PARAMETER)["closing_weight"] == ICOSAHEDRON_PARAMETER
