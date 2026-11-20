"""The rewrite relation and the stochastic attempt kernel answer different questions."""

from __future__ import annotations

import math

import numpy as np
import pytest

from oph_fep_synergetics.core.native_dynamics import (
    binary_attempt_counts,
    binary_dirichlet_census,
    binary_ensemble_census,
    first_entry_census,
    native_rewrite_census,
    repair_contrast_lumping_census,
    repair_distance_lumping_census,
)


def test_native_and_gated_rewrite_relations_have_opposite_obstructions():
    result = native_rewrite_census()
    assert result["ungated_normal_forms"] == 0
    assert result["ungated_has_cycle"]
    assert result["ungated_confluent"]
    assert result["binary_is_closed_connected_class"]
    assert result["all_reach_binary_class"]
    assert result["gated_strict_descent"]
    assert result["gated_normal_forms"] == 66
    assert result["gated_normal_forms_are_binary"]
    assert result["gated_normal_forms_per_doubled_state"] == 5
    assert result["gated_normal_form_count_is_uniform"]


def test_binary_attempt_kernel_keeps_waits_and_preserves_uniform_law():
    states, counts = binary_attempt_counts()
    assert counts.shape == (66, 66)
    assert np.all(counts.sum(axis=1) == 60)
    assert np.array_equal(counts, counts.T)
    i = states.index((1, 1, *([0] * 10)))
    assert counts[i, i] == 52
    assert sorted(counts[i][counts[i] > 0]) == [1] * 8 + [52]
    spectrum = np.linalg.eigvalsh(counts / 60)
    assert spectrum[-1] == pytest.approx(1.0)
    assert 1 - spectrum[-2] == pytest.approx((5 - math.sqrt(5)) / 60)
    result = binary_ensemble_census()
    assert result["binary_chain_stationary_dimension"] == 1
    assert result["binary_chain_mean_operator_intertwines"]
    assert result["binary_chain_sample_kl_monotone"]
    assert result["binary_chain_sample_chi_squared_bound"]


def test_binary_dirichlet_form_equals_mean_conditional_variance():
    assert binary_dirichlet_census() == {
        "binary_states": 66,
        "seams_checked": 30,
        "two_state_fibres": 300,
        "coefficient_matrices_agree": True,
        "adjacent_indicator_energy": "1/495",
    }


def test_native_first_entry_has_geometric_time_and_incident_pair_law():
    assert first_entry_census() == {
        "doubled_states": 12,
        "binary_states": 66,
        "attempts_per_state": 60,
        "entry_blocks_match": True,
        "survival_per_attempt": "5/6",
        "hazard_per_attempt": "1/6",
        "mean_attempts": "6",
        "variance_attempts": "30",
        "entry_pairs_per_doubled_state": 5,
        "uniform_doubled_entry_pairs": 30,
        "uniform_doubled_entry_pair_mass": "1/30",
        "stationary_pair_mass": "1/66",
        "entry_to_stationary_support_ratio": "11/5",
        "entry_to_stationary_total_variation": "6/11",
        "entry_one_port_occupancy": "1/6",
        "stationary_one_port_occupancy": "1/6",
        "entry_one_port_marginals_match": True,
    }


def test_repair_distance_lumping_is_exact_and_has_a_distinct_spectrum():
    result = repair_distance_lumping_census()
    assert result == {
        "class_sizes": (12, 30, 30, 6),
        "attempts_per_state": 60,
        "attempt_count_matrix": ((50, 10, 0, 0), (0, 56, 4, 0), (0, 4, 54, 2), (0, 0, 10, 50)),
        "all_rows_match": True,
        "class_uniform_lift": True,
        "stationary_class_masses": ("0", "5/11", "5/11", "1/11"),
        "stationarity_verified": True,
        "binary_eigenvalue_sum": "5/3",
        "binary_eigenvalue_product": "311/450",
        "binary_eigenvalue_discriminant": "1/75",
    }
    quotient = np.asarray(result["attempt_count_matrix"]) / result["attempts_per_state"]
    rates = np.sort(np.linalg.eigvals(quotient))
    np.testing.assert_allclose(rates, [(25 - math.sqrt(3)) / 30, 5 / 6, (25 + math.sqrt(3)) / 30, 1])
    assert 1 - rates[-2] > (5 - math.sqrt(5)) / 60


def test_carrier_contrast_variance_recovers_the_binary_distance_quotient():
    carrier = repair_contrast_lumping_census()["carrier"]
    distance_quotient = repair_distance_lumping_census()
    assert carrier["attempts_per_state"] == 60
    assert carrier["classes"] == (
        (1, "4", "11/30", 30),
        (2, "6 + 2 sqrt 5", "7/15", 30),
        (3, "10 + 2 sqrt 5", "1/2", 6),
    )
    assert carrier["attempt_count_matrix"] == tuple(
        row[1:] for row in distance_quotient["attempt_count_matrix"][1:]
    )
    assert carrier["stationary_class_masses"] == ("5/11", "5/11", "1/11")
    assert carrier["distance_lumpable"] is True
    assert carrier["contrast_lumpable"] is True
    assert carrier["chord_partition_matches"] is True
    assert carrier["all_rows_match_token_moves"] is True
    assert carrier["class_uniform_lift"] is True
    assert carrier["stationarity_verified"] is True
    assert carrier["distance_two_adjacency_counts_match"] is True


def test_fcc_repair_has_a_contrast_quotient_but_no_distance_quotient():
    cuboctahedron = repair_contrast_lumping_census()["cuboctahedron"]
    assert cuboctahedron["attempts_per_state"] == 48
    assert cuboctahedron["classes"] == (
        (1, "4", "11/24", 24),
        (2, "8", "7/12", 12),
        (2, "12", "5/8", 24),
        (3, "16", "2/3", 6),
    )
    assert cuboctahedron["attempt_count_matrix"] == (
        (44, 2, 2, 0),
        (4, 40, 4, 0),
        (2, 2, 42, 2),
        (0, 0, 8, 40),
    )
    assert cuboctahedron["stationary_class_masses"] == ("4/11", "2/11", "4/11", "1/11")
    assert cuboctahedron["distance_lumpable"] is False
    assert cuboctahedron["contrast_lumpable"] is True
    assert cuboctahedron["chord_partition_matches"] is True
    assert cuboctahedron["all_rows_match_token_moves"] is True
    assert cuboctahedron["class_uniform_lift"] is True
    assert cuboctahedron["stationarity_verified"] is True
    assert cuboctahedron["distance_two_adjacency_counts_match"] is True
    assert tuple(row[1:] for row in cuboctahedron["distance_two_witnesses"]) == (
        ("8", "7/12", 2, "1/12"),
        ("12", "5/8", 1, "1/24"),
    )
