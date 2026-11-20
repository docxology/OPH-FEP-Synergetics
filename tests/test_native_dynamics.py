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
