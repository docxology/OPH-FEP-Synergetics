"""The rewrite relation and the stochastic attempt kernel answer different questions."""

from __future__ import annotations

import math

import numpy as np
import pytest

from oph_fep_synergetics.core.native_dynamics import (
    binary_attempt_counts,
    binary_ensemble_census,
    native_rewrite_census,
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
