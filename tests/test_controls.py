"""SC3: the obstruction and dynamics controls from 'Agreement and Surprise'.

Every asserted number is a value the source states literally (the rational
affine-path contrast, the Routh first column, the Lyapunov determinant, the
mediator threshold) or a structural invariant the source states qualitatively
(five transition multiplicities, the four-cycle trap, the solenoidal drift).
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
import pytest

from oph_fep_synergetics.core import (
    D_STABILITY_A,
    affine_path_contrast_exact,
    affine_path_contrast_float,
    blanket_dynamics_check,
    contraction_displacement_bounds,
    d_stability_charpoly,
    d_stability_eigenvalues,
    d_stability_expected_charpoly,
    d_stability_free_energy_form,
    d_stability_matrix,
    d_stability_unstable,
    directed_gradient_control,
    equilibrium_propagation_identity,
    four_cycle_local_minimum,
    gated_descent_census,
    logistic_embedding_checks,
    logistic_embedding_jacobian,
    lyapunov_form,
    mediator_expected,
    mediator_pairwise_decomposition,
    mediator_pairwise_margin,
    mediator_schur,
    mismatch_one_coin,
    mismatch_one_excursion,
    mismatch_one_moves,
    path_multiplicity_census,
    path_multiplicity_values,
    reduction_control,
    replicator_potential_check,
    routh_first_column,
    solenoidal_check,
    strict_descent_on_broken_seams,
    surprise_sandwich,
    two_node_map,
    two_node_orbit_min_radius,
    two_node_relaxation_jacobian,
    two_node_relaxation_orbit,
    two_node_sync_orbit,
)

ONE_AT_ZERO = (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


# ---------------------------------------------------------------------------
# 1. Descent on the gated domain, drift at mismatch one
# ---------------------------------------------------------------------------


def test_gated_census_classifies_every_transition():
    census = gated_descent_census()
    assert census == {"strict": 120, "zero": 600, "noop": 3960, "states": 78}
    assert census["strict"] + census["zero"] + census["noop"] == 78 * 30 * 2 == 4680


def test_every_broken_seam_repair_lowers_the_energy_strictly():
    assert strict_descent_on_broken_seams() is True


def test_mismatch_one_offers_two_moves_and_neither_descends():
    waits, swaps = mismatch_one_moves(ONE_AT_ZERO, 0, 1)
    assert waits.current == 1
    assert swaps.current == 0
    assert waits.energy_decrease == 0.0
    assert swaps.energy_decrease == 0.0
    assert waits.after != swaps.after
    assert waits.after == (0, 1) + (0,) * 10


def test_mismatch_one_coin_moves_the_unit_to_the_emptier_port():
    assert mismatch_one_coin(ONE_AT_ZERO, 0, 1) == 0
    assert mismatch_one_coin((0, 1) + (0,) * 10, 0, 1) == 1


def test_mismatch_one_excursion_travels_without_lowering_energy():
    walk = mismatch_one_excursion(steps=6)
    assert len(walk) == 7
    assert len(set(walk)) == 7
    for state in walk:
        assert sum(state) == 1
        assert sum(load * load for load in state) / 2.0 == 0.5
    for before, after in zip(walk, walk[1:]):
        assert before != after


def test_mismatch_one_excursion_needs_a_positive_step_count():
    with pytest.raises(ValueError):
        mismatch_one_excursion(steps=0)


# ---------------------------------------------------------------------------
# 2. The D-stability family
# ---------------------------------------------------------------------------


def test_d_stability_matrix_is_the_scaled_generator():
    assert np.allclose(
        D_STABILITY_A, np.array([[-1.0, 2.0, 0.0], [2.0, -1.0, 3.0], [0.0, -3.0, -1.0]])
    )
    for k in (0.5, 1.0, 2.0, 4.0, 7.0):
        assert np.allclose(d_stability_matrix(k), np.diag([k, 1.0, 1.0]) @ D_STABILITY_A)


def test_characteristic_polynomial_matches_its_closed_form():
    for k in (0.5, 1.0, 2.0, 4.0, 7.0):
        expected = np.array([1.0, k + 2.0, 10.0 - 2.0 * k, 6.0 * k])
        assert np.allclose(d_stability_charpoly(k), expected)
        assert np.allclose(d_stability_charpoly(k), d_stability_expected_charpoly(k))


def test_routh_first_column_matches_the_closed_form():
    assert routh_first_column(4.0) == (1.0, 6.0, -2.0, 24.0)
    for k in (0.5, 1.0, 2.0, 4.0, 7.0):
        expected = (1.0, k + 2.0, (10.0 - 2.0 * k) - 6.0 * k / (k + 2.0), 6.0 * k)
        assert np.allclose(routh_first_column(k), expected)


def test_d_stability_is_unstable_exactly_when_a_root_is_in_the_right_half_plane():
    for k in (0.5, 1.0, 2.0, 4.0, 7.0):
        eigenvalues = d_stability_eigenvalues(k)
        assert len(eigenvalues) == 3
        assert d_stability_unstable(k) == bool(np.any(eigenvalues.real > 0.0))
    assert d_stability_unstable(1.0) is False
    assert d_stability_unstable(4.0) is True
    assert sum(1 for value in d_stability_eigenvalues(4.0) if value.real > 0.0) == 2


def test_logistic_embedding_is_the_quadratic_field_and_its_jacobian():
    assert logistic_embedding_checks() == {"f_zero_norm": 0.0, "jacobian_error": 0.0}
    assert np.allclose(
        logistic_embedding_jacobian(),
        np.array([[-1.0, 2.0, 0.0], [2.0, -1.0, 3.0], [0.0, -3.0, -1.0]]),
    )


# ---------------------------------------------------------------------------
# 3. Solenoidal drift and the stationary blanket
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("q,entropy", [(0.5, 0.5), (1.5, 4.5)])
def test_solenoidal_field_is_stationary_and_produces_entropy(q, entropy):
    values = solenoidal_check(q)
    assert values["stationary_residual"] == 0.0
    assert abs(values["spectral_real"] + 1.0) < 1e-12
    assert abs(abs(values["spectral_imag"]) - q) < 1e-9
    assert abs(values["entropy_production"] - 2.0 * q * q) < 1e-12
    assert abs(values["entropy_production"] - entropy) < 1e-12


def test_stationary_blanket_rows_are_the_declared_drifts():
    values = blanket_dynamics_check(0.75)
    assert values["conditional_independence"] is True
    assert abs(values["mu_drift_external"] - 1.5) < 1e-12
    assert abs(values["mu_drift_blanket"] + 0.25) < 1e-12
    assert abs(values["mu_drift_internal"] + 2.0) < 1e-12
    assert abs(values["cross_drift_difference"] - 1.5) < 1e-12


# ---------------------------------------------------------------------------
# 4. The reduction criterion and the mediator
# ---------------------------------------------------------------------------


def test_reduction_control_reports_the_observation_split():
    values = reduction_control()
    assert values["same_observation"] == 1.0
    assert values["first_after"] == 1.0
    assert values["second_after"] == -1.0
    assert values["response_distance"] == 2.0
    assert values["worst_case_error_lower_bound"] == 1.0


def test_affine_path_contrast_matches_the_source_rational():
    exact = affine_path_contrast_exact(steps=12)
    assert isinstance(exact, Fraction)
    assert exact == Fraction(921633999, 40000000000)
    assert abs(float(exact) - 0.023040849975) < 1e-18
    assert abs(affine_path_contrast_float(steps=12) - float(exact)) < 1e-15
    shorter = affine_path_contrast_exact(steps=4)
    assert isinstance(shorter, Fraction)
    assert abs(float(shorter) - affine_path_contrast_float(steps=4)) < 1e-15


def test_mediator_schur_complement_matches_its_closed_form():
    for epsilon in (0.2, 1.0 / 3.0, 0.5, 1.0):
        assert np.allclose(mediator_schur(epsilon), mediator_expected(epsilon))


def test_mediator_exists_exactly_at_a_third():
    assert abs(mediator_pairwise_margin(1.0 / 3.0)) < 1e-15
    assert mediator_pairwise_margin(0.5) == 0.5
    assert abs(mediator_pairwise_margin(0.2) + 0.4) < 1e-15
    assert mediator_pairwise_decomposition(0.5) is not None
    assert np.allclose(mediator_pairwise_decomposition(0.5), mediator_schur(0.5))
    assert mediator_pairwise_decomposition(0.2) is None


# ---------------------------------------------------------------------------
# 5. The two-node orbit: existence without reachability
# ---------------------------------------------------------------------------


def test_clipped_map_is_the_scaled_rotation_with_a_box():
    assert np.allclose(two_node_map([0.0, 0.0]), [0.0, 0.0])
    assert np.allclose(two_node_map([0.1, 0.2]), [-0.2, 0.6])
    assert np.allclose(two_node_map([0.25, -0.5]), [1.0, -0.5])
    assert np.allclose(two_node_map([1.0, 0.0]), [1.0, 1.0])
    assert np.allclose(two_node_map([4.0, -4.0]), [1.0, 0.0])


def test_synchronous_orbit_never_reaches_the_origin():
    assert two_node_orbit_min_radius(steps=500, skip=100) >= 0.99
    assert two_node_sync_orbit(steps=5).shape == (6, 2)


def test_relaxation_jacobian_has_the_stated_eigenvalues():
    assert np.allclose(two_node_relaxation_jacobian(), np.array([[1.0, -2.0], [2.0, 1.0]]))
    eigenvalues = np.linalg.eigvals(two_node_relaxation_jacobian())
    assert sorted(round(value.real, 9) for value in eigenvalues) == [1.0, 1.0]
    assert sorted(round(abs(value.imag), 9) for value in eigenvalues) == [2.0, 2.0]


def test_relaxation_orbit_has_the_reported_radius_and_period():
    values = two_node_relaxation_orbit(steps=60_000, skip=20_000)
    assert values["min_radius"] > 0.84
    assert values["max_radius"] > values["min_radius"]
    assert abs(values["period"] - 7.58) < 0.1


def test_relaxation_orbit_rejects_invalid_integration_parameters():
    with pytest.raises(ValueError):
        two_node_relaxation_orbit(steps=10, skip=20)


# ---------------------------------------------------------------------------
# 6. Carrier invariants of the native law
# ---------------------------------------------------------------------------


def test_path_multiplicities_are_the_stated_values():
    census = path_multiplicity_census()
    assert sorted(census) == [1, 2, 50, 52]
    assert path_multiplicity_values() == (1, 2, 50, 52)
    assert sum(multiplicity * count for multiplicity, count in census.items()) == 4680


def test_four_cycle_trap_is_a_local_minimum():
    values = four_cycle_local_minimum()
    assert values["compatible_states"] == ((0, 0, 0, 0), (1, 1, 1, 1))
    assert values["compatible_violations"] == (0, 0)
    assert values["trap"] == (0, 0, 1, 1)
    assert values["trap_violations"] == 2
    assert values["flip_violations"] == (2, 2, 2, 2)
    assert values["is_local_minimum"] is True


# ---------------------------------------------------------------------------
# 7. Free-energy forms
# ---------------------------------------------------------------------------


def test_lyapunov_form_matches_the_stated_matrices():
    form = lyapunov_form(D_STABILITY_A)
    assert set(form) == {"P", "Q", "H"}
    expected_p = np.array([[11.0, 4.0, -8.0], [4.0, 5.0, -2.0], [-8.0, -2.0, 9.0]]) / 3.0
    expected_q = np.array([[0.0, 6.0, 4.0], [-6.0, 0.0, 13.0], [-4.0, -13.0, 0.0]]) / 3.0
    assert np.allclose(form["P"], expected_p)
    assert np.allclose(form["Q"], expected_q)
    assert np.allclose(form["P"], form["P"].T)
    assert np.allclose(form["Q"], -form["Q"].T)


def test_d_stability_free_energy_form_checks_out():
    values = d_stability_free_energy_form()
    assert values["P_matches_source"] is True
    assert values["Q_matches_source"] is True
    assert values["Q_is_antisymmetric"] is True
    assert values["scaled_is_unstable"] is True
    assert abs(values["determinant"] - 115.0 / 27.0) < 1e-12
    assert values["residual"] < 1e-12
    assert all(minor > 0.0 for minor in values["leading_minors"])
    assert abs(values["leading_minors"][-1] - values["determinant"]) < 1e-12


def test_directed_gradient_is_not_the_forward_response():
    values = directed_gradient_control()
    assert np.allclose(values["jacobian"], np.array([[-1.0, 2.0], [0.0, -1.0]]))
    assert np.allclose(values["direct"], [1.0, 0.0])
    assert np.allclose(values["adjoint"], [1.0, 2.0])
    assert not np.allclose(values["direct"], values["adjoint"])


def test_equilibrium_propagation_recovers_the_gradient_in_the_limit():
    values = equilibrium_propagation_identity()
    assert values["exact_gradient"] == -1.25
    assert values["closed_form_error"] < 1e-9
    assert abs(values["distance_to_gradient"] - values["expected_distance_to_gradient"]) < 1e-12
    assert 0.0 < values["distance_to_gradient"] < 1e-5


def test_surprise_sandwich_brackets_the_bayesian_surprise(chain_model):
    interior = chain_model.blocks["L_xx"]
    gradient = chain_model.blocks["L_xo"] @ np.array([1.0, 2.0])
    bayesian = 0.5 * float(gradient @ np.linalg.solve(interior, gradient))
    values = surprise_sandwich(interior, gradient, bayesian)
    assert set(values) == {"bayesian_surprise", "lower_bound", "upper_bound", "inside"}
    eigenvalues = np.linalg.eigvalsh(interior)
    norm_squared = float(gradient @ gradient)
    assert abs(values["lower_bound"] - norm_squared / (2.0 * eigenvalues[-1])) < 1e-12
    assert abs(values["upper_bound"] - norm_squared / (2.0 * eigenvalues[0])) < 1e-12
    assert values["lower_bound"] <= bayesian <= values["upper_bound"]
    assert values["inside"] is True


def test_contraction_displacement_bounds_bracket_the_displacement():
    for rate, disagreement in ((0.5, 0.25), (0.25, 0.4), (0.75, 1.0)):
        values = contraction_displacement_bounds(rate, disagreement)
        assert abs(values["lower"] - disagreement / (1.0 + rate)) < 1e-12
        assert abs(values["upper"] - disagreement / (1.0 - rate)) < 1e-12
        assert values["lower"] <= values["upper"]
    at_zero = contraction_displacement_bounds(0.0, 0.3)
    assert abs(at_zero["lower"] - 0.3) < 1e-12
    assert abs(at_zero["upper"] - 0.3) < 1e-12


def test_replicator_field_is_a_mobility_gradient():
    values = replicator_potential_check()
    assert values["field_equals_mobility_gradient"] is True
    assert values["rates_agree"] is True
    assert values["field_norm"] > 0.0
    assert math.isclose(
        values["potential_rate_direct"], values["potential_rate_via_mobility"], rel_tol=1e-12
    )
