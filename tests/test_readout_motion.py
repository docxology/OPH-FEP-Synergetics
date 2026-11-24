"""Exact moving-readout controls on a square with changing weights."""

from fractions import Fraction

from oph_fep_synergetics.core.readout_motion import moving_readout_census


def test_moving_readouts_admit_either_residual_sign_and_orthogonal_cancellation():
    result = moving_readout_census()
    assert result["state_count"] == 6
    assert result["covariance_scale"] == Fraction(4, 3)
    assert result["defects"] == (Fraction(5, 9), Fraction(-5, 9), 0)
    assert result["hidden_component_norms"] == (
        (Fraction(25, 36), Fraction(4, 9)),
        (Fraction(25, 36), Fraction(4, 9)),
        (Fraction(25, 36), Fraction(1, 4)),
    )
    assert result["native_matches_transport"] and result["split_factorisation_matches"]
    assert result["uniform_stationary"]


def test_consecutive_triples_miss_a_longer_gap_obstruction():
    result = moving_readout_census()
    assert result["adjacent_covariances"] == (Fraction(2, 3), Fraction(5, 6), Fraction(2, 3))
    assert result["adjacent_defects"] == (0, 0)
    assert result["endpoint_covariance"] == Fraction(5, 6)
    assert result["endpoint_defects"] == (Fraction(25, 54), Fraction(25, 54))


def test_transport_and_common_invariance_give_two_closure_controls():
    result = moving_readout_census()
    assert result["changing_basis_closes"]
    assert result["transported_moving_closes"]
    assert result["minimum_common_dimensions"] == (2, 3)
