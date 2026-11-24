"""Exact controls for the Jitterbug moving-coordinate covariance."""

from fractions import Fraction

import pytest

from oph_fep_synergetics.core.jitterbug_motion import jitterbug_motion_census, jitterbug_motion_defect
from oph_fep_synergetics.geometry.exact import Surd5
from oph_fep_synergetics.geometry.jitterbug import ICOSAHEDRON_PARAMETER


def test_contracting_coordinate_covariance_does_not_close():
    census = jitterbug_motion_census()
    assert census["microstate_count"] == 66
    assert census["readout_covariance_scale"] == Fraction(20, 33)
    row = census["samples"]["contracting"]
    assert row["poses"] == (Fraction(3, 5), Fraction(1, 2), Fraction(1, 3))
    assert row["raw_coordinate_covariances"] == (
        Fraction(164575, 248292), Fraction(1885, 2508), Fraction(245, 363), Fraction(25, 33)
    )
    assert row["raw_coordinate_lag_defect"] == Fraction(-841, 124146)
    assert row["coordinate_lag_defect"] == Fraction(-841, 75240)
    assert jitterbug_motion_defect(*row["poses"]) == row["coordinate_lag_defect"]
    assert row["coordinate_isotropic"] and row["formula_matches"] and row["uniform_stationary"]


def test_six_odd_modes_close_through_holds_contraction_and_reversal():
    census = jitterbug_motion_census()
    assert census["odd_readout_dimension"] == 6
    assert census["centred_fixed_odd_covariance"]
    assert census["bar_and_matching_preserve_odd_basis"] == (True, True)
    assert census["samples"]["frozen"]["coordinate_lag_defect"] == 0
    assert census["samples"]["unchanged_first"]["coordinate_lag_defect"] == 0
    assert census["samples"]["unchanged_second"]["coordinate_lag_defect"] == 0
    assert census["samples"]["expanding"]["coordinate_lag_defect"] == Fraction(-667, 59400)
    assert census["samples"]["reversing"]["coordinate_lag_defect"] == Fraction(841, 125400)
    for row in census["samples"].values():
        assert row["uniform_stationary"] and row["formula_matches"]
        assert row["augmented_lag_closes"] and row["augmented_mean_invariant"]


def test_defect_requires_distinct_joint_poses_and_accepts_golden_ratio():
    with pytest.raises(ValueError, match="noncoincident"):
        jitterbug_motion_defect(1, Fraction(1, 2), 0)
    with pytest.raises(ValueError, match="noncoincident"):
        jitterbug_motion_defect(2, Fraction(1, 2), Fraction(1, 3))
    residual = jitterbug_motion_defect(Fraction(3, 4), ICOSAHEDRON_PARAMETER, Fraction(1, 2))
    assert Surd5.of(residual).sign() == -1
