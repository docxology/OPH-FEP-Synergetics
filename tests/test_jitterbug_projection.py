"""Static marginal and dynamic lag witnesses for the same hidden modes."""

from fractions import Fraction as F

import pytest

from oph_fep_synergetics.core.jitterbug_projection import (
    jitterbug_pair_projection,
    jitterbug_pair_projection_census,
)


def test_half_pose_pair_projection_exact_witness():
    row = jitterbug_pair_projection(F(3, 2))
    assert row["attempt_mass"] == 33 and row["difference_precision"] == 7
    assert (row["pair_mean_odd_precision"], row["pair_mean_even_precision"]) == (F(48, 7), 12)
    assert row["pair_mean_opposite_precision"] == F(4, 7)
    assert row["pair_mean_adjacent_variance"] == F(11, 48)
    assert row["pair_mean_opposite_variance"] == F(7, 24)
    assert row["odd_precision_softening"] == F(8, 7)
    assert (row["odd_lag_one_retention"], row["even_lag_one_retention"]) == (F(31, 33), F(10, 11))
    assert row["odd_lag_two_defect"] == F(1, 1089)
    assert row["precision_softening_to_lag_defect"] == F(8712, 7)


def test_pair_projection_matches_direct_gaussian_and_native_routes():
    rows = jitterbug_pair_projection_census()
    assert tuple(row["closing_weight"] for row in rows) == (0, F(3, 2), F(15, 4))
    assert tuple(row["pair_mean_opposite_precision"] for row in rows) == (1, F(4, 7), F(8, 23))
    assert tuple(row["odd_lag_two_defect"] for row in rows) == (F(1, 576), F(1, 1089), F(4, 8649))
    for row in rows:
        assert row["lag_defect_rank"] == 3
        assert all(row[key] for key in (
            "matrix_pair_blocks_match", "matrix_gaussian_marginal_matches",
            "native_pair_lag_matches", "matrix_common_coupling_matches",
        ))
    rows[0]["lag_defect_rank"] = 0
    assert jitterbug_pair_projection_census()[0]["lag_defect_rank"] == 3


def test_pair_projection_weight_domain_and_large_weight_control():
    with pytest.raises(ValueError, match="must be nonnegative"):
        jitterbug_pair_projection(-1)
    row = jitterbug_pair_projection(10_000)
    assert 0 < row["pair_mean_opposite_precision"] < F(1, 5000)
    assert 0 < 8 - row["pair_mean_odd_precision"] < F(1, 2500)
    assert 0 < row["odd_lag_two_defect"] < F(1, 60_000**2)
