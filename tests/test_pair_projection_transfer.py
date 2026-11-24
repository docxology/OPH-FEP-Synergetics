"""Pair-projection witnesses on small graphs other than the Jitterbug."""

from fractions import Fraction as F

from oph_fep_synergetics.core.pair_projection_transfer import (
    pair_projection_transfer_census,
    switching_pair_census,
)


def test_pair_projection_sign_can_disagree_with_lag_closure():
    rows = pair_projection_transfer_census()
    four, six = rows["path4"], rows["path6"]
    assert four["endpoint_precision"] == F(-2, 3) and four["site_precision_is_laplacian"]
    assert four["lag_two_defect"] == ((F(1, 72), F(-1, 72)), (F(-1, 72), F(1, 72)))
    assert four["lag_two_defect_trace"] == F(1, 36)
    assert six["endpoint_precision"] == F(4, 35) and not six["site_precision_is_laplacian"]
    assert six["literal_site_precision"] == (
        (F(24, 35), F(-4, 5), F(4, 35)),
        (F(-4, 5), F(8, 5), F(-4, 5)),
        (F(4, 35), F(-4, 5), F(24, 35)),
    )
    assert six["lag_two_defect_trace"] == F(3, 200)
    for row in rows.values():
        assert row["gaussian_marginal_matches"] and row["native_lags_match"]
        assert row["balance_formula_matches"]
        assert row["pair_balanced"] == row["static_correction_vanishes"] == (row["lag_two_defect_trace"] == 0)


def test_external_balance_and_count_lumpability_differ():
    rows = pair_projection_transfer_census()
    cycle, complete = rows["cycle4"], rows["complete4"]
    assert cycle["pair_balanced"] and complete["pair_balanced"]
    assert cycle["endpoint_precision"] == -2 and complete["endpoint_precision"] == -4
    assert cycle["return_probabilities"] == (F(1, 8), F(1, 16))
    assert complete["return_probabilities"] == (F(1, 12), F(1, 12))
    for row in rows.values():
        assert row["uniform_cross_blocks"] == row["strongly_lumpable"]
    assert not cycle["strongly_lumpable"] and complete["strongly_lumpable"]


def test_switching_lumpable_pair_records_creates_count_memory_and_lag_defect():
    result = switching_pair_census()
    assert result["adjacent_covariances"] == (0, 0)
    assert result["endpoint_covariance"] == result["lag_two_residual"] == F(4, 9)
    assert result["return_probabilities"] == (F(11, 20), F(5, 24))
    assert result["all_fixed_partitions_lumpable"] and result["all_fixed_partitions_balanced"]
