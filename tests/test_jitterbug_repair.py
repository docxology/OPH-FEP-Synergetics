"""Exact witnesses for the collapsed binary repair model."""

from fractions import Fraction

from oph_fep_synergetics.core.jitterbug_repair import jitterbug_pair_repair_census


def test_collapsed_counts_have_memory_at_finite_weight():
    row = jitterbug_pair_repair_census()
    assert row["matching_conserves_pair_loads"] and row["repair_reversible"]
    assert row["witness_coin_counts"] == (0, 2)
    assert row["witness_probabilities_half_pose"] == (0, Fraction(1, 33))
    assert row["stationary_return_probabilities_half_pose"] == (Fraction(1, 132), Fraction(1, 66))
    assert row["stationary_return_bias_factor"] == 2


def test_fast_fibres_and_averaged_repair_keep_the_octahedral_mean():
    row = jitterbug_pair_repair_census()
    assert row["microstate_count"] == 66 and row["pair_load_count"] == 21
    assert row["fibre_sizes"] == (1, 4) and row["fibre_sizes_match"]
    assert row["positive_rates"] == (Fraction(1, 4), Fraction(1, 2), Fraction(1))
    assert row["fast_fibres_connected"] and row["averaged_rates_match"] and row["detailed_balance"]
    assert row["mean_matches_oct_laplacian"] and row["isotropic_centred_count_covariance"]
    assert row["centred_count_covariance"] == Fraction(10, 33)
