"""The repair centroid has exact Gaussian lag covariance and stationary memory."""

from oph_fep_synergetics.core.centroid_dynamics import centroid_gaussian_census


def test_gaussian_lag_covariance_does_not_make_the_stationary_centroid_markov():
    census = centroid_gaussian_census()
    expected = {
        "carrier": {
            "attempts_per_state": 60,
            "coordinate_eigenvalue": "5 - 1 sqrt 5",
            "retention": "11/12 + 1/60 sqrt 5",
            "retention_parts": ("11/12", "1/60"),
            "one_axis_variance": "50/33 + 10/33 sqrt 5",
            "target_centroid": (0, 2, 0),
            "zero_fibre_jump_counts": (0, 0, 0, 0, 1, 1),
            "next_target_given_zero": "1/180",
            "next_target_given_previous_target_and_zero": "1/60",
        },
        "cuboctahedron": {
            "attempts_per_state": 48,
            "coordinate_eigenvalue": "2",
            "retention": "23/24",
            "retention_parts": ("23/24", "0"),
            "one_axis_variance": "40/33",
            "target_centroid": (0, 1, -1),
            "zero_fibre_jump_counts": (0, 0, 1, 1, 1, 1),
            "next_target_given_zero": "1/72",
            "next_target_given_previous_target_and_zero": "1/48",
        },
    }
    for geometry, values in expected.items():
        result = census[geometry]
        assert {key: result[key] for key in values} == values
        assert result["zero_fibre_size"] == 6
        assert result["graph_matches_reference"] is True
        assert result["uniform_stationary"] is True
        assert result["coordinate_mode_verified"] is True
        assert result["conditional_mean_verified"] is True
        assert result["stationary_variance_verified"] is True
        assert result["witness_matches_edge_moves"] is True
        assert result["witness_matches_reversible_formula"] is True
        assert result["stationary_centroid_non_markov"] is True
