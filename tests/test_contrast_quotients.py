"""Independent exact checks of positive and obstructed repair quotients."""

from oph_fep_synergetics.core.contrast_quotients import (
    bipartite_contrast_census,
    distance_regular_repair_census,
    periodic_fcc_contrast_census,
)


def test_four_cube_distance_classes_give_the_contrast_quotient():
    assert distance_regular_repair_census() == {
        "vertices": 16,
        "degree": 4,
        "class_sizes": (32, 48, 32, 8),
        "attempts_per_state": 64,
        "attempt_count_matrix": ((58, 6, 0, 0), (4, 56, 4, 0), (0, 6, 56, 2), (0, 0, 8, 56)),
        "resistances": ("15/32", "7/12", "61/96", "2/3"),
        "uniform_stationary": True,
        "all_rows_match_hamming_moves": True,
        "distance_matrix_coefficients_match": True,
        "biggs_matches_grounded": True,
        "resistance_strictly_increases": True,
    }


def test_complete_bipartite_family_has_resistance_collisions_with_different_moves():
    smallest = bipartite_contrast_census(2)
    assert smallest == {
        "side_sizes": (2, 3),
        "class_sizes": (1, 6, 3),
        "attempts_per_state": 12,
        "attempt_count_matrix": ((6, 6, 0), (1, 9, 2), (0, 4, 8)),
        "resistances": ("2/3", "2/3", "1"),
        "collision": True,
        "uniform_stationary": True,
        "all_rows_match_bipartite_moves": True,
        "grounded_matches_formula": True,
        "next_larger_side_probabilities": ("0", "1/6"),
        "stationary_next_target_given_collision": "1/7",
        "stationary_next_target_given_previous_target_and_collision": "1/6",
        "stationary_probabilities_match_formula": True,
    }
    for size in (3, 4):
        result = bipartite_contrast_census(size)
        assert result["collision"]
        assert result["uniform_stationary"]
        assert result["all_rows_match_bipartite_moves"]
        assert result["grounded_matches_formula"]
        assert result["stationary_probabilities_match_formula"]
        assert result["next_larger_side_probabilities"] == ("0", f"1/{2 * (size + 1)}")


def test_periodic_fcc_collides_at_the_same_distance_degrees_and_common_neighbour_count():
    result = periodic_fcc_contrast_census()
    assert result == {
        "period": 4,
        "vertices": 32,
        "degree": 12,
        "attempts_per_state": 384,
        "target_variance": "11/64",
        "witnesses": (
            ((0, 2, 2), 2, 4, "17/96", "0"),
            ((1, 1, 2), 2, 4, "17/96", "1/192"),
        ),
        "fourier_matches_grounded": True,
        "native_matches_geometry": True,
        "collision_displacement_attempt_counts": ((0, 3), (2, 12)),
        "stationary_next_target_given_collision": "1/240",
        "stationary_next_target_given_previous_target_and_collision": "1/192",
    }
    result["degree"] = 0
    assert periodic_fcc_contrast_census()["degree"] == 12
