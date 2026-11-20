"""Exact and independently reconstructed covariance remainders on a torus."""

from fractions import Fraction

import numpy as np

from oph_fep_synergetics.core.surface_duality import toroidal_duality_census


def test_toroidal_reciprocal_covariance_leaves_two_harmonic_periods():
    census = toroidal_duality_census()
    assert census["cells"] == (9, 18, 9)
    assert census["genus"] == 1
    assert census["changed_edge"] == (0, 1)
    assert census["changed_weight"] == 4
    assert census["entries_checked"] == 324
    assert census["boundary_of_boundary_zero"] is True

    unit, reciprocal = census["unit"], census["reciprocal"]
    assert unit["uniform_direction_blocks_match"] is True
    assert (unit["primal_diagonal"], unit["dual_diagonal"], unit["remainder_diagonal"]) == (
        Fraction(4, 9), Fraction(4, 9), Fraction(1, 9),
    )
    assert (reciprocal["primal_diagonal"], reciprocal["dual_diagonal"],
            reciprocal["remainder_diagonal"]) == (Fraction(16, 21), Fraction(1, 6), Fraction(1, 14))
    assert unit["period_gram"] == ((1, 0), (0, 1))
    assert reciprocal["period_gram"] == ((Fraction(8, 7), 0), (0, 1))
    assert census["rescaled_reciprocal"] == reciprocal
    for result in (unit, reciprocal):
        assert result["harmonic_periods_match"] is True
        assert result["is_projector"] is True
        assert result["cross_pairing_rank"] == 0
        assert result["remainder_trace"] == result["remainder_rank"] == 2


def test_toroidal_nonreciprocal_defect_has_fixed_trace_and_extra_rank():
    census = toroidal_duality_census()
    copied = census["copied"]
    assert copied["primal_diagonal"] == copied["dual_diagonal"] == Fraction(16, 21)
    assert copied["remainder_diagonal"] == Fraction(-11, 21)
    assert copied["remainder_trace"] == 2
    assert copied["remainder_rank"] == 4
    assert copied["cross_pairing_rank"] == 1
    assert copied["is_projector"] is False

    # Build aligned columns directly from shifts of the periodic coordinates.
    # This avoids the production face-boundary assembly and grounded inverse.
    grid = np.arange(9).reshape(3, 3)
    anchors = np.eye(9)
    east = anchors[np.roll(grid, -1, axis=1).ravel()]
    north = anchors[np.roll(grid, -1, axis=0).ravel()]
    west = anchors[np.roll(grid, 1, axis=1).ravel()]
    south = anchors[np.roll(grid, 1, axis=0).ravel()]
    primal = np.column_stack(((east - anchors).T, (north - anchors).T))
    dual = np.column_stack(((anchors - south).T, (west - anchors).T))
    weights = np.ones(18)
    weights[0] = 2
    for name, left, right in (
        ("unit", primal, dual),
        ("reciprocal", primal * weights, dual / weights),
        ("copied", primal * weights, dual * weights),
    ):
        defect = np.eye(18) - np.linalg.pinv(left) @ left - np.linalg.pinv(right) @ right
        row = census[name]
        np.testing.assert_allclose(np.trace(defect), float(row["remainder_trace"]))
        np.testing.assert_allclose(defect[0, 0], float(row["remainder_diagonal"]))
        assert np.linalg.matrix_rank(defect, tol=1e-10) == row["remainder_rank"]
        spectrum = np.linalg.eigvalsh(defect)
        assert np.count_nonzero(spectrum < -1e-10) == row["cross_pairing_rank"]
        np.testing.assert_allclose(spectrum[-2:], (1.0, 1.0), atol=1e-10)
        if name == "unit":
            np.testing.assert_allclose(defect, np.kron(np.eye(2), np.full((9, 9), 1 / 9)), atol=1e-10)
