"""Exact invariance and a separate one-edge witness on a larger root shell."""

from fractions import Fraction as F

from oph_fep_synergetics.core.root_shells import root_shell_lag_census, root_shell_quadratics


def test_root_shell_quadratics_change_order_across_the_24_cell():
    cases = root_shell_quadratics()
    assert tuple(case[0] for case in cases) == (3, 4, 5, 6, 7)
    for dimension, count, edges, frame, precision, variance, matches in cases:
        assert count == 2 * dimension * (dimension - 1)
        assert edges == count * 2 * (dimension - 2)
        assert frame == (dimension - 2, 2)
        assert precision == (2 * dimension * (dimension - 2), 8 * (dimension - 2))
        assert variance == tuple(F(1, value) for value in precision)
        assert matches
    assert cases[1][5] == (F(1, 16), F(1, 16))


def test_one_weighted_24_cell_seam_breaks_single_step_gaussian_covariance_closure():
    values = root_shell_lag_census()
    assert values["unit_coordinate_modes"]
    assert values["unit"]["routes_match"] and values["weighted"]["routes_match"]
    assert values["unit"]["lag_one"] == F(47, 48)
    assert values["unit"]["defect_trace"] == 0
    assert values["weighted"]["lag_one"] == F(91, 96)
    assert values["weighted"]["lag_two"] == F(2819, 3072)
    assert values["weighted"]["rank_one_chord"]
    # The endpoint difference has squared projection 1/6 onto the
    # coordinate space and complementary squared norm 11/6.
    one_edge_formula = F(1, 6) * F(11, 6) * F(1, 4)**2
    assert values["weighted"]["defect_trace"] == one_edge_formula == F(11, 576)
    assert values["weighted"]["lag_two"] - values["weighted"]["lag_one"]**2 == one_edge_formula
