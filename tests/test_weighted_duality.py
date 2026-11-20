"""Reciprocal and copied conductance witnesses on both spherical pairs."""

from fractions import Fraction

from oph_fep_synergetics.core.weighted_duality import weighted_planar_duality_census


def test_reciprocal_conductances_complement_after_normalisation():
    census = weighted_planar_duality_census()
    carrier = census["carrier"]
    ivm = census["ivm"]
    assert (carrier["edge_count"], carrier["entries_checked"]) == (30, 900)
    assert (ivm["edge_count"], ivm["entries_checked"]) == (24, 576)
    assert (carrier["primal_variance"], carrier["reciprocal_dual_variance"]) == (
        Fraction(44, 63), Fraction(19, 63),
    )
    assert (ivm["primal_variance"], ivm["reciprocal_dual_variance"]) == (
        Fraction(44, 57), Fraction(13, 57),
    )
    assert carrier["reciprocal_raw_variance_sum"] == Fraction(29, 21)
    assert ivm["reciprocal_raw_variance_sum"] == Fraction(21, 19)
    assert carrier["weighted_tree_counts"] == (10886400, 2721600)
    assert ivm["weighted_tree_counts"] == (787968, 196992)
    for row in (carrier, ivm):
        assert row["primal_weight"] == 4
        assert row["reciprocal_defect_entries"] == row["rescaled_reciprocal_defect_entries"] == 0
        assert row["rank_one_checks_agree"] is True


def test_copying_the_primal_conductance_to_the_dual_breaks_complementarity():
    census = weighted_planar_duality_census()
    carrier = census["carrier"]
    ivm = census["ivm"]
    assert carrier["copied_dual_variance"] == Fraction(76, 87)
    assert ivm["copied_dual_variance"] == Fraction(52, 63)
    assert carrier["copied_variance_sum"] == Fraction(2872, 1827)
    assert ivm["copied_variance_sum"] == Fraction(1912, 1197)
    assert (carrier["copied_defect_entries"], ivm["copied_defect_entries"]) == (676, 484)
    assert carrier["copied_cut_pairing_nonzeros"] == ivm["copied_cut_pairing_nonzeros"] == 4
