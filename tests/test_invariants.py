"""SC2: spectral and probabilistic invariants of the carrier's Laplacian."""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

from oph_fep_synergetics.core import (
    blanket_sizes,
    gmrf_precision,
    icosahedral_carrier,
    independence_matches_graph,
    independence_pattern,
    marginal_variance,
    markov_blanket,
    pairwise_precision,
    proper_precision,
    reciprocal_spectrum_sum,
    reduced_laplacian_determinant,
    spanning_tree_number,
    tree_number_from_spectrum,
)
from oph_fep_synergetics.core import carrier as carrier_module


def test_spanning_tree_number_is_exact():
    assert spanning_tree_number() == 5184000


def test_reduced_laplacian_determinant_is_the_tree_number():
    assert reduced_laplacian_determinant() == 5184000


def test_tree_number_from_spectrum_matches_the_exact_count():
    assert tree_number_from_spectrum() == pytest.approx(5184000, rel=1e-6)


def test_reciprocal_spectrum_sum_is_seven_thirds():
    assert reciprocal_spectrum_sum() == Fraction(7, 3)


def test_marginal_variance_is_seven_thirty_sixths():
    assert marginal_variance() == Fraction(7, 36)


def test_marginal_variance_matches_the_pseudoinverse():
    carrier = icosahedral_carrier()
    pseudoinverse = np.linalg.pinv(carrier.laplacian())
    assert float(marginal_variance()) == pytest.approx(float(pseudoinverse[0][0]), abs=1e-9)
    assert float(np.trace(pseudoinverse)) == pytest.approx(7.0 / 3.0, abs=1e-9)


def test_gmrf_precision_is_the_carrier_laplacian():
    carrier = icosahedral_carrier()
    assert np.allclose(gmrf_precision(carrier), carrier.laplacian())


def test_proper_precision_shifts_only_the_diagonal():
    carrier = icosahedral_carrier()
    laplacian = carrier.laplacian()
    precision = proper_precision(carrier)
    assert np.allclose(precision, laplacian + (1.0 / 12.0) * np.eye(12), atol=1e-12)
    off_diagonal = ~np.eye(12, dtype=bool)
    assert np.allclose(precision[off_diagonal], laplacian[off_diagonal], atol=1e-12)
    assert np.all(np.linalg.eigvalsh(precision) > 0.0)


def test_independence_matches_the_graph():
    assert independence_matches_graph() is True


def test_independence_pattern_is_symmetric_with_sixty_seams():
    pattern = independence_pattern()
    assert len(pattern) == 12
    assert all(len(row) == 12 for row in pattern)
    assert all(pattern[i][j] == pattern[j][i] for i in range(12) for j in range(12))
    assert all(pattern[i][i] is False for i in range(12))
    ordered = sum(1 for i in range(12) for j in range(12) if i != j and pattern[i][j])
    assert ordered == 60


def test_pairwise_precision_is_minus_one_on_every_seam():
    precision = pairwise_precision()
    assert len(precision) == 60
    assert all(value == -1 for value in precision.values())


def test_blanket_sizes_are_five_everywhere():
    assert blanket_sizes() == (5,) * 12


def test_markov_blanket_is_the_neighbourhood():
    carrier = icosahedral_carrier()
    for port in range(12):
        assert markov_blanket(carrier, port) == tuple(carrier.neighbours(port))


def test_carrier_module_is_ophs_labelled_carrier():
    carrier = icosahedral_carrier()
    assert [carrier.neighbours(port) for port in range(12)] == list(carrier_module.OPH_NEIGHBOURS)
