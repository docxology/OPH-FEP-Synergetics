"""Checks of quadratic precision in a laboratory frame and a rotated frame."""

from fractions import Fraction
from itertools import combinations, product

import numpy as np

from oph_fep_synergetics.core.quadratic_modes import quadratic_mode_census


def test_quadratic_precision_distinguishes_rotated_variance_on_one_shell():
    frame, precision, variances = ("frame", "precision", "variance")
    census = quadratic_mode_census()
    assert census["carrier"][frame] == (Fraction(8, 5),) * 5
    assert census["carrier"][precision] == (Fraction(48, 5),) * 5
    assert census["fcc"][frame] == (1, 1, 2, 2, 2)
    assert census["fcc"][precision] == (6, 6, 8, 8, 8)
    for name, ratio in (("carrier", 1), ("fcc", Fraction(4, 3))):
        diagonal, _, off_diagonal, *_ = census[name][variances]
        assert diagonal / off_diagonal == ratio

    for name, scale in (("carrier", (1 + np.sqrt(5)) / 2), ("fcc", 1)):
        # Build the geometry independently of the production vertex order.
        vertices = np.array([
            np.roll((0, a, scale * b), axis) for axis in range(3)
            for a, b in product((-1, 1), repeat=2)
        ])
        pairs = list(combinations(vertices, 2))
        edge_squared = min(np.sum((a - b) ** 2) for a, b in pairs)
        seams = [(a, b) for a, b in pairs if np.isclose(np.sum((a - b) ** 2), edge_squared)]
        radius_squared = 1 + scale**2
        # Rotating diag(1,-1,0) by 45 degrees about z produces the xy form.
        def fields(point):
            x, y, _ = point
            return np.array([x**2 - y**2, 2*x*y]) / (radius_squared * np.sqrt(2))

        gradient = np.asarray([fields(a) - fields(b) for a, b in seams])
        covariance = np.linalg.inv(gradient.T @ gradient)
        expected = [float(census[name][variances][index]) for index in (0, 2)]
        np.testing.assert_allclose(covariance, np.diag(expected), atol=1e-12)
