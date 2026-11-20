"""Central-plane record sets and Gaussian hemisphere independence."""

import math
from itertools import product

import numpy as np
import pytest

from oph_fep_synergetics.core.carrier import geometric_carrier
from oph_fep_synergetics.core.geometric_screening import central_plane_blanket_census


def test_exact_central_plane_and_minimum_separator_censuses():
    result = central_plane_blanket_census()
    assert result["oph_isomorphism"] is True
    ico = result["icosahedron"]
    cubo = result["cuboctahedron"]
    assert ico["plane_classes"] == ((4, 4, 4, 2, "10/9", 15),)
    assert cubo["plane_classes"] == ((4, 4, 4, 0, "1", 3), (6, 3, 3, 0, "1", 4))
    assert (ico["vertex_spanned_planes"], ico["screening_planes"]) == (15, 0)
    assert (cubo["vertex_spanned_planes"], cubo["screening_planes"]) == (7, 7)
    assert (ico["minimum_separator_size"], cubo["minimum_separator_size"]) == (5, 4)
    assert ico["minimum_separator_component_counts"] == (((1, 6), 12),)
    assert cubo["minimum_separator_component_counts"] == (((1, 7), 12), ((4, 4), 3))
    assert ico["inverse_matches_seams"] is True
    assert cubo["inverse_matches_seams"] is True
    assert tuple(map(len, ico["axis_plane_partition"])) == (4, 4, 4)
    assert tuple(map(len, cubo["axis_plane_partition"])) == (4, 4, 4)
    assert len(ico["axis_plane_crossing_seams"]) == 2
    assert cubo["axis_plane_crossing_seams"] == ()
    assert ico["axis_plane_precision_determinants"] == (480, 480, 207360)
    assert cubo["axis_plane_precision_determinants"] == (192, 192, 36864)


def test_hemisphere_information_from_the_proper_joint_covariance():
    icosahedron = geometric_carrier()
    ico_points = np.asarray(icosahedron.vertices)
    # Independently enumerate the nearest FCC sites in Cartesian coordinates.
    cubo_points = np.array(
        [(0, a, b) for a, b in product((-1, 1), repeat=2)]
        + [(a, b, 0) for a, b in product((-1, 1), repeat=2)]
        + [(a, 0, b) for a, b in product((-1, 1), repeat=2)]
    )
    differences = cubo_points[:, None, :] - cubo_points[None, :, :]
    adjacent = np.sum(differences**2, axis=2) == 2
    cubo_laplacian = np.diag(adjacent.sum(axis=1)) - adjacent.astype(int)
    for laplacian, plane_values, expected in (
        (icosahedron.laplacian(), ico_points[:, 0], math.log(10 / 9) / 2),
        (cubo_laplacian, cubo_points[:, 0], 0.0),
        (cubo_laplacian, cubo_points.sum(axis=1), 0.0),
    ):
        positive = np.flatnonzero(plane_values > 0)
        negative = np.flatnonzero(plane_values < 0)
        records = np.flatnonzero(plane_values == 0)
        free = np.concatenate((positive, negative))
        # Ground the joint field by adding mass only at the records, then
        # condition via a Schur complement of the full covariance rather
        # than taking a principal precision minor.
        precision = laplacian.copy()
        precision[records, records] += 1
        joint = np.linalg.inv(precision)
        cross = joint[np.ix_(free, records)]
        covariance = joint[np.ix_(free, free)] - cross @ np.linalg.solve(
            joint[np.ix_(records, records)], cross.T
        )
        split = len(positive)
        signs_and_logs = [
            np.linalg.slogdet(matrix)
            for matrix in (covariance[:split, :split], covariance[split:, split:], covariance)
        ]
        assert all(sign == 1 for sign, _ in signs_and_logs)
        first, second, joint_log = (log for _, log in signs_and_logs)
        information = (first + second - joint_log) / 2
        assert information == pytest.approx(expected, abs=1e-12)
        block = covariance[:split, split:]
        if expected:
            assert np.linalg.matrix_rank(block, tol=1e-12) == 2
        else:
            assert np.allclose(block, 0.0, rtol=0.0, atol=1e-12)
