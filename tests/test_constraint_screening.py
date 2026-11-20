"""Exact sum constraints and geometric separators on FCC boundaries."""

import math
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from oph_fep_synergetics.core.constraint_screening import (
    fcc_sum_constraint_census,
    separated_sum_constraint,
)
from oph_fep_synergetics.geometry.lattice import generate_shell
from oph_fep_synergetics.geometry.quadray import _to_xyz_int


def _d3_shell(frequency):
    return [
        point for point in product(range(-frequency, frequency + 1), repeat=3)
        if sum(point) % 2 == 0
        and max(max(map(abs, point)), sum(map(abs, point)) // 2) == frequency
    ]


def test_separating_plane_families_on_larger_fcc_boundaries():
    moves = _d3_shell(1)
    normals = ((1, 0, 0), (0, 1, 0), (0, 0, 1),
               *((1, y, z) for y, z in product((-1, 1), repeat=2)))
    for frequency in range(1, 5):
        sites = set(_d3_shell(frequency))
        assert sites == {
            tuple(value // 2 for value in _to_xyz_int(point)) for point in generate_shell(frequency)
        }
        for normal in normals:
            sides = {point: sum(x * y for x, y in zip(point, normal, strict=True)) for point in sites}
            on_plane = 4 * frequency if sum(map(abs, normal)) == 1 else 6 * frequency
            assert sum(value == 0 for value in sides.values()) == on_plane
            assert sum(value > 0 for value in sides.values()) == (len(sites) - on_plane) // 2
            for point, value in sides.items():
                for move in moves:
                    neighbour = tuple(a + b for a, b in zip(point, move, strict=True))
                    assert value * sides.get(neighbour, 0) >= 0


@pytest.mark.parametrize(("frequency", "plane", "partition", "side_variance", "ratio", "hard"), (
    (1, "axis", (4, 4, 4), F(2), F(9, 5), F(-1, 16)),
    (1, "diagonal", (6, 3, 3), F(3, 2), F(25, 16), F(-1, 12)),
    (2, "axis", (8, 17, 17), F(1123, 48), F(1371241, 110112), F(-867, 35936)),
    (2, "diagonal", (12, 15, 15), F(216, 11), F(51529, 4873), F(-49, 2112)),
))
def test_exact_fcc_witness_against_independent_joint_gaussians(frequency, plane, partition, side_variance, ratio, hard):
    key = f"frequency_{frequency}_{plane}"
    witness = fcc_sum_constraint_census().pop(key)
    assert witness["partition_sizes"] == partition
    assert witness["site_count"] == 10 * frequency**2 + 2
    assert witness["seam_count"] == {1: 24, 2: 96}[frequency]
    assert witness["crossing_seams"] == 0
    assert witness["side_sum_variances"] == (side_variance, side_variance)
    assert witness["soft_information_ratio"] == ratio
    assert witness["centered_antipodal_cross_covariance"] == hard

    points = np.array(_d3_shell(frequency))
    delta = points[:, None, :] - points[None, :, :]
    adjacency = (delta**2).sum(axis=2) == 2
    laplacian = np.diag(adjacency.sum(axis=1)) - adjacency
    values = points[:, 0] if plane == "axis" else points.sum(axis=1)
    positive = tuple(map(int, np.flatnonzero(values > 0)))
    negative = tuple(map(int, np.flatnonzero(values < 0)))
    records = np.flatnonzero(values == 0)
    free = positive + negative
    split = len(positive)
    analysis = separated_sum_constraint(laplacian.tolist(), positive, negative)

    # Condition the full Q^+ model using its covariance on the observed sites.
    # This never uses the row-sum or rank-one update from the implementation.
    projector = np.ones(laplacian.shape) / len(points)
    full_centered = np.linalg.inv(laplacian + projector) - projector
    cross_records = full_centered[np.ix_(free, records)]
    hard_conditioned = full_centered[np.ix_(free, free)] - cross_records @ np.linalg.solve(
        full_centered[np.ix_(records, records)], cross_records.T,
    )
    cross = hard_conditioned[:split, split:]
    assert np.allclose(cross, np.array(analysis.centered_cross_covariance, dtype=float), atol=1e-13)
    assert np.all(cross < 0)
    assert np.linalg.matrix_rank(cross, tol=1e-12) == 1
    assert np.linalg.matrix_rank(hard_conditioned, tol=1e-12) == len(free) - 1
    for marginal in (hard_conditioned[:split, :split], hard_conditioned[split:, split:]):
        assert np.linalg.eigvalsh(marginal).min() > 0
    assert np.allclose(hard_conditioned.sum(axis=1), 0, atol=1e-12)
    first = tuple(points[positive[0]])
    antipode = tuple(-coord for coord in first)
    column = next(col for col, i in enumerate(negative) if tuple(points[i]) == antipode)
    assert witness["witness_site"] == first
    assert cross[0, column] == pytest.approx(float(hard))

    # The noisy sum adds an all-ones matrix to the grounded precision at
    # unit observation variance. Three covariance determinants compute I.
    soft = np.linalg.inv(laplacian[np.ix_(free, free)] + np.ones((len(free), len(free))))
    assert np.allclose(soft[:split, split:], np.array(analysis.soft_cross_covariance, dtype=float), atol=1e-13)
    parts = (soft[:split, :split], soft[split:, split:], soft)
    first, second, joint = (2 * np.log(np.linalg.cholesky(part).diagonal()).sum() for part in parts)
    assert first + second - joint == pytest.approx(math.log(ratio), abs=1e-12)
    witness["partition_sizes"] = ()
    assert fcc_sum_constraint_census()[key]["partition_sizes"] == partition


def test_weighted_path_generalization_and_rejections():
    # Fix the middle of a path with seam conductances 2 and 3/2.
    path = ((2, -2, 0), (-2, F(7, 2), F(-3, 2)), (0, F(-3, 2), F(3, 2)))
    analysis = separated_sum_constraint(path, (0,), (2,), noise_variance=F(5, 6))
    assert analysis.side_sum_variances == (F(1, 2), F(2, 3))
    assert analysis.centered_cross_covariance == ((F(-2, 7),),)
    assert analysis.soft_cross_covariance == ((F(-1, 6),),)
    assert analysis.soft_canonical_correlation_squared == F(1, 6)
    assert analysis.soft_information_ratio == F(6, 5)

    bad_inputs = (
        (((0, 0), (0,)), (0,), (1,), 1, "square"),
        (path, (0,), (0,), 1, "partition"),
        (path, (), (2,), 1, "partition"),
        (path, (0,), (2,), 0, "positive noise"),
        (((1, -1, 0), (1, -1, 0), (0, 0, 0)), (0,), (2,), 1, "symmetric"),
        (path, (0,), (1,), 1, "do not separate"),
        (((0, 0, 0),) * 3, (0,), (2,), 1, "connected"),
    )
    for laplacian, pos, neg, noise, error in bad_inputs:
        with pytest.raises(ValueError, match=error):
            separated_sum_constraint(laplacian, pos, neg, noise_variance=noise)
