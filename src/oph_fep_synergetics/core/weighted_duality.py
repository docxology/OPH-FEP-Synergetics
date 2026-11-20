"""Reciprocal conductances and normalised Gaussian edge covariances."""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from functools import cache
from math import lcm

from .invariants import _exact_inverse, laplacian_determinant
from .planar_duality import Edge, Matrix, PlanarDuality, carrier_planar_duality, ivm_planar_duality


def _normalised_covariance(
    edges: tuple[Edge, ...], rank: int, roots: tuple[Fraction, ...],
) -> tuple[Matrix, Fraction]:
    """Build a weighted precision and cofactor, then measure edge contrasts."""
    zero = Fraction(0)
    laplacian = [[zero] * (rank + 1) for _ in range(rank + 1)]
    for (tail, head), root in zip(edges, roots, strict=True):
        conductance = root**2
        laplacian[tail][tail] += conductance
        laplacian[head][head] += conductance
        laplacian[tail][head] -= conductance
        laplacian[head][tail] -= conductance
    inverse = _exact_inverse([row[1:] for row in laplacian[1:]])
    green = [[zero] * (rank + 1)] + [[zero, *row] for row in inverse]
    covariance = tuple(
        tuple(
            left * right * (green[a][c] - green[a][d] - green[b][c] + green[b][d])
            for (c, d), right in zip(edges, roots, strict=True)
        )
        for (a, b), left in zip(edges, roots, strict=True)
    )
    denominator = lcm(*(root.denominator**2 for root in roots))
    integral_laplacian = [[int(denominator * value) for value in row] for row in laplacian]
    cofactor = Fraction(laplacian_determinant(integral_laplacian), denominator**rank)
    return covariance, cofactor


def _rank_one_update(unit: Matrix, root: Fraction) -> Matrix:
    """Change the first edge weight by applying the grounded rank-one identity."""
    change = root**2 - 1
    scale = change / (1 + change * unit[0][0])
    roots = (root,) + (Fraction(1),) * (len(unit) - 1)
    return tuple(
        tuple(roots[i] * roots[j] * (entry - scale * unit[i][0] * unit[0][j])
              for j, entry in enumerate(row))
        for i, row in enumerate(unit)
    )


def _sum_defect(primal: Matrix, dual: Matrix) -> int:
    return sum(left + right != int(i == j)
               for i, (primal_row, dual_row) in enumerate(zip(primal, dual, strict=True))
               for j, (left, right) in enumerate(zip(primal_row, dual_row, strict=True)))


def _nonzero_cut_pairings(
    pair: PlanarDuality, primal: tuple[Fraction, ...], dual: tuple[Fraction, ...],
) -> int:
    """Count nonzero entries of B diag(sqrt(g h)) C transpose directly."""
    entries: dict[Edge, Fraction] = defaultdict(Fraction)
    for (a, b), (c, d), left, right in zip(
        pair.primal_edges, pair.dual_edges, primal, dual, strict=True,
    ):
        gain = left * right
        entries[a, c] += gain
        entries[a, d] -= gain
        entries[b, c] -= gain
        entries[b, d] += gain
    return sum(value != 0 for value in entries.values())


@cache
def _weighted_pair_census(pair: PlanarDuality) -> tuple[tuple[str, object], ...]:
    count = len(pair.primal_edges)
    primal_rank, dual_rank = pair.cut_ranks
    roots = (Fraction(2),) + (Fraction(1),) * (count - 1)
    reciprocal_roots = tuple(1 / root for root in roots)
    rescaled_roots = tuple(3 * root for root in reciprocal_roots)
    primal, primal_trees = _normalised_covariance(pair.primal_edges, primal_rank, roots)
    reciprocal, dual_trees = _normalised_covariance(pair.dual_edges, dual_rank, reciprocal_roots)
    copied, _ = _normalised_covariance(pair.dual_edges, dual_rank, roots)
    rescaled, _ = _normalised_covariance(pair.dual_edges, dual_rank, rescaled_roots)
    weight = roots[0]**2
    return (
        ("edge_count", count),
        ("entries_checked", count * count),
        ("edge", pair.primal_edges[0]),
        ("primal_weight", weight),
        ("primal_variance", primal[0][0]),
        ("reciprocal_dual_variance", reciprocal[0][0]),
        ("copied_dual_variance", copied[0][0]),
        ("copied_variance_sum", primal[0][0] + copied[0][0]),
        ("reciprocal_raw_variance_sum", primal[0][0] / weight + reciprocal[0][0] * weight),
        ("reciprocal_defect_entries", _sum_defect(primal, reciprocal)),
        ("rescaled_reciprocal_defect_entries", _sum_defect(primal, rescaled)),
        ("copied_defect_entries", _sum_defect(primal, copied)),
        ("copied_cut_pairing_nonzeros", _nonzero_cut_pairings(pair, roots, roots)),
        ("rank_one_checks_agree", (
            _rank_one_update(pair.primal_covariance, roots[0]) == primal
            and _rank_one_update(pair.dual_covariance, reciprocal_roots[0]) == reciprocal
            and _rank_one_update(pair.dual_covariance, roots[0]) == copied
        )),
        ("weighted_tree_counts", (primal_trees, dual_trees)),
    )


def weighted_planar_duality_census() -> dict[str, dict[str, object]]:
    """Compare reciprocal, rescaled reciprocal, and copied seam conductances.

    Each graph uses unit conductance except that the first primal seam has
    conductance four. Reciprocal conductances on the dual may be multiplied
    by a common positive factor without changing its normalised covariance.
    """
    return {
        "carrier": dict(_weighted_pair_census(carrier_planar_duality())),
        "ivm": dict(_weighted_pair_census(ivm_planar_duality())),
    }
