"""Loss of Gaussian graph screening under a global sum constraint."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import cache
from typing import Sequence

from ..geometry.lattice import generate_shell
from ..geometry.quadray import _to_xyz_int
from .invariants import _exact_inverse
from .network import laplacian_matrix


@dataclass(frozen=True)
class SumConstraintScreening:
    """Exact covariance and information certificates after fixing the records."""

    side_sum_variances: tuple[Fraction, Fraction]
    centered_cross_covariance: tuple[tuple[Fraction, ...], ...]
    soft_cross_covariance: tuple[tuple[Fraction, ...], ...]
    soft_canonical_correlation_squared: Fraction
    soft_information_ratio: Fraction


def separated_sum_constraint(
    laplacian: Sequence[Sequence[int | Fraction]],
    positive: Sequence[int],
    negative: Sequence[int],
    *,
    noise_variance: int | Fraction = 1,
) -> SumConstraintScreening:
    """Condition a grounded graph Gaussian on a hard or noisy global sum.

    The two input sides have no seam between them. All other vertices are the
    grounded records. The hard sum is the conditional law of the zero-sum
    centred vertex field given those records. For the soft law the grounded
    field's total sum is observed with independent Gaussian noise of the
    given strictly positive variance. Covariances are independent of the
    observed values. The soft information ratio is ``exp(2 * I)`` for the
    mutual information between sides in the conditional law.
    """
    size = len(laplacian)
    if any(len(row) != size for row in laplacian):
        raise ValueError("the Laplacian must be square")
    free = set(positive) | set(negative)
    if (min(len(positive), len(negative)) == 0 or len(free) != len(positive) + len(negative)
            or not free < set(range(size))):
        raise ValueError("the two sides and the grounded records must form a nonempty partition")
    noise = Fraction(noise_variance)
    if noise <= 0:
        raise ValueError("the soft constraint requires positive noise variance")

    precision = tuple(tuple(Fraction(value) for value in row) for row in laplacian)
    if any(sum(row) != 0 or any(
        value != precision[j][i] or (i != j and value > 0)
        for j, value in enumerate(row)
    ) for i, row in enumerate(precision)):
        raise ValueError("expected a symmetric seam Laplacian")
    if any(precision[i][j] for i in positive for j in negative):
        raise ValueError("grounded records do not separate the two sides")

    seen = {0}
    frontier = [0]
    while frontier:
        vertex = frontier.pop()
        reached = {i for i, value in enumerate(precision[vertex]) if value < 0} - seen
        seen.update(reached)
        frontier.extend(reached)
    if len(seen) != size:
        raise ValueError("the original graph must be connected")

    response = []
    for side in (positive, negative):
        covariance = _exact_inverse([[precision[i][j] for j in side] for i in side])
        # Covariance of each site with its own side sum before conditioning.
        response.append(tuple(sum(row, Fraction(0)) for row in covariance))
    first, second = response
    a, b = sum(first), sum(second)
    total = a + b
    cross_numerator = tuple(tuple(-u * v for v in second) for u in first)
    centered_cross = tuple(tuple(value / total for value in row) for row in cross_numerator)
    soft_cross = tuple(tuple(value / (total + noise) for value in row) for row in cross_numerator)
    rho_squared = a * b / ((a + noise) * (b + noise))
    return SumConstraintScreening(
        side_sum_variances=(a, b),
        centered_cross_covariance=centered_cross,
        soft_cross_covariance=soft_cross,
        soft_canonical_correlation_squared=rho_squared,
        soft_information_ratio=1 / (1 - rho_squared),
    )


@cache
def _fcc_sum_constraint_rows() -> dict[str, dict[str, object]]:
    result = {}
    for frequency in (1, 2):
        sites = tuple(sorted(tuple(value // 2 for value in _to_xyz_int(site))
                             for site in generate_shell(frequency)))
        adjacency = tuple(tuple(
            int(sum((a - b)**2 for a, b in zip(point, other, strict=True)) == 2)
            for other in sites
        ) for point in sites)
        laplacian = laplacian_matrix(adjacency)
        for name, values in (("axis", tuple(site[0] for site in sites)),
                             ("diagonal", tuple(map(sum, sites)))):
            positive = tuple(i for i, value in enumerate(values) if value > 0)
            negative = tuple(i for i, value in enumerate(values) if value < 0)
            analysis = separated_sum_constraint(laplacian, positive, negative)
            witness = sites[positive[0]]
            antipode = negative.index(sites.index(tuple(-value for value in witness)))
            result[f"frequency_{frequency}_{name}"] = {
                "site_count": len(sites),
                "seam_count": sum(map(sum, adjacency)) // 2,
                "partition_sizes": (values.count(0), len(positive), len(negative)),
                "crossing_seams": sum(adjacency[i][j] for i in positive for j in negative),
                "side_sum_variances": analysis.side_sum_variances,
                "soft_information_ratio": analysis.soft_information_ratio,
                "soft_canonical_correlation_squared": analysis.soft_canonical_correlation_squared,
                "centered_antipodal_cross_covariance": analysis.centered_cross_covariance[0][antipode],
                "witness_site": witness,
            }
    return result


def fcc_sum_constraint_census() -> dict[str, dict[str, object]]:
    """Compare two separating plane classes on the first two FCC boundaries."""
    return {name: row.copy() for name, row in _fcc_sum_constraint_rows().items()}
