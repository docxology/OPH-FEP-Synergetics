"""Exact tests of binary repair quotients beyond the two twelve-port shells."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from functools import cache, lru_cache
from itertools import combinations, product
from typing import Sequence

import numpy as np

from .native_dynamics import _attempt_counts
from .network import Adjacency, graph_distances, pair_resistances
from .repair import native_repair

Pair = tuple[int, int]
Site = tuple[int, int, int]


def _categorical_attempts(
    adjacency: Adjacency, categories: Sequence[int]
) -> tuple[np.ndarray, tuple[int, ...], int, bool]:
    """Aggregate native draws by destination category, with one row per microstate."""
    ports = len(adjacency)
    pairs = tuple(combinations(range(ports), 2))
    seams = tuple(pair for pair in pairs if adjacency[pair[0]][pair[1]])
    states = tuple(tuple(int(port in pair) for port in range(ports)) for pair in pairs)
    counts = _attempt_counts(states, seams)
    membership = np.eye(1 + max(categories), dtype=np.int64)[list(categories)]
    sizes = tuple(membership.sum(axis=0).tolist())
    draws = 2 * len(seams)
    uniform_stationary = bool(np.array_equal(counts, counts.T) and np.all(counts.sum(axis=1) == draws))
    return counts @ membership, sizes, draws, uniform_stationary


def distance_regular_repair_census() -> dict[str, object]:
    """Check the distance-regular criterion on the four-dimensional cube."""
    return _distance_regular_repair_census().copy()


@cache
def _distance_regular_repair_census() -> dict[str, object]:
    dimension = 4
    order = 2**dimension
    hamming = np.asarray([
        [(left ^ right).bit_count() for right in range(order)] for left in range(order)
    ])
    adjacency = (hamming == 1).astype(np.int64)
    layers = tuple((hamming == step).astype(np.int64) for step in range(dimension + 1))
    inward = hamming
    outward = dimension - hamming
    distance_matrix_coefficients_match = all(
        np.array_equal(
            adjacency @ layer,
            inward * (hamming - 1 == step) + outward * (hamming + 1 == step),
        )
        for step, layer in enumerate(layers)
    )

    distances = tuple(hamming[left, right] - 1 for left, right in combinations(range(order), 2))
    projected, sizes, draws, uniform = _categorical_attempts(adjacency, distances)
    quotient = np.zeros((dimension, dimension), dtype=np.int64)
    for step in range(1, dimension + 1):
        if step > 1:
            quotient[step - 1, step - 2] = 2 * step
        if step < dimension:
            quotient[step - 1, step] = 2 * (dimension - step)
    np.fill_diagonal(quotient, draws - quotient.sum(axis=1))

    # Biggs potentials from the independently known Hamming intersection array.
    potentials = [Fraction(order - 1)]
    for step in range(1, dimension):
        potentials.append((step * potentials[-1] - dimension) / (dimension - step))
    biggs = tuple(2 * sum(potentials[:step]) / draws for step in range(1, dimension + 1))
    grounded = pair_resistances(adjacency.tolist())
    return {
        "vertices": order,
        "degree": dimension,
        "class_sizes": sizes,
        "attempts_per_state": draws,
        "attempt_count_matrix": tuple(tuple(row) for row in quotient.tolist()),
        "resistances": tuple(map(str, biggs)),
        "uniform_stationary": uniform,
        "all_rows_match_hamming_moves": bool(np.array_equal(projected, quotient[list(distances)])),
        "distance_matrix_coefficients_match": distance_matrix_coefficients_match,
        "biggs_matches_grounded": all(resistance == biggs[hamming[left, right] - 1]
                                     for left, right, resistance in grounded),
        "resistance_strictly_increases": all(left < right for left, right in zip(biggs, biggs[1:])),
    }


def _stationary_probabilities(successes: Sequence[int], draws: int) -> tuple[Fraction, Fraction]:
    """Compare next-class probabilities under uniform reversible stationarity."""
    total = sum(successes)
    # Reversibility weights the shared fibre by the same probabilities of
    # leaving it when the preceding step was observed in the target class.
    return (
        Fraction(total, draws * len(successes)),
        Fraction(sum(count**2 for count in successes), draws * total),
    )


def bipartite_contrast_census(small_side: int) -> dict[str, object]:
    """Check the analytic resistance collision family K_(m,m+1) at a given m >= 2."""
    return _bipartite_contrast_census(small_side).copy()


@lru_cache(maxsize=8)
def _bipartite_contrast_census(small_side: int) -> dict[str, object]:
    if small_side < 2:
        raise ValueError("the smaller side must contain at least two vertices")
    m = small_side
    n = m + 1
    ports = m + n
    adjacency = tuple(
        tuple(int((left < m) != (right < m)) for right in range(ports))
        for left in range(ports)
    )
    # Pair categories count the tokens on the larger side.
    types = tuple(int(left >= m) + int(right >= m) for left, right in combinations(range(ports), 2))
    projected, sizes, draws, uniform = _categorical_attempts(adjacency, types)
    quotient = (
        (draws - 2 * n, 2 * n, 0),
        (m - 1, draws - (2 * m - 1), m),
        (0, 2 * m, draws - 2 * m),
    )
    resistances = (Fraction(2, n), Fraction(m + n - 1, m * n), Fraction(2, m))
    grounded = tuple(resistance for _, _, resistance in pair_resistances(adjacency))
    shared_successes = tuple(int(row[2]) for row, kind in zip(projected, types, strict=True) if kind < 2)
    stationary, history = _stationary_probabilities(shared_successes, draws)
    return {
        "side_sizes": (m, n),
        "class_sizes": sizes,
        "attempts_per_state": draws,
        "attempt_count_matrix": quotient,
        "resistances": tuple(map(str, resistances)),
        "collision": resistances[0] == resistances[1] != resistances[2],
        "uniform_stationary": uniform,
        "all_rows_match_bipartite_moves": bool(np.array_equal(projected, np.asarray(quotient)[list(types)])),
        "grounded_matches_formula": grounded == tuple(resistances[kind] for kind in types),
        "next_larger_side_probabilities": tuple(
            str(Fraction(int(projected[types.index(kind), 2]), draws)) for kind in (0, 1)
        ),
        "stationary_next_target_given_collision": str(stationary),
        "stationary_next_target_given_previous_target_and_collision": str(history),
        "stationary_probabilities_match_formula": (stationary, history) == (
            Fraction(1, 3 * m + 1), Fraction(1, 2 * n),
        ),
    }


def _fcc_displacement(first: Site, second: Site) -> Site:
    x, y, z = ((a - b) % 4 for a, b in zip(first, second, strict=True))
    return x, y, z


def _fcc_fourier_resistances(
    sites: tuple[Site, ...], steps: tuple[Site, ...]
) -> dict[Site, Fraction]:
    """Diagonalise the periodic FCC graph using the characters of (Z/4Z)^3."""
    cosine = (1, 0, -1, 0)

    def real_character(character: Site, displacement: Site) -> int:
        return cosine[sum(a * b for a, b in zip(character, displacement, strict=True)) % 4]

    characters = tuple(product(range(4), repeat=3))
    eigenvalues = {
        k: len(steps) - sum(real_character(k, step) for step in steps) for k in characters
    }
    # Characters k and k+(2,2,2) agree on the even-parity subgroup. Using both
    # accounts for the factor two in |1 - chi(d)|² = 2(1 - Re chi(d)).
    return {
        displacement: sum(
            Fraction(1 - real_character(k, displacement), eigenvalue)
            for k, eigenvalue in eigenvalues.items() if eigenvalue
        ) / len(sites)
        for displacement in sites
    }


def periodic_fcc_contrast_census() -> dict[str, object]:
    """Find a resistance collision on the degree-twelve periodic FCC lattice mod four."""
    return _periodic_fcc_contrast_census().copy()


@cache
def _periodic_fcc_contrast_census() -> dict[str, object]:
    period = 4
    origin = (0, 0, 0)
    sites = tuple(site for site in product(range(period), repeat=3) if sum(site) % 2 == 0)
    steps = tuple(step for step in product(range(-1, 2), repeat=3) if sum(x * x for x in step) == 2)
    modular_steps = {_fcc_displacement(step, origin) for step in steps}
    adjacency = tuple(
        tuple(int(_fcc_displacement(left, right) in modular_steps) for right in sites)
        for left in sites
    )
    distances = graph_distances(adjacency)
    grounded = {(left, right): value for left, right, value in pair_resistances(adjacency)}
    seams = tuple(pair for pair in grounded if adjacency[pair[0]][pair[1]])
    draws = 2 * len(seams)
    fourier = _fcc_fourier_resistances(sites, steps)
    shared_variance = fourier[(0, 2, 2)]
    target_variance = fourier[(0, 0, 2)]

    def geometric_histogram(displacement: Site) -> dict[Fraction, int]:
        # In displacement space the two tokens supply the same move counts
        # by symmetry under step reversal. A move onto the held token waits.
        targets = Counter(
            fourier[new]
            for step in steps
            if (new := _fcc_displacement(displacement, step)) != origin
        )
        histogram = {variance: 2 * count for variance, count in targets.items()}
        variance = fourier[displacement]
        histogram[variance] = histogram.get(variance, 0) + draws - sum(histogram.values())
        return histogram

    witnesses = []
    native_matches_geometry = True
    for displacement in ((0, 2, 2), (1, 1, 2)):
        endpoint = sites.index(displacement)
        pair: Pair = (sites.index(origin), endpoint)
        state = tuple(int(port in pair) for port in range(len(sites)))
        counts: Counter[Fraction] = Counter()
        for seam in seams:
            for coin in (0, 1):
                moved = native_repair(state, seam, coin).after
                left, right = (i for i, token in enumerate(moved) if token)
                counts[grounded[left, right]] += 1
        native_matches_geometry &= counts == geometric_histogram(displacement)
        common_neighbours = sum(
            left * right for left, right in zip(adjacency[pair[0]], adjacency[pair[1]], strict=True)
        )
        witnesses.append((
            displacement, distances[pair[0]][pair[1]], common_neighbours, str(grounded[pair]),
            str(Fraction(counts[target_variance], draws)),
        ))

    # Translation invariance means listing each directed displacement from
    # the origin yields the stationary collision-fibre mixture. Each
    # unordered pair occurs twice in the directed count.
    shared_successes = tuple(
        geometric_histogram(displacement).get(target_variance, 0)
        for displacement, variance in fourier.items() if variance == shared_variance
    )
    stationary, history = _stationary_probabilities(shared_successes, draws)
    return {
        "period": period,
        "vertices": len(sites),
        "degree": len(modular_steps),
        "attempts_per_state": draws,
        "target_variance": str(target_variance),
        "witnesses": tuple(witnesses),
        "fourier_matches_grounded": all(
            variance == fourier[_fcc_displacement(sites[left], sites[right])]
            for (left, right), variance in grounded.items()
        ),
        "native_matches_geometry": native_matches_geometry,
        "collision_displacement_attempt_counts": tuple(sorted(Counter(shared_successes).items())),
        "stationary_next_target_given_collision": str(stationary),
        "stationary_next_target_given_previous_target_and_collision": str(history),
    }
