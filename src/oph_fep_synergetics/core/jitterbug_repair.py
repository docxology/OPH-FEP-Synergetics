"""Pair totals and fast within-pair repair at the Jitterbug limit."""

from __future__ import annotations

from fractions import Fraction
from functools import cache
from itertools import combinations
from math import comb, prod

import numpy as np

from ..geometry.jitterbug import jitterbug_matching_edges, jitterbug_rigid_edges, jitterbug_snapshot
from .native_dynamics import _attempt_counts


@cache
def _pair_repair_census() -> dict:
    limit = jitterbug_snapshot(0)
    mapping = limit["joint_to_vertex"]
    bars = jitterbug_rigid_edges()
    matching = jitterbug_matching_edges()
    number_of_joints = len(mapping)
    number_of_sites = limit["vertex_count"]
    tokens = tuple(combinations(range(number_of_joints), 2))
    states = tuple(tuple(int(i in pair) for i in range(number_of_joints)) for pair in tokens)
    loads_by_state = tuple(tuple(sum(mapping[i] == site for i in pair) for site in range(number_of_sites))
                           for pair in tokens)
    loads = tuple(sorted(set(loads_by_state)))
    indices = {row: i for i, row in enumerate(loads)}
    fibres = np.eye(len(loads), dtype=int)[[indices[row] for row in loads_by_state]]
    multiplicities = fibres.sum(axis=0)

    # A bar attempt carries a fair coin. The projected count includes stays.
    # Matching attempts always stay in the current fibre.
    bar_counts = _attempt_counts(states, bars)
    bars_to_fibres = bar_counts @ fibres
    match_counts = _attempt_counts(states, matching)
    matching_conserves = np.array_equal(match_counts @ fibres, 2 * len(matching) * fibres)
    # Check that the matching chain connects every pair of microstates in a
    # fibre. There are at most two occupied dimers to flip independently.
    fast_support = (match_counts > 0).astype(int)
    fast_fibres_connected = np.array_equal(fast_support @ fast_support > 0, fibres @ fibres.T > 0)
    same_loads = [bars_to_fibres[fibres[:, i] == 1] for i in range(len(loads))]
    witness = next((int(rows[:, j].min()), int(rows[:, j].max()))
                   for i, rows in enumerate(same_loads) for j in range(len(loads))
                   if i != j and rows[:, j].min() == 0 and rows[:, j].max() == 2)
    split_to_double = next(tuple(map(int, rows[:, j]))
                           for rows in same_loads for j, target in enumerate(loads)
                           if 2 in target and len(rows) == 4 and sorted(rows[:, j]) == [0, 0, 1, 1])
    # Under uniform microstate stationarity and reversibility, the chance of
    # arriving from the double fibre equals the chance of leaving for it.
    # Conditioning on that arrival biases the split fibre by its exit chance.
    denominator = 2 * (len(bars) + Fraction(3, 2) * len(matching))
    blind_coin_mean = Fraction(sum(split_to_double), len(split_to_double))
    returned_coin_mean = Fraction(sum(n**2 for n in split_to_double), sum(split_to_double))

    totals = fibres.T @ bars_to_fibres - np.diag(2 * len(bars) * multiplicities)
    averaged = np.asarray([[Fraction(int(count), 2 * int(size)) for count in row]
                           for row, size in zip(totals, multiplicities, strict=True)], dtype=object)
    expected = np.full(averaged.shape, Fraction(0), dtype=object)
    laplacian = np.zeros((number_of_sites, number_of_sites), dtype=int)
    directed_edges = tuple(oriented for first, second in limit["bar_edges"]
                           for oriented in ((first, second), (second, first)))
    for i, j in directed_edges:
        laplacian[i, i] += 1
        laplacian[i, j] -= 1
    for index, row in enumerate(loads):
        for i, j in directed_edges:
            rate = Fraction(row[i] * (2 - row[j]), 4)
            if rate:
                destination = tuple(value + (site == j) - (site == i) for site, value in enumerate(row))
                expected[index, indices[destination]] = rate
                expected[index, index] -= rate

    load_matrix = np.asarray(loads, dtype=object)
    stationary = np.asarray([Fraction(int(size), len(states)) for size in multiplicities], dtype=object)
    flux = averaged * stationary[:, None]
    centred = load_matrix - Fraction(2, number_of_sites)
    covariance = centred.T @ (centred * stationary[:, None])
    projected_coeff = covariance[0, 0] - covariance[0, 1]
    return {
        "microstate_count": len(states),
        "pair_load_count": len(loads),
        "fibre_sizes": tuple(sorted(set(map(int, multiplicities)))),
        "positive_rates": tuple(sorted(set(averaged[averaged > 0]))),
        "fibre_sizes_match": all(size == prod(comb(2, n) for n in row)
                                 for row, size in zip(loads, multiplicities, strict=True)),
        "witness_coin_counts": witness,
        "witness_probabilities_half_pose": tuple(n / denominator for n in witness),
        "stationary_return_probabilities_half_pose": (blind_coin_mean / denominator, returned_coin_mean / denominator),
        "stationary_return_bias_factor": returned_coin_mean / blind_coin_mean,
        "repair_reversible": np.array_equal(bar_counts, bar_counts.T) and np.array_equal(match_counts, match_counts.T),
        "matching_conserves_pair_loads": matching_conserves,
        "fast_fibres_connected": fast_fibres_connected,
        "averaged_rates_match": np.array_equal(averaged, expected),
        "detailed_balance": np.array_equal(flux, flux.T),
        "mean_matches_oct_laplacian": np.array_equal(averaged @ load_matrix, -load_matrix @ laplacian * Fraction(1, 2)),
        "centred_count_covariance": projected_coeff,
        "isotropic_centred_count_covariance": np.array_equal(
            covariance, projected_coeff * (np.eye(number_of_sites, dtype=object) - Fraction(1, number_of_sites)),
        ),
    }


def jitterbug_pair_repair_census() -> dict:
    """Check the finite obstruction and the limiting two-capacity generator."""
    return _pair_repair_census().copy()
