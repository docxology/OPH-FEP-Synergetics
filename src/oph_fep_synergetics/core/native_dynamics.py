"""Finite rewrite and ensemble dynamics of the carrier's total-load-two shell."""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache

import numpy as np

from .carrier import icosahedral_carrier
from .controls import repair_datum
from .repair import energy, is_broken, native_repair, shell_gap

State = tuple[int, ...]


def _reachable(start: State, relation: dict[State, set[State]]) -> set[State]:
    reached = {start}
    todo = [start]
    while todo:
        for successor in relation[todo.pop()] - reached:
            reached.add(successor)
            todo.append(successor)
    return reached


@lru_cache(maxsize=1)
def _rewrite_census() -> dict[str, int | bool]:
    carrier = icosahedral_carrier()
    states = carrier.shell_configurations()
    binary = set(carrier.binary_shell_configurations())
    ungated: dict[State, set[State]] = {state: set() for state in states}
    gated: dict[State, set[State]] = {state: set() for state in states}
    for state in states:
        for seam in carrier.edges:
            for coin in (0, 1):
                target = native_repair(state, seam, coin).after
                if target != state:
                    ungated[state].add(target)
                    if is_broken(state, *seam):
                        gated[state].add(target)
    reach = {state: _reachable(state, ungated) for state in states}
    normal_forms = {state for state in states if not gated[state]}
    gated_counts = {
        len(_reachable(state, gated) & normal_forms) for state in states if state not in binary
    }
    return {
        "ungated_normal_forms": sum(not targets for targets in ungated.values()),
        "ungated_has_cycle": any(state in ungated[target] for state in states for target in ungated[state]),
        "binary_is_closed_connected_class": all(reach[state] == binary for state in binary),
        "all_reach_binary_class": all(binary <= seen for seen in reach.values()),
        "ungated_confluent": all(
            bool(reach[left] & reach[right])
            for seen in reach.values()
            for left in seen
            for right in seen
        ),
        "gated_strict_descent": all(
            energy(target) < energy(state) for state, targets in gated.items() for target in targets
        ),
        "gated_normal_forms": len(normal_forms),
        "gated_normal_forms_are_binary": normal_forms == binary,
        "gated_normal_forms_per_doubled_state": min(gated_counts),
        "gated_normal_form_count_is_uniform": len(gated_counts) == 1,
    }


def native_rewrite_census() -> dict[str, int | bool]:
    """Exhaust the actual state-changing moves with and without the broken-seam gate.

    Reachability includes the state itself. Waits are attempts but are excluded
    from the rewrite relation, so they neither create cycles nor hide normal forms.
    """
    return dict(_rewrite_census())


def binary_attempt_counts() -> tuple[tuple[State, ...], np.ndarray]:
    """Return states and integer transition counts for uniform seam and fair coin."""
    carrier = icosahedral_carrier()
    states = carrier.binary_shell_configurations()
    indices = {state: i for i, state in enumerate(states)}
    counts = np.zeros((len(states), len(states)), dtype=np.int64)
    for i, state in enumerate(states):
        targets = Counter(native_repair(state, seam, coin).after
                          for seam in carrier.edges for coin in (0, 1))
        for target, multiplicity in targets.items():
            counts[i, indices[target]] = multiplicity
    return states, counts


def binary_dirichlet_census() -> dict[str, int | bool | str]:
    """Compare exact coefficients of the native Dirichlet and fibre variance forms.

    Under the uniform binary-shell law, a two-state fibre contributes
    ``(f(x) - f(y))**2 / (2 * len(states))`` to its seam's expected
    conditional variance. A singleton contributes zero. Averaging over seams
    gives both forms the common denominator ``2 * seams * len(states)``.
    """
    states, counts = binary_attempt_counts()
    carrier = icosahedral_carrier()
    draws = 2 * carrier.seam_count
    dirichlet = draws * np.eye(len(states), dtype=np.int64) - counts
    variance = np.zeros_like(dirichlet)
    two_state_fibres = 0
    for seam in carrier.edges:
        fibres: dict[tuple[int, State], list[int]] = defaultdict(list)
        for i, state in enumerate(states):
            fibres[repair_datum(state, seam)].append(i)
        for fibre in fibres.values():
            if len(fibre) == 1:
                continue
            left, right = fibre
            two_state_fibres += 1
            variance[left, left] += 1
            variance[right, right] += 1
            variance[left, right] -= 1
            variance[right, left] -= 1

    adjacent = carrier.edges[0]
    indicator = states.index(tuple(int(i in adjacent) for i in range(carrier.port_count)))
    indicator_energy = Fraction(int(dirichlet[indicator, indicator]), draws * len(states))
    return {
        "binary_states": len(states),
        "seams_checked": carrier.seam_count,
        "two_state_fibres": two_state_fibres,
        "coefficient_matrices_agree": bool(np.array_equal(dirichlet, variance)),
        "adjacent_indicator_energy": str(indicator_energy),
    }


@lru_cache(maxsize=1)
def _binary_ensemble_census() -> dict[str, int | bool]:
    """Check uniform stationarity, the spectral gap and ensemble entropy descent.

    Double stochasticity implies KL monotonicity for every starting law by data
    processing. A point mass gives a separate numerical trajectory check. The
    chi-square bound uses the squared nontrivial eigenvalue; it is not asserted
    to be a sharp rate of KL decay.
    """
    states, counts = binary_attempt_counts()
    carrier = icosahedral_carrier()
    draws = 2 * carrier.seam_count
    transition = counts / draws
    eigenvalues = np.linalg.eigvalsh(transition)
    beta = 1 - shell_gap()
    loads = np.asarray(states)
    mean_step = draws * np.eye(carrier.port_count, dtype=np.int64) - carrier.laplacian()
    p = np.eye(1, len(states))[0]
    previous = np.log(len(states))
    chi_squared = len(states) - 1
    kl_monotone = chi_squared_bounded = True
    for _ in range(32):
        p = p @ transition
        support = p[p > 0]
        kl = float(support @ np.log(support * len(states)))
        chi_squared *= beta**2
        kl_monotone &= kl <= previous + 1e-12
        chi_squared_bounded &= len(states) * float(np.sum((p - 1 / len(states))**2)) <= chi_squared + 1e-12
        previous = kl
    return {
        "binary_chain_states": len(states),
        "binary_chain_draws": draws,
        "binary_chain_symmetric": bool(np.array_equal(counts, counts.T)),
        "binary_chain_doubly_stochastic": bool(
            np.all(counts.sum(axis=1) == draws) and np.all(counts.sum(axis=0) == draws)
        ),
        "binary_chain_stationary_dimension": int(np.count_nonzero(abs(eigenvalues - 1) < 1e-12)),
        "binary_chain_gap_matches_port_gap": bool(abs(eigenvalues[-2] - beta) < 1e-12),
        "binary_chain_nonnegative_spectrum": bool(eigenvalues[0] >= -1e-12),
        "binary_chain_mean_operator_intertwines": bool(np.array_equal(counts @ loads, loads @ mean_step)),
        "binary_chain_sample_kl_monotone": bool(kl_monotone),
        "binary_chain_sample_chi_squared_bound": bool(chi_squared_bounded),
    }


def binary_ensemble_census() -> dict[str, int | bool]:
    """Report finite-chain checks without exposing the cached values for mutation."""
    return dict(_binary_ensemble_census())
