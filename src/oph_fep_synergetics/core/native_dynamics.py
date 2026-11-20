"""Finite rewrite and ensemble dynamics of the carrier's total-load-two shell."""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from typing import Callable

import numpy as np

from ..geometry.exact import Surd5, exact_icosahedron_vertices
from ..geometry.lattice import cuboctahedron_adjacency, cuboctahedron_neighbors
from ..geometry.quadray import squared_distance
from .carrier import icosahedral_carrier
from .controls import repair_datum
from .network import (
    Adjacency,
    graph_distances,
    graphs_are_isomorphic,
    intersection_array,
    pair_resistances,
)
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


def _attempt_counts(states: tuple[State, ...], seams: tuple[tuple[int, int], ...]) -> np.ndarray:
    indices = {state: i for i, state in enumerate(states)}
    counts = np.zeros((len(states), len(states)), dtype=np.int64)
    for i, state in enumerate(states):
        targets = Counter(native_repair(state, seam, coin).after
                          for seam in seams for coin in (0, 1))
        for target, multiplicity in targets.items():
            counts[i, indices[target]] = multiplicity
    return counts


def binary_attempt_counts() -> tuple[tuple[State, ...], np.ndarray]:
    """Return states and integer transition counts for uniform seam and fair coin."""
    carrier = icosahedral_carrier()
    states = carrier.binary_shell_configurations()
    return states, _attempt_counts(states, carrier.edges)


@lru_cache(maxsize=1)
def _full_shell_attempt_counts() -> tuple[tuple[State, ...], np.ndarray]:
    carrier = icosahedral_carrier()
    states = carrier.shell_configurations()
    return states, _attempt_counts(states, carrier.edges)


def first_entry_census() -> dict[str, int | bool | str]:
    """Check the transient block and first binary entry against incidence.

    A doubled state waits in place until one of its incident seams is drawn.
    That draw yields a binary adjacent pair for either coin. The geometric
    waiting time is independent of the entered pair for each initial port.
    """
    carrier = icosahedral_carrier()
    draws = 2 * carrier.seam_count
    states, counts = _full_shell_attempt_counts()
    doubled = [i for i, state in enumerate(states) if 2 in state]
    binary = [i for i, state in enumerate(states) if 2 not in state]
    ports = [states[i].index(2) for i in doubled]
    occupied = [frozenset(i for i, load in enumerate(states[index]) if load) for index in binary]
    seams = {frozenset(edge) for edge in carrier.edges}

    transient = counts[np.ix_(doubled, doubled)]
    entrance = counts[np.ix_(doubled, binary)]
    exit_counts = entrance.sum(axis=1)
    expected_waits = np.diag([draws - 2 * carrier.degree(port) for port in ports])
    expected_entrance = np.asarray(
        [[2 * int(port in pair and pair in seams) for pair in occupied] for port in ports],
        dtype=np.int64,
    )
    blocks_match = (
        np.array_equal(transient, expected_waits)
        and np.array_equal(entrance, expected_entrance)
        and not np.any(counts[np.ix_(binary, doubled)])
        and np.all(exit_counts == exit_counts[0])
    )

    hazard = Fraction(int(exit_counts[0]), draws)
    survival = Fraction(int(transient[0, 0]), draws)
    # Normalize each exit row before averaging over uniform initial doubled ports.
    exit_laws = [
        tuple(Fraction(int(mass), int(total)) for mass in row)
        for row, total in zip(entrance, exit_counts, strict=True)
    ]
    entry = tuple(sum(column) / len(doubled) for column in zip(*exit_laws, strict=True))
    uniform_pair_mass = Fraction(1, len(binary))
    support = [mass for mass in entry if mass]
    entry_marginals = tuple(
        sum(mass for pair, mass in zip(occupied, entry, strict=True) if port in pair)
        for port in ports
    )
    stationary_marginals = tuple(
        sum(port in pair for pair in occupied) * uniform_pair_mass for port in ports
    )
    return {
        "doubled_states": len(doubled),
        "binary_states": len(binary),
        "attempts_per_state": draws,
        "entry_blocks_match": bool(blocks_match),
        "survival_per_attempt": str(survival),
        "hazard_per_attempt": str(hazard),
        "mean_attempts": str(1 / hazard),
        "variance_attempts": str(survival / hazard**2),
        "entry_pairs_per_doubled_state": int(np.count_nonzero(entrance[0])),
        "uniform_doubled_entry_pairs": len(support),
        "uniform_doubled_entry_pair_mass": str(support[0]),
        "stationary_pair_mass": str(uniform_pair_mass),
        "entry_to_stationary_support_ratio": str(support[0] / uniform_pair_mass),
        "entry_to_stationary_total_variation": str(sum(abs(p - uniform_pair_mass) for p in entry) / 2),
        "entry_one_port_occupancy": str(entry_marginals[0]),
        "stationary_one_port_occupancy": str(stationary_marginals[0]),
        "entry_one_port_marginals_match": entry_marginals == stationary_marginals,
    }


def repair_distance_lumping_census() -> dict[str, object]:
    """Check the distance quotient, its uniform lift, and its binary spectrum.

    The two token positions coincide in a doubled state, whose distance is
    zero. For positive distances, one successful swap moves one token along
    a seam. The intersection array counts its inward and outward options.
    A move into the other occupied port leaves the state unchanged.
    """
    carrier = icosahedral_carrier()
    adjacency = carrier.adjacency_matrix().astype(int)
    distances = graph_distances(adjacency)
    outward, inward = intersection_array(adjacency)
    states, counts = _full_shell_attempt_counts()
    positions = [tuple(port for port, load in enumerate(state) for _ in range(load)) for state in states]
    classes = np.asarray([distances[a][b] for a, b in positions])
    class_count = len(outward) + 1
    membership = np.eye(class_count, dtype=np.int64)[classes]
    sizes = membership.sum(axis=0)
    draws = 2 * carrier.seam_count

    quotient = np.zeros((class_count, class_count), dtype=np.int64)
    for distance, count in enumerate(outward):
        # At distance zero both coins split the double on the chosen seam.
        # At positive distance one coin swaps and either end can move.
        quotient[distance, distance + 1] = 2 * count
    for distance, count in enumerate(inward[1:], start=2):
        quotient[distance, distance - 1] = 2 * count
    np.fill_diagonal(quotient, draws - quotient.sum(axis=1))

    # The two matrix identities test projection for all initial laws and lift
    # for class-uniform laws. Cross multiplication keeps the lift exact.
    all_rows_match = np.array_equal(counts @ membership, membership @ quotient)
    class_uniform_lift = np.array_equal(
        (membership.T @ counts) * sizes[classes],
        sizes[:, None] * quotient[:, classes],
    )
    stationary_counts = sizes.copy()
    stationary_counts[0] = 0
    stationary_masses = tuple(
        str(Fraction(int(count), int(stationary_counts.sum()))) for count in stationary_counts
    )
    stationarity_verified = np.array_equal(stationary_counts @ quotient, draws * stationary_counts)

    # The binary block has one root at one after normalisation. Its trace and
    # sum of principal minors give the other roots' sum and product exactly.
    binary_counts = quotient[1:, 1:]
    trace = int(np.trace(binary_counts))
    second_coefficient = (trace**2 - int(np.trace(binary_counts @ binary_counts))) // 2
    eigenvalue_sum = Fraction(trace - draws, draws)
    eigenvalue_product = Fraction(second_coefficient, draws**2) - eigenvalue_sum
    return {
        "class_sizes": tuple(sizes.tolist()),
        "attempts_per_state": draws,
        "attempt_count_matrix": tuple(tuple(row) for row in quotient.tolist()),
        "all_rows_match": all_rows_match,
        "class_uniform_lift": class_uniform_lift,
        "stationary_class_masses": stationary_masses,
        "stationarity_verified": stationarity_verified,
        "binary_eigenvalue_sum": str(eigenvalue_sum),
        "binary_eigenvalue_product": str(eigenvalue_product),
        "binary_eigenvalue_discriminant": str(eigenvalue_sum**2 - 4 * eigenvalue_product),
    }


def _contrast_lumping_case(
    adjacency: Adjacency,
    squared_chord: Callable[[int, int], Fraction | Surd5],
) -> dict[str, object]:
    """Compare the distance and Gaussian contrast-variance partitions."""
    pairs_and_variances = pair_resistances(adjacency)
    pairs = tuple((left, right) for left, right, _ in pairs_and_variances)
    variances = tuple(variance for _, _, variance in pairs_and_variances)
    chords = tuple(squared_chord(left, right) for left, right in pairs)
    distances = graph_distances(adjacency)
    pair_distances = tuple(distances[left][right] for left, right in pairs)
    seams = tuple((left, right) for left, right in pairs if adjacency[left][right])
    states = tuple(tuple(int(port in pair) for port in range(len(adjacency))) for pair in pairs)
    counts = _attempt_counts(states, seams)
    draws = 2 * len(seams)

    variance_classes = sorted(set(variances))
    class_indices = np.asarray([variance_classes.index(value) for value in variances])
    membership = np.eye(len(variance_classes), dtype=np.int64)[class_indices]
    sizes = membership.sum(axis=0)
    projected = counts @ membership
    representatives = [variances.index(value) for value in variance_classes]
    quotient = projected[representatives]
    distance_indices = np.asarray(pair_distances) - 1
    distance_membership = np.eye(max(pair_distances), dtype=np.int64)[distance_indices]
    distance_counts = counts @ distance_membership
    distance_lumpable = all(
        np.all(rows == rows[0])
        for rows in (distance_counts[distance_indices == group]
                     for group in range(distance_membership.shape[1]))
    )

    # Independently count occupied-to-empty moves. Each seam with one token
    # contributes one move under one of its two coins. All other draws wait.
    variance_index = dict(zip(pairs, class_indices, strict=True))
    move_counts = np.zeros_like(projected)
    for row, (left, right) in enumerate(pairs):
        for mover, fixed in ((left, right), (right, left)):
            for destination, edge in enumerate(adjacency[mover]):
                if edge and destination != fixed:
                    target = tuple(sorted((fixed, destination)))
                    move_counts[row, variance_index[target]] += 1
        move_counts[row, class_indices[row]] += draws - int(move_counts[row].sum())

    witnesses = {}
    adjacency_counts_match = True
    for row, (left, right) in enumerate(pairs):
        if pair_distances[row] != 2:
            continue
        common_neighbors = sum(
            bool(first and second)
            for first, second in zip(adjacency[left], adjacency[right], strict=True)
        )
        adjacent_attempts = int(distance_counts[row, 0])
        adjacency_counts_match &= adjacent_attempts == 2 * common_neighbors
        witnesses.setdefault(variances[row], (
            (left, right), str(chords[row]), str(variances[row]), common_neighbors,
            str(Fraction(adjacent_attempts, draws)),
        ))

    return {
        "attempts_per_state": draws,
        "classes": tuple(
            (pair_distances[row], str(chords[row]), str(variances[row]), int(size))
            for row, size in zip(representatives, sizes, strict=True)
        ),
        "attempt_count_matrix": tuple(tuple(row) for row in quotient.tolist()),
        "contrast_lumpable": bool(np.array_equal(projected, membership @ quotient)),
        "distance_lumpable": distance_lumpable,
        "chord_partition_matches": (
            len(set(zip(chords, variances, strict=True))) == len(set(chords)) == len(variance_classes)
        ),
        "all_rows_match_token_moves": bool(np.array_equal(projected, move_counts)),
        "class_uniform_lift": bool(np.array_equal(
            (membership.T @ counts) * sizes[class_indices],
            sizes[:, None] * quotient[:, class_indices],
        )),
        "stationary_class_masses": tuple(str(Fraction(int(size), len(states))) for size in sizes),
        "stationarity_verified": bool(np.array_equal(sizes @ quotient, draws * sizes)),
        "distance_two_adjacency_counts_match": adjacency_counts_match,
        "distance_two_witnesses": tuple(witnesses[variance] for variance in sorted(witnesses)),
    }


def repair_contrast_lumping_census() -> dict[str, dict[str, object]]:
    """Test exact binary repair quotients on the carrier and the FCC first shell.

    The same fair, uniform local rule is run on each graph. The shell is the
    cuboctahedral graph of nearest neighbours of the FCC origin, with the
    origin omitted. The Gaussian is the centred intrinsic Gaussian with the
    unit graph Laplacian as precision, so contrast variances are resistances.
    Chords use edge length two. Witness indices follow the respective exact
    vertex constructors.
    """
    vertices = exact_icosahedron_vertices()

    def carrier_chord(left: int, right: int) -> Surd5:
        difference = (x - y for x, y in zip(vertices[left], vertices[right], strict=True))
        return Surd5.of(sum(delta * delta for delta in difference))

    carrier_adjacency = tuple(
        tuple(int(carrier_chord(left, right) == 4) for right in range(len(vertices)))
        for left in range(len(vertices))
    )
    oph_adjacency = icosahedral_carrier().adjacency_matrix().astype(int).tolist()
    if not graphs_are_isomorphic(carrier_adjacency, oph_adjacency):
        raise AssertionError("the exact icosahedron must realise the OPH carrier")

    sites = cuboctahedron_neighbors()
    return {
        "carrier": _contrast_lumping_case(carrier_adjacency, carrier_chord),
        "cuboctahedron": _contrast_lumping_case(
            cuboctahedron_adjacency(),
            lambda left, right: Fraction(squared_distance(sites[left], sites[right]), 2),
        ),
    }


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
