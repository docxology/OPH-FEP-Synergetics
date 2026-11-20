"""The carrier read as a network: resistance, distance shells and subgraph counts.

The Laplacian already carries two readings in this package: the spanning-tree
number, from its determinant, and the Gaussian Markov random field, from its
inverse. This module adds the readings that need no programme of their own but
that both vocabularies can be checked against:

* the **effective resistance** between two ports, which is the same matrix seen
  as an electrical network and, in the Gaussian reading, the variance of the
  potential difference. It is computed exactly, over the rationals, by
  grounding one port and solving ``L x = e_i - e_j``.
* the **distance partition**: the shells around a port and the intersection
  array that says how a neighbour of a shell member is distributed.
* the **subgraph census**: the triangles, the induced (chordless) five-cycles,
  the perfect matchings and the Hamiltonian cycles (as undirected edge sets), all
  counted by exhaustion rather than quoted.
* the **spectral chord comparison**: the carrier's regular icosahedral
  coordinates give the contribution of the first positive Laplacian eigenspace
  to each pair variance, up to a common scale.

The network calculations take only an adjacency matrix. The chord comparison
also uses the exact coordinates of the regular icosahedron.
"""

from __future__ import annotations

import itertools
from collections import Counter
from fractions import Fraction
from typing import Callable, Iterator, Mapping, Sequence

import numpy as np

__all__ = [
    "carrier_resistance_geometry",
    "cuboctahedral_resistance_geometry",
    "distance_distribution",
    "effective_resistance",
    "graph_distances",
    "graphs_are_isomorphic",
    "intersection_array",
    "kirchhoff_index",
    "laplacian_matrix",
    "resistance_ladder",
    "subgraph_census",
]

Adjacency = Sequence[Sequence[int]]


def laplacian_matrix(adjacency: Adjacency) -> list[list[int]]:
    """Return the graph Laplacian of an adjacency matrix, in exact integers."""
    size = len(adjacency)
    return [
        [(sum(adjacency[i]) if i == j else -int(adjacency[i][j])) for j in range(size)]
        for i in range(size)
    ]


def graph_distances(adjacency: Adjacency) -> tuple[tuple[int, ...], ...]:
    """Return the shortest-path distance matrix, by Floyd-Warshall."""
    size = len(adjacency)
    unreachable = size + 1
    distances = [
        [0 if i == j else (1 if adjacency[i][j] else unreachable) for j in range(size)]
        for i in range(size)
    ]
    for through in range(size):
        for i in range(size):
            for j in range(size):
                step = distances[i][through] + distances[through][j]
                if step < distances[i][j]:
                    distances[i][j] = step
    return tuple(tuple(row) for row in distances)


def distance_distribution(adjacency: Adjacency, source: int = 0) -> tuple[int, ...]:
    """Return the number of ports at each distance from ``source``, nearest first."""
    distances = graph_distances(adjacency)[source]
    return tuple(sum(1 for value in distances if value == step) for step in range(max(distances) + 1))


def intersection_array(adjacency: Adjacency, source: int = 0) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Return the intersection array ``(b, c)`` in its standard form.

    For a graph of diameter ``d`` the array is ``{b_0, ..., b_(d-1); c_1, ...,
    c_d}``: ``b[k]`` counts the neighbours of a vertex at distance ``k`` from
    ``source`` that lie at distance ``k + 1``, and ``c[k - 1]`` counts the
    neighbours of a vertex at distance ``k`` that lie at distance ``k - 1``.
    ``b_0`` is the degree and ``b_d = 0`` is omitted, so the icosahedron's array
    is ``{5, 2, 1; 1, 2, 5}``. One representative per shell is enough exactly
    when the graph is distance-regular, which the caller asserts by comparing
    representatives.
    """
    distances = graph_distances(adjacency)[source]
    reach = max(distances)
    representatives = {step: distances.index(step) for step in range(reach + 1)}
    back = []
    onward = []
    for step in sorted(representatives):
        port = representatives[step]
        neighbours = [other for other in range(len(adjacency)) if adjacency[port][other]]
        if step < reach:
            onward.append(sum(1 for other in neighbours if distances[other] == step + 1))
        if step > 0:
            back.append(sum(1 for other in neighbours if distances[other] == step - 1))
    return tuple(onward), tuple(back)


def _grounded_inverse(adjacency: Adjacency) -> list[list[Fraction]]:
    """Invert the Laplacian exactly, with port 0 grounded (row and column dropped).

    One elimination serves every pair: the grounded inverse determines each
    effective resistance directly, so the pair formulas below share it.
    """
    size = len(adjacency) - 1
    laplacian = laplacian_matrix(adjacency)
    matrix = [[Fraction(laplacian[i + 1][j + 1]) for j in range(size)] for i in range(size)]
    inverse = [[Fraction(1 if row == column else 0) for column in range(size)] for row in range(size)]
    for column in range(size):
        pivot = next((row for row in range(column, size) if matrix[row][column] != 0), None)
        if pivot is None:
            raise AssertionError("a grounded Laplacian must be invertible")
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            inverse[column], inverse[pivot] = inverse[pivot], inverse[column]
        diagonal = matrix[column][column]
        matrix[column] = [value / diagonal for value in matrix[column]]
        inverse[column] = [value / diagonal for value in inverse[column]]
        for row in range(size):
            factor = matrix[row][column]
            if row != column and factor:
                matrix[row] = [
                    value - factor * other for value, other in zip(matrix[row], matrix[column])
                ]
                inverse[row] = [
                    value - factor * other for value, other in zip(inverse[row], inverse[column])
                ]
    return inverse


def _resistance_from_inverse(
    inverse: Sequence[Sequence[Fraction]], source: int, target: int
) -> Fraction:
    """Read the effective resistance off a grounded Laplacian inverse.

    With port 0 grounded, the pair solve ``(e_source - e_target) . L_reduced^-1 .
    (e_source - e_target)`` expands to ``g_ii + g_jj - 2 g_ij`` over the reduced
    indices, the grounded port contributing zero.
    """
    def diagonal(vertex: int) -> Fraction:
        return Fraction(0) if vertex == 0 else inverse[vertex - 1][vertex - 1]

    off_diagonal = Fraction(0)
    if source and target:
        off_diagonal = inverse[source - 1][target - 1]
    return diagonal(source) + diagonal(target) - 2 * off_diagonal


def _pair_resistances(
    adjacency: Adjacency, inverse: Sequence[Sequence[Fraction]]
) -> Iterator[tuple[int, int, Fraction]]:
    """Yield ``(source, target, resistance)`` for every unordered pair."""
    for source, target in itertools.combinations(range(len(adjacency)), 2):
        yield source, target, _resistance_from_inverse(inverse, source, target)


def effective_resistance(adjacency: Adjacency, source: int, target: int) -> Fraction:
    """Return the effective resistance between two ports, exactly.

    Grounding port 0 and inverting the reduced Laplacian once puts every pair's
    resistance on the inverse's own entries; the value is gauge-invariant and
    rational for an integer graph.
    """
    return _resistance_from_inverse(_grounded_inverse(adjacency), source, target)


def resistance_ladder(adjacency: Adjacency) -> tuple[tuple[int, Fraction, int], ...]:
    """Return ``(distance, resistance, pairs)`` for every pair class, nearest first."""
    inverse = _grounded_inverse(adjacency)
    distances = graph_distances(adjacency)
    classes: dict[int, set[Fraction]] = {}
    counts: dict[int, int] = {}
    for source, target, resistance in _pair_resistances(adjacency, inverse):
        step = distances[source][target]
        classes.setdefault(step, set()).add(resistance)
        counts[step] = counts.get(step, 0) + 1
    ladder = []
    for step in sorted(classes):
        if len(classes[step]) != 1:
            raise AssertionError("this graph has more than one resistance in a distance class")
        ladder.append((step, classes[step].pop(), counts[step]))
    return tuple(ladder)


def kirchhoff_index(adjacency: Adjacency) -> Fraction:
    """Return the sum of the effective resistances over all unordered pairs."""
    inverse = _grounded_inverse(adjacency)
    total = Fraction(0)
    for _, _, resistance in _pair_resistances(adjacency, inverse):
        total += resistance
    return total


def carrier_resistance_geometry() -> dict[str, object]:
    """Compare edge-two icosahedral chords with intrinsic Gaussian variances.

    The chord configuration spans the first positive eigenspace of the unit
    carrier Laplacian. Verify this in Q(sqrt(5)) and derive the scale from
    the squared norm of its coordinate columns. The full Gaussian covariance
    embeds the square root of resistance in dimension n-1 because the
    grounded inverse is the positive definite Gram matrix relative to port 0.
    """
    from ..geometry.exact import Surd5, exact_icosahedron_vertices
    from .carrier import icosahedral_carrier

    vertices = exact_icosahedron_vertices()
    chords = [
        [sum(delta * delta for delta in (a - b for a, b in zip(left, right))) for right in vertices]
        for left in vertices
    ]
    adjacency = [[int(chord == 4) for chord in row] for row in chords]
    oph_adjacency = icosahedral_carrier().adjacency_matrix().astype(int).tolist()
    if not graphs_are_isomorphic(adjacency, oph_adjacency):
        raise AssertionError("the exact chord adjacency must realize the OPH carrier")

    inverse = _grounded_inverse(adjacency)
    distances = graph_distances(adjacency)
    classes: dict[int, list[tuple[Fraction, Surd5]]] = {}
    for source, target, resistance in _pair_resistances(adjacency, inverse):
        step = distances[source][target]
        classes.setdefault(step, []).append((resistance, Surd5.of(chords[source][target])))
    shells = []
    for step, pairs in sorted(classes.items()):
        (resistance, chord), = set(pairs)
        shells.append((step, resistance, chord, len(pairs)))

    laplacian = np.asarray(laplacian_matrix(adjacency), dtype=np.int64)
    identity = np.eye(len(vertices), dtype=np.int64)
    shifted = laplacian - 5 * identity
    annihilator = laplacian @ (laplacian - 6 * identity) @ (shifted @ shifted - 5 * identity)
    if np.any(annihilator):
        raise AssertionError("the carrier Laplacian must satisfy its spectral polynomial")
    # An integer trace forces the two irrational roots to have equal
    # multiplicity. Invertibility after grounding leaves one zero eigenvalue.
    positive_rank = len(inverse)
    lowest_multiplicity = (6 * positive_rank - int(np.trace(laplacian))) // 2

    coordinates = np.asarray(vertices, dtype=object)
    lowest_eigenvalue = Surd5(5, -1)
    if not np.array_equal(laplacian @ coordinates, coordinates * lowest_eigenvalue):
        raise AssertionError("the exact coordinate columns must be lowest positive eigenvectors")
    column_gram = coordinates.T @ coordinates
    chord_rank = coordinates.shape[1]
    coordinate_norm_squared = column_gram[0, 0]
    if (
        coordinate_norm_squared == 0
        or not np.array_equal(column_gram, np.eye(chord_rank, dtype=object) * coordinate_norm_squared)
        or chord_rank != lowest_multiplicity
    ):
        raise AssertionError("the scaled coordinate columns must be a basis of the eigenspace")

    factor = coordinate_norm_squared * lowest_eigenvalue
    _, edge_resistance, edge_chord, _ = shells[0]
    return {
        "shells": tuple(shells),
        "ranks": (positive_rank, chord_rank),
        "lowest_eigenvalue": lowest_eigenvalue,
        "lowest_multiplicity": lowest_multiplicity,
        "coordinate_norm_squared": coordinate_norm_squared,
        "chord_to_lowest_variance_factor": factor,
        "lowest_mode_variances": tuple(chord / factor for _, _, chord, _ in shells),
        "full_variance_proportional": all(
            chord * edge_resistance == edge_chord * resistance
            for _, resistance, chord, _ in shells
        ),
    }


def cuboctahedral_resistance_geometry() -> dict[str, object]:
    """Compare chords and intrinsic contrast variances on the first IVM shell.

    Use the exact integer embedding, with squared edge eight. The reported
    squared chords use edge two, so their squares and coordinate Gram are
    halved. Check the Laplacian projector formulas in integer arithmetic and
    compare their full contrast variance with a separate grounded solve.
    """
    from ..geometry.lattice import cuboctahedron_adjacency, cuboctahedron_neighbors
    from ..geometry.quadray import EMBEDDING

    adjacency = cuboctahedron_adjacency()
    xyz = np.asarray([site.as_tuple() for site in cuboctahedron_neighbors()]) @ np.asarray(EMBEDDING).T
    laplacian = np.asarray(laplacian_matrix(adjacency), dtype=np.int64)
    identity = np.eye(len(adjacency), dtype=np.int64)
    eigenvalues = (0, 2, 4, 6)
    spectrum = []
    spectral_inverse = np.zeros_like(laplacian, dtype=object)
    for eigenvalue in eigenvalues:
        other_roots = [root for root in eigenvalues if root != eigenvalue]
        numerator = identity
        denominator = 1
        for root in other_roots:
            numerator = numerator @ (laplacian - root * identity)
            denominator *= eigenvalue - root
        if not np.array_equal(laplacian @ numerator, eigenvalue * numerator):
            raise AssertionError("the cuboctahedral Laplacian must satisfy its spectral polynomial")
        multiplicity = Fraction(int(np.trace(numerator)), denominator)
        spectrum.append((eigenvalue, int(multiplicity)))
        if eigenvalue:
            spectral_inverse += numerator.astype(object) * Fraction(1, denominator * eigenvalue)

    smallest = eigenvalues[1]
    gram = xyz.T @ xyz
    coordinate_norm_squared = Fraction(int(gram[0, 0]), 2)
    if (
        not np.array_equal(laplacian @ xyz, smallest * xyz)
        or not np.array_equal(gram, gram[0, 0] * np.eye(xyz.shape[1], dtype=np.int64))
        or coordinate_norm_squared <= 0
        or spectrum[1][1] != xyz.shape[1]
    ):
        raise AssertionError("the shell coordinates must form a basis of the first positive eigenspace")
    factor = smallest * coordinate_norm_squared

    grounded = _grounded_inverse(adjacency)
    graph_distance = graph_distances(adjacency)
    classes = Counter()
    for left, right, resistance in _pair_resistances(adjacency, grounded):
        delta = xyz[left] - xyz[right]
        squared_chord = Fraction(int(delta @ delta), 2)
        from_spectrum = (spectral_inverse[left, left] + spectral_inverse[right, right]
                         - 2 * spectral_inverse[left, right])
        if from_spectrum != resistance:
            raise AssertionError("spectral and grounded variances must agree on every pair")
        classes[graph_distance[left][right], squared_chord, resistance, squared_chord / factor] += 1

    return {
        "spectrum": tuple(spectrum),
        "chord_to_lowest_variance_factor": factor,
        "ranks": (len(grounded), xyz.shape[1]),
        "pair_classes": tuple((*group, count) for group, count in sorted(classes.items())),
        "kirchhoff_index": sum(group[2] * count for group, count in classes.items()),
    }


def subgraph_census(adjacency: Adjacency) -> dict[str, int]:
    """Count the triangles, induced five-cycles, perfect matchings and Hamiltonian cycles.

    A five-cycle is counted when its five vertices induce exactly a cycle, so the
    carrier's count is the twelve vertex links (the dual dodecahedron's faces);
    the graph has more five-cycles with chords.
    """
    size = len(adjacency)
    triangles = sum(
        all(adjacency[i][j] for i, j in itertools.combinations(corner, 2))
        for corner in itertools.combinations(range(size), 3)
    )
    pentagons = 0
    for candidate in itertools.combinations(range(size), 5):
        degrees = [sum(adjacency[vertex][other] for other in candidate) for vertex in candidate]
        if all(degree == 2 for degree in degrees):
            pentagons += 1
    return {
        "triangles": triangles,
        "pentagons": pentagons,
        "perfect_matchings": _perfect_matchings(adjacency),
        "hamiltonian_cycles": _hamiltonian_cycles(adjacency),
    }


def _perfect_matchings(adjacency: Adjacency) -> int:
    size = len(adjacency)
    memo: dict[int, int] = {}

    def count(remaining: int) -> int:
        """Count the matchings of the vertices still unpaired.

        The cursor walks the candidates while `rest` stays fixed: recursing on
        the cursor instead would drop the partners already tried.
        """
        if remaining == 0:
            return 1
        if remaining in memo:
            return memo[remaining]
        first = (remaining & -remaining).bit_length() - 1
        rest = remaining ^ (1 << first)
        cursor = rest
        total = 0
        while cursor:
            partner = (cursor & -cursor).bit_length() - 1
            if adjacency[first][partner]:
                total += count(rest ^ (1 << partner))
            cursor ^= 1 << partner
        memo[remaining] = total
        return total

    return count((1 << size) - 1)


def _hamiltonian_cycles(adjacency: Adjacency) -> int:
    size = len(adjacency)
    memo: dict[tuple[int, int], int] = {}

    def count(mask: int, last: int) -> int:
        if mask == (1 << size) - 1:
            return 1 if adjacency[last][0] else 0
        key = (mask, last)
        if key in memo:
            return memo[key]
        total = 0
        for following in range(size):
            if not mask >> following & 1 and adjacency[last][following]:
                total += count(mask | (1 << following), following)
        memo[key] = total
        return total

    return count(1, 0) // 2

def _isomorphism_backtrack(
    size: int,
    order: Sequence[int],
    signature_ok: Callable[[int, int], bool],
    mapped_ok: Callable[[int, int, Mapping[int, int]], bool],
) -> dict[int, int] | None:
    """Backtrack a left-to-right vertex mapping, or return ``None``.

    Places the left vertices in ``order``, trying each unused right vertex of
    the same ``size``: ``signature_ok`` prunes on static vertex invariants and
    ``mapped_ok`` prunes against the placements already made. Callers keep
    their own pruning rules; the search itself is shared.
    """
    mapping: dict[int, int] = {}
    used: set[int] = set()

    def place(position: int) -> bool:
        if position == len(order):
            return True
        vertex = order[position]
        for candidate in range(size):
            if candidate in used or not signature_ok(vertex, candidate):
                continue
            if not mapped_ok(vertex, candidate, mapping):
                continue
            mapping[vertex] = candidate
            used.add(candidate)
            if place(position + 1):
                return True
            used.discard(candidate)
            del mapping[vertex]
        return False

    if not place(0):
        return None
    return mapping


def graphs_are_isomorphic(left: Adjacency, right: Adjacency) -> bool:
    """Whether two adjacency matrices are related by a vertex relabelling.

    Backtracking on the most constrained vertex first, with the degree sequence
    as the pruning rule. It is exhaustive, and it is small enough for the
    twelve-vertex graphs this package compares.
    """
    size = len(left)
    if size != len(right) or sorted(sum(row) for row in left) != sorted(sum(row) for row in right):
        return False
    order = sorted(range(size), key=lambda vertex: -sum(left[vertex]))

    def signature_ok(vertex: int, candidate: int) -> bool:
        return sum(right[candidate]) == sum(left[vertex])

    def mapped_ok(vertex: int, candidate: int, mapping: Mapping[int, int]) -> bool:
        return all(left[vertex][other] == right[candidate][mapping[other]] for other in mapping)

    return _isomorphism_backtrack(size, order, signature_ok, mapped_ok) is not None
