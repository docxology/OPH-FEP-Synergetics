"""The carrier read as a network: resistance, distance shells and subgraph counts.

The Laplacian already carries two readings in this package: the spanning-tree
number, from its determinant, and the Gaussian Markov random field, from its
inverse. This module adds the readings that need no programme of their own but
that both vocabularies can be checked against:

* the **effective resistance** between two ports, which is the same matrix seen
  as an electrical network and, in the Gaussian reading, the marginal variance
  of the potential difference. It is computed exactly, over the rationals, by
  grounding one port and solving ``L x = e_i - e_j``.
* the **distance partition**: the shells around a port and the intersection
  array that says how a neighbour of a shell member is distributed.
* the **subgraph census**: the triangles, the five-cycles, the perfect matchings
  and the Hamiltonian cycles, all counted by exhaustion rather than quoted.

Nothing here imports a geometry module: an adjacency matrix is all these
functions need, so the same code counts the carrier, the cuboctahedron and any
line graph put beside them.
"""

from __future__ import annotations

import itertools
from fractions import Fraction
from typing import Sequence

__all__ = [
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
    """Return ``(b, c)``: neighbours outwards and back, shell by shell.

    ``b[k - 1]`` counts the neighbours of a vertex in shell ``k`` that lie in
    shell ``k + 1`` and ``c[k - 1]`` those in shell ``k - 1``. One representative
    per shell is enough exactly when the graph is distance-regular, which the
    caller asserts by comparing representatives.
    """
    distances = graph_distances(adjacency)[source]
    reach = max(distances)
    representatives = {step: distances.index(step) for step in range(1, reach + 1)}
    back = []
    onward = []
    for step in sorted(representatives):
        port = representatives[step]
        neighbours = [other for other in range(len(adjacency)) if adjacency[port][other]]
        back.append(sum(1 for other in neighbours if distances[other] == step - 1))
        onward.append(sum(1 for other in neighbours if distances[other] == step + 1))
    return tuple(onward), tuple(back)


def _solve_grounded(laplacian: Sequence[Sequence[int]], rhs: Sequence[int]) -> list[Fraction]:
    """Solve ``L x = rhs`` exactly with port 0 grounded (row and column dropped)."""
    size = len(laplacian) - 1
    matrix = [[Fraction(laplacian[i + 1][j + 1]) for j in range(size)] for i in range(size)]
    vector = [Fraction(value) for value in rhs]
    for column in range(size):
        pivot = next((row for row in range(column, size) if matrix[row][column] != 0), None)
        if pivot is None:
            raise AssertionError("a grounded Laplacian must be invertible")
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        vector[column], vector[pivot] = vector[pivot], vector[column]
        diagonal = matrix[column][column]
        for row in range(column + 1, size):
            factor = matrix[row][column] / diagonal
            if factor:
                for col in range(column, size):
                    matrix[row][col] -= factor * matrix[column][col]
                vector[row] -= factor * vector[column]
    solution = [Fraction(0)] * size
    for column in range(size - 1, -1, -1):
        total = vector[column] - sum(matrix[column][col] * solution[col] for col in range(column + 1, size))
        solution[column] = total / matrix[column][column]
    return solution


def effective_resistance(adjacency: Adjacency, source: int, target: int) -> Fraction:
    """Return the effective resistance between two ports, exactly.

    Grounding port 0 and solving ``L x = e_source - e_target`` gives the
    potentials; the resistance is ``(e_source - e_target) . x``, which is
    gauge-invariant and rational for an integer graph.
    """
    laplacian = laplacian_matrix(adjacency)
    size = len(adjacency) - 1
    rhs = [0] * size
    if source:
        rhs[source - 1] += 1
    if target:
        rhs[target - 1] -= 1
    solution = _solve_grounded(laplacian, rhs)
    return sum(Fraction(rhs[index]) * solution[index] for index in range(size))


def resistance_ladder(adjacency: Adjacency) -> tuple[tuple[int, Fraction, int], ...]:
    """Return ``(distance, resistance, pairs)`` for every pair class, nearest first."""
    distances = graph_distances(adjacency)
    classes: dict[int, set[Fraction]] = {}
    counts: dict[int, int] = {}
    for source, target in itertools.combinations(range(len(adjacency)), 2):
        step = distances[source][target]
        classes.setdefault(step, set()).add(effective_resistance(adjacency, source, target))
        counts[step] = counts.get(step, 0) + 1
    ladder = []
    for step in sorted(classes):
        if len(classes[step]) != 1:
            raise AssertionError("an edge-transitive carrier has one resistance per distance class")
        ladder.append((step, classes[step].pop(), counts[step]))
    return tuple(ladder)


def kirchhoff_index(adjacency: Adjacency) -> Fraction:
    """Return the sum of the effective resistances over all unordered pairs."""
    total = Fraction(0)
    for source, target in itertools.combinations(range(len(adjacency)), 2):
        total += effective_resistance(adjacency, source, target)
    return total


def subgraph_census(adjacency: Adjacency) -> dict[str, int]:
    """Count the triangles, five-cycles, perfect matchings and Hamiltonian cycles."""
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
    mapping: dict[int, int] = {}
    used: set[int] = set()

    def place(position: int) -> bool:
        if position == size:
            return True
        vertex = order[position]
        for candidate in range(size):
            if candidate in used or sum(right[candidate]) != sum(left[vertex]):
                continue
            if any(left[vertex][other] != right[candidate][mapping[other]] for other in mapping):
                continue
            mapping[vertex] = candidate
            used.add(candidate)
            if place(position + 1):
                return True
            used.discard(candidate)
            del mapping[vertex]
        return False

    return place(0)
