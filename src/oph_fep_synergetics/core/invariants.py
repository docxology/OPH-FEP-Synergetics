"""Spectral and probabilistic invariants of a carrier's Laplacian.

Three views of one matrix, each belonging to a different programme:

* the **spanning-tree number**, by the matrix-tree theorem: the number of ways the
  seam structure can be reduced to a connected acyclic seam set.  For the
  twelve-port carrier it is exactly ``5184000``, and the same integer is
  ``(1/n) * prod(lambda != 0)`` over the closed-form spectrum.
* the **intrinsic Gaussian Markov random field** whose precision is the
  Laplacian: its conditional-independence structure is exactly the carrier's
  non-adjacency, so every port's Markov blanket is its five neighbours.
* the **exact marginal variance** ``7/36`` of a port, obtained from the
  reciprocal spectrum in closed form rather than from a numerical inverse.

The reciprocal-spectrum identity is exact in the rationals: the three copies of
``5 - sqrt(5)`` and the three of ``5 + sqrt(5)`` contribute
``3 * (5 - sqrt 5 + 5 + sqrt 5) / 20 = 3/2``, and the five copies of ``6``
contribute ``5/6``, so ``trace(L^+) = 7/3`` and every vertex-transitive diagonal
entry is ``7/36``.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Sequence

import numpy as np

from .carrier import IcosahedralCarrier, icosahedral_carrier
from .network import laplacian_matrix

__all__ = [
    "blanket_sizes",
    "gmrf_precision",
    "independence_matches_graph",
    "independence_pattern",
    "laplacian_determinant",
    "marginal_variance",
    "markov_blanket",
    "pairwise_precision",
    "proper_precision",
    "reciprocal_spectrum_sum",
    "reduced_laplacian_determinant",
    "spanning_tree_number",
    "spanning_tree_number_of_adjacency",
    "tree_number_from_spectrum",
]


def laplacian_determinant(laplacian: Sequence[Sequence[int]]) -> int:
    """Return the determinant of a Laplacian's grounded minor, exactly.

    The minor is taken over ``Fraction`` arithmetic, so the result is an exact
    integer rather than a floating-point shadow of one. Every graph whose
    adjacency matrix is available can be counted with this function, which is
    why the carrier and the cuboctahedron are compared by the same code path.
    """
    size = len(laplacian) - 1
    matrix = [
        [Fraction(int(laplacian[i][j])) for j in range(1, size + 1)]
        for i in range(1, size + 1)
    ]
    determinant = Fraction(1)
    for column in range(size):
        pivot = next((row for row in range(column, size) if matrix[row][column] != 0), None)
        if pivot is None:
            return 0
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            determinant = -determinant
        determinant *= matrix[column][column]
        pivot_value = matrix[column][column]
        for row in range(column + 1, size):
            factor = matrix[row][column] / pivot_value
            if factor:
                for col in range(column, size):
                    matrix[row][col] -= factor * matrix[column][col]
    if determinant.denominator != 1:
        raise AssertionError("a reduced Laplacian determinant must be an integer")
    return int(determinant)


def reduced_laplacian_determinant(carrier: IcosahedralCarrier | None = None) -> int:
    """Return the exact determinant of the reduced Laplacian as an integer.

    The determinant is taken over ``Fraction`` arithmetic, so the result is an
    exact integer rather than a floating-point shadow of one.
    """
    carrier = carrier or icosahedral_carrier()
    laplacian = np.rint(carrier.laplacian()).astype(int)
    return laplacian_determinant([[int(value) for value in row] for row in laplacian])


def spanning_tree_number_of_adjacency(adjacency: Sequence[Sequence[int]]) -> int:
    """Return the number of spanning trees of any integer graph, exactly."""
    return laplacian_determinant(laplacian_matrix(adjacency))


def spanning_tree_number(carrier: IcosahedralCarrier | None = None) -> int:
    """Return the number of spanning trees of the carrier, exactly."""
    return reduced_laplacian_determinant(carrier)


def tree_number_from_spectrum(carrier: IcosahedralCarrier | None = None) -> float:
    """Return ``(1/n) * prod(lambda != 0)`` from the computed spectrum."""
    carrier = carrier or icosahedral_carrier()
    spectrum = carrier.laplacian_spectrum()
    nonzero = [float(value) for value in spectrum if value > 1e-9]
    product = 1.0
    for value in nonzero:
        product *= value
    return product / carrier.port_count


def reciprocal_spectrum_sum(carrier: IcosahedralCarrier | None = None) -> Fraction:
    """Return ``sum(1/lambda)`` over the nonzero spectrum, exactly.

    The surd pairs cancel: ``1/(5 - sqrt 5) + 1/(5 + sqrt 5) = 1/2``, so three
    copies of each pair contribute ``3/2`` and the five copies of ``6``
    contribute ``5/6``.
    """
    carrier = carrier or icosahedral_carrier()
    spectrum = sorted(float(value) for value in carrier.laplacian_spectrum())
    nonzero = [value for value in spectrum if value > 1e-9]
    surd_pairs = [value for value in nonzero if abs(value - 6.0) > 1e-9]
    sixes = len(nonzero) - len(surd_pairs)
    if len(surd_pairs) % 2 != 0:
        raise AssertionError("the surd eigenvalues must come in conjugate pairs")
    pair_reciprocal = Fraction(1, 1) / Fraction(2, 1)  # 1/(5-sqrt5) + 1/(5+sqrt5) = 1/2
    return Fraction(len(surd_pairs) // 2, 1) * pair_reciprocal + Fraction(sixes, 6)


def marginal_variance(carrier: IcosahedralCarrier | None = None) -> Fraction:
    """Return the exact marginal variance of a port under the intrinsic GMRF.

    The carrier is vertex-transitive, so every diagonal entry of ``L^+`` is
    ``trace(L^+) / n = (7/3) / 12``.
    """
    carrier = carrier or icosahedral_carrier()
    return reciprocal_spectrum_sum(carrier) / carrier.port_count


def gmrf_precision(carrier: IcosahedralCarrier | None = None) -> np.ndarray:
    """Return the intrinsic GMRF precision: the carrier Laplacian itself."""
    carrier = carrier or icosahedral_carrier()
    return carrier.laplacian()


def proper_precision(carrier: IcosahedralCarrier | None = None) -> np.ndarray:
    """Return ``L + (1/n) I``, a positive-definite version of the precision.

    Adding the smallest possible shift makes the field proper without changing
    any off-diagonal entry, so the conditional-independence structure is
    unaltered.
    """
    carrier = carrier or icosahedral_carrier()
    laplacian = gmrf_precision(carrier)
    shift = 1.0 / carrier.port_count
    return laplacian + shift * np.eye(carrier.port_count)


def independence_pattern(carrier: IcosahedralCarrier | None = None) -> tuple[tuple[bool, ...], ...]:
    """Return the off-diagonal nonzeros of the precision as a boolean pattern."""
    carrier = carrier or icosahedral_carrier()
    precision = gmrf_precision(carrier)
    return tuple(
        tuple(i != j and abs(float(precision[i][j])) > 1e-12 for j in range(carrier.port_count))
        for i in range(carrier.port_count)
    )


def independence_matches_graph(carrier: IcosahedralCarrier | None = None) -> bool:
    """Whether conditional dependence holds exactly on the seams and nowhere else."""
    carrier = carrier or icosahedral_carrier()
    pattern = independence_pattern(carrier)
    adjacency = np.rint(carrier.adjacency_matrix()).astype(int)
    return all(pattern[i][j] == bool(adjacency[i][j]) for i in range(carrier.port_count) for j in range(carrier.port_count))


def pairwise_precision(carrier: IcosahedralCarrier | None = None) -> dict[tuple[int, int], int]:
    """Return the exact pairwise precision, ``-1`` on every seam."""
    carrier = carrier or icosahedral_carrier()
    precision = np.rint(carrier.laplacian()).astype(int)
    return {
        (i, j): int(precision[i][j])
        for i in range(carrier.port_count)
        for j in range(carrier.port_count)
        if i != j and precision[i][j] != 0
    }


def markov_blanket(carrier: IcosahedralCarrier | None, port: int) -> tuple[int, ...]:
    """Return the Markov blanket of a port: exactly its neighbours."""
    carrier = carrier or icosahedral_carrier()
    return tuple(carrier.neighbours(port))


def blanket_sizes(carrier: IcosahedralCarrier | None = None) -> tuple[int, ...]:
    """Return the blanket size of every port, in port order."""
    carrier = carrier or icosahedral_carrier()
    return tuple(len(markov_blanket(carrier, port)) for port in range(carrier.port_count))
