"""Variational free energy as precision-weighted interface disagreement.

For a finite linear-Gaussian cellular sheaf with energy

.. math:: E(x, o) = \\tfrac12 z^T L z, \\qquad z = (x, o),

the free-energy dictionary of the accompanying review is exact:

* perception with records held fixed is sheaf diffusion and converges to the
  harmonic extension ``x_hat(o) = -L_xx^{-1} L_xo o``;
* the record precision is the Schur complement
  ``S = L_oo - L_ox L_xx^{-1} L_xo`` plus the record prior precision.
  A proper Gaussian record density gives an absolute surprisal ``-ln p(o)``;
* the Bayesian surprise of new records is ``\\tfrac12 g^T L_xx^{-1} g`` at the
  old rest, where ``g`` is the update residual there;
* the variational free energy is ``surprisal + KL[q || posterior]`` and is
  therefore minimised exactly at the posterior.

The module also exposes the three named surprises (numerical, Bayesian, and
surprisal) so the manuscript can distinguish what repair removes from what it
cannot.  The nonlinear, non-Gaussian and multi-modal cases are outside this
implementation and are stated as scope limits in the manuscript.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from typing import Sequence

import numpy as np

from .carrier import icosahedral_carrier
from .invariants import laplacian_determinant
from .network import graph_distances, laplacian_matrix, resistance_ladder
from .sheaf import CellularSheaf, linear_chain_sheaf


def _as_vector(value: Sequence[float], length: int, name: str) -> np.ndarray:
    vector = np.asarray(value, dtype=float)
    if vector.shape != (length,):
        raise ValueError(f"{name} must have shape ({length},)")
    if not np.isfinite(vector).all():
        raise ValueError(f"{name} must be finite")
    return vector


def _symmetrise(matrix: np.ndarray) -> np.ndarray:
    return 0.5 * (matrix + matrix.T)


def _covariance_and_logdet(
    value: np.ndarray, dimension: int, name: str
) -> tuple[np.ndarray, float]:
    covariance = np.asarray(value, dtype=float)
    if covariance.shape != (dimension, dimension):
        raise ValueError(f"{name} must have shape ({dimension}, {dimension})")
    if not np.isfinite(covariance).all() or not np.allclose(
        covariance, covariance.T, rtol=1e-12, atol=0.0
    ):
        raise ValueError(f"{name} must be finite and symmetric")
    # Use the original entries when already symmetric to preserve their scale.
    if not np.array_equal(covariance, covariance.T):
        covariance = _symmetrise(covariance)
    try:
        factor = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError(f"{name} must be positive definite") from exc
    return covariance, 2.0 * float(np.log(np.diag(factor)).sum())


def gaussian_kl(
    mu_0: Sequence[float],
    sigma_0: np.ndarray,
    mu_1: Sequence[float],
    sigma_1: np.ndarray,
) -> float:
    """The KL divergence ``KL[N(mu_0, sigma_0) || N(mu_1, sigma_1)]``."""
    first = np.asarray(mu_0, dtype=float)
    if first.ndim != 1:
        raise ValueError("Gaussian means must be vectors")
    dimension = first.size
    first = _as_vector(first, dimension, "first mean")
    second = _as_vector(mu_1, dimension, "second mean")
    covariance_0, logdet_0 = _covariance_and_logdet(sigma_0, dimension, "first covariance")
    covariance_1, logdet_1 = _covariance_and_logdet(sigma_1, dimension, "second covariance")
    difference = second - first
    solved = np.linalg.solve(covariance_1, covariance_0)
    trace_term = float(np.trace(solved))
    quadratic = float(difference @ np.linalg.solve(covariance_1, difference))
    return 0.5 * (trace_term + quadratic - dimension + logdet_1 - logdet_0)


@dataclass(frozen=True)
class SurpriseTriple:
    """Three quantities called surprise, kept numerically separate."""

    numerical: float
    bayesian: float
    surprisal: float


class GaussianSheafModel:
    """A sheaf model with an optional isotropic Gaussian record prior."""

    def __init__(
        self,
        sheaf: CellularSheaf,
        record_nodes: Sequence[int],
        *,
        record_prior_precision: float | None = None,
    ) -> None:
        if record_prior_precision is not None and (
            not math.isfinite(record_prior_precision) or record_prior_precision <= 0.0
        ):
            raise ValueError("record prior precision must be finite and positive")
        self._record_prior_precision = record_prior_precision
        self.sheaf = sheaf
        self.record_nodes = tuple(sorted(int(node) for node in record_nodes))
        if len(set(self.record_nodes)) != len(self.record_nodes):
            raise ValueError("record nodes must be distinct")
        if any(not 0 <= node < sheaf.node_count for node in self.record_nodes):
            raise IndexError("record node outside the sheaf")
        self.interior_nodes = tuple(
            node for node in range(sheaf.node_count) if node not in self.record_nodes
        )
        offsets = sheaf.offsets()
        record_indices = [
            index
            for node in self.record_nodes
            for index in range(offsets[node], offsets[node + 1])
        ]
        interior_indices = [
            index
            for node in self.interior_nodes
            for index in range(offsets[node], offsets[node + 1])
        ]
        self._record_indices = tuple(record_indices)
        self._interior_indices = tuple(interior_indices)
        laplacian = sheaf.laplacian()
        self._l_oo = laplacian[np.ix_(record_indices, record_indices)]
        self._l_xx = laplacian[np.ix_(interior_indices, interior_indices)]
        self._l_xo = laplacian[np.ix_(interior_indices, record_indices)]
        self._l_ox = self._l_xo.T
        if record_prior_precision is not None:
            joint_precision = laplacian.copy()
            joint_precision[record_indices, record_indices] += record_prior_precision
            try:
                factor = np.linalg.cholesky(joint_precision)
            except np.linalg.LinAlgError as exc:
                raise ValueError("the record prior does not ground the joint precision") from exc
            self._joint_log_normaliser = (
                0.5 * laplacian.shape[0] * math.log(2.0 * math.pi)
                - float(np.log(np.diag(factor)).sum())
            )
        # RHS-independent quantities derived from the fixed ``_l_xx`` are
        # computed once and cached.  ``np.linalg`` exposes no reusable LU
        # factors, so solve-based results stay on ``np.linalg.solve``.
        self._schur: np.ndarray | None = None
        self._l_xx_inv: np.ndarray | None = None
        self._largest_eigenvalue: float | None = None

    # -- dimensions ---------------------------------------------------------
    @property
    def record_dim(self) -> int:
        return len(self._record_indices)

    @property
    def interior_dim(self) -> int:
        return len(self._interior_indices)

    @property
    def blocks(self) -> dict[str, np.ndarray]:
        return {
            "L_xx": self._l_xx.copy(),
            "L_xo": self._l_xo.copy(),
            "L_ox": self._l_ox.copy(),
            "L_oo": self._l_oo.copy(),
        }

    def _records(self, records: Sequence[float]) -> np.ndarray:
        return _as_vector(records, self.record_dim, "records")

    def _residual(self, state: np.ndarray, records: np.ndarray) -> np.ndarray:
        """The energy gradient at ``state``: ``L_xx @ state + L_xo @ records``."""
        return self._l_xx @ state + self._l_xo @ records

    # -- inference ----------------------------------------------------------
    def harmonic_extension(self, records: Sequence[float]) -> np.ndarray:
        """Return the unique interior minimiser for the given records."""
        observed = self._records(records)
        if self.interior_dim == 0:
            return np.zeros(0)
        return np.linalg.solve(self._l_xx, -(self._l_xo @ observed))

    def schur_complement(self) -> np.ndarray:
        """Return the record Schur complement before adding the prior."""
        if self.interior_dim == 0:
            return self._l_oo.copy()
        if self._schur is None:
            solved = np.linalg.solve(self._l_xx, self._l_xo)
            self._schur = _symmetrise(self._l_oo - self._l_ox @ solved)
        return self._schur.copy()

    def posterior(self, records: Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
        """Return the Gaussian posterior ``(mean, covariance)`` over interiors."""
        mean = self.harmonic_extension(records)
        if self.interior_dim == 0:
            return mean, np.zeros((0, 0))
        if self._l_xx_inv is None:
            self._l_xx_inv = np.linalg.inv(self._l_xx)
        return mean, _symmetrise(self._l_xx_inv)

    def minimum_interface_energy(self, records: Sequence[float]) -> float:
        """Return the minimum sheaf energy ``o^T S o / 2`` at fixed records."""
        observed = self._records(records)
        return 0.5 * float(observed @ self.schur_complement() @ observed)

    def unnormalised_record_potential(self, records: Sequence[float]) -> float:
        """Return ``-ln integral exp(-E(x, o)) dx``, before the record prior."""
        value = self.minimum_interface_energy(records)
        if self.interior_dim:
            _, logdet = np.linalg.slogdet(self._l_xx)
            value += 0.5 * logdet - 0.5 * self.interior_dim * math.log(2.0 * math.pi)
        return value

    def record_surprisal(self, records: Sequence[float]) -> float:
        """Return the negative log density of the proper Gaussian records."""
        if self._record_prior_precision is None:
            raise ValueError("record surprisal requires a record prior")
        observed = self._records(records)
        precision = self.schur_complement() + self._record_prior_precision * np.eye(self.record_dim)
        sign, logdet = np.linalg.slogdet(precision)
        if sign <= 0 or not math.isfinite(logdet):
            raise ValueError("the marginal record precision does not define a Gaussian density")
        normaliser = 0.5 * (self.record_dim * math.log(2.0 * math.pi) - logdet)
        return 0.5 * float(observed @ precision @ observed) + normaliser

    def surprisal(self, records: Sequence[float], *, prior_quadratic: float = 0.0) -> float:
        """Return record surprisal with a prior, or the legacy potential otherwise.

        The legacy ``prior_quadratic`` term shifts an unnormalised potential.
        Use :meth:`record_surprisal` to explicitly require a probability density.
        """
        if self._record_prior_precision is not None:
            if prior_quadratic != 0.0:
                raise ValueError("use the configured record prior for normalised surprisal")
            return self.record_surprisal(records)
        return self.unnormalised_record_potential(records) + float(prior_quadratic)

    def excess_free_energy(
        self,
        mean: Sequence[float],
        covariance: np.ndarray,
        records: Sequence[float],
    ) -> float:
        """Return ``KL[q || posterior]``, the part of free energy that repair removes."""
        posterior_mean, posterior_covariance = self.posterior(records)
        return gaussian_kl(mean, covariance, posterior_mean, posterior_covariance)

    def variational_free_energy(
        self,
        mean: Sequence[float],
        covariance: np.ndarray,
        records: Sequence[float],
        *,
        prior_quadratic: float = 0.0,
    ) -> float:
        """Return ``-ln p(o) + KL[q || posterior]`` with a record prior.

        Without a prior, the legacy bound uses the unnormalised record
        potential in place of ``-ln p(o)``. Both attain equality at the posterior.
        """
        return self.surprisal(records, prior_quadratic=prior_quadratic) + self.excess_free_energy(
            mean, covariance, records
        )

    def energy_minus_entropy(
        self,
        mean: Sequence[float],
        covariance: np.ndarray,
        records: Sequence[float],
    ) -> float:
        """Compute free energy directly for ``q = N(mean, covariance)``.

        This is the first form of the variational free energy, energy minus
        entropy, evaluated without passing through the surprisal or the
        posterior: ``E_q[E] = E(mean, o) + (1/2) tr(L_xx covariance)`` because the
        energy is quadratic, and ``H[q] = (1/2) logdet(2 pi e covariance)``.  It
        is an independent route to the number :meth:`variational_free_energy`
        computes as surprisal plus divergence, so comparing the two tests the
        identity of the free-energy decomposition rather than restating it.
        With a record prior, the energy includes that prior and the joint
        log normaliser computed from the full precision.
        """
        mean = _as_vector(mean, self.interior_dim, "mean")
        observed = self._records(records)
        covariance, logdet = _covariance_and_logdet(
            covariance, self.interior_dim, "recognition covariance"
        )
        energy = 0.5 * float(
            mean @ self._l_xx @ mean
            + 2.0 * mean @ self._l_xo @ observed
            + observed @ self._l_oo @ observed
        )
        if self._record_prior_precision is not None:
            energy += (
                0.5 * self._record_prior_precision * float(observed @ observed)
                + self._joint_log_normaliser
            )
        if self.interior_dim == 0:
            return energy
        expected_energy = energy + 0.5 * float(np.trace(self._l_xx @ covariance))
        entropy = 0.5 * (self.interior_dim * (1.0 + math.log(2.0 * math.pi)) + logdet)
        return expected_energy - entropy

    def constrained_minimiser(self, records: Sequence[float]) -> np.ndarray:
        """Minimise ``(1/2) z^T L z`` subject to ``z_o = records`` through its KKT system.

        For Gaussians of a fixed precision ``L`` the relative entropy of a mean
        ``z`` to the reference at the origin is ``(1/2) z^T L z``, so this is the
        constrained relative-entropy minimisation read in the mean.  It is solved
        as one saddle-point system over all coordinates and a Lagrange
        multiplier per record -- a different linear system from the interior
        solve :meth:`harmonic_extension` uses -- and the interior block of its
        solution is returned.
        """
        observed = self._records(records)
        laplacian = self.sheaf.laplacian()
        total = laplacian.shape[0]
        selector = np.zeros((self.record_dim, total))
        for row, index in enumerate(self._record_indices):
            selector[row, index] = 1.0
        kkt = np.block(
            [
                [laplacian, selector.T],
                [selector, np.zeros((self.record_dim, self.record_dim))],
            ]
        )
        rhs = np.concatenate([np.zeros(total), observed])
        solution = np.linalg.solve(kkt, rhs)
        return solution[list(self._interior_indices)]

    def perceive(
        self,
        records: Sequence[float],
        *,
        initial: Sequence[float] | None = None,
        tolerance: float = 1e-10,
        max_steps: int = 10_000,
    ) -> tuple[np.ndarray, float, int]:
        """Run clamped gradient flow on the sheaf energy.

        The step size is ``1 / lambda_max(L_xx)``, which is stable for the
        gradient system.  Returns ``(interior, residual_norm, iterations)``.
        """
        observed = self._records(records)
        if self.interior_dim == 0:
            return np.zeros(0), 0.0, 0
        if initial is None:
            state = np.zeros(self.interior_dim)
        else:
            state = _as_vector(initial, self.interior_dim, "initial")
        if self._largest_eigenvalue is None:
            self._largest_eigenvalue = float(np.linalg.eigvalsh(self._l_xx)[-1])
        step = 1.0 / self._largest_eigenvalue if self._largest_eigenvalue > 0 else 1.0
        residual = self._residual(state, observed)
        iterations = 0
        while float(np.linalg.norm(residual)) > tolerance and iterations < max_steps:
            state = state - step * residual
            residual = self._residual(state, observed)
            iterations += 1
        return state, float(np.linalg.norm(residual)), iterations

    # -- surprises ----------------------------------------------------------
    def numerical_surprise(self, old_records: Sequence[float], new_records: Sequence[float]) -> float:
        """The norm of the update residual at the old rest under new records."""
        old_rest = self.harmonic_extension(old_records)
        gradient = self._residual(old_rest, self._records(new_records))
        return float(np.linalg.norm(gradient))

    def bayesian_surprise(self, old_records: Sequence[float], new_records: Sequence[float]) -> float:
        """The KL between successive posteriors, ``0.5 g^T L_xx^{-1} g``."""
        if self.interior_dim == 0:
            return 0.0
        old_rest = self.harmonic_extension(old_records)
        gradient = self._residual(old_rest, self._records(new_records))
        return 0.5 * float(gradient @ np.linalg.solve(self._l_xx, gradient))

    def record_injectivity(self) -> bool:
        """Whether the records separate global sections."""
        return self.sheaf.record_injectivity(self.record_nodes)


def three_surprises(
    model: GaussianSheafModel,
    old_records: Sequence[float],
    new_records: Sequence[float],
) -> SurpriseTriple:
    """Return numerical surprise, Bayesian surprise and surprisal.

    The first two are read at the old rest after the records change; the third
    is a property of the new records alone.
    """
    return SurpriseTriple(
        numerical=model.numerical_surprise(old_records, new_records),
        bayesian=model.bayesian_surprise(old_records, new_records),
        surprisal=model.surprisal(new_records),
    )


def demo_chain_model(
    node_count: int = 4, *, record_prior_precision: float | None = None
) -> GaussianSheafModel:
    """A scalar chain with both endpoints observed; useful for examples/tests."""
    if node_count < 3:
        raise ValueError("a demo chain needs at least three nodes")
    sheaf = linear_chain_sheaf(node_count)
    return GaussianSheafModel(
        sheaf, (0, node_count - 1), record_prior_precision=record_prior_precision
    )


def asymmetric_seam_energy(disagreement: np.ndarray) -> np.ndarray:
    """An oriented spring with twice the stiffness for positive disagreement."""
    return np.where(disagreement <= 0.0, 0.5, 1.0) * disagreement**2


@lru_cache(maxsize=1)
def _carrier_grounded_forest_census() -> dict[str, object]:
    carrier = icosahedral_carrier()
    adjacency = carrier.adjacency_matrix().astype(int).tolist()
    laplacian = laplacian_matrix(adjacency)
    distances = graph_distances(adjacency)
    resistance_by_distance = {distance: value for distance, value, _ in resistance_ladder(adjacency)}
    # Compute the tree factor from the closed carrier spectrum, independently
    # of the principal-determinant routine. The conjugate eigenvalues pair to
    # 5² - 5, with three copies of each eigenvalue and five copies of 6.
    trees = Fraction((5**2 - 5) ** 3 * 6**5, carrier.port_count)
    pair_determinants: dict[tuple[int, ...], int] = {}
    triple_determinants: dict[tuple[int, ...], int] = {}
    layouts = {}
    matches_contrast_covariance = True

    for size, determinants in ((2, pair_determinants), (3, triple_determinants)):
        grouped: dict[tuple[int, ...], list[tuple[int, ...]]] = defaultdict(list)
        for records in combinations(range(carrier.port_count), size):
            steps = tuple(distances[a][b] for a, b in combinations(records, 2))
            grouped[tuple(sorted(steps))].append(records)
            determinant = laplacian_determinant(laplacian, grounded=records)
            determinants[records] = determinant

            # Jacobi's complementary-minor identity uses the contrast
            # covariance on the additional records after grounding the first.
            contrast_determinant = resistance_by_distance[steps[0]]
            if size == 3:
                r_ab, r_ac, r_bc = (resistance_by_distance[step] for step in steps)
                covariance = (r_ab + r_ac - r_bc) / 2
                contrast_determinant = r_ab * r_ac - covariance**2
            matches_contrast_covariance &= determinant == trees * contrast_determinant

        rows = []
        for pattern, record_sets in sorted(grouped.items()):
            (determinant,) = {determinants[records] for records in record_sets}
            rows.append((pattern, len(record_sets), determinant, record_sets[0]))
        layouts[size] = tuple(rows)

    best_pair = max(pair_determinants.values())
    best_pairs = {
        records for records, determinant in pair_determinants.items() if determinant == best_pair
    }
    greedy_extensions = {
        records: determinant
        for records, determinant in triple_determinants.items()
        if any(pair in best_pairs for pair in combinations(records, 2))
    }
    best_greedy = max(greedy_extensions.values())
    best_triple = max(triple_determinants.values())
    worst_triple = min(triple_determinants.values())
    return {
        "pair_layouts": layouts[2],
        "triple_layouts": layouts[3],
        "pairs_checked": len(pair_determinants),
        "triples_checked": len(triple_determinants),
        "matches_contrast_covariance": matches_contrast_covariance,
        "best_pair_example": min(best_pairs),
        "best_pair_determinant": best_pair,
        "best_triple_example": next(
            records for records, determinant in triple_determinants.items() if determinant == best_triple
        ),
        "best_triple_determinant": best_triple,
        "worst_triple_example": next(
            records for records, determinant in triple_determinants.items() if determinant == worst_triple
        ),
        "worst_triple_determinant": worst_triple,
        "greedy_triple_example": next(
            records for records, determinant in greedy_extensions.items() if determinant == best_greedy
        ),
        "greedy_triple_determinants": tuple(sorted(Counter(greedy_extensions.values()).items())),
        "optimum_to_greedy_ratio": str(Fraction(best_triple, best_greedy)),
    }


def carrier_grounded_forest_census() -> dict[str, object]:
    """Compare two- and three-record precision minors on the unit carrier.

    A layout row gives the sorted pairwise graph distances, the number of
    record sets, the grounded determinant, and the first set in port order.
    The contrast-covariance route uses exact effective resistances from a
    separate grounded inverse and the tree factor from the closed spectrum.
    The greedy extensions start from any pair that maximises the determinant.
    """
    return _carrier_grounded_forest_census().copy()
