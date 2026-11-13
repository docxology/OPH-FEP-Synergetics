"""Variational free energy as precision-weighted interface disagreement.

For a finite linear-Gaussian cellular sheaf with energy

.. math:: E(x, o) = \\tfrac12 z^T L z, \\qquad z = (x, o),

the free-energy dictionary of the accompanying review is exact:

* perception with records held fixed is sheaf diffusion and converges to the
  harmonic extension ``x_hat(o) = -L_xx^{-1} L_xo o``;
* the surprisal ``-ln p(o)`` is ``\\tfrac12 o^T S o`` plus a Gaussian
  normaliser, where ``S = L_oo - L_ox L_xx^{-1} L_xo`` is the Schur complement;
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
from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .sheaf import CellularSheaf, linear_chain_sheaf


def _as_vector(value: Sequence[float], length: int, name: str) -> np.ndarray:
    vector = np.asarray(value, dtype=float)
    if vector.shape != (length,):
        raise ValueError(f"{name} must have shape ({length},)")
    return vector


def _symmetrise(matrix: np.ndarray) -> np.ndarray:
    return 0.5 * (matrix + matrix.T)


def gaussian_kl(
    mu_0: Sequence[float],
    sigma_0: np.ndarray,
    mu_1: Sequence[float],
    sigma_1: np.ndarray,
) -> float:
    """The KL divergence ``KL[N(mu_0, sigma_0) || N(mu_1, sigma_1)]``."""
    mu_0 = np.asarray(mu_0, dtype=float)
    mu_1 = np.asarray(mu_1, dtype=float)
    sigma_0 = np.asarray(sigma_0, dtype=float)
    sigma_1 = np.asarray(sigma_1, dtype=float)
    if mu_0.shape != mu_1.shape:
        raise ValueError("Gaussian means must share a shape")
    dimension = mu_0.shape[0]
    if sigma_0.shape != (dimension, dimension) or sigma_1.shape != (dimension, dimension):
        raise ValueError("Gaussian covariances must be square and match the means")
    sigma_1 = _symmetrise(sigma_1)
    sigma_0 = _symmetrise(sigma_0)
    sign_1, logdet_1 = np.linalg.slogdet(sigma_1)
    sign_0, logdet_0 = np.linalg.slogdet(sigma_0)
    if sign_1 <= 0 or sign_0 <= 0:
        raise ValueError("Gaussian covariances must be positive definite")
    difference = mu_1 - mu_0
    solved = np.linalg.solve(sigma_1, sigma_0)
    trace_term = float(np.trace(solved))
    quadratic = float(difference @ np.linalg.solve(sigma_1, difference))
    return 0.5 * (trace_term + quadratic - dimension + logdet_1 - logdet_0)


@dataclass(frozen=True)
class SurpriseTriple:
    """Three quantities called surprise, kept numerically separate."""

    numerical: float
    bayesian: float
    surprisal: float


class GaussianSheafModel:
    """A linear-Gaussian sheaf model with record (boundary) and interior nodes."""

    def __init__(self, sheaf: CellularSheaf, record_nodes: Sequence[int]) -> None:
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
        """Return the marginal precision of the records."""
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

    def surprisal(self, records: Sequence[float], *, prior_quadratic: float = 0.0) -> float:
        """Return ``-ln p(o)`` for the Gaussian sheaf model, up to a constant in ``o``.

        Theorem 7 item 3 of the source gives
        ``-ln p(o) = (1/2) o^T S o + rho(o) + const`` for
        ``p(x, o) ∝ exp(-E(x, o) - rho(o))``, with ``S`` the Schur complement.
        The returned value is the whole ``o``-dependent part,
        ``(1/2) o^T S o + rho(o)``, together with the Gaussian-integral
        prefactor ``(1/2) logdet(L_xx) - (d_x / 2) log(2 pi)``, which depends on
        the model and not on ``o``.  The absolute level is therefore
        gauge-dependent and may be negative.  Differences are exact surprisal
        differences: ``surprisal(o) - surprisal(o') = -ln p(o) + ln p(o')``.
        """
        observed = self._records(records)
        marginal = self.schur_complement()
        value = 0.5 * float(observed @ marginal @ observed) + float(prior_quadratic)
        if self.interior_dim:
            _, logdet = np.linalg.slogdet(self._l_xx)
            value += 0.5 * logdet - 0.5 * self.interior_dim * math.log(2.0 * math.pi)
        return value

    def excess_free_energy(
        self,
        mean: Sequence[float],
        covariance: np.ndarray,
        records: Sequence[float],
    ) -> float:
        """Return ``KL[q || posterior]``, the part of free energy that repair removes."""
        mean = _as_vector(mean, self.interior_dim, "mean")
        covariance = np.asarray(covariance, dtype=float)
        if covariance.shape != (self.interior_dim, self.interior_dim):
            raise ValueError("covariance has the wrong shape")
        if self.interior_dim == 0:
            return 0.0
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
        """Return ``-ln p(o) + KL[q || posterior]``.

        This is the FEP bound written for the sheaf model: it is an upper bound
        on surprisal with equality exactly at the posterior.
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
        """Return ``E_q[E(x, o)] - H[q]`` computed directly, for ``q = N(mean, covariance)``.

        This is the first form of the variational free energy, energy minus
        entropy, evaluated without passing through the surprisal or the
        posterior: ``E_q[E] = E(mean, o) + (1/2) tr(L_xx covariance)`` because the
        energy is quadratic, and ``H[q] = (1/2) logdet(2 pi e covariance)``.  It
        is an independent route to the number :meth:`variational_free_energy`
        computes as surprisal plus divergence, so comparing the two tests the
        identity of the free-energy decomposition rather than restating it.
        """
        mean = _as_vector(mean, self.interior_dim, "mean")
        observed = self._records(records)
        covariance = _symmetrise(np.asarray(covariance, dtype=float))
        energy = 0.5 * float(
            mean @ self._l_xx @ mean
            + 2.0 * mean @ self._l_xo @ observed
            + observed @ self._l_oo @ observed
        )
        if self.interior_dim == 0:
            return energy
        sign, logdet = np.linalg.slogdet(covariance)
        if sign <= 0:
            raise ValueError("the recognition covariance must be positive definite")
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
        target = self.harmonic_extension(records)
        if self.interior_dim == 0:
            return target, 0.0, 0
        if initial is None:
            state = np.zeros(self.interior_dim)
        else:
            state = _as_vector(initial, self.interior_dim, "initial")
        if self._largest_eigenvalue is None:
            self._largest_eigenvalue = float(np.linalg.eigvalsh(self._l_xx)[-1])
        step = 1.0 / self._largest_eigenvalue if self._largest_eigenvalue > 0 else 1.0
        residual = self._residual(state, self._records(records))
        iterations = 0
        while float(np.linalg.norm(residual)) > tolerance and iterations < max_steps:
            state = state - step * residual
            residual = self._residual(state, self._records(records))
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


def demo_chain_model(node_count: int = 4) -> GaussianSheafModel:
    """A scalar chain with both endpoints observed; useful for examples/tests."""
    if node_count < 3:
        raise ValueError("a demo chain needs at least three nodes")
    sheaf = linear_chain_sheaf(node_count)
    return GaussianSheafModel(sheaf, record_nodes=(0, node_count - 1))
