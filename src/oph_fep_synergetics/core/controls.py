"""Exact controls that locate isomorphisms, disagreements and complements.

Each function is a small, self-contained mathematical experiment drawn from
the formal comparison of OPH, the FEP, and Synergetic geometry.  The controls
are deliberately exact where the source result is exact (integer arithmetic,
rational arithmetic, explicit matrices) and numerically checked where the
source statement is analytic.  They are the project's negative controls: a
bridge that cannot survive them is not a bridge.

The controls fall into seven groups:

1. **Descent versus drift** -- the native repair law strictly lowers the
   quadratic energy on broken seams, but admits zero-descent moves at
   mismatch one.  This is OPH's PR-01 obstruction, not a defect of this code.
2. **Law under-determination** -- compatible states do not select a law;
   D-stability and the solenoidal family exhibit the freedom.
3. **State sufficiency** -- a projection that drops a record fails the
   intervention test by a computable margin.
4. **Interface reduction** -- eliminating a mediator can create a genuinely
   three-way boundary term that no pairwise model represents below a stated
   precision.
5. **Existence without reachability** -- the clipped two-node map has exactly one
   equilibrium and a periodic orbit; the relaxation orbit's minimum radius and
   period are reported numerically, not proved.
6. **Carrier invariants of the native law** -- the transition multiplicities of
   the seventy-eight-state shell, and the four-cycle local minimum of an
   independently defined discrepancy.
7. **Free-energy form and its boundaries** -- the Lyapunov form of a stable
   law, the directed-gradient boundary, the two-settle identity, the surprise
   sandwich, and the contraction displacement bounds.
"""

from __future__ import annotations

import itertools
import math
from fractions import Fraction
from typing import Sequence

import numpy as np

from .repair import (
    energy,
    mismatch,
    native_repair,
    predict_energy_decrease,
    shell_gap,
)


# ---------------------------------------------------------------------------
# 1. Descent on the gated domain, drift at mismatch one
# ---------------------------------------------------------------------------


def gated_descent_census() -> dict[str, int]:
    """Exhaustively classify every native move on the 78-state shell.

    Returns counts of strict-descent moves (broken seams), zero-descent moves
    (mismatch one), and no-op moves.  The quadratic decrease is checked against
    the closed form ``(d^2 - (|d| mod 2)) / 4``.
    """
    from .carrier import icosahedral_carrier

    carrier = icosahedral_carrier()
    states = carrier.shell_configurations(total=2)
    strict = zero = noop = 0
    for state in states:
        for a, b in carrier.edges:
            difference = mismatch(state, a, b)
            for coin in (0, 1):
                step = native_repair(state, (a, b), coin)
                if step.after == step.before:
                    noop += 1
                    continue
                predicted = predict_energy_decrease(state, a, b)
                observed = step.energy_decrease
                if not math.isclose(observed, predicted, rel_tol=0.0, abs_tol=1e-12):
                    raise AssertionError(
                        f"energy decrease mismatch at {state} seam {(a, b)} coin {coin}"
                    )
                if observed > 0:
                    strict += 1
                elif abs(difference) == 1:
                    zero += 1
                else:  # pragma: no cover - only |d| in {0,1,2} on this shell
                    raise AssertionError("unexpected zero-descent move")
    return {"strict": strict, "zero": zero, "noop": noop, "states": len(states)}


def mismatch_one_moves(
    loads: Sequence[int],
    a: int,
    b: int,
) -> tuple[RepairStep, RepairStep]:
    """The two outcomes of the native law at mismatch one, for coins ``0`` and ``1``.

    At ``|d| = 1`` the law either moves the odd unit (``|j| = 1``) or waits
    (``j = 0``), and the fair coin decides which.  Neither outcome changes the
    quadratic energy, because ``(d^2 - (|d| mod 2)) / 4`` vanishes at
    ``|d| = 1``.  This pair is the executable content of the "waits or swaps at
    mismatch one" reading of OPH's PR-01 obstruction.
    """
    if abs(mismatch(loads, a, b)) != 1:
        raise ValueError("this control is defined at mismatch one")
    return (native_repair(loads, (a, b), 0), native_repair(loads, (a, b), 1))


def mismatch_one_coin(loads: Sequence[int], a: int, b: int) -> int:
    """The coin that places the odd unit at the emptier of the two ports.

    The law's current is ``j = (d - (2c - 1)(|d| mod 2)) // 2``, so at ``|d| = 1``
    the coin that moves the unit from the fuller port to the emptier one depends
    on the seam's orientation: ``0`` when ``d > 0`` and ``1`` when ``d < 0``.
    A coin drawn independently of the orientation therefore waits half the time.
    """
    difference = mismatch(loads, a, b)
    if abs(difference) != 1:
        raise ValueError("this control is defined at mismatch one")
    return 0 if difference > 0 else 1


def mismatch_one_excursion(steps: int = 6) -> tuple[tuple[int, ...], ...]:
    """Walk one unit around the carrier with zero quadratic-energy change.

    Starts from the one-unit shell and repairs consecutive seams at mismatch
    one, always taking the coin that moves the odd unit to the emptier port,
    never stepping straight back, and preferring a port the walk has not yet
    visited.  Every visited state carries energy ``1/2``, so a schedule that
    applies the local rule without the broken-seam gate keeps moving while the
    descent measure stays flat.  This is the executable content of OPH's PR-01
    obstruction: accepted repairs need not lower a measure.
    """
    if steps < 1:
        raise ValueError("steps must be positive")
    from .carrier import icosahedral_carrier

    carrier = icosahedral_carrier()
    state = [0] * carrier.port_count
    state[0] = 1
    visited = [tuple(state)]
    port = 0
    previous: int | None = None
    seen = {0}
    for _ in range(steps):
        fresh = [
            candidate
            for candidate in carrier.neighbours(port)
            if candidate not in seen and candidate != previous
        ]
        neighbour = (
            fresh[0]
            if fresh
            else next(
                candidate for candidate in carrier.neighbours(port) if candidate != previous
            )
        )
        step = native_repair(
            tuple(state), (port, neighbour), mismatch_one_coin(state, port, neighbour)
        )
        state = list(step.after)
        if state == list(step.before):  # pragma: no cover - the adapted coin moves
            raise AssertionError("the orientation-adapted coin must move the odd unit")
        visited.append(tuple(state))
        if energy(step.after) != energy(step.before):  # pragma: no cover - invariant
            raise AssertionError("mismatch-one move changed the quadratic energy")
        seen.add(neighbour)
        previous, port = port, neighbour
    return tuple(visited)


def strict_descent_on_broken_seams() -> bool:
    """Whether every broken-seam repair strictly lowers the quadratic energy."""
    census = gated_descent_census()
    return census["strict"] > 0 and census["zero"] > 0


# ---------------------------------------------------------------------------
# 2. Compatible states do not select a law
# ---------------------------------------------------------------------------


D_STABILITY_A = np.array(
    [
        [-1.0, 2.0, 0.0],
        [2.0, -1.0, 3.0],
        [0.0, -3.0, -1.0],
    ]
)


def d_stability_matrix(k: float) -> np.ndarray:
    """Return ``D_k A`` for ``D_k = diag(k, 1, 1)``."""
    if k <= 0:
        raise ValueError("local gain k must be positive")
    return np.diag([k, 1.0, 1.0]) @ D_STABILITY_A


def d_stability_charpoly(k: float) -> np.ndarray:
    """Return the coefficients of ``det(z I - D_k A)``, highest power first."""
    return np.poly(d_stability_matrix(k))


def d_stability_expected_charpoly(k: float) -> tuple[float, ...]:
    """The closed form ``z^3 + (k+2) z^2 + (10-2k) z + 6k``."""
    return (1.0, k + 2.0, 10.0 - 2.0 * k, 6.0 * k)


def routh_first_column(k: float) -> tuple[float, float, float, float]:
    """The first column of the Routh array for the cubic characteristic polynomial."""
    a2 = k + 2.0
    a1 = 10.0 - 2.0 * k
    a0 = 6.0 * k
    if a2 == 0:
        raise ValueError("Routh array is singular at this k")
    return (1.0, a2, (a2 * a1 - a0) / a2, a0)


def d_stability_eigenvalues(k: float) -> np.ndarray:
    return np.linalg.eigvals(d_stability_matrix(k))


def d_stability_unstable(k: float, *, tolerance: float = 1e-9) -> bool:
    """Whether the positive diagonal scaling destabilises the common rest."""
    return bool(np.any(np.real(d_stability_eigenvalues(k)) > tolerance))


def logistic_embedding_jacobian() -> np.ndarray:
    """The Jacobian at zero of the bounded-activation embedding of ``A``.

    With ``W = 4(A + I)`` and ``f(v) = -v + W sigma(v) + b`` for the logistic
    ``sigma``, the bias ``b = -W (1/2, 1/2, 1/2)^T`` makes ``f(0) = 0``, and
    ``sigma'(0) = 1/4`` gives ``Df(0) = -I + W / 4 = A``.  The D-stability
    control therefore occurs inside a bounded-activation neural family, with
    zero self-weights.
    """
    weight = 4.0 * (D_STABILITY_A + np.eye(3))
    return -np.eye(3) + 0.25 * weight


def logistic_embedding_checks() -> dict[str, float]:
    weight = 4.0 * (D_STABILITY_A + np.eye(3))
    bias = -weight @ np.full(3, 0.5)
    zero = np.zeros(3)
    sigma_zero = 0.5
    value = -zero + weight @ np.full(3, sigma_zero) + bias
    return {
        "f_zero_norm": float(np.linalg.norm(value)),
        "jacobian_error": float(np.linalg.norm(logistic_embedding_jacobian() - D_STABILITY_A)),
    }


def solenoidal_check(q: float, *, dimension: int = 2) -> dict[str, float]:
    """One stationary density, a family of laws indexed by the solenoidal part.

    For ``H = Gamma = I`` and ``Q = [[0, q], [-q, 0]]``, the diffusion
    ``dx = -(Gamma - Q) H x dt + sqrt(2 Gamma) dW`` has stationary law
    ``N(0, H^{-1})`` for every ``q``, spectrum ``-1 +/- i q``, and entropy
    production ``2 q^2``.
    """
    if dimension != 2:
        raise ValueError("this control is two-dimensional")
    h = np.eye(2)
    gamma = np.eye(2)
    q_matrix = np.array([[0.0, q], [-q, 0.0]])
    drift = -(gamma - q_matrix) @ h
    covariance = np.linalg.inv(h)
    residual = drift @ covariance + covariance @ drift.T + 2.0 * gamma
    eigenvalues = np.linalg.eigvals(drift)
    entropy_production = float(np.trace(q_matrix.T @ np.linalg.inv(gamma) @ q_matrix @ h))
    return {
        "stationary_residual": float(np.linalg.norm(residual)),
        "spectral_real": float(np.real(eigenvalues[0])),
        "spectral_imag": float(abs(np.imag(eigenvalues[0]))),
        "entropy_production": entropy_production,
    }


def blanket_dynamics_check(q: float) -> dict[str, float | bool]:
    """A stationary blanket that is not a dynamical interface.

    Orders the states as external ``eta``, blanket ``b`` and internal ``mu``.
    The precision has a zero ``(mu, eta)`` entry, so the stationary Gaussian
    makes ``mu`` and ``eta`` conditionally independent given ``b``.  The drift
    of ``mu`` is ``2 q eta + (q - 1) b - 2 mu``, so two states with equal
    ``(b, mu)`` and different ``eta`` have different internal motion whenever
    ``q != 0``.
    """
    h = np.array(
        [
            [2.0, 1.0, 0.0],
            [1.0, 2.0, 1.0],
            [0.0, 1.0, 2.0],
        ]
    )
    q_matrix = np.zeros((3, 3))
    q_matrix[0, 2] = -q
    q_matrix[2, 0] = q
    gamma = np.eye(3)
    drift = -(gamma - q_matrix) @ h
    row = drift[2]
    conditional = abs(h[0, 2]) <= 1e-12
    difference = 2.0 * abs(q)
    return {
        "conditional_independence": bool(conditional),
        "mu_drift_external": float(row[0]),
        "mu_drift_blanket": float(row[1]),
        "mu_drift_internal": float(row[2]),
        "cross_drift_difference": float(difference),
    }


# ---------------------------------------------------------------------------
# 3. State sufficiency and the intervention test
# ---------------------------------------------------------------------------


def reduction_control() -> dict[str, float]:
    """Two histories with the same observed state and different interventions.

    ``x`` is observed, ``h`` is a retained hidden value, and
    ``F_u(x, h) = ((1-u)x + u h, h)`` for ``u in {0, 1}``.  At ``u=0`` both
    ``(0, 1)`` and ``(0, -1)`` are fixed with the same observed state; the
    intervention ``u=1`` sends their observed values to ``1`` and ``-1``.
    """
    def transition(state: tuple[float, float], u: int) -> tuple[float, float]:
        x, h = state
        return ((1.0 - u) * x + u * h, h)

    first = (0.0, 1.0)
    second = (0.0, -1.0)
    same_observation = first[0] == second[0]
    observed_after = (
        transition(first, 1)[0],
        transition(second, 1)[0],
    )
    distance = abs(observed_after[0] - observed_after[1])
    return {
        "same_observation": 1.0 if same_observation else 0.0,
        "first_after": observed_after[0],
        "second_after": observed_after[1],
        "response_distance": distance,
        "worst_case_error_lower_bound": distance / 2.0,
    }


def affine_path_contrast_exact(*, steps: int = 12) -> Fraction:
    """The exact schedule-induced contrast difference.

    ``R_beta(x) = 0.4 x + 0.3 + beta`` is iterated from the free root ``1/2`` of
    ``R_0``.  Full iteration and half-damped iteration have the same globally
    attracting root for every ``beta``, so they preserve every answer, yet they
    change the centered contrast that a two-settle learner reads.  After twelve
    steps from the free root the difference is exactly
    ``921633999 / 40000000000`` (= 0.023040849975).
    """
    if steps < 1:
        raise ValueError("steps must be positive")
    beta = Fraction(1, 10)

    def endpoint(sign: int, damped: bool) -> Fraction:
        x = Fraction(1, 2)
        for _ in range(steps):
            target = Fraction(2, 5) * x + Fraction(3, 10) + sign * beta
            x = (x + target) / 2 if damped else target
        return x

    def contrast(damped: bool) -> Fraction:
        return (endpoint(1, damped) - endpoint(-1, damped)) / (2 * beta)

    return abs(contrast(True) - contrast(False))


def affine_path_contrast_float(*, steps: int = 12) -> float:
    return float(affine_path_contrast_exact(steps=steps))


# ---------------------------------------------------------------------------
# 4. Interface reduction
# ---------------------------------------------------------------------------


def mediator_schur(epsilon: float) -> np.ndarray:
    """The marginal record precision after eliminating a two-dimensional mediator.

    The mediator ``w`` meets three scalar records through ``x - w1``,
    ``y - w2`` and ``z - w1 - w2``; the records carry prior precision
    ``epsilon``.  The result is ``(1/3) v v^T + epsilon I`` with
    ``v = (1, 1, -1)``.
    """
    if epsilon < 0:
        raise ValueError("epsilon must be nonnegative")
    l_xx = np.array([[2.0, 1.0], [1.0, 2.0]])
    l_xo = np.array([[-1.0, 0.0, -1.0], [0.0, -1.0, -1.0]])
    l_oo = np.eye(3) * (1.0 + epsilon)
    return l_oo - l_xo.T @ np.linalg.solve(l_xx, l_xo)


def mediator_expected(epsilon: float) -> np.ndarray:
    vector = np.array([1.0, 1.0, -1.0])
    return (1.0 / 3.0) * np.outer(vector, vector) + epsilon * np.eye(3)


def mediator_pairwise_margin(epsilon: float) -> float:
    """The slack in the trace inequality ``1 + 3 epsilon >= 2``, i.e. ``3 epsilon - 1``.

    A pairwise representation gives each of the three pairs a positive
    semidefinite block with off-diagonal magnitude ``1/3``, so its two diagonal
    entries satisfy ``a_i a_j >= 1/9`` and therefore ``a_i + a_j >= 2/3``.
    Summing the three pairs puts at least two on the pairwise trace, against the
    mediator's ``1 + 3 epsilon``.  The trace inequality holds exactly when
    ``epsilon >= 1/3``, where this slack vanishes.
    """
    return 3.0 * epsilon - 1.0


def mediator_pairwise_decomposition(epsilon: float) -> np.ndarray | None:
    """Return a pairwise energy matrix reconstructing the mediator, if one exists."""
    if epsilon < 1.0 / 3.0 - 1e-12:
        return None
    diagonal = epsilon - 1.0 / 3.0
    result = np.zeros((3, 3))
    for i, j in itertools.combinations(range(3), 2):
        sign = -1.0 if {i, j} == {0, 2} or {i, j} == {1, 2} else 1.0
        block = np.array([[1.0 / 3.0, sign / 3.0], [sign / 3.0, 1.0 / 3.0]])
        result[np.ix_([i, j], [i, j])] += block
    result += diagonal * np.eye(3)
    return result


# ---------------------------------------------------------------------------
# 5. Existence without reachability
# ---------------------------------------------------------------------------


def _clip(value: np.ndarray) -> np.ndarray:
    return np.clip(value, -1.0, 1.0)


def two_node_map(state: Sequence[float]) -> np.ndarray:
    """The clipped two-node repair map from the review's Proposition 15."""
    x, y = (float(v) for v in state)
    return _clip(np.array([2.0 * x - 2.0 * y, 2.0 * x + 2.0 * y]))


def two_node_sync_orbit(
    *,
    steps: int = 500,
    start: tuple[float, float] = (0.1, 0.2),
) -> np.ndarray:
    """Iterate the synchronous clipped map and return the orbit."""
    if steps < 1:
        raise ValueError("steps must be positive")
    state = np.array(start, dtype=float)
    orbit = [state.copy()]
    for _ in range(steps):
        state = two_node_map(state)
        orbit.append(state.copy())
    return np.asarray(orbit)


def two_node_relaxation_jacobian() -> np.ndarray:
    """The Jacobian at the origin of ``xdot = R(x) - x`` (eigenvalues ``1 +/- 2i``)."""
    return np.array([[1.0, -2.0], [2.0, 1.0]])


def two_node_orbit_min_radius(*, steps: int = 500, skip: int = 100) -> float:
    orbit = two_node_sync_orbit(steps=steps)
    return float(np.min(np.linalg.norm(orbit[skip:], axis=1)))


def two_node_relaxation_orbit(
    *,
    start: tuple[float, float] = (0.1, 0.2),
    dt: float = 1e-3,
    steps: int = 100_000,
    skip: int = 20_000,
) -> dict[str, float]:
    """The periodic orbit of the two-node relaxation ``xdot = R(x) - x``.

    Proposition 15's clipped map has the origin as its only equilibrium, and the
    relaxation Jacobian there has eigenvalues ``1 +/- 2i``, so the origin is a
    source and every nonconstant trajectory approaches a periodic orbit.  This
    measures that orbit: the source reports a minimum radius of at least ``0.84``
    and a period of ``7.58`` relaxation time units for this control.  The values
    are numerical, not proved.

    The synchronous clipped map is a different object; :func:`two_node_sync_orbit`
    and :func:`two_node_orbit_min_radius` measure that one.
    """
    if dt <= 0.0 or steps <= 1 or skip < 0 or skip >= steps:
        raise ValueError("invalid integration parameters")

    def field(a: float, b: float) -> tuple[float, float]:
        first = 2.0 * a - 2.0 * b
        second = 2.0 * a + 2.0 * b
        first = -1.0 if first < -1.0 else (1.0 if first > 1.0 else first)
        second = -1.0 if second < -1.0 else (1.0 if second > 1.0 else second)
        return first - a, second - b

    x, y = float(start[0]), float(start[1])
    stride = max(1, steps // 20_000)
    samples: list[tuple[float, float]] = []
    for index in range(steps):
        k1 = field(x, y)
        k2 = field(x + 0.5 * dt * k1[0], y + 0.5 * dt * k1[1])
        k3 = field(x + 0.5 * dt * k2[0], y + 0.5 * dt * k2[1])
        k4 = field(x + dt * k3[0], y + dt * k3[1])
        x += (dt / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
        y += (dt / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
        if index >= skip and index % stride == 0:
            samples.append((x, y))
    radii = [math.hypot(a, b) for a, b in samples]
    anchor = samples[0]
    # Skip one relaxation time unit before looking for the first return, so the
    # search reports the fundamental period rather than a later multiple of it.
    warmup = min(max(1, int(1.0 / (stride * dt))), len(samples) - 1)
    period = float("nan")
    for index in range(warmup, len(samples)):
        a, b = samples[index]
        if math.hypot(a - anchor[0], b - anchor[1]) < 1e-3:
            period = float(index * stride * dt)
            break
    return {
        "min_radius": min(radii),
        "max_radius": max(radii),
        "period": period,
    }


# ---------------------------------------------------------------------------
# 6. Carrier invariants of the native law
# ---------------------------------------------------------------------------


def path_multiplicity_census() -> dict[int, int]:
    """Count native transitions by multiplicity over the seventy-eight states.

    Thirty seams and two coins give sixty transitions out of every
    configuration, so the uniform schedule reaches a target with probability
    ``m / 60`` where ``m`` is that target's multiplicity.  Over the whole shell
    the multiplicities are exactly ``{1, 2, 50, 52}``: ``1`` for the eight
    neighbour moves out of a binary state, ``2`` for the ten moves out of a
    singly-loaded state (both coins place the unit on the same seam), ``52`` for
    the no-op self-loop of a binary state and ``50`` for that of a singly-loaded
    state.  The conditional surprisal of a path is therefore a sum of
    ``log(60 / m)`` terms drawn from this set.
    """
    from .carrier import icosahedral_carrier

    carrier = icosahedral_carrier()
    census: dict[tuple[tuple[int, ...], tuple[int, ...]], int] = {}
    for state in carrier.shell_configurations(total=2):
        for a, b in carrier.edges:
            for coin in (0, 1):
                step = native_repair(state, (a, b), coin)
                key = (step.before, step.after)
                census[key] = census.get(key, 0) + 1
    counts: dict[int, int] = {}
    for multiplicity in census.values():
        counts[multiplicity] = counts.get(multiplicity, 0) + 1
    return dict(sorted(counts.items()))


def path_multiplicity_values() -> tuple[int, ...]:
    """The distinct native transition multiplicities, ``(1, 2, 50, 52)``."""
    return tuple(path_multiplicity_census())


def four_cycle_local_minimum() -> dict[str, object]:
    """A compatible state and a non-compatible local minimum of a discrepancy.

    On a four-cycle with equality constraints the states ``0000`` and ``1111``
    are compatible, while ``0011`` violates two edges and every single-bit flip
    also violates two, so a rule accepting only strict decreases of the
    violation count stops there.  An energy-descent argument therefore needs a
    completeness argument before termination can imply agreement.
    """
    cycle = ((0, 1), (1, 2), (2, 3), (3, 0))

    def violations(state: tuple[int, ...]) -> int:
        return sum(1 for a, b in cycle if state[a] != state[b])

    trap = (0, 0, 1, 1)
    flips = tuple(
        tuple(1 - value if position == index else value for position, value in enumerate(trap))
        for index in range(4)
    )
    flip_counts = tuple(violations(state) for state in flips)
    return {
        "compatible_states": ((0, 0, 0, 0), (1, 1, 1, 1)),
        "compatible_violations": (violations((0, 0, 0, 0)), violations((1, 1, 1, 1))),
        "trap": trap,
        "trap_violations": violations(trap),
        "flip_violations": flip_counts,
        "is_local_minimum": all(count >= violations(trap) for count in flip_counts),
    }


# ---------------------------------------------------------------------------
# 7. Free-energy form and its boundaries
# ---------------------------------------------------------------------------


def lyapunov_form(
    matrix: np.ndarray,
    *,
    gamma: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Solve ``A P + P A^T = -2 Gamma`` and return ``P``, ``Q`` and ``H``.

    Proposition 10 writes a stable linear law as ``A = -(Gamma - Q) H`` with
    ``Gamma`` and ``H`` positive definite and ``Q`` antisymmetric; the
    admissible ``H^{-1}`` are exactly the ``P`` for which ``A P + P A^T`` is
    negative definite, an open convex cone.
    """
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("A must be square")
    dimension = matrix.shape[0]
    gamma = np.eye(dimension) if gamma is None else np.asarray(gamma, dtype=float)
    if gamma.shape != (dimension, dimension):
        raise ValueError("Gamma must match A")
    generator = np.kron(np.eye(dimension), matrix) + np.kron(matrix, np.eye(dimension))
    p = np.linalg.solve(generator, (-2.0 * gamma).reshape(-1)).reshape(dimension, dimension)
    p = 0.5 * (p + p.T)
    q = 0.5 * (matrix @ p - p @ matrix.T)
    return {"P": p, "Q": q, "H": np.linalg.inv(p)}


def d_stability_free_energy_form() -> dict[str, object]:
    """Proposition 10's displayed ``P`` and ``Q`` for the three-node control.

    With ``Gamma = I`` the stable ``A`` admits the Helmholtz form with
    ``P = (1/3)[[11, 4, -8], [4, 5, -2], [-8, -2, 9]]``, whose leading minors are
    ``11/3``, ``13/3`` and ``115/27``, and
    ``Q = (1/3)[[0, 6, 4], [-6, 0, 13], [-4, -13, 0]]``.  The scaled matrix
    ``D_4 A`` is unstable and admits no such form at all.
    """
    form = lyapunov_form(D_STABILITY_A)
    residual = float(
        np.linalg.norm(D_STABILITY_A @ form["P"] + form["P"] @ D_STABILITY_A.T + 2.0 * np.eye(3))
    )
    expected_p = np.array([[11.0, 4.0, -8.0], [4.0, 5.0, -2.0], [-8.0, -2.0, 9.0]]) / 3.0
    expected_q = np.array([[0.0, 6.0, 4.0], [-6.0, 0.0, 13.0], [-4.0, -13.0, 0.0]]) / 3.0
    return {
        "P": form["P"],
        "Q": form["Q"],
        "residual": residual,
        "P_matches_source": bool(np.allclose(form["P"], expected_p, atol=1e-12)),
        "Q_matches_source": bool(np.allclose(form["Q"], expected_q, atol=1e-12)),
        "leading_minors": tuple(
            float(np.linalg.det(form["P"][:size, :size])) for size in range(1, 4)
        ),
        "determinant": float(np.linalg.det(form["P"])),
        "Q_is_antisymmetric": bool(np.allclose(form["Q"], -form["Q"].T, atol=1e-12)),
        "scaled_is_unstable": d_stability_unstable(4.0),
    }


def directed_gradient_control() -> dict[str, np.ndarray]:
    """The forward response and the exact gradient of a directed network.

    For the stable Jacobian ``J = [[-1, 2], [0, -1]]`` and loss derivative
    ``c = (1, 0)``, the direct nudging response is ``-J^{-1} c = (1, 0)``, while
    the adjoint that carries generic parameter sensitivity is
    ``-J^{-T} c = (1, 2)``.  A directed network's forward response therefore need
    not carry its exact gradient, which is where the two-settle rule leaves the
    energy-based family.
    """
    jacobian = np.array([[-1.0, 2.0], [0.0, -1.0]])
    derivative = np.array([1.0, 0.0])
    return {
        "jacobian": jacobian,
        "direct": -np.linalg.solve(jacobian, derivative),
        "adjoint": -np.linalg.solve(jacobian.T, derivative),
    }


def equilibrium_propagation_identity(
    *,
    theta: float = 0.25,
    target: float = 1.5,
    nudge: float = 1e-6,
    steps: int = 200,
) -> dict[str, float]:
    """Two settled states give the loss gradient on a scalar energy-based model.

    For ``F(theta, beta, x) = E(theta, x) + beta C(x)`` with a differentiable
    stationary branch ``x^beta_theta``, the loss ``J(theta) = C(x^0_theta)``
    satisfies ``dJ/dtheta = lim (1/beta)(dE/dtheta(theta, x^beta) -
    dE/dtheta(theta, x^0))``.  The model here is ``E = (x - theta)^2 / 2`` with
    ``C = (x - target)^2 / 2``, whose stationary branch is
    ``x^beta = (theta + beta * target) / (1 + beta)``; the contrast is exactly
    ``(theta - target) / (1 + beta)`` and the exact gradient is ``theta - target``.
    Both the closed form and an actual relaxation to each settled state are
    reported.
    """
    if nudge <= 0.0:
        raise ValueError("the nudge must be positive")
    if steps < 1:
        raise ValueError("steps must be positive")

    def stationary(beta: float) -> float:
        return (theta + beta * target) / (1.0 + beta)

    def settle(beta: float) -> float:
        rate = 0.9 / (1.0 + beta)
        state = 0.0
        for _ in range(steps):
            gradient = (state - theta) + beta * (state - target)
            state -= rate * gradient
            if abs(gradient) < 1e-15:
                break
        return state

    def energy_derivative(state: float) -> float:
        return -(state - theta)

    exact = theta - target
    closed_form = exact / (1.0 + nudge)
    relaxed = (energy_derivative(settle(nudge)) - energy_derivative(settle(0.0))) / nudge
    return {
        "closed_form_contrast": closed_form,
        "relaxed_contrast": relaxed,
        "exact_gradient": exact,
        "closed_form_error": abs(closed_form - relaxed),
        "distance_to_gradient": abs(exact - closed_form),
        "expected_distance_to_gradient": abs(exact) * nudge / (1.0 + nudge),
    }


def surprise_sandwich(
    interior_block: np.ndarray,
    residual: np.ndarray,
    bayesian_surprise: float,
) -> dict[str, float | bool]:
    """Theorem 7 item 5: the Bayesian surprise is pinned by the update residual.

    With ``g`` the update residual at the old rest and ``L_xx`` the interior
    block, ``(1/2) g^T L_xx^{-1} g`` lies between ``||g||^2 / (2 lambda_max)``
    and ``||g||^2 / (2 lambda_min)`` for the extreme eigenvalues of ``L_xx``.
    """
    interior_block = np.asarray(interior_block, dtype=float)
    residual = np.asarray(residual, dtype=float)
    eigenvalues = np.linalg.eigvalsh(interior_block)
    smallest = float(eigenvalues[0])
    largest = float(eigenvalues[-1])
    norm_squared = float(residual @ residual)
    lower = norm_squared / (2.0 * largest)
    upper = norm_squared / (2.0 * smallest)
    value = float(bayesian_surprise)
    return {
        "bayesian_surprise": value,
        "lower_bound": lower,
        "upper_bound": upper,
        "inside": bool(lower <= value <= upper),
    }


def contraction_displacement_bounds(rate: float, disagreement: float) -> dict[str, float]:
    """Theorem 20: ``D / (1 + q) <= ||x* - x|| <= D / (1 - q)``.

    ``D`` is the disagreement ``||R(x, e') - x||`` at any state ``x`` under the
    new evidence, and ``q < 1`` is the contraction rate.  These bound the
    equilibrium displacement; they bound no trajectory length or propagation
    time.
    """
    if not 0.0 <= rate < 1.0:
        raise ValueError("the contraction rate must lie in [0, 1)")
    if disagreement < 0.0:
        raise ValueError("the disagreement must be nonnegative")
    return {
        "lower": disagreement / (1.0 + rate),
        "upper": disagreement / (1.0 - rate),
    }


def replicator_potential_check(
    payoff: np.ndarray | None = None,
    frequencies: tuple[float, float, float] = (0.5, 0.3, 0.2),
) -> dict[str, float | bool]:
    """Selection as a potential flow on the frequency simplex.

    For ``xdot_i = x_i ((A x)_i - x^T A x)`` with symmetric payoffs and
    ``F(x) = (1/2) x^T A x``, the field is ``M(x) grad F`` with
    ``M = diag(x) - x x^T``, so ``Fdot = sum_i x_i ((A x)_i - x^T A x)^2 >= 0``:
    the specified selection model climbs a potential.  Absent types stay absent
    without mutation, and a boundary rest point can still be invadable, which is
    why this is a comparison class rather than a shared mechanism.
    """
    if payoff is None:
        payoff = np.array([[0.0, 1.0, -1.0], [1.0, 0.0, 2.0], [-1.0, 2.0, 0.0]])
    payoff = np.asarray(payoff, dtype=float)
    if payoff.shape != (3, 3) or not np.allclose(payoff, payoff.T):
        raise ValueError("the control uses a symmetric three-type payoff")
    state = np.asarray(frequencies, dtype=float)
    if abs(float(state.sum()) - 1.0) > 1e-12 or np.any(state < 0.0):
        raise ValueError("frequencies must be a probability vector")
    mean = float(state @ payoff @ state)
    gradient = payoff @ state
    field = state * (gradient - mean)
    mobility = np.diag(state) - np.outer(state, state)
    direct_rate = float(sum(state[i] * (gradient[i] - mean) ** 2 for i in range(3)))
    return {
        "field_equals_mobility_gradient": bool(np.allclose(field, mobility @ gradient)),
        "potential_rate_direct": direct_rate,
        "potential_rate_via_mobility": float(gradient @ mobility @ gradient),
        "rates_agree": bool(
            abs(direct_rate - float(gradient @ mobility @ gradient)) < 1e-15
        ),
        "field_norm": float(np.linalg.norm(field)),
    }
