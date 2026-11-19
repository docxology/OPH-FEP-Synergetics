"""The native seam-repair law on an integer-load carrier.

Observer Patch Holography's executable observer stores a nonnegative integer
load on each of twelve ports.  A seam ``(a, b)`` is *broken* when the two loads
differ by at least two.  Repairing a seam moves whole units across it so that
the two loads differ by at most one; an odd unit is placed by a fair coin.  The
law conserves total load and decreases the quadratic energy ``sum(x_i^2)/2`` by
``(d^2 - (|d| mod 2)) / 4`` for a seam whose load difference is ``d``.

This module implements the law exactly on Python integers, exposes the finite
shell of configurations with total load two, and records the spectral gap of
the carrier's classical relaxation.  It does not implement OPH's quantum event
algebra or its continuum limit; those are cited, not reimplemented.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from fractions import Fraction
from operator import index
from typing import Iterable, Sequence

from .carrier import IcosahedralCarrier, icosahedral_carrier


@dataclass(frozen=True)
class RepairStep:
    """One accepted native repair.

    ``current`` is the signed integer current ``j`` in OPH's notation:
    ``y = x + j(e_b - e_a)``.  A positive current moves load from port ``a``
    to port ``b``.
    """

    seam: tuple[int, int]
    coin: int
    before: tuple[int, ...]
    after: tuple[int, ...]
    current: int
    energy_before: Fraction
    energy_after: Fraction

    @property
    def energy_decrease(self) -> Fraction:
        return self.energy_before - self.energy_after


def energy(loads: Sequence[int]) -> Fraction:
    """The exact quadratic load energy ``sum(x_i^2) / 2``."""
    return Fraction(sum(index(load) ** 2 for load in loads), 2)


def _validated_loads(loads: Sequence[int]) -> tuple[int, ...]:
    if any(int(load) != load or load < 0 for load in loads):
        raise ValueError("loads must be nonnegative integers")
    return tuple(int(load) for load in loads)


def mismatch(loads: Sequence[int], a: int, b: int) -> int:
    """The signed load difference ``x_a - x_b`` across a seam."""
    return int(loads[a]) - int(loads[b])


def is_broken(loads: Sequence[int], a: int, b: int) -> bool:
    """Whether a seam has mismatch magnitude at least two."""
    return abs(mismatch(loads, a, b)) >= 2


def native_repair(
    loads: Sequence[int],
    seam: tuple[int, int],
    coin: int,
) -> RepairStep:
    """Apply the native repair law to one seam.

    Parameters
    ----------
    loads:
        Nonnegative integer loads, one per port.
    seam:
        Two distinct port indices ``(a, b)``.
    coin:
        ``0`` or ``1``; used only when the mismatch is odd.
    """
    if coin not in (0, 1):
        raise ValueError("coin must be 0 or 1")
    coin = int(coin)
    a, b = seam
    if a == b:
        raise ValueError("a seam joins two distinct ports")
    if a < 0 or b < 0 or a >= len(loads) or b >= len(loads):
        raise IndexError("seam port outside the load vector")
    before = _validated_loads(loads)
    difference = before[a] - before[b]
    numerator = difference - (2 * coin - 1) * (abs(difference) % 2)
    current = numerator // 2
    after = list(before)
    after[a] -= current
    after[b] += current
    if min(after) < 0:  # pragma: no cover - follows from the load range
        raise AssertionError("native repair produced a negative load")
    return RepairStep(
        seam=(a, b),
        coin=coin,
        before=before,
        after=tuple(after),
        current=current,
        energy_before=energy(before),
        energy_after=energy(after),
    )


def predict_energy_decrease(loads: Sequence[int], a: int, b: int) -> int:
    """The exact quadratic-energy decrease of repairing seam ``(a, b)``."""
    difference = mismatch(loads, a, b)
    return (difference * difference - (abs(difference) % 2)) // 4


def classify_native_outcome(
    loads: Sequence[int],
    seam: tuple[int, int],
    coin: int,
) -> dict[str, object]:
    """Classify the native repair of one seam as strict, zero, or noop.

    The shared primitive behind the census controls: the outcome is
    ``"noop"`` when the law leaves the configuration alone, ``"strict"``
    when the quadratic energy strictly decreases, and ``"zero"`` when the
    configuration changes but the energy stays flat (mismatch one).
    """
    step = native_repair(loads, seam, coin)
    if step.after == step.before:
        kind = "noop"
    elif step.energy_decrease > 0:
        kind = "strict"
    else:
        kind = "zero"
    return {"kind": kind, "step": step}


def broken_seams(loads: Sequence[int], seams: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    """Return the broken seams in deterministic order."""
    return tuple((a, b) for a, b in seams if is_broken(loads, a, b))


def simulate_repairs(
    loads: Sequence[int],
    seams: Sequence[tuple[int, int]],
    *,
    steps: int,
    coins: Sequence[int] | None = None,
    rng: random.Random | None = None,
    stop_at_rest: bool = True,
    min_mismatch: int = 2,
) -> tuple[list[RepairStep], tuple[int, ...]]:
    """Run native repairs under a declared seam schedule.

    The schedule is the cyclic list ``seams``.  A seam is repaired when its
    mismatch magnitude is at least ``min_mismatch``; the default two is OPH's
    declared gate for a *broken* seam.  Setting ``min_mismatch=1`` exposes the
    local law's mismatch-one moves, which can swap a unit without lowering the
    quadratic energy; that is the executable content of PR-01's obstruction.
    If ``stop_at_rest`` is true the run ends when no seam passes the gate.
    ``coins`` supplies a deterministic coin sequence; when it is omitted a local
    ``random.Random`` instance is used (seeded by the caller for reproducibility).
    """
    if min_mismatch < 1:
        raise ValueError("min_mismatch must be at least one")
    state = _validated_loads(loads)
    history: list[RepairStep] = []
    if coins is None:
        generator = rng if rng is not None else random.Random(0)
        coin_stream = (generator.randrange(2) for _ in range(steps * max(1, len(seams))))
    else:
        coin_stream = iter(coins)
    for _ in range(steps):
        if stop_at_rest and not any(
            abs(mismatch(state, a, b)) >= min_mismatch for a, b in seams
        ):
            break
        for a, b in seams:
            if abs(mismatch(state, a, b)) < min_mismatch:
                continue
            step = native_repair(state, (a, b), next(coin_stream))
            state = step.after
            history.append(step)
    return history, state


def shell_gap() -> float:
    """The OPH carrier gap ``(5 - sqrt(5)) / 60``."""
    return (5.0 - math.sqrt(5.0)) / 60.0


def mean_native_step(
    loads: Sequence[int],
    seams: Sequence[tuple[int, int]],
) -> tuple[Fraction, ...]:
    """The exact mean outcome of one native repair under the uniform schedule.

    One attempt draws a seam uniformly from ``seams`` and a fair coin, so the
    mean is the average of :func:`native_repair` over the ``2 * len(seams)``
    equally likely draws.  The result is returned in exact rationals.  On the
    twelve-port carrier it equals ``x - L x / 60`` for the carrier Laplacian
    ``L``: the mean current across a seam is half its mismatch for either
    parity, so the mean of the integer law is one gradient step of the
    seam-disagreement energy ``(1/2) x^T L x`` with step ``1 / (2 |E|)``.
    """
    if not seams:
        raise ValueError("at least one seam is needed")
    draws = 2 * len(seams)
    totals = [0] * len(loads)
    for seam in seams:
        for coin in (0, 1):
            for port, load in enumerate(native_repair(loads, seam, coin).after):
                totals[port] += load
    return tuple(Fraction(total, draws) for total in totals)


def binary_shell_configurations(
    carrier: IcosahedralCarrier | None = None,
) -> tuple[tuple[int, ...], ...]:
    """The 66 binary load vectors with two unit loads on distinct ports.

    This is the module-level convenience wrapper around
    :meth:`IcosahedralCarrier.binary_shell_configurations`, matching OPH's
    ``model/observer/native.py`` binary shell.
    """
    graph = carrier if carrier is not None else icosahedral_carrier()
    return graph.binary_shell_configurations()


def all_shell_configurations(
    carrier: IcosahedralCarrier | None = None,
    total: int = 2,
) -> tuple[tuple[int, ...], ...]:
    """All 78 nonnegative integer load vectors of the given total."""
    graph = carrier if carrier is not None else icosahedral_carrier()
    return graph.shell_configurations(total=total)


def carrier_spectral_gap(carrier: IcosahedralCarrier | None = None) -> float:
    """The icosahedron Laplacian gap ``5 - sqrt(5)``."""
    graph = carrier if carrier is not None else icosahedral_carrier()
    return graph.spectral_gap()


def relaxation_gap_from_laplacian(carrier: IcosahedralCarrier | None = None) -> float:
    """The seam-average gap ``lambda_1 / (2W)`` for total seam weight ``W``.

    For the unweighted icosahedron ``lambda_1 = 5 - sqrt(5)`` and ``W = 30``,
    so this equals :func:`shell_gap`.
    """
    graph = carrier if carrier is not None else icosahedral_carrier()
    total_weight = float(graph.seam_count)
    return graph.spectral_gap() / (2.0 * total_weight)
