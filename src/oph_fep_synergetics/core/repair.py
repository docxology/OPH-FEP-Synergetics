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
    energy_before: float
    energy_after: float

    @property
    def energy_decrease(self) -> float:
        return self.energy_before - self.energy_after


def energy(loads: Sequence[int]) -> float:
    """The quadratic seam energy ``sum(x_i^2) / 2``."""
    return sum(x * x for x in loads) / 2.0


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
    a, b = seam
    if a == b:
        raise ValueError("a seam joins two distinct ports")
    if a < 0 or b < 0 or a >= len(loads) or b >= len(loads):
        raise IndexError("seam port outside the load vector")
    if any(int(x) != x or x < 0 for x in loads):
        raise ValueError("loads must be nonnegative integers")
    before = tuple(int(x) for x in loads)
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


def predict_energy_decrease(loads: Sequence[int], a: int, b: int) -> float:
    """The exact quadratic-energy decrease of repairing seam ``(a, b)``."""
    difference = mismatch(loads, a, b)
    return (difference * difference - (abs(difference) % 2)) / 4.0


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
    state = tuple(int(x) for x in loads)
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
