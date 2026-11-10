"""SC3: the native repair law, its conservation law and its energy identity."""

from __future__ import annotations

import math
import random

import pytest

from oph_fep_synergetics.core import (
    all_shell_configurations,
    broken_seams,
    carrier_spectral_gap,
    energy,
    icosahedral_carrier,
    is_broken,
    mismatch,
    native_repair,
    predict_energy_decrease,
    relaxation_gap_from_laplacian,
    shell_gap,
    simulate_repairs,
)

ONE_AT_ZERO = (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


def test_energy_is_half_the_sum_of_squares():
    assert energy((0, 0)) == 0.0
    assert energy((2, 0)) == 2.0
    assert energy((3, 4)) == 12.5
    assert energy(ONE_AT_ZERO) == 0.5


def test_mismatch_is_signed():
    assert mismatch((3, 1), 0, 1) == 2
    assert mismatch((1, 3), 0, 1) == -2


@pytest.mark.parametrize(
    "loads,expected",
    [
        ((0, 0), False),
        ((1, 0), False),
        ((0, 1), False),
        ((2, 0), True),
        ((0, 2), True),
        ((5, 0), True),
    ],
)
def test_is_broken_needs_mismatch_two(loads, expected):
    assert is_broken(loads, 0, 1) is expected


def test_even_mismatch_splits_the_difference_exactly():
    step = native_repair((4, 0), (0, 1), coin=0)
    assert step.seam == (0, 1)
    assert step.coin == 0
    assert step.before == (4, 0)
    assert step.current == 2
    assert step.after == (2, 2)
    assert step.energy_before == 8.0
    assert step.energy_after == 4.0
    assert step.energy_decrease == 4.0


def test_odd_mismatch_places_the_odd_unit_by_coin():
    low = native_repair((3, 0), (0, 1), coin=0)
    high = native_repair((3, 0), (0, 1), coin=1)
    assert (low.current, high.current) == (2, 1)
    assert (low.after, high.after) == ((1, 2), (2, 1))
    assert abs(mismatch(low.after, 0, 1)) == 1
    assert abs(mismatch(high.after, 0, 1)) == 1
    assert low.energy_after == high.energy_after == 2.5


def test_an_unbroken_seam_is_never_moved_by_the_coin():
    for coin in (0, 1):
        step = native_repair((1, 0), (0, 1), coin=coin)
        if coin == 0:
            assert step.after == (0, 1)
        else:
            assert step.after == (1, 0)


def test_repair_conserves_load_and_never_increases_energy(carrier):
    """The whole 78-state shell, both coins, every seam, checked exhaustively."""
    checked = 0
    for state in all_shell_configurations():
        for a, b in carrier.edges:
            for coin in (0, 1):
                step = native_repair(state, (a, b), coin)
                assert sum(step.after) == sum(state)
                assert min(step.after) >= 0
                assert step.energy_after <= step.energy_before
                assert abs(mismatch(step.after, a, b)) <= 1
                assert step.energy_decrease == predict_energy_decrease(state, a, b)
                checked += 1
    assert checked == 78 * 30 * 2 == 4680


def test_predict_energy_decrease_is_the_closed_form():
    assert predict_energy_decrease((0, 0), 0, 1) == 0.0
    assert predict_energy_decrease((1, 0), 0, 1) == 0.0
    assert predict_energy_decrease((2, 0), 0, 1) == 1.0
    assert predict_energy_decrease((3, 0), 0, 1) == 2.0
    assert predict_energy_decrease((6, 0), 0, 1) == 9.0
    assert predict_energy_decrease((7, 0), 0, 1) == (49 - 1) / 4.0


def test_repair_rejects_bad_arguments():
    with pytest.raises(ValueError):
        native_repair((1, 0), (0, 1), coin=2)
    with pytest.raises(ValueError):
        native_repair((1, 0), (0, 0), coin=0)
    with pytest.raises(IndexError):
        native_repair((1, 0), (0, 5), coin=0)
    with pytest.raises(ValueError):
        native_repair((1, -1), (0, 1), coin=0)
    with pytest.raises(ValueError):
        native_repair((1.5, 0), (0, 1), coin=0)


def test_broken_seams_keeps_carrier_order(carrier):
    state = [0] * 12
    state[0] = 3
    expected = tuple((a, b) for a, b in carrier.edges if abs(mismatch(state, a, b)) >= 2)
    assert (0, 1) in expected
    assert broken_seams(state, carrier.edges) == expected


def test_gated_schedule_reaches_rest_and_lowers_energy(carrier):
    start = [0] * 12
    start[0] = 6
    history, final = simulate_repairs(start, carrier.edges, steps=200)
    assert history
    assert sum(final) == 6
    assert broken_seams(final, carrier.edges) == ()
    assert energy(final) < energy(start)
    assert all(step.energy_decrease >= 1.0 for step in history)


def test_ungated_schedule_moves_at_mismatch_one_without_descent(carrier):
    """PR-01's obstruction: an accepted repair that lowers no measure at all."""
    history, final = simulate_repairs(
        ONE_AT_ZERO, [(0, 1)], steps=1, coins=(0,), min_mismatch=1
    )
    assert len(history) == 1
    assert history[0].current == 1
    assert history[0].energy_decrease == 0.0
    assert final == (0, 1) + (0,) * 10


def test_schedule_rejects_a_zero_gate():
    with pytest.raises(ValueError):
        simulate_repairs((1, 0), [(0, 1)], steps=1, min_mismatch=0)


def test_schedule_is_deterministic_for_a_seeded_generator(carrier):
    start = [0] * 12
    start[2] = 1
    start[7] = 1
    first = simulate_repairs(start, carrier.edges, steps=40, rng=random.Random(11))
    second = simulate_repairs(start, carrier.edges, steps=40, rng=random.Random(11))
    assert first == second
    other = simulate_repairs(start, carrier.edges, steps=40, rng=random.Random(12))
    assert other[1] == first[1]


def test_schedule_accepts_an_explicit_coin_stream(carrier):
    start = [0] * 12
    start[0] = 3
    pinned, _ = simulate_repairs(start, carrier.edges, steps=20, coins=[0] * 600)
    again, _ = simulate_repairs(start, carrier.edges, steps=20, coins=[0] * 600)
    assert pinned == again
    assert all(step.coin == 0 for step in pinned)


def test_gap_helpers_agree_with_the_carrier(carrier):
    exact = 5.0 - math.sqrt(5.0)
    assert abs(carrier_spectral_gap() - exact) < 1e-12
    assert abs(carrier_spectral_gap(carrier) - carrier.spectral_gap()) < 1e-15
    assert abs(shell_gap() - carrier_spectral_gap() / 60.0) < 1e-15
    assert abs(relaxation_gap_from_laplacian(carrier) - shell_gap()) < 1e-15
