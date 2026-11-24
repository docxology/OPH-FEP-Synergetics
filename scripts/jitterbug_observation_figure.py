"""Show the covariance and precision consequences of Jitterbug readouts."""

from __future__ import annotations

from fractions import Fraction

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

from oph_fep_synergetics.core.jitterbug_motion import jitterbug_motion_defect
from oph_fep_synergetics.core.jitterbug_projection import jitterbug_pair_projection

BLUE = "#0072B2"
VERMILLION = "#D55E00"
GREEN = "#009E73"
GREY = "#6B7280"
INK = "#1B2A3A"


def draw_jitterbug_observation_limits() -> tuple[plt.Figure, str]:
    """Return the manuscript figure and its accessible description."""
    figure = plt.figure(figsize=(10.5, 4.6))
    panels = figure.add_gridspec(1, 2)
    motion = figure.add_subplot(panels[0])
    pair_grid = panels[1].subgridspec(2, 1, hspace=0.16)
    precision = figure.add_subplot(pair_grid[0])
    lag = figure.add_subplot(pair_grid[1], sharex=precision)

    a, b = Fraction(3, 5), Fraction(1, 2)
    final_poses = [Fraction(step, 100) for step in range(25, 62)]
    final_defects = [float(jitterbug_motion_defect(a, b, c)) for c in final_poses]
    motion.axhline(0, color=GREY, linewidth=1.0, linestyle=":", zorder=1)
    motion.plot(list(map(float, final_poses)), final_defects, color=BLUE, linewidth=2.35)
    controls = (
        (Fraction(1, 3), "Contraction", VERMILLION, "D", (8, -16), "left", "top"),
        (b, "Hold", GREEN, "o", (-5, 17), "right", "bottom"),
        (a, "Reversal", VERMILLION, "s", (-4, 10), "right", "bottom"),
    )
    exact_defects = []
    for pose, label, colour, marker, offset, align, valign in controls:
        defect = jitterbug_motion_defect(a, b, pose)
        exact_defects.append(defect)
        point = (float(pose), float(defect))
        motion.scatter(*point, color=colour, marker=marker, s=54,
                       edgecolor="white", linewidth=0.85, zorder=6)
        motion.annotate(
            f"{label}\n${defect}$", point, xytext=offset, textcoords="offset points",
            ha=align, va=valign, fontsize=10.1, linespacing=1.25, color=INK,
        )
    motion.set(xlim=(0.235, 0.624), ylim=(-0.020, 0.014),
               xlabel="Final pose $c$", ylabel="Normalised covariance residual")
    motion.set_xticks((0.25, 1 / 3, 0.5, 0.6), labels=("$1/4$", "$1/3$", "$1/2$", "$3/5$"))
    motion.set_yticks((-0.02, -0.01, 0, 0.01))
    motion.yaxis.set_major_formatter(StrMethodFormatter("{x:g}"))
    motion.set_title("A   Moving geometric coordinates", loc="left", fontweight="bold")
    motion.text(0.035, 0.945, f"Departure poses $a={a}$, $b={b}$",
                transform=motion.transAxes, va="top", fontsize=10.3, color=INK)

    weight_grid = [Fraction(step, 8) for step in range(97)]
    samples = [jitterbug_pair_projection(weight) for weight in weight_grid]
    weights = [float(row["closing_weight"]) for row in samples]
    reference = jitterbug_pair_projection(Fraction(3, 2))
    reference_weight = reference["closing_weight"]
    opposite = reference["pair_mean_opposite_precision"]
    lag_defect = reference["odd_lag_two_defect"]
    for axis, key, colour, marker, scale in (
        (precision, "pair_mean_opposite_precision", BLUE, "s", 1),
        (lag, "odd_lag_two_defect", GREEN, "o", 1000),
    ):
        axis.plot(weights, [float(row[key] * scale) for row in samples], color=colour, linewidth=2.25)
        axis.axvline(float(reference_weight), color=GREY, linestyle=":", linewidth=1.0, zorder=1)
        value = reference[key]
        point = (float(reference_weight), float(value * scale))
        axis.scatter(*point, color=colour, marker=marker, s=55, edgecolor="white", linewidth=0.85, zorder=5)
        axis.annotate(f"${value}$", point, xytext=(13, 8), textcoords="offset points",
                      color=INK, fontsize=10.3, ha="left", va="bottom")
        axis.set_xlim(-0.2, 12.2)
        axis.grid(False, axis="x")
    precision.set_title("B   Finite pair projections", loc="left", fontweight="bold")
    precision.set(ylabel="Opposite-site\nprecision entry", ylim=(0, 1.10))
    precision.set_yticks((0, 0.5, 1))
    precision.tick_params(axis="x", bottom=False, labelbottom=False)
    lag.set(ylabel="Normalised pair-total\nlag-two defect ($\\times 10^{-3}$)",
            xlabel="Matching conductance $w$ (unit bars)", ylim=(0, 1.96))
    lag.set_xticks((0, 1.5, 4, 8, 12), labels=("$0$", "$3/2$", "$4$", "$8$", "$12$"))
    lag.set_yticks((0, 0.5, 1, 1.5))
    for axis in (precision, lag):
        axis.yaxis.set_major_formatter(StrMethodFormatter("{x:g}"))
    motion.grid(False, axis="x")

    contraction, hold, reversal = exact_defects
    alt = (
        f"Panel A varies the final Jitterbug pose after one tuned repair attempt at each of "
        f"the departure poses {a} and {b}. The normalised Gaussian Markov covariance residual "
        f"is {contraction} for contraction to {controls[0][0]}, {hold} for holding the second "
        f"pose, and {reversal} for reversal to {a}. Panel B uses separate scales for the "
        "positive opposite-site entry in the pair-mean precision and the normalised odd "
        "lag-two covariance defect of pair totals. Both quantities decrease with matching "
        f"conductance. At weight {reference_weight} they are {opposite} and {lag_defect}, "
        "respectively."
    )
    return figure, alt
