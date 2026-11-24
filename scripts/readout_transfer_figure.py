"""Display finite readout transfer controls from the executable census."""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from oph_fep_synergetics.core.pair_projection_transfer import pair_projection_transfer_census, switching_pair_census
from oph_fep_synergetics.core.readout_motion import moving_readout_census

INK = "#1b2a3a"
MUTED = "#526371"
BLUE = "#0072B2"
ORANGE = "#BC4B0A"
GREEN = "#007C66"
GRID = "#D7E2E9"


def _label(axis, x, y, message, *, color=INK, weight="normal", size=10, align="left"):
    axis.text(x, y, message, color=color, fontweight=weight, fontsize=size,
              ha=align, va="center", transform=axis.transAxes)


def _pair_path(axis, row, y):
    pairs = row["pairs"]
    coordinates = {vertex: 0.045 + 0.068 * k for k, vertex in enumerate(vertex for pair in pairs for vertex in pair)}
    for first, second in pairs:
        axis.add_patch(FancyBboxPatch((coordinates[first] - 0.022, y - 0.063), 0.112, 0.126,
                                     boxstyle="round,pad=0,rounding_size=0.022", transform=axis.transAxes,
                                     facecolor="#E3F0F5", edgecolor="none"))
    for first, second in row["edges"]:
        axis.plot((coordinates[first], coordinates[second]), (y, y), transform=axis.transAxes,
                  color=MUTED, linewidth=1.8)
    axis.scatter(list(coordinates.values()), [y] * len(coordinates), s=23,
                 c=INK, linewidths=0, transform=axis.transAxes, zorder=4)


def draw_readout_transfer():
    """Return the transfer figure and its full text alternative."""
    motion = moving_readout_census()
    graphs = pair_projection_transfer_census()
    switched = switching_pair_census()
    figure, panels = plt.subplots(2, 2, figsize=(10.6, 5.65))
    ax, bx, cx, dx = panels.flat
    figure.subplots_adjust(left=0.035, right=0.975, top=0.93, bottom=0.035, wspace=0.17, hspace=0.35)
    for axis, title in zip(panels.flat, (
        "A   A longer lag survives on the weighted square",
        "B   Path precision signs and lag defects",
        "C   Fixed pairs can preserve lags and retain memory",
        "D   Changing pairings on the tetrahedron",
    )):
        axis.set_axis_off()
        axis.set_title(title, loc="left", fontweight="bold", fontsize=12.3, color=INK)

    xs = (0.09, 0.36, 0.63, 0.90)
    arc = FancyArrowPatch((xs[0], 0.65), (xs[-1], 0.65), transform=ax.transAxes,
                          connectionstyle="arc3,rad=-0.24", arrowstyle="-|>",
                          color=ORANGE, linewidth=1.8, mutation_scale=12)
    ax.add_patch(arc)
    _label(ax, 0.50, 0.89, f"Endpoint residual  {motion['endpoint_defects'][0]}", color=ORANGE,
           weight="bold", align="center", size=11.3)
    ax.plot(xs, [0.56] * len(xs), color=BLUE, linewidth=2, transform=ax.transAxes)
    ax.scatter(xs, [0.56] * len(xs), color=BLUE, s=35, edgecolors="white", linewidth=0.8,
               transform=ax.transAxes, zorder=5)
    for x, frame in zip(xs, ("$x+y$", "$x$", "$x$", "$x+y$")):
        _label(ax, x, 0.44, frame, align="center", size=11)
    for left, right, covariance in zip(xs, xs[1:], motion["adjacent_covariances"]):
        _label(ax, (left + right) / 2, 0.63, str(covariance), color=BLUE, align="center", size=9)
    _label(ax, 0.025, 0.24, "Blue values give adjacent covariances.")
    _label(ax, 0.025, 0.10, "Both consecutive three-time residuals are zero.", weight="bold")

    _label(bx, 0.03, 0.83, "Consecutive pairs", weight="bold", size=9.5)
    _label(bx, 0.50, 0.83, "Endpoint\nprecision", weight="bold", size=9.5)
    _label(bx, 0.79, 0.83, "Lag-defect\ntrace", weight="bold", size=9.5)
    for case, y in (("path4", 0.57), ("path6", 0.29)):
        row = graphs[case]
        _pair_path(bx, row, y)
        _label(bx, 0.50, y, f"{row['endpoint_precision']:+}", color=ORANGE if row["endpoint_precision"] > 0 else BLUE,
               weight="bold", size=11.5)
        _label(bx, 0.79, y, str(row["lag_two_defect_trace"]), weight="bold", size=11.5)
    _label(bx, 0.03, 0.09, "The shorter path has ordinary site conductances.", color=MUTED)

    labels = ("Unit square", "Complete graph")
    _label(cx, 0.405, 0.81, "Present alone", weight="bold", size=9.5)
    _label(cx, 0.72, 0.81, "Past + present", weight="bold", size=9.5)
    for label, key, y in zip(labels, ("cycle4", "complete4"), (0.57, 0.34)):
        with_past, present_only = graphs[key]["return_probabilities"]
        cx.plot((0.02, 0.97), (y + 0.11, y + 0.11), color=GRID, linewidth=0.8, transform=cx.transAxes)
        _label(cx, 0.03, y, label, size=10.4)
        _label(cx, 0.49, y, str(present_only), color=BLUE, align="center", size=11.5)
        _label(cx, 0.82, y, str(with_past), color=ORANGE if with_past != present_only else GREEN,
               weight="bold", align="center", size=11.5)
    _label(cx, 0.03, 0.11, "Return probabilities after a split middle count.", color=MUTED)

    for x, letter in zip((0.13, 0.50, 0.87), ("A", "B", "A")):
        color = BLUE if letter == "A" else GREEN
        dx.add_patch(FancyBboxPatch((x - 0.059, 0.60), 0.118, 0.19, transform=dx.transAxes,
                                    boxstyle="round,pad=0,rounding_size=0.022",
                                    facecolor=color, edgecolor="none"))
        _label(dx, x, 0.695, letter, color="white", weight="bold", align="center", size=12)
    for x in (0.313, 0.683):
        _label(dx, x, 0.695, "→", align="center", color=MUTED, size=17)
    _label(dx, 0.50, 0.44, f"Endpoint covariance residual  {switched['lag_two_residual']}",
           align="center", color=ORANGE, weight="bold", size=11)
    history, present = switched["return_probabilities"]
    _label(dx, 0.50, 0.25, f"Return probability  {present} alone  ·  {history} with past",
           align="center", size=10.3)
    _label(dx, 0.50, 0.09, "Each pairing gives an exact count quotient when held fixed.",
           color=MUTED, align="center", size=9.3)

    short, long = graphs["path4"], graphs["path6"]
    square_history, square_present = graphs["cycle4"]["return_probabilities"]
    complete_history, complete_present = graphs["complete4"]["return_probabilities"]
    alt = (
        "Four panels compare readout transfer. A, on the weighted square the readouts are x plus y, x, x, and x plus y. "
        f"The two consecutive three-time residuals are zero and the endpoint residual is {motion['endpoint_defects'][0]}. "
        f"B, for consecutive pairs on the shorter unit path, endpoint precision is {short['endpoint_precision']} "
        f"and lag-defect trace is {short['lag_two_defect_trace']}. On the longer unit path these are "
        f"positive {long['endpoint_precision']} and {long['lag_two_defect_trace']}. "
        f"C, conditioned on a split middle count, the unit square's return probability is {square_present} "
        f"and becomes {square_history} given a prior return target. Both complete-graph rates are {complete_present} "
        f"and {complete_history}, respectively. D, on the tetrahedron the readouts alternate A, B, A, "
        f"giving covariance residual {switched['lag_two_residual']}. Return probability is {present} "
        f"conditioned on the present and {history} conditioned on both past and present."
    )
    return figure, alt
