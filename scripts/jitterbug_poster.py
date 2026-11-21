"""Draw the Jitterbug graphical abstract from the verified geometry and model."""

from __future__ import annotations

from fractions import Fraction

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle

from oph_fep_synergetics.core.jitterbug_inference import (
    jitterbug_coordinate_inference,
    weighted_matching_inference,
)
from oph_fep_synergetics.geometry import jitterbug

PAPER = "#FAF9F5"
INK = "#152B3A"
MUTED = "#465763"
BLUE = "#176981"
PALE_BLUE = "#B9DCE2"
ORANGE = "#A03D16"
TEAL = "#17644C"
RULE = "#CBD5D4"
REPOSITORY = "https://github.com/docxology/OPH-FEP-Synergetics"


def _pose(axis, parameter, rotation, bars, diagonals, highlight=()):
    """Draw the bars and hull diagonals in the same projection and scale."""
    points = np.asarray(jitterbug.scaled_jitterbug_vertices(parameter))
    view = points @ rotation.T
    if highlight:
        axis.add_patch(Polygon(view[list(highlight), :2], facecolor=PALE_BLUE,
                               edgecolor="none", alpha=0.85, zorder=1))
    else:
        # At the limit every edge carries two bars. Draw its geometry once.
        points, index = np.unique(points, axis=0, return_inverse=True)
        bars = sorted({tuple(sorted((index[a], index[b]))) for a, b in bars})
        view = points @ rotation.T
    segments = [
        (view[[a, b], 2].mean(), a, b, colour, width, style)
        for edges, colour, width, style in (
            (bars, BLUE, 1.65, "-"), (diagonals, ORANGE, 1.8, (0, (3, 2))),
        )
        for a, b in edges
    ]
    for depth, a, b, colour, width, style in sorted(segments):
        axis.plot(*view[[a, b], :2].T, color=colour, linewidth=width, linestyle=style,
                  alpha=0.95 if depth > 0 else 0.43, solid_capstyle="round", zorder=2)
    axis.scatter(*view[:, :2].T, s=12, color=INK, edgecolors=PAPER, linewidths=0.45, zorder=3)
    axis.set(xlim=(-2.07, 2.07), ylim=(-2.05, 2.05), aspect="equal")
    axis.set_axis_off()


def _seam(axis, position, *, added):
    """Compare a single matching contrast at fixed positions and nearby bars."""
    x = (position - 0.96, position, position + 0.96)
    # The four bars depict the two two-edge routes around the former square.
    for apex in (0.07, 0.93):
        axis.plot(x, (0.5, apex, 0.5), color=BLUE, alpha=0.55, linewidth=1.3)
    axis.scatter([position, position], [0.07, 0.93], s=12, color=BLUE, zorder=3)
    if added:
        axis.plot((x[0], x[2]), (0.5, 0.5), color=ORANGE, linestyle=(0, (4, 2)),
                  linewidth=2)
    axis.scatter((x[0], x[2]), (0.5, 0.5), marker="o", s=34, color=INK,
                 edgecolors="white", linewidths=0.55, zorder=3)


def draw_jitterbug_poster(rotation):
    """Return the poster figure and a text alternative for the manuscript."""
    bars = jitterbug.jitterbug_rigid_edges()
    matching = jitterbug.jitterbug_matching_edges()
    triangles = jitterbug.jitterbug_rigid_triangles()
    golden = jitterbug.ICOSAHEDRON_PARAMETER
    middle = Fraction(4, 5)
    open_pose = jitterbug.jitterbug_snapshot(1)
    carrier_pose = jitterbug.jitterbug_snapshot(golden)
    limit_pose = jitterbug.jitterbug_snapshot(0)
    joints = open_pose["vertex_count"]
    limit_joints = limit_pose["vertex_count"]
    carrier_edges = carrier_pose["hull_edge_count"]
    before = weighted_matching_inference(0)["closing_pair_variance"]
    after = weighted_matching_inference(1)["closing_pair_variance"]
    open_weight = jitterbug_coordinate_inference(1)["tuned_weight"]
    carrier_weight = jitterbug_coordinate_inference(golden)["tuned_weight"]

    figure = plt.figure(figsize=(8.5, 11), facecolor=PAPER)

    def text(x, y, content, size, *, colour=INK, weight="normal", align="left", **kwargs):
        return figure.text(x, y, content, fontsize=size, color=colour, fontweight=weight,
                           ha=align, va="top", linespacing=1.28, **kwargs)

    def rule(y):
        figure.add_artist(Line2D([0.049, 0.951], [y, y], transform=figure.transFigure,
                                 color=RULE, linewidth=0.9))

    figure.add_artist(Rectangle((0, 0.993), 1, 0.007, transform=figure.transFigure,
                                facecolor=BLUE, edgecolor="none"))
    text(0.050, 0.960, "SYNERGETICS  /  OBSERVER PATCH HOLOGRAPHY  /  FREE ENERGY PRINCIPLE", 9.8,
         colour=BLUE, weight="bold")
    text(0.049, 0.924, "Fuller’s Jitterbug", 35, weight="bold")
    text(0.050, 0.868, "From an FCC shell to the OPH carrier", 20, weight="bold")
    text(0.050, 0.822,
         f"{len(triangles)} rigid triangles hinge at their corners. As they rotate, the outline of their\n"
         f"{joints} joints changes shape while all {len(bars)} bars keep the same length.", 12, colour=MUTED)

    rule(0.767)
    text(0.050, 0.749, "Follow the same linkage as it contracts", 15, weight="bold")
    legend = figure.add_axes((0.051, 0.692, 0.90, 0.030), facecolor=PAPER)
    legend.set_axis_off()
    legend.legend(
        handles=[
            Line2D([], [], color=BLUE, linewidth=2, label="Fixed bar"),
            Rectangle((0, 0), 1, 1, facecolor=PALE_BLUE, label="Same rigid triangle"),
            Line2D([], [], color=ORANGE, linewidth=2, linestyle=(0, (3, 2)),
                   label="Hull diagonal (no bar)"),
        ], loc="center left", ncol=3, bbox_to_anchor=(-0.012, 0.45), fontsize=10.2,
        frameon=False, columnspacing=2.7, handlelength=2.3,
    )

    landmarks = (
        (1, "Cuboctahedron", f"FCC coordination shell\n{len(bars)} edges"),
        (middle, "Contracting pose", f"Square openings turn skew\n{len(matching)} hull diagonals appear"),
        (golden, "Regular icosahedron", f"OPH carrier geometry\n{carrier_edges} equal hull edges"),
        (0, "Octahedral limit", f"Joints meet in pairs\n{limit_joints} distinct positions"),
    )
    open_points = np.asarray(jitterbug.scaled_jitterbug_vertices(1)) @ rotation.T
    highlight = max(triangles, key=lambda face: float(open_points[list(face), 2].mean()))
    for index, (parameter, title, label) in enumerate(landmarks):
        center = 0.156 + index * 0.230
        if parameter == golden:
            figure.add_artist(FancyBboxPatch(
                (center - 0.107, 0.484), 0.214, 0.206, boxstyle="round,pad=0.009,rounding_size=0.011",
                transform=figure.transFigure, facecolor="#EEF3EF", edgecolor="none", zorder=-1,
            ))
        axis = figure.add_axes((center - 0.104, 0.532, 0.208, 0.160), facecolor="none")
        _pose(axis, parameter, rotation, bars, matching if 0 < float(parameter) < 1 else (),
              () if not parameter else highlight)
        text(center, 0.531, title, 10.8, align="center", weight="bold")
        text(center, 0.508, label, 9.4, align="center", colour=MUTED)
        if index < len(landmarks) - 1:
            text(center + 0.115, 0.631, "→", 11, align="center", colour=MUTED)

    rule(0.461)
    text(0.050, 0.442, "What changes when the hull diagonals count as seams?", 15, weight="bold")
    text(0.050, 0.408,
         f"At the same icosahedral pose, add all {len(matching)} diagonals. The diagram shows one.",
         10.8, colour=MUTED)
    comparison = figure.add_axes((0.15, 0.304, 0.70, 0.078), facecolor=PAPER)
    comparison.set(xlim=(-3.9, 3.9), ylim=(0, 1))
    comparison.set_axis_off()
    _seam(comparison, -2.05, added=False)
    _seam(comparison, 2.05, added=True)
    comparison.text(0, 0.5, "→", color=TEAL, fontsize=16, va="center", ha="center")
    text(0.314, 0.308, f"{len(bars)} unit bars", 11, align="center", weight="bold")
    text(0.684, 0.308, f"{carrier_edges} unit seams", 11, align="center", weight="bold")
    text(0.314, 0.285, f"Full-graph contrast variance {before}", 11, align="center", colour=BLUE)
    text(0.684, 0.285, f"Full-graph contrast variance {after}", 11, align="center", colour=ORANGE)

    rule(0.250)
    text(0.050, 0.231, "Match the readouts at each fixed pose", 14, weight="bold")
    text(0.050, 0.199, "OPH side: sample fair binary repairs on seams\nin proportion to their weights.",
         10.4, colour=MUTED)
    text(0.050, 0.158, "FEP side: use the weighted graph Laplacian\nas the centred Gaussian precision.",
         10.4, colour=MUTED)
    text(0.050, 0.117, "At this weight, repair’s coordinate lag covariances\nmatch those of a Gaussian Markov process.",
         10.2, colour=INK, weight="bold")
    text(0.050, 0.079, "With unit bars alone, precision stays constant.", 9.2, colour=MUTED)

    weight_axis = figure.add_axes((0.655, 0.120, 0.271, 0.098), facecolor=PAPER)
    parameters = [1 - Fraction(step, 72) * (1 - golden) for step in range(73)]
    weights = [jitterbug_coordinate_inference(t)["tuned_weight"] for t in parameters]
    weight_axis.plot(list(map(float, parameters)), list(map(float, weights)), color=TEAL, linewidth=2)
    weight_axis.scatter([1, float(golden)], [float(open_weight), float(carrier_weight)], color=TEAL,
                        edgecolor=PAPER, linewidth=0.8, zorder=3, s=38)
    weight_axis.set(xlim=(1.025, float(golden) - 0.024), ylim=(-0.14, 1.14))
    weight_axis.set_xticks([1, float(golden)], ["Cuboctahedron", "Regular ico."])
    weight_axis.set_yticks([float(open_weight), float(carrier_weight)], [str(open_weight), str(carrier_weight)])
    weight_axis.tick_params(length=0, pad=4, labelsize=9)
    weight_axis.grid(visible=True, axis="y", color=RULE, alpha=0.8)
    weight_axis.grid(visible=False, axis="x")
    weight_axis.spines[["top", "right"]].set_visible(False)
    text(0.642, 0.235, r"Tuned diagonal weight $w=1/t-t$", 10, colour=TEAL, weight="bold")
    text(0.645, 0.091, "$t$ is the ratio of the two nonzero joint\ncoordinate magnitudes (small / large).", 8.4, colour=MUTED)

    rule(0.061)
    text(0.050, 0.047, "Daniel Ari Friedman · Bernhard Mueller  |  Observers, Blankets, and Screens", 9, weight="bold")
    text(0.050, 0.027, "Fuller, Synergetics, fig. 460.08", 8, colour=BLUE,
         url="https://rwgrayprojects.com/synergetics/s04/figs/f6008.html")
    text(0.267, 0.027, "Stachel, The HEUREKA-Polyhedron", 8, colour=BLUE,
         url="https://www.geometrie.tuwien.ac.at/stachel/060_Heureka.pdf")
    text(0.951, 0.027, "github.com/docxology/OPH-FEP-Synergetics", 8, colour=BLUE, align="right", url=REPOSITORY)

    alt = (
        f"Fuller's Jitterbug is a linkage of {len(triangles)} rigid triangles and {len(bars)} fixed bars. "
        "Four poses at a common scale show the cuboctahedral FCC coordination shell, a contracting "
        f"configuration with {len(matching)} dashed hull diagonals, the regular icosahedral hull with "
        f"{joints} joints and {carrier_edges} equal hull edges matching the OPH carrier, and the "
        f"octahedral limit with {limit_joints} joint positions. The same rigid triangle is shaded until "
        "the limit. At the regular icosahedral pose, the comparison keeps all coordinates fixed and adds "
        f"unit weight to the diagonals. The intrinsic Gaussian contrast variance across each diagonal drops "
        f"from {before} to {after}. The coordinate lag covariances of weighted fair binary seam repair "
        "match those of a Gaussian Markov process at each fixed pose from the cuboctahedron to the regular "
        "icosahedron only for matching weight one over t minus t. The plotted weight rises from zero to one. "
        "With just the fixed bars, precision is constant along the motion."
    )
    return figure, alt
