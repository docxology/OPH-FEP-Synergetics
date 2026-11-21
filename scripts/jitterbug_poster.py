"""Draw the Jitterbug graphical abstract from the verified geometry and model."""

from __future__ import annotations

from fractions import Fraction

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle

from oph_fep_synergetics.core.jitterbug_inference import jitterbug_contraction_census
from oph_fep_synergetics.core.jitterbug_repair import jitterbug_pair_repair_census
from oph_fep_synergetics.geometry import jitterbug

PAPER = "#FAF9F5"
INK = "#152B3A"
MUTED = "#465763"
BLUE = "#176981"
PALE_BLUE = "#B9DCE2"
ORANGE = "#A03D16"
TEAL = "#17644C"
RULE = "#CBD5D4"
POST_ICO = "#F5EBE3"
REPOSITORY = "https://github.com/docxology/OPH-FEP-Synergetics"


def _pose(axis, parameter, rotation, bars, diagonals, *, highlight=()):
    """Draw the selected seams in a shared projection and scale."""
    points = np.asarray(jitterbug.scaled_jitterbug_vertices(parameter))
    if parameter == 0:
        # At the limit every edge carries two bars. Draw its geometry once.
        points, index = np.unique(points, axis=0, return_inverse=True)
        bars = sorted({tuple(sorted((index[a], index[b]))) for a, b in bars})
    view = points @ rotation.T
    if highlight:
        axis.add_patch(Polygon(view[list(highlight), :2], facecolor=PALE_BLUE,
                               edgecolor="none", alpha=0.85, zorder=1))
    segments = [
        (view[[a, b], 2].mean(), a, b, colour, width, style)
        for edges, colour, width, style in (
            (bars, BLUE, 1.65, "-"), (diagonals, ORANGE, 1.8, (0, (3, 2))),
        )
        for a, b in edges
    ]
    for depth, a, b, colour, width, style in sorted(segments):
        axis.plot(*view[[a, b], :2].T, color=colour, linewidth=width, linestyle=style,
                  alpha=0.95 if depth > 0 or not bars else 0.43, solid_capstyle="round", zorder=2)
    axis.scatter(*view[:, :2].T, s=12, color=INK, edgecolors=PAPER, linewidths=0.45, zorder=3)
    axis.set(xlim=(-2.07, 2.07), ylim=(-2.05, 2.05), aspect="equal")
    axis.set_axis_off()


def draw_jitterbug_poster(rotation):
    """Return the poster figure and a text alternative for the manuscript."""
    bars = jitterbug.jitterbug_rigid_edges()
    matching = jitterbug.jitterbug_matching_edges()
    triangles = jitterbug.jitterbug_rigid_triangles()
    golden = jitterbug.ICOSAHEDRON_PARAMETER
    half = Fraction(1, 2)
    carrier = jitterbug.jitterbug_snapshot(golden)
    post_ico = jitterbug.jitterbug_snapshot(half)
    limit = jitterbug.jitterbug_snapshot(0)
    inference = jitterbug_contraction_census()
    comparisons = {row["case"]: row for row in inference["continuation_samples"]}
    half_tuned = comparisons["half-tuned"]
    third_inverse_square = comparisons["third-inverse-square"]
    repair = jitterbug_pair_repair_census()
    rank = post_ico["vertex_count"] - half_tuned["nearest_only_nullity"]
    nullity = half_tuned["nearest_only_centered_nullity"]
    decay = inference["tuned_limit"]["mass_clock_covariance_decay_rate"]
    shared_bar_weight, = set(limit["bar_multiplicities"])

    figure = plt.figure(figsize=(8.5, 11), facecolor=PAPER)

    def text(x, y, content, size, *, colour=INK, weight="normal", align="left", **kwargs):
        return figure.text(x, y, content, fontsize=size, color=colour, fontweight=weight,
                           ha=align, va="top", linespacing=1.24, **kwargs)

    def rule(y):
        figure.add_artist(Line2D([0.049, 0.951], [y, y], transform=figure.transFigure,
                                 color=RULE, linewidth=0.9))

    figure.add_artist(Rectangle((0, 0.993), 1, 0.007, transform=figure.transFigure,
                                facecolor=BLUE, edgecolor="none"))
    text(0.050, 0.960, "SYNERGETICS  /  OBSERVER PATCH HOLOGRAPHY  /  FREE ENERGY PRINCIPLE", 9.8,
         colour=BLUE, weight="bold")
    text(0.049, 0.929, "Fuller’s Jitterbug", 35, weight="bold")
    text(0.050, 0.872, "From an FCC shell to paired octahedral sites", 20, weight="bold")
    text(0.050, 0.826,
         f"{len(triangles)} rigid triangles turn about their joined corners while the {len(bars)} bars keep their length.\n"
         "The regular icosahedral hull supplies the seam graph for observer patch holography (OPH).\n"
         r"The ratio $t$ compares each joint’s smaller in-plane coordinate magnitude to its larger one.",
         11.0, colour=MUTED)

    rule(0.765)
    text(0.050, 0.747, "Four poses at the same bar length and scale", 15, weight="bold")
    legend = figure.add_axes((0.051, 0.683, 0.90, 0.030), facecolor=PAPER)
    legend.set_axis_off()
    legend.legend(
        handles=[
            Line2D([], [], color=BLUE, linewidth=2, label="Fixed bar"),
            Rectangle((0, 0), 1, 1, facecolor=PALE_BLUE, label="Same rigid triangle"),
            Line2D([], [], color=ORANGE, linewidth=2, linestyle=(0, (3, 2)), label="Hull diagonal"),
        ], loc="center left", ncol=3, bbox_to_anchor=(-0.012, 0.45), fontsize=9.8,
        frameon=False, columnspacing=2.7, handlelength=2.3,
    )

    landmarks = (
        (1, "Cuboctahedron", r"$t=1$", "FCC coordination shell"),
        (golden, "Regular icosahedron", r"$t=(\sqrt{5}-1)/2$",
         f"{carrier['hull_edge_count']} equal hull edges"),
        (half, "After the icosahedron", r"$t=1/2$", "Diagonals shorter than bars"),
        (0, "Octahedral limit", r"$t=0$", f"{limit['vertex_count']} sites, two joints each"),
    )
    open_points = np.asarray(jitterbug.scaled_jitterbug_vertices(1)) @ rotation.T
    highlight = max(triangles, key=lambda face: float(open_points[list(face), 2].mean()))
    for index, (parameter, title, label, detail) in enumerate(landmarks):
        center = 0.156 + index * 0.230
        if parameter == half:
            figure.add_artist(FancyBboxPatch(
                (center - 0.104, 0.476), 0.208, 0.194, boxstyle="round,pad=0.007,rounding_size=0.010",
                transform=figure.transFigure, facecolor=POST_ICO, edgecolor="none", zorder=-1,
            ))
        axis = figure.add_axes((center - 0.104, 0.544, 0.208, 0.138), facecolor="none")
        _pose(axis, parameter, rotation, bars, matching if 0 < float(parameter) < 1 else (),
              highlight=highlight if parameter else ())
        text(center, 0.541, title, 10.4, align="center", weight="bold",
             colour=ORANGE if parameter == half else INK)
        text(center, 0.519, label, 9.8, align="center", colour=MUTED)
        text(center, 0.496, detail, 9.1, align="center", colour=MUTED)
        if index < len(landmarks) - 1:
            text(center + 0.115, 0.619, "→", 10.5, align="center", colour=MUTED)

    rule(0.454)
    text(0.050, 0.439, "After the icosahedron, the graph and the weights matter", 14, weight="bold")
    text(0.050, 0.401, r"Compare two seam graphs at $t=1/2$", 11, weight="bold")
    comparison = figure.add_axes((0.046, 0.278, 0.161, 0.103), facecolor="none")
    _pose(comparison, half, rotation, (), matching)
    text(0.125, 0.281, "Nearest seams only", 8.7, align="center", colour=ORANGE)
    text(0.216, 0.369,
         f"Full hull  ·  {post_ico['hull_edge_count']} seams\nOne connected graph",
         9.9, colour=MUTED)
    text(0.216, 0.322,
         f"Nearest only  ·  {half_tuned['nearest_only_seams']} disjoint dimers\n"
         f"Precision rank {rank}, with {nullity} extra\nnull modes after centring",
         9.9, colour=MUTED)

    text(0.535, 0.401, r"Tune the diagonal weight $w=1/t-t$", 11.0, weight="bold", colour=TEAL)
    text(0.535, 0.368,
         "With unit bars, tuning makes fair binary seam\n"
         "repair’s coordinate lag covariances match a\n"
         "Gaussian Markov process.",
         9.9, colour=MUTED)
    text(0.535, 0.308,
         rf"At $t=1/2$, tuning needs $w={half_tuned['tuned_weight']}$ instead of $1$." "\n"
         rf"At $t=1/3$, inverse-square gives ${third_inverse_square['closing_weight']}$."
         rf" Tuning needs ${third_inverse_square['tuned_weight']}$.",
         9.6, colour=INK)

    rule(0.256)
    text(0.050, 0.240, "The paired limit has two different readouts", 14, weight="bold")
    text(0.050, 0.206, "Gaussian field", 11.3, weight="bold", colour=BLUE)
    text(0.050, 0.174, rf"Centred shared values have precision ${shared_bar_weight}L_{{\mathrm{{oct}}}}$.", 10.8)
    text(0.050, 0.142,
         "As the diagonal weight diverges, each pair\n"
         "shares one value. Two bars per octahedral\n"
         "edge give its doubled conductance.", 10.0, colour=MUTED)

    text(0.493, 0.206, "Fair binary repair", 11.3, weight="bold", colour=TEAL)
    text(0.493, 0.174,
         "Pair totals freeze on the ordinary attempt clock.\n"
         rf"On the faster $\lfloor ({len(bars)}+{len(matching)}w)s\rfloor$ clock, the averaged" "\n"
         rf"capacity-two process has coordinate lag factor $e^{{-{decay}s}}$.",
         10.0, colour=MUTED)
    text(0.493, 0.112,
         f"At finite weights, the {repair['pair_load_count']} possible pair-load states\n"
         "retain memory of which joints hold the loads.", 9.8, colour=INK)

    rule(0.061)
    text(0.050, 0.047, "Daniel Ari Friedman · Bernhard Mueller  |  Observers, Blankets, and Screens", 9, weight="bold")
    text(0.050, 0.027, "Fuller, Synergetics, fig. 460.08", 8, colour=BLUE,
         url="https://rwgrayprojects.com/synergetics/s04/figs/f6008.html")
    text(0.267, 0.027, "Stachel, The HEUREKA-Polyhedron", 8, colour=BLUE,
         url="https://www.geometrie.tuwien.ac.at/stachel/060_Heureka.pdf")
    text(0.951, 0.027, "github.com/docxology/OPH-FEP-Synergetics", 8, colour=BLUE, align="right", url=REPOSITORY)

    alt = (
        f"Fuller's Jitterbug is a linkage of {len(triangles)} rigid triangles and {len(bars)} fixed bars. "
        "At a common scale, the FCC shell contracts through the regular icosahedral OPH carrier "
        f"to a highlighted pose at t=1/2 and then to {limit['vertex_count']} paired octahedral sites. "
        f"At t=1/2 the hull has {post_ico['hull_edge_count']} seams. The nearest-only drawing shows "
        f"{half_tuned['nearest_only_seams']} disjoint dimers, precision rank {rank}, and {nullity} extra centred null modes. "
        "Matching repair and Gaussian coordinate lag covariance requires diagonal weight 1/t minus t, "
        f"including {half_tuned['tuned_weight']} at t=1/2. Inverse-square weighting fails at t=1/3. "
        "The Gaussian limit on shared site values has precision "
        f"{shared_bar_weight} times the octahedral Laplacian. Pair totals freeze on the ordinary attempt clock. "
        f"The accelerated, averaged capacity-two repair has coordinate lag exp(-{decay}s). "
        "At finite weights the pair totals retain memory of which joints carry the loads."
    )
    return figure, alt
