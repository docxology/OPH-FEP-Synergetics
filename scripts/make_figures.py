#!/usr/bin/env python3
"""Render the manuscript figures and their accessibility registry.

Thin orchestrator: every plotted value comes from the package's own functions,
so a figure cannot disagree with the prose or with the probes.  The script also
writes ``output/figures/figure_registry.json``, which the rendering
infrastructure reads to attach alternative text to each body figure.

Exit codes:
    0   every figure and the registry were written
    1   a figure could not be produced
"""

from __future__ import annotations

import json
import math
import os
import sys
from collections import Counter
from collections.abc import Iterable
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from oph_fep_synergetics import geometry  # noqa: E402
from oph_fep_synergetics.analysis.registry import (  # noqa: E402
    CORRESPONDENCES,
    FEP_SYNERGETICS,
    KINDS,
    OPH_FEP,
    PAIRS,
    SYNERGETICS_OPH,
    THREE_WAY,
    by_key,
    by_pair,
)
from oph_fep_synergetics.core import (  # noqa: E402
    centroid_dynamics,
    contrast_quotients,
    constraint_screening,
    controls,
    free_energy,
    geometric_screening,
    invariants,
    jitterbug_inference,
    native_dynamics,
    network,
    planar_duality,
    quadratic_modes,
    repair,
    root_shells,
    sheaf,
    surface_duality,
    symmetry,
    weighted_duality,
)
from oph_fep_synergetics.core.carrier import (  # noqa: E402
    OPH_ORIENTED_BASE_FACES,
    IcosahedralCarrier,
    icosahedral_carrier,
)
from oph_fep_synergetics.geometry import jitterbug, lattice  # noqa: E402
from oph_fep_synergetics.manuscript_variables import as_dict, format_value  # noqa: E402

FIGURES: Path = ROOT / "output" / "figures"

# Okabe-Ito colour-blind-safe palette, used for every figure except the dark cover.
BLUE = "#0072B2"
VERMILLION = "#D55E00"
GREEN = "#009E73"
ORANGE = "#E69F00"
SKY = "#56B4E9"
PURPLE = "#CC79A7"
YELLOW = "#F0E442"
BLACK = "#000000"
INK = "#1B2A3A"  # text and edge colour: a near-black, not a data colour
GREY = "#6B7280"  # neutral for context, never a status

KIND_COLOURS = {
    "isomorphism": BLUE,
    "disagreement": VERMILLION,
    "complementarity": GREEN,
}
STATUS_MARKERS = {"exact": "o", "exact-on-domain": "s", "obstructed": "X"}
ALT_TEXT: dict[str, str] = {}


def _counted(count: int, singular: str, plural: str) -> str:
    """Format a count with the singular or plural noun."""
    return f"{count} {singular if count == 1 else plural}"


def _enumerated(values: Iterable[object]) -> str:
    """Format a short bound sequence as "a, b and c" for alt text."""
    items = [format_value(item) for item in values]
    if len(items) < 3:
        return " and ".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


def _style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 200,
            "savefig.dpi": 250,
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "text.color": INK,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "axes.labelcolor": INK,
            "axes.edgecolor": "#718192",
            "axes.titlepad": 9,
            "axes.axisbelow": True,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "xtick.color": INK,
            "ytick.color": INK,
            "legend.fontsize": 10,
            "mathtext.default": "regular",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": "#A2AFBD",
            "grid.alpha": 0.32,
            "grid.linewidth": 0.6,
            "legend.frameon": False,
        }
    )


def _save(
    figure: plt.Figure, name: str, alt: str, *, facecolor: str = "white"
) -> Path:
    FIGURES.mkdir(parents=True, exist_ok=True)
    target = FIGURES / name
    figure.tight_layout()
    figure.savefig(target, bbox_inches="tight", facecolor=facecolor)
    plt.close(figure)
    ALT_TEXT[name] = alt
    return target


def carrier_figure() -> Path:
    """Connect the carrier's labelled incidence with its Laplacian modes."""
    carrier, view = _carrier_view()
    flat = view[:, :2]
    figure, (left, right) = plt.subplots(1, 2, figsize=(10.2, 4.4))
    front_edges: set[tuple[int, int]] = set()
    for face in OPH_ORIENTED_BASE_FACES:
        polygon = view[list(face)]
        normal = np.cross(polygon[1] - polygon[0], polygon[2] - polygon[0])
        if normal @ polygon.mean(axis=0) < 0:
            normal = -normal
        if normal[2] <= 0:
            continue
        depth = normal[2] / np.linalg.norm(normal)
        left.fill(polygon[:, 0], polygon[:, 1], color=SKY, alpha=0.15 + 0.3 * depth,
                  linewidth=0)
        front_edges.update(tuple(sorted((vertex, face[i - 1]))) for i, vertex in enumerate(face))
    for edge in carrier.edges:
        front = edge in front_edges
        x, y = flat[list(edge)].T
        left.plot(x, y, color=INK if front else GREY, lw=1.35 if front else 0.95,
                  ls="-" if front else "--", zorder=3 if front else 2)
    left.scatter(flat[:, 0], flat[:, 1], s=245, color="white", edgecolor=INK,
                 linewidth=1.25, zorder=4)
    for port, (x, y) in enumerate(flat):
        left.text(x, y, str(port), ha="center", va="center", fontsize=10, weight="bold",
                  color=INK, zorder=5)
    left.set_aspect("equal", adjustable="datalim")
    left.margins(0.10)
    left.set_axis_off()
    left.set_title("Icosahedral carrier", loc="left")
    left.set_xlabel("")
    left.text(0.5, -0.035, f"{carrier.port_count} ports  /  {carrier.seam_count} seams  /  "
              f"{carrier.face_count} faces", transform=left.transAxes, va="top", ha="center",
              color=INK)

    marks, counts = np.unique(np.round(carrier.laplacian_spectrum(), 10), return_counts=True)
    lam1_value = carrier.spectral_gap()
    normalised = lam1_value / (2.0 * carrier.seam_count)
    for mark, count in zip(marks, counts, strict=True):
        first_mode = math.isclose(mark, lam1_value)
        color = VERMILLION if first_mode else BLUE
        right.vlines(mark, 0, count, color=color, linewidth=2.4)
        right.scatter(mark, count, marker="D" if first_mode else "o", color=color,
                      edgecolor="white", linewidth=0.8, zorder=3, s=62)
        right.annotate(f"{count} mode{'s' if count != 1 else ''}", (mark, count),
                       textcoords="offset points", xytext=(0, 7), ha="center",
                       fontsize=9.5, color=INK)
    right.annotate("", (0.0, 0.45), (lam1_value, 0.45),
                   arrowprops={"arrowstyle": "<->", "color": VERMILLION, "lw": 1.25})
    right.text(lam1_value / 2, 0.62, "first gap", color=VERMILLION, ha="center", fontsize=9.5)
    right.text(0.03, 0.97, f"Normalised gap  {normalised:.6f}", transform=right.transAxes,
               va="top", fontsize=10.5, color=INK)
    right.set_xticks(marks, ("0", r"$5-\sqrt{5}$", "6", r"$5+\sqrt{5}$"))
    right.set_yticks(range(int(counts.max()) + 1))
    right.set_xlabel("Laplacian eigenvalue")
    right.set_ylabel("multiplicity")
    right.set_ylim(0, float(counts.max()) + 1.8)
    right.set_xlim(-0.85, float(marks.max()) + 0.85)
    right.set_title("Spectrum and first gap", loc="left")

    return _save(
        figure,
        "fig_carrier.png",
        f"Left: the {carrier.port_count} labelled ports in the carrier's spectral embedding, "
        f"using foreground shading and dashed rear seams to show depth. The {carrier.seam_count} "
        f"seams bound {carrier.face_count} triangles. Each port meets {carrier.degree(0)} others. "
        "Right: the Laplacian spectrum, with the zero eigenvalue once, five minus root five "
        "three times, six five times and five plus root five three times. The first positive "
        f"mode is highlighted, and the normalised gap is {normalised:.6f}.",
    )


def _rotation(theta: float, phi: float) -> np.ndarray:
    about_x = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, math.cos(theta), -math.sin(theta)],
            [0.0, math.sin(theta), math.cos(theta)],
        ]
    )
    about_y = np.array(
        [
            [math.cos(phi), 0.0, math.sin(phi)],
            [0.0, 1.0, 0.0],
            [-math.sin(phi), 0.0, math.cos(phi)],
        ]
    )
    return about_y @ about_x


def _edge_clearance(flat: np.ndarray, edges: np.ndarray) -> float:
    """The smallest distance from a port to an edge it does not belong to.

    Two edges that nearly coincide in the projection, or a port that sits on an
    edge, both drive this to zero, so maximising it picks a viewing angle with
    no doubled edges.
    """
    start = flat[edges[:, 0]]
    stop = flat[edges[:, 1]]
    direction = stop - start
    length_squared = np.einsum("ij,ij->i", direction, direction)
    if np.any(length_squared <= np.finfo(float).eps):
        return 0.0
    offset = flat[:, None, :] - start[None, :, :]
    fraction = np.clip(np.einsum("vej,ej->ve", offset, direction) / length_squared, 0.0, 1.0)
    nearest = start[None, :, :] + fraction[:, :, None] * direction[None, :, :]
    distance = np.linalg.norm(flat[:, None, :] - nearest, axis=2)
    own = (np.arange(flat.shape[0])[:, None] == edges[None, :, 0]) | (
        np.arange(flat.shape[0])[:, None] == edges[None, :, 1]
    )
    return float(np.min(np.where(own, np.inf, distance)))


def _clearest_view(points: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Rotate a solid to separate projected vertices from unrelated edges."""
    best_score = -1.0
    best_view = points
    for theta in np.linspace(0.0, math.pi, 36, endpoint=False):
        for phi in np.linspace(0.0, math.pi, 18, endpoint=False):
            view = points @ _rotation(theta, phi).T
            score = _edge_clearance(view[:, :2], edges)
            if score > best_score:
                best_score = score
                best_view = view
    return best_view


def _carrier_view() -> tuple[IcosahedralCarrier, np.ndarray]:
    """Return the carrier and its rotated spectral embedding, including depth."""
    carrier = icosahedral_carrier()
    eigenvalues, eigenvectors = np.linalg.eigh(carrier.laplacian())
    lam1 = 5.0 - math.sqrt(5.0)
    columns = [index for index, value in enumerate(eigenvalues) if abs(value - lam1) < 1e-9]
    if len(columns) != 3:
        raise RuntimeError(f"expected a three-dimensional lambda_1 eigenspace, got {len(columns)}")
    points = eigenvectors[:, columns] * math.sqrt(carrier.port_count)
    return carrier, _clearest_view(points, np.array(carrier.edges))


def _projected_carrier() -> tuple[IcosahedralCarrier, np.ndarray]:
    """Return the carrier and its planar spectral projection."""
    carrier, view = _carrier_view()
    return carrier, view[:, :2]


def registry_figure() -> Path:
    """List all registry keys with kind, contiguous status bands and status markers."""
    kinds = KINDS
    columns = {
        kind: sorted(
            (entry for entry in CORRESPONDENCES if entry.kind == kind),
            key=lambda entry: (entry.status, entry.key),
        )
        for kind in kinds
    }
    rows = max(len(entries) for entries in columns.values())
    figure, axes = plt.subplots(1, 3, figsize=(10.8, 5.7), width_ratios=(1, 1.12, 1))
    for axis, kind in zip(axes, kinds):
        entries = columns[kind]
        color = KIND_COLOURS[kind]
        status_groups = {status: index for index, status in enumerate(dict.fromkeys(
            entry.status for entry in entries
        ))}
        axis.axhspan(-1.65, -0.56, color=color, alpha=0.14, linewidth=0)
        axis.text(0.05, -1.105, kind.capitalize(), weight="bold", va="center", color=INK,
                  fontsize=12)
        axis.text(0.95, -1.105, str(len(entries)), weight="bold", va="center", ha="right",
                  color=INK, fontsize=13)
        for row, entry in enumerate(entries):
            shade = 0.09 if status_groups[entry.status] % 2 == 0 else 0.045
            axis.axhspan(row - 0.5, row + 0.5, color=color, alpha=shade, linewidth=0)
            if row > 0 and entry.status != entries[row - 1].status:
                axis.axhline(row - 0.5, color=color, linewidth=0.9, alpha=0.65)
            axis.scatter(
                [0.065],
                [row],
                marker=STATUS_MARKERS.get(entry.status, "D"),
                s=44,
                color=color,
                edgecolor=INK,
                linewidth=0.45,
                zorder=3,
            )
            axis.text(0.13, row, entry.key, va="center", fontsize=9.4, family="monospace",
                      color=INK)
        axis.set_xlim(0, 1)
        axis.set_ylim(rows + 1.55, -1.72)
        axis.set_axis_off()
    present = sorted({entry.status for entry in CORRESPONDENCES})
    handles = [
        plt.Line2D(
            [],
            [],
            marker=STATUS_MARKERS.get(status, "D"),
            linestyle="",
            markerfacecolor="white",
            markeredgecolor=INK,
            markeredgewidth=1.2,
            label=status,
            markersize=8,
        )
        for status in present
    ]
    figure.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.018),
        ncol=len(handles),
        fontsize=10.5,
    )
    figure.suptitle(f"{len(CORRESPONDENCES)} executable comparisons", fontsize=15, weight="bold")
    kind_counts = {kind: len(columns[kind]) for kind in kinds}
    return _save(
        figure,
        "fig_registry.png",
        f"The {len(CORRESPONDENCES)} correspondences arranged in three columns by kind. "
        f"{_counted(kind_counts['isomorphism'], 'isomorphism', 'isomorphisms')} on the left, "
        f"{_counted(kind_counts['disagreement'], 'disagreement', 'disagreements')} in the middle, "
        f"{_counted(kind_counts['complementarity'], 'complementarity', 'complementarities')} on "
        "the right. Each labelled registry key sits in a shaded group and carries a marker for its "
        "status: a circle for exact, a square for exact-on-domain, and a cross for obstructed. "
        "Column headings give the totals by kind.",
    )


def descent_figure() -> Path:
    """Show how often energy changes and where one unit can move at fixed energy."""
    census = controls.gated_descent_census()
    excursion = controls.mismatch_one_excursion(steps=6)
    figure, (left, right) = plt.subplots(1, 2, figsize=(10.7, 3.8))
    labels = ("Strict descent", "Zero descent", "No-op")
    values = [census["strict"], census["zero"], census["noop"]]
    total = sum(values)
    bars = left.barh(labels, values, color=[GREEN, ORANGE, GREY], height=0.54)
    for bar, value in zip(bars, values):
        left.annotate(
            f"{value:,}  ({value / total:.1%})",
            (value, bar.get_y() + bar.get_height() / 2),
            textcoords="offset points",
            xytext=(7, 0),
            va="center",
            fontsize=10,
            color=INK,
        )
    left.invert_yaxis()
    left.set_axisbelow(True)
    left.grid(False, axis="y")
    left.set_xlim(0, max(values) * 1.42)
    left.set_xlabel("Native attempts (count)")
    left.set_title(f"{total:,} attempts on the load shell", loc="left")

    ports = [state.index(1) for state in excursion]
    steps = np.arange(len(ports))
    energy = 0.5 * sum(value * value for value in excursion[0])
    right.plot(steps, ports, color=BLUE, linestyle=":", linewidth=1.55, zorder=2)
    right.scatter(steps, ports, s=58, color=BLUE, edgecolor="white", linewidth=1.1, zorder=3)
    right.set_xticks(steps)
    right.set_xlim(-0.35, len(steps) - 0.65)
    right.set_yticks(sorted(set(ports)))
    right.set_ylim(min(ports) - 0.55, max(ports) + 0.55)
    right.set_xlabel("Repair step")
    right.set_ylabel("Occupied port")
    right.set_title(f"One unit moves at fixed energy {energy:g}", loc="left")
    return _save(
        figure,
        "fig_descent.png",
        f"Left: counts and percentages for {values[0]} strict descents, {values[1]} zero "
        f"descents and {values[2]} no-ops out of {total} attempts. Right: discrete steps of "
        f"a one-unit excursion visit ports {_enumerated(ports)} in order. Dotted segments "
        f"guide the eye between consecutive steps. The load energy remains {energy:g}.",
    )


def dstability_figure() -> Path:
    """Eigenvalue loci and the sign of the third Routh entry."""
    values = as_dict()
    samples = controls.D_STABILITY_SAMPLES
    onset = float(values["d_stability.routh_onset"])
    k4_third = float(controls.routh_first_column(4.0)[2])
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.1, 4.1), width_ratios=(1.15, 1))
    spectra = {k: controls.d_stability_eigenvalues(k) for k in samples}
    reals = np.concatenate([spectra[k].real for k in samples])
    imag_max = max(float(np.abs(spectrum.imag).max()) for spectrum in spectra.values())
    x_limit = max(0.5, 2 * float(reals.max()))
    left.axvspan(0, x_limit, color=VERMILLION, alpha=0.08, linewidth=0)
    left.text(
        0.96, 0.96, "Unstable side", transform=left.transAxes, ha="right", va="top",
        fontsize=9.5,
        color=VERMILLION,
    )
    cmap = plt.get_cmap("viridis")
    handles = []
    for k, position in zip(samples, np.linspace(0.14, 0.84, len(samples)), strict=True):
        colour = cmap(position)
        eigenvalues = spectra[k]
        stable = eigenvalues.real <= 1e-9
        left.scatter(
            eigenvalues.real[stable],
            eigenvalues.imag[stable],
            s=62,
            color=colour,
            edgecolor=INK,
            linewidth=0.7,
            zorder=3,
        )
        left.scatter(
            eigenvalues.real[~stable],
            eigenvalues.imag[~stable],
            s=85,
            marker="X",
            color=colour,
            edgecolor=INK,
            linewidth=0.9,
            zorder=3,
        )
        handles.append(
            plt.Line2D([], [], marker="o", linestyle="", color=colour, markeredgecolor=INK,
                       markersize=7, label=f"{k:g}")
        )
    left.axvline(0.0, color=BLACK, linewidth=1.0, linestyle="--")
    left.set_xscale("symlog", linthresh=0.25)
    left.set_xlim(float(reals.min()) * 1.35, x_limit)
    left.set_ylim(-imag_max * 1.14, imag_max * 1.30)
    left.set_xticks([-6, -2, -0.5, 0, 0.25, 0.5])
    left.set_xticklabels(["−6", "−2", "−0.5", "0", "0.25", "0.5"])
    left.set_xlabel("Real part (symlog scale)")
    left.set_ylabel("Imaginary part")
    left.set_title(r"Eigenvalues of $D_kA$", loc="left")
    left.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.24),
        ncol=len(samples),
        title=r"Sampled couplings $k$",
        fontsize=9,
        handletextpad=0.25,
        columnspacing=0.65,
    )

    grid = np.linspace(0.05, max(samples), 400)
    third = np.array([controls.routh_first_column(k)[2] for k in grid])
    right.fill_between(grid, 0, third, where=third >= 0, color=BLUE, alpha=0.12)
    right.fill_between(grid, 0, third, where=third <= 0, color=VERMILLION, alpha=0.12)
    right.plot(grid, third, color=INK, linewidth=2.0)
    right.axhline(0.0, color=BLACK, linewidth=0.9)
    right.axvline(onset, color=VERMILLION, linestyle="--", linewidth=1.1)
    right.annotate(
        r"$k = \sqrt{10}$",
        xy=(onset, 0.0),
        xytext=(onset + 0.3, 4),
        fontsize=10.5,
        color=VERMILLION,
        arrowprops={"arrowstyle": "->", "color": VERMILLION, "lw": 1.0},
    )
    right.scatter([4.0], [k4_third], color=VERMILLION, zorder=3, s=40)
    right.annotate(
        f"$k=4$ gives ${k4_third:g}$".replace("-", "−"),
        xy=(4.0, k4_third),
        xytext=(4.3, 1.35),
        fontsize=10,
        color=VERMILLION,
        arrowprops={"arrowstyle": "->", "color": VERMILLION, "lw": 1.0},
    )
    right.set_xlabel("Coupling $k$")
    right.set_ylabel("Third Routh entry")
    right.set_ylim(float(third.min()) * 1.12, float(third.max()) * 1.10)
    right.set_title("Stability is lost when the entry crosses zero", loc="left")
    return _save(
        figure,
        "fig_dstability.png",
        "Two panels. Left: the eigenvalues of D sub k times A in the complex plane for "
        f"{len(samples)} couplings. The real axis uses a symmetric log scale to open the region "
        "near zero. Coupling is distinguished by colour and the positive-real-part samples use "
        "cross markers on the shaded side. "
        "Right: the third Routh entry, twice the quantity ten minus k squared divided by "
        "the quantity k plus two, against the "
        f"coupling, crossing zero at {format_value(values['d_stability.routh_onset'])} and equal "
        f"to {format_value(values['d_stability.routh_k4_third'])} at k equal to four.",
    )


def mediator_figure() -> Path:
    """Compare the mediator feasibility threshold with an independent reduction control."""
    bound = as_dict()
    threshold = float(bound["mediator.threshold"])
    below = float(bound["mediator.epsilon_below"])
    margin_below = float(bound["mediator.margin_below"])
    figure, (left, right) = plt.subplots(1, 2, figsize=(10.6, 3.7))
    grid = np.linspace(0.0, 1.0, 200)
    margins = [controls.mediator_pairwise_margin(value) for value in grid]
    left.axvspan(0.0, threshold, color=VERMILLION, alpha=0.10, lw=0)
    left.axvspan(threshold, float(grid[-1]), color=GREEN, alpha=0.085, lw=0)
    left.text(0.025, 0.95, "Pairwise form\nfails", transform=left.transAxes, va="top",
              fontsize=9.5, color=VERMILLION)
    left.text(0.975, 0.08, "Pairwise form\nexists", transform=left.transAxes, ha="right",
              fontsize=9.5, color=INK)
    left.plot(grid, margins, color=BLUE, linewidth=2.0)
    left.axhline(0.0, color=BLACK, linewidth=0.9)
    left.axvline(threshold, color=VERMILLION, linestyle="--", linewidth=1.0)
    left.scatter([below, threshold], [margin_below, 0.0], color=VERMILLION, edgecolor="white",
                 linewidth=0.9, s=65, zorder=3)
    fraction = Fraction(threshold).limit_denominator(100)
    left.annotate(
        f"$\\varepsilon={fraction}$",
        xy=(threshold, 0.0),
        xytext=(threshold + 0.1, -0.64),
        fontsize=10.5,
        color=VERMILLION,
        arrowprops={"arrowstyle": "->", "color": VERMILLION, "lw": 1.0},
    )
    left.annotate(
        f"$\\varepsilon={below:g}$\nmargin {margin_below:.1f}".replace("-", "−"),
        xy=(below, margin_below),
        xytext=(0.035, 0.75),
        fontsize=10,
        color=VERMILLION,
        arrowprops={"arrowstyle": "->", "color": VERMILLION, "lw": 1.0},
    )
    left.set_xlim(0, float(grid[-1]))
    left.set_xlabel(r"Mediator coupling $\varepsilon$")
    left.set_ylabel(r"Pairwise margin $3\varepsilon-1$")
    left.set_title("Mediator feasibility", loc="left")

    reduction = controls.reduction_control()
    first = reduction["first_after"]
    second = reduction["second_after"]
    floor = reduction["worst_case_error_lower_bound"]
    midpoint = (first + second) / 2
    right.axhline(0, color=GREY, linewidth=1)
    for value, history, color in ((first, "History A", BLUE), (second, "History B", VERMILLION)):
        right.scatter(value, 0, s=110, marker="o", color=color, edgecolor="white",
                      linewidth=1.0, zorder=3)
        right.text(value, -0.18, history, color=color, ha="center", fontsize=10, weight="bold")
        right.annotate("", (value, 0.35), (midpoint, 0.35),
                       arrowprops={"arrowstyle": "|-|", "color": color, "lw": 1.2})
        right.text((midpoint + value) / 2, 0.43, f"Error {floor:g}", ha="center",
                   fontsize=10, color=color)
    right.scatter(midpoint, 0, s=110, marker="D", facecolor="white", edgecolor=INK,
                  linewidth=1.3, zorder=3)
    right.text(midpoint, -0.18, "Shared prediction", ha="center", color=INK, fontsize=10)
    right.text(midpoint, 0.75, f"True responses are {reduction['response_distance']:g} apart",
               ha="center", color=INK, fontsize=10.5)
    margin = 0.22 * reduction["response_distance"]
    right.set_xlim(min(first, second) - margin, max(first, second) + margin)
    right.set_ylim(-0.4, 1)
    right.set_xticks([second, midpoint, first])
    right.set_yticks([])
    right.spines["left"].set_visible(False)
    right.grid(False)
    right.set_xlabel("Response to the same intervention")
    right.set_title("Two histories share one observation", loc="left")
    return _save(
        figure,
        "fig_mediator.png",
        "Two panels. Left: the pairwise margin three epsilon minus one against epsilon, crossing "
        f"zero at the threshold {fraction}, with the point at epsilon equal to {below:g} marked at "
        f"{margin_below:.1f}. Right: a number line for the separate reduction control. Two "
        f"histories share an observation and produce responses {first:g} and {second:g} to the "
        f"same intervention. The best shared prediction lies at their midpoint and still "
        f"incurs absolute error {floor:g} on both histories.",
    )


def blanket_figure() -> Path:
    """The drift that a stationary law does not fix, and free entropy production."""
    bound = as_dict()
    q = float(bound["blanket.q"])
    drift = controls.blanket_dynamics_check(q)
    figure, (left, right) = plt.subplots(1, 2, figsize=(10.7, 3.85))
    names = ("Exterior", "Blanket", "Internal")
    values = [
        float(drift["mu_drift_external"]),
        float(drift["mu_drift_blanket"]),
        float(drift["mu_drift_internal"]),
    ]
    for index, value in enumerate(values):
        exterior = index == 0
        color = VERMILLION if exterior else BLUE
        left.hlines(index, 0, value, color=color, linewidth=2.4)
        left.scatter(value, index, s=76, color=color, marker="D" if exterior else "o",
                     edgecolor="white", linewidth=0.8, zorder=3)
        left.annotate(
            f"{value:+g}".replace("-", "−"),
            (value, index),
            textcoords="offset points",
            xytext=(7 if value >= 0 else -7, 0),
            ha="left" if value >= 0 else "right",
            va="center",
            fontsize=10.5,
            color=INK,
        )
    left.axvline(0.0, color=BLACK, linewidth=0.9)
    span = max(values) - min(values)
    left.set_xlim(min(values) - 0.30 * span, max(values) + 0.30 * span)
    left.set_yticks(range(len(names)), names)
    left.set_ylim(2.55, -0.65)
    left.grid(False, axis="y")
    left.set_xlabel(f"Coefficient in internal drift at $q={q:g}$")
    left.set_title("The exterior drives the interior", loc="left")

    marks = [float(q) for q in bound["solenoidal.q_grid"]]
    grid = np.linspace(min(marks), max(marks), 200)
    production = [controls.solenoidal_check(value)["entropy_production"] for value in grid]
    right.fill_between(grid, 0, production, color=BLUE, alpha=0.12, linewidth=0)
    right.plot(grid, production, color=BLUE, linewidth=2.0)
    marked = [controls.solenoidal_check(value)["entropy_production"] for value in marks]
    right.scatter(marks, marked, marker="D", color=VERMILLION, edgecolor="white",
                  linewidth=0.8, zorder=3, s=60)
    for value, height in zip(marks, marked, strict=True):
        right.annotate(
            f"{height:.3g}",
            (value, height),
            textcoords="offset points",
            xytext=(-7, 8),
            ha="right",
            fontsize=10.5,
            color=VERMILLION,
        )
    right.set_xlim(-0.075 * max(marks), 1.14 * max(marks))
    right.set_ylim(-0.04 * max(marked), 1.13 * max(marked))
    right.set_xlabel("Solenoidal parameter $q$")
    right.set_ylabel(r"Entropy production $2q^2$")
    right.set_title("Entropy production varies at fixed density", loc="left")
    return _save(
        figure,
        "fig_blanket.png",
        f"Left: three signed coefficients of the internal drift at q equal to "
        f"{format_value(bound['blanket.q'])}, namely {_enumerated(bound['blanket.drift'])}, "
        "in exterior, blanket and internal order. A diamond marks the direct exterior "
        "coefficient while the stationary law gives conditional independence of exterior "
        "and interior. Right: entropy production against the solenoidal parameter q, "
        "rising as two q squared, with the marked values "
        f"{_enumerated(bound['solenoidal.entropy_production'])} at q equal to "
        f"{_enumerated(bound['solenoidal.q_grid'])}.",
    )


def orbit_figure() -> Path:
    """The relaxation orbit and its radius against time."""
    samples, sample_dt = controls.two_node_relaxation_trajectory()
    radius = np.linalg.norm(samples, axis=1)
    measured = controls.two_node_relaxation_orbit()
    r_min = float(measured["min_radius"])
    r_max = float(measured["max_radius"])
    period = float(measured["period"])
    window = 5.0 * period
    sync = controls.two_node_sync_orbit(steps=math.ceil(window))
    sync_radius = np.linalg.norm(sync, axis=1)
    sync_floor = controls.two_node_orbit_min_radius()

    figure = plt.figure(figsize=(10.5, 5.0))
    grid = figure.add_gridspec(2, 2, width_ratios=(1.0, 1.35))
    left = figure.add_subplot(grid[:, 0])
    top = figure.add_subplot(grid[0, 1])
    bottom = figure.add_subplot(grid[1, 1])

    theta = np.linspace(0.0, 2.0 * math.pi, 400)
    for bound, fill, style in (
        (r_max, matplotlib.colors.to_rgba(SKY, 0.18), "--"),
        (r_min, left.get_facecolor(), "-"),
    ):
        left.fill(bound * np.cos(theta), bound * np.sin(theta), color=fill, lw=0,
                  zorder=1)
        left.plot(
            bound * np.cos(theta),
            bound * np.sin(theta),
            color=GREY,
            linestyle=style,
            linewidth=1.1,
            zorder=2,
        )
    left.plot(samples[:, 0], samples[:, 1], color=BLUE, linewidth=1.8, zorder=3)
    arrow_step = round(0.025 * period / sample_dt)
    left.annotate("", xy=samples[2 * arrow_step], xytext=samples[arrow_step],
                  arrowprops={"arrowstyle": "-|>", "color": BLUE, "lw": 1.6},
                  zorder=4)
    left.scatter([0.0], [0.0], color=VERMILLION, s=67, marker="X", zorder=4)
    left.annotate(
        "Unstable\nfixed point",
        xy=(0.0, 0.0), xytext=(0.0, -0.33), ha="center", va="top",
        fontsize=10.5, color=VERMILLION,
        arrowprops={"arrowstyle": "-", "color": VERMILLION, "lw": 1.1, "shrinkB": 5},
    )
    left.set_aspect("equal")
    left.set_xlim(-1.17 * r_max, 1.17 * r_max)
    left.set_ylim(-1.17 * r_max, 1.17 * r_max)
    left.set_xlabel("x")
    left.set_ylabel("y")
    left.set_title("Relaxation orbit in the phase plane", loc="left")

    times = np.arange(radius.size) * sample_dt
    keep = times <= window
    top.axhspan(r_min, r_max, color=SKY, alpha=0.18, lw=0)
    top.plot(times[keep], radius[keep], color=BLUE, linewidth=1.65)
    top.axhline(r_min, color=GREY, linewidth=0.9)
    top.axhline(r_max, color=GREY, linestyle="--", linewidth=0.9)
    arrow_y = r_min - 0.085
    top.annotate(
        "",
        xy=(period, arrow_y),
        xytext=(0.0, arrow_y),
        arrowprops={"arrowstyle": "|-|", "color": INK, "lw": 1.1, "shrinkA": 0, "shrinkB": 0},
    )
    top.text(period * 1.12, arrow_y, f"Full phase-plane return ≈ {period:.3f}", va="center",
             fontsize=10, color=INK)
    top.set_xlim(0.0, window)
    top.set_ylim(r_min - 0.16, r_max + 0.08)
    top.set_yticks([r_min, r_max], [f"{r_min:.3f}", f"{r_max:.3f}"])
    top.set_xlabel("Continuous time $t$")
    top.set_ylabel("Radius")
    top.set_title("Continuous relaxation radius", loc="left")

    steps = np.arange(sync_radius.size)
    show = steps <= window
    bottom.plot(steps[show], sync_radius[show], color=VERMILLION, marker="s", markersize=4.4,
                linestyle=":", linewidth=1.0, markeredgecolor=INK, markeredgewidth=0.4)
    bottom.axhline(sync_floor, color=GREY, linestyle="--", linewidth=1.0)
    bottom.set_xlim(-0.5, window + 0.5)
    bottom.set_xlabel("Discrete iteration $n$")
    bottom.set_ylabel("Radius")
    bottom.set_title("Synchronous map radius", loc="left")
    bottom.text(0.98, 0.21, f"Post-transient floor {sync_floor:g}", transform=bottom.transAxes,
                ha="right", color=INK, fontsize=10)
    return _save(
        figure,
        "fig_orbit.png",
        "Three panels. Left: the relaxation orbit in the plane, a closed curve that neither "
        "escapes nor reaches the marked unstable fixed point at the origin. An arrow shows the "
        "direction of time. The shaded annulus "
        f"ranges from radius {r_min:.4f} to {r_max:.4f}. Top right: the continuous relaxation "
        f"radius over five full phase-plane periods, with a return-time estimate of {period:.3f} "
        "marked. The radius oscillates several times within each return. Bottom "
        f"right: samples of the synchronous clipped map by discrete iteration, showing a radius "
        f"floor of {sync_floor:g} after the transient.",
    )


def _plot_tetravolumes(
    axis: plt.Axes, volumes: dict[str, Fraction], *, icosahedron: float | None = None
) -> None:
    """Plot tetravolumes by position on a positive logarithmic axis."""
    entries = [
        (name.replace("_", " "), float(volume), str(volume), "o", BLUE)
        for name, volume in volumes.items()
    ]
    if icosahedron is not None:
        entries.append(("carrier icosahedron", icosahedron, f"≈ {icosahedron:.4f}", "D", VERMILLION))
    entries.sort(key=lambda entry: entry[1])
    for row, (_, value, label, marker, colour) in enumerate(entries):
        axis.axhline(row, color=GREY, linewidth=0.55, alpha=0.22, zorder=0)
        axis.scatter(value, row, s=60, marker=marker, color=colour, edgecolor="white",
                     linewidth=0.65, zorder=3)
        axis.annotate(
            label, (value, row), textcoords="offset points", xytext=(8, 0),
            va="center", color=INK, fontsize=10,
        )
    axis.set_xscale("log")
    axis.set_xlim(0.75, max(entry[1] for entry in entries) * 1.8)
    axis.set_yticks(range(len(entries)), [entry[0] for entry in entries], fontsize=10)
    axis.set_ylim(len(entries) - 0.35, -0.65)
    axis.set_xticks((1, 3, 10, 30, 100), ("1", "3", "10", "30", "100"))
    axis.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    axis.grid(False)
    axis.tick_params(axis="y", length=0)
    axis.spines["left"].set_visible(False)
    axis.set_xlabel("Tetravolume (log scale)")


def ivm_figure() -> Path:
    """IVM shell sizes and tetravolume ratios."""
    shells = list(range(0, 9))
    counts = [lattice.shell_count(index) for index in shells]
    cumulative = [lattice.cumulative_count(index) for index in shells]

    figure, (left, right) = plt.subplots(
        1, 2, figsize=(10.5, 4.15), gridspec_kw={"width_ratios": (1.35, 1.0)}
    )
    colours = [BLUE] * len(shells)
    colours[1] = VERMILLION
    bars = left.bar(shells, counts, color=colours, width=0.68, zorder=2)
    for bar, value in zip(bars, counts):
        left.annotate(
            str(value),
            (bar.get_x() + bar.get_width() / 2, value),
            textcoords="offset points",
            xytext=(0, 3),
            ha="center",
            fontsize=10,
            color=INK,
        )
    left.set_xticks(shells, [f"{k}\n{total:,}" for k, total in zip(shells, cumulative)])
    left.set_xlabel("Shell $k$ and cumulative sites through $k$")
    left.set_ylabel("Sites in shell")
    left.set_title(r"Coordination shells  ·  $10k^2+2$ for $k\geq 1$", loc="left")
    left.set_ylim(0, max(counts) * 1.16)
    left.grid(False, axis="x")
    left.annotate(
        f"{counts[1]} nearest neighbours", xy=(1, counts[1]), xytext=(0.03, 0.44),
        textcoords="axes fraction", color=INK, fontsize=10,
        arrowprops={"arrowstyle": "-", "color": VERMILLION, "linewidth": 0.9,
                    "shrinkB": 14},
    )

    _plot_tetravolumes(right, geometry.derived_volume_ratios())
    right.set_title("Volumes relative to the IVM unit tetrahedron", loc="left")
    return _save(
        figure,
        "fig_ivm.png",
        "Two panels. Left: the site count of each IVM shell, ten k squared plus two for shells "
        "one onward, as bars for shells zero to eight, with the cumulative count printed under "
        "each shell number. The first shell is highlighted as the twelve nearest neighbours. Right: tetravolume "
        "ratios as dots on a logarithmic scale, with the tetrahedron at one, the cube at three, the "
        "octahedron at four, the rhombic dodecahedron at six, the cuboctahedron at twenty and "
        "the truncated octahedron at ninety-six.",
    )


def subdivision_figure() -> Path:
    """The same quadratic sequence in three places."""
    frequencies = tuple(range(1, 7))
    vertices = [lattice.geodesic_vertex_count(value) for value in frequencies]
    faces = [lattice.geodesic_face_count(value) for value in frequencies]
    if vertices != [lattice.goldberg_face_count(value) for value in frequencies] or vertices != [
        lattice.shell_count(value) for value in frequencies
    ] or faces != [lattice.goldberg_vertex_count(value) for value in frequencies]:
        raise ValueError("the plotted subdivision counts must coincide within each family")

    figure, axis = plt.subplots(figsize=(9.0, 4.15))
    for counts, marker, colour, style, offset, label in (
        (faces, "s", VERMILLION, "--", 9, r"$20f^2$   Geodesic faces = Goldberg vertices"),
        (vertices, "o", BLUE, "-", -13, r"$10f^2+2$   Geodesic vertices = Goldberg faces = IVM shell sites"),
    ):
        # The labels identify both first endpoints where the markers would overlap.
        axis.plot(frequencies, counts, style, marker=marker, markersize=7, color=colour,
                  markevery=slice(1, None), linewidth=1.9, markeredgecolor="white",
                  label=label, zorder=3)
        for value, count in zip(frequencies, counts):
            axis.annotate(str(count), (value, count), textcoords="offset points", xytext=(0, offset),
                          ha="center", fontsize=10, color=INK)
    axis.set_xticks(frequencies)
    axis.set_xlim(0.63, 6.37)
    axis.set_ylim(-42, max(faces) * 1.12)
    axis.set_yticks(np.arange(0, max(faces) * 1.12, 150))
    axis.set_xlabel("Subdivision frequency or IVM shell index $f$")
    axis.set_ylabel("Count")
    axis.set_title("Geodesic and Goldberg counts exchange vertices and faces", loc="left")
    axis.grid(False, axis="x")
    axis.legend(loc="upper left", fontsize=10)
    return _save(
        figure,
        "fig_subdivision.png",
        "A line chart against subdivision frequency from one to six. Geodesic vertex count, "
        "Goldberg face count and the IVM shell count coincide on ten f squared plus two, shown "
        "by solid circles. Geodesic face count and Goldberg vertex count coincide on twenty f "
        "squared, shown by dashed squares. Every plotted count is labelled.",
    )


def _cover_solid(
    axis: plt.Axes,
    view: np.ndarray,
    faces: tuple[tuple[int, ...], ...],
    edges: tuple[tuple[int, int], ...],
    center: tuple[float, float],
    face_color: str,
    line_color: str,
    *,
    numbered: bool = False,
) -> None:
    """Shade convex faces and distinguish foreground and rear edges."""
    view = view + np.array([*center, 0.0])
    flat = view[:, :2]
    rgb = np.array(matplotlib.colors.to_rgb(face_color))
    light = np.array([-0.47, 0.56, 0.68])
    foreground_edges = set()
    foreground_ports = set()
    for face in faces:
        polygon = view[list(face)]
        normal = np.cross(polygon[1] - polygon[0], polygon[2] - polygon[0])
        normal /= np.linalg.norm(normal)
        if normal @ (polygon.mean(axis=0) - np.array([*center, 0.0])) < 0:
            normal = -normal
        if normal[2] <= 0:
            continue
        brightness = 0.32 + 0.62 * max(0.0, float(normal @ light))
        axis.fill(
            polygon[:, 0], polygon[:, 1],
            color=rgb * brightness, linewidth=0, zorder=2,
        )
        foreground_edges.update(
            tuple(sorted((port, face[index - 1]))) for index, port in enumerate(face)
        )
        foreground_ports.update(face)

    for edge in edges:
        if edge in foreground_edges:
            continue
        start, end = flat[list(edge)].T
        axis.plot(start, end, color=line_color, alpha=0.40, linewidth=1.05,
                  linestyle=(0, (2, 3)), zorder=3)
    for edge in sorted(foreground_edges):
        start, end = flat[list(edge)].T
        axis.plot(
            start, end, color=line_color, alpha=0.20, linewidth=5.4,
            solid_capstyle="round", zorder=4,
        )
        axis.plot(start, end, color=line_color, linewidth=1.55,
                  solid_capstyle="round", zorder=5)
    for port, (x, y) in enumerate(flat):
        foreground = port in foreground_ports
        axis.scatter(
            x, y, s=168 if numbered else 46, color="#0D202C",
            edgecolor=line_color if foreground else "#A8B8BF", linewidth=1.25, zorder=6,
        )
        if numbered:
            axis.text(x, y, str(port), ha="center", va="center", fontsize=9.5,
                      color="#F3FBFA", fontweight="bold", zorder=7)


def _cover_contrasts(figure: plt.Figure, width: float, height: float,
                     rows: Iterable[tuple], blue: str, gold: str) -> str:
    """Draw the two carrier distance ladders at a common relative scale."""
    steps, resistance, squared_chords, _ = zip(*rows, strict=True)
    chart = figure.add_axes((5.21 / width, 4.56 / height, 3.10 / width, 2.06 / height))
    chart.patch.set_alpha(0)
    for spine in ("top", "right"):
        chart.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        chart.spines[spine].set_color("#8AA0B4")
    chart.tick_params(colors="#DCE8F3", labelsize=12, length=3)
    chart.grid(axis="y", linestyle=":", color="#8AA0B4", alpha=0.44)
    chart.grid(False, axis="x")
    chart.set_xticks(steps)
    chart.set_yticks((1, 2, 3, 4))
    chart.set(xlim=(0.85, 3.65), ylim=(0.75, 4.04))
    chart.set_title("Adjacent-pair value = 1", color="#DCE8F3", fontsize=14, pad=9)
    chart.set_xlabel("Graph distance (seams)", color="#DCE8F3", fontsize=13, labelpad=6)
    for numbers, color, marker, style, fill, size in (
        (squared_chords, gold, "D", "-", "none", 8.5),
        (resistance, blue, "o", "--", blue, 5.2),
    ):
        relative = [float(number / numbers[0]) for number in numbers]
        chart.plot(steps, relative, linewidth=2.1, color=color, linestyle=style,
                   marker=marker, markersize=size, markerfacecolor=fill,
                   markeredgecolor=color, markeredgewidth=1.2)
        chart.annotate(f"{relative[-1]:.2f}×", (steps[-1], relative[-1]),
                       textcoords="offset points", xytext=(7, 0),
                       va="center", color=color, fontsize=11.5)
    return _enumerated(resistance)


def cover_figure() -> Path:
    """Summarise the finite models, quantitative comparisons and registry."""
    values = as_dict()
    kinds = Counter(entry.kind for entry in CORRESPONDENCES)
    carrier, view = _carrier_view()
    sites = lattice.cuboctahedron_neighbors()
    shell_edges = lattice.cuboctahedron_edges()
    shell_faces = lattice.cuboctahedron_faces()
    shell = _clearest_view(np.array([geometry.to_xyz(site) for site in sites]),
                          np.array(shell_edges))
    carrier_center = (2.48, 5.33)
    shell_center = (11.02, 5.33)
    width, height = 13.5, 9.3
    paper, secondary = "#F5F7F9", "#C0D1E1"
    teal, blue, gold = "#A6F5EB", "#A7CFFF", "#FFD2A1"
    figure, axis = plt.subplots(figsize=(width, height))
    figure.subplots_adjust(left=0, right=1, top=1, bottom=0)
    axis.set_axis_off()
    axis.set(xlim=(0, width), ylim=(0, height), aspect="equal")

    # Soft illumination behind the two solids. All strokes on the solids trace
    # actual edges of the carrier or of the close-packed coordination shell.
    x, y = np.meshgrid(np.linspace(0, width, 810), np.linspace(0, height, 558))
    field = np.broadcast_to(np.array([7, 16, 28]), (*x.shape, 3)).astype(float).copy()
    for (cx, cy), spread, color in (
        (carrier_center, 4.8, (10, 48, 58)),
        ((6.75, 5.55), 7.4, (16, 29, 51)),
        (shell_center, 4.8, (51, 26, 20)),
    ):
        glow = np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / spread)
        field += glow[..., None] * color
    axis.imshow(field.astype(np.uint8), origin="lower", extent=(0, width, 0, height),
                interpolation="bilinear", zorder=0)
    axis.text(0.52, 8.65, "Geometry, repair and Gaussian inference", va="center", color=paper,
              fontsize=28, fontweight="bold")
    axis.text(0.54, 8.17, "Two shells anchor the models. Transfer tests compare results across graphs.",
              va="center", color=secondary, fontsize=14.5)
    for level in (7.86, 3.11, 1.53):
        axis.hlines(level, 0.52, 12.98, color="#89A9BB", linewidth=0.85, alpha=0.48)
    axis.vlines((4.59, 8.91), 1.76, 7.56, color="#89A9BB", linewidth=0.7, alpha=0.38)

    _cover_solid(axis, view * (1.45 / np.linalg.norm(view[0])), carrier.faces, carrier.edges,
                 carrier_center, "#286D7D", teal, numbered=True)
    _cover_solid(axis, shell * (1.43 / np.linalg.norm(shell[0])),
                 shell_faces["triangles"] + shell_faces["squares"], shell_edges,
                 shell_center, "#A75D39", gold)
    contrast_rows = network.carrier_resistance_geometry()["shells"]
    resistance_alt = _cover_contrasts(figure, width, height, contrast_rows, blue, gold)
    edge_length = math.sqrt(float(contrast_rows[0][2]))
    for cx, color, label, title, finding, lines in (
        (carrier_center[0], teal, "OPH", "Icosahedral carrier",
         "Fair binary repair",
         ("Uniform resampling on a seam", "given its total and outside loads.",
          "Distance fixes the next-class law.")),
        (6.75, blue, "FREE ENERGY PRINCIPLE", "Gaussian port contrasts",
         "Gaussian contrast",
         (f"Chord² = {values['chord.variance_factor']} × lowest-mode variance",
          f"At carrier edge length {edge_length:g}.",
          f"Full covariance rank: {values['chord.full_rank']}.")),
        (shell_center[0], gold, "SYNERGETICS", "Cuboctahedral FCC shell",
         "FCC plane screening",
         ("Axis and diagonal planes separate halves.",
          "Fixing on-plane sites screens them.",
          "Conditioning on a global sum couples them.")),
    ):
        axis.text(cx, 7.53, label, ha="center", va="center", color=color,
                  fontsize=13, fontweight="bold")
        axis.text(cx, 7.11, title, ha="center", va="center", color=paper,
                  fontsize=18.5, fontweight="bold")
        axis.text(cx, 2.83, finding, ha="center", va="center", color=color,
                  fontsize=16, fontweight="bold")
        for y_value, line in zip((2.44, 2.11, 1.78), lines, strict=True):
            axis.text(cx, y_value, line, ha="center", va="center", color=paper, fontsize=13.5)

    for cx, y_value, color, label in (
        (carrier_center[0], 3.77, paper, f"{carrier.port_count} ports  ·  {carrier.seam_count} seams"),
        (carrier_center[0], 3.43, secondary,
         f"Load two  ·  {values['shell.total_two']} states, {values['shell.transitions']:,} attempts"),
        (6.75, 3.78, gold, "◇  Squared chord"),
        (6.75, 3.43, blue, "●  Full Gaussian contrast variance"),
        (shell_center[0], 3.77, paper, f"{len(sites)} FCC sites  ·  {len(shell_edges)} edges"),
        (shell_center[0], 3.43, secondary,
         f"{len(shell_faces['triangles'])} triangles  ·  {len(shell_faces['squares'])} squares"),
    ):
        axis.text(cx, y_value, label, ha="center", va="center", color=color, fontsize=14)

    total = sum(kinds.values())
    axis.text(0.54, 1.19, f"{total} executable claims", va="center", color=paper,
              fontsize=16.5, fontweight="bold")
    axis.text(12.96, 1.19, f"including {len(by_pair()[THREE_WAY])} three-way claims", ha="right",
              va="center", color=secondary, fontsize=14)
    start = 0.54
    for kind, color, label in (
        ("isomorphism", "#8EC9F9", "isomorphisms"),
        ("disagreement", "#FFD2A1", "disagreements"),
        ("complementarity", "#93E6CF", "complementarities"),
    ):
        count = kinds[kind]
        span = 12.42 * count / total
        axis.broken_barh([(start, span - 0.04)], (0.42, 0.15), facecolors=color, edgecolors="none")
        x, align = (12.96, "right") if kind == "complementarity" else (start, "left")
        axis.text(x, 0.80, f"{count} {label}", ha=align, va="center", color=color, fontsize=13)
        start += span

    filename = "fig_cover.png"
    FIGURES.mkdir(parents=True, exist_ok=True)
    target = FIGURES / filename
    figure.savefig(target, dpi=260, facecolor="#07101C")
    plt.close(figure)
    ALT_TEXT[filename] = (
        f"Three readings of the finite models. The OPH icosahedron has {carrier.port_count} "
        f"numbered ports and {carrier.seam_count} seams. Its total-two shell has "
        f"{values['shell.total_two']} states and {values['shell.transitions']:,} repair attempts. "
        f"The cuboctahedral FCC shell has {len(sites)} sites and {len(shell_edges)} edges. "
        "Solid and dashed lines show foreground and rear edges. A centre chart compares squared "
        "icosahedral chords with full Gaussian contrast variance by graph distance, both divided "
        f"by their adjacent-pair value. Squared chords rise faster. Full intrinsic variances are "
        f"{resistance_alt}. Squared chords "
        f"equal {values['chord.variance_factor']} times the lowest-mode variance at edge length "
        f"{edge_length:g}. On the binary carrier shell, fair repair resamples the selected seam "
        "uniformly given its total and all outside loads. With uniform seam sampling, occupied-port "
        "distance fixes the next distance distribution. Coordinate and signed body-diagonal planes "
        "separate the strict FCC halves at every positive shell frequency. Fixing the plane sites "
        "screens the grounded field. A global sum couples the halves. "
        f"The {total} claims include {kinds['isomorphism']} isomorphisms, {kinds['disagreement']} "
        f"disagreements and {kinds['complementarity']} complementarities. "
        f"{len(by_pair()[THREE_WAY])} claims use all three frameworks."
    )
    return target


def _kind_phrases(counts: Counter) -> list[str]:
    """The non-zero kind counts of a ledger, one phrase each."""
    names = (
        ("isomorphism", "isomorphism", "isomorphisms"),
        ("disagreement", "disagreement", "disagreements"),
        ("complementarity", "complementarity", "complementarities"),
    )
    return [_counted(counts[kind], one, many) for kind, one, many in names if counts[kind]]


def ga_triangle_figure() -> Path:
    """Map where claims sit, then name a representative result for each pair."""
    figure, (diagram, ledger) = plt.subplots(
        1, 2, figsize=(11.5, 4.9), gridspec_kw={"width_ratios": (0.79, 1)}
    )
    grouped = by_pair()
    entries_by_key = by_key()
    vertices = {"OPH": (0.50, 0.78), "FEP": (0.17, 0.24), "Synergetics": (0.83, 0.24)}
    examples = (
        (OPH_FEP, "OPH", "FEP", "spectral-gap", "Mean repair and intrinsic precision share a Laplacian."),
        (SYNERGETICS_OPH, "OPH", "Synergetics", "jitterbug-rigid-motion", "The Jitterbug reaches the regular carrier hull."),
        (FEP_SYNERGETICS, "FEP", "Synergetics", "cuboctahedral-distance-split", "First-shell pairs two steps apart can differ in variance."),
    )
    for axis in (diagram, ledger):
        axis.set(xlim=(0.0, 1.0), ylim=(0.0, 1.0))
        axis.set_axis_off()
    diagram.text(0.5, 0.99, "Scope of each check", ha="center", va="top",
                 fontsize=15, color=INK, fontweight="bold")
    edge_descriptions = []
    for index, (pair, first, second, key, result) in enumerate(examples):
        if entries_by_key[key].pair != pair:
            raise ValueError(f"{key} belongs to {entries_by_key[key].pair}, expected {pair}")
        entries = grouped[pair]
        phrases = _kind_phrases(Counter(entry.kind for entry in entries))
        start, end = np.array([vertices[first], vertices[second]])
        diagram.plot(
            (start[0], end[0]), (start[1], end[1]), color=BLUE, linewidth=2.0, zorder=1,
        )
        midpoint = (start + end) / 2
        diagram.text(
            *midpoint, str(len(entries)), ha="center", va="center", color=BLUE,
            fontsize=14, fontweight="bold", zorder=3,
            bbox={"facecolor": "white", "edgecolor": BLUE, "boxstyle": "round,pad=0.21"},
        )
        y = 0.97 - 0.245 * index
        ledger.text(0.02, y, f"{first} and {second}", va="top", color=INK,
                    fontsize=14, fontweight="bold")
        ledger.text(0.02, y - 0.09, result, va="top", color=BLUE, fontsize=12)
        ledger.text(0.02, y - 0.16, " · ".join(phrases), va="top", color=INK, fontsize=10.5)
        ledger.hlines(y - 0.21, 0.02, 0.98, color="#CBD5E1", linewidth=0.8)
        edge_descriptions.append(
            f"{first} and {second} have {', '.join(phrases)}. {result} Probe: {key}."
        )

    for name, (x, y) in vertices.items():
        diagram.scatter(x, y, s=60, color=INK, zorder=4)
        offset = 0.055 if name == "OPH" else -0.062
        diagram.text(x, y + offset, name, ha="center", va="center", fontsize=14,
                     fontweight="bold", color=INK)
    three = grouped[THREE_WAY]
    diagram.text(
        0.50, 0.42, f"{len(three)} three-way", ha="center", va="center",
        fontsize=13, color=INK, fontweight="bold",
    )
    diagram.text(
        0.5, 0.085, f"{len(CORRESPONDENCES)} executable claims", ha="center",
        va="center", color=INK, fontsize=12, fontweight="bold",
    )
    diagram.text(0.5, 0.015, "Specify the model  →  Derive  →  Check", ha="center",
                 va="center", color=BLUE, fontsize=11)
    ledger.text(0.02, 0.23, "Three-way claims", va="top", color=INK, fontsize=14, fontweight="bold")
    ledger.text(0.02, 0.14, "Equal variance can hide distinct repair rates on periodic FCC.",
                va="top", color=BLUE, fontsize=12)
    ledger.text(0.02, 0.07, "Transfer tests on periodic lattices, root shells and a torus",
                va="top", color=INK, fontsize=10.5)
    return _save(
        figure,
        "fig_ga_triangle.png",
        "A triangle connects Observer Patch Holography (OPH), the Free Energy Principle (FEP) "
        f"and Fuller's Synergetics. Its edge totals and centre account for all "
        f"{len(CORRESPONDENCES)} executable claims. " + " ".join(edge_descriptions)
        + f" The centre counts {len(three)} three-way claims. For instance, on a periodic FCC graph, "
        "equal variance can conceal different repair probabilities. The workflow specifies the model, "
        "derives a comparison and checks it against the probe.",
    )


def ga_carrier_figure() -> Path:
    """Show one Markov blanket and the carrier's invariant ledger."""
    carrier, view = _carrier_view()
    flat = view[:, :2]
    target = 0
    neighbours = invariants.markov_blanket(carrier, target)
    figure, (left, right) = plt.subplots(
        1, 2, figsize=(11.0, 5.2), gridspec_kw={"width_ratios": (0.92, 1)}
    )
    front_edges = set()
    for face in carrier.faces:
        triangle = view[list(face)]
        a, b, c = triangle
        normal = np.cross(b - a, c - a)
        if normal[2] * (normal @ triangle.mean(axis=0)) > 0:
            left.fill(*triangle[:, :2].T, color=SKY, alpha=0.18, linewidth=0, zorder=1)
            front_edges.update(tuple(sorted((index, face[i - 1]))) for i, index in enumerate(face))

    for edge in carrier.edges:
        foreground = edge in front_edges
        blanket_edge = target in edge
        left.plot(
            *flat[list(edge)].T, color=BLUE if blanket_edge else GREY,
            linewidth=2.3 if blanket_edge else 1.1, linestyle="-" if foreground else "--",
            alpha=1.0 if blanket_edge else 0.62, zorder=2.5 if blanket_edge else 2,
        )
    groups = (
        ((target,), "D", INK, INK, "Target port"),
        (neighbours, "o", "#D8EEFA", BLUE, "Markov blanket"),
        (tuple(port for port in range(carrier.port_count) if port not in (target, *neighbours)),
         "o", "white", GREY, "Other ports"),
    )
    handles = []
    for members, marker, face, edge, label in groups:
        handles.append(left.scatter(
            *flat[list(members)].T, s=290 if marker == "D" else 335,
            marker=marker, color=face, edgecolor=edge, linewidth=1.45, label=label, zorder=3,
        ))
    for port, (x, y) in enumerate(flat):
        left.text(
            x, y, str(port),
            ha="center",
            va="center",
            fontsize=11,
            fontweight="bold",
            color="white" if port == target else INK,
            zorder=4,
        )
    left.set_aspect("equal")
    left.set_axis_off()
    left.set_title(f"Carrier blanket of port {target}", fontsize=15, fontweight="bold")
    left.text(
        0.5, -0.03, f"{carrier.port_count} ports  ·  {carrier.seam_count} seams",
        transform=left.transAxes, ha="center", va="top", color=INK, fontsize=12,
    )
    left.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.09),
                ncol=3, fontsize=9.5, markerscale=0.55, handletextpad=0.4, columnspacing=1.25)
    left.margins(0.12)

    cohomology = sheaf.cohomology_dimensions(sheaf.carrier_sheaf(carrier))
    group = symmetry.carrier_automorphisms(carrier)
    rotations = symmetry.rotation_subgroup(carrier)
    orientation_reversing = symmetry.reflection_subgroup(carrier)
    marginal = invariants.marginal_variance(carrier)
    blanket = invariants.blanket_sizes(carrier)
    precision = invariants.pairwise_precision(carrier)
    tree_number = invariants.spanning_tree_number(carrier)

    blanket_values = sorted(set(blanket))
    blanket_text = (
        str(blanket_values[0]) if len(blanket_values) == 1 else ", ".join(map(str, blanket_values))
    )
    precision_values = sorted(set(precision.values()))
    precision_text = ", ".join(str(value) for value in precision_values)

    sections = (
        (0.98, "Graph and agreement sheaf", (
            ("Spanning trees", f"{tree_number:,}"),
            (r"Dimensions $(H^0, H^1)$", f"({cohomology['h0']}, {cohomology['h1']})"),
        )),
        (0.70, "Symmetry", (
            ("Automorphism group order", str(len(group))),
            ("Rotations / orientation-reversing", f"{len(rotations)} / {len(orientation_reversing)}"),
        )),
        (0.42, "Gaussian field", (
            (r"Centred intrinsic variance $(L^+)_{ii}$", str(marginal)),
            ("Proper Gaussian blanket size", blanket_text),
            ("Precision on a seam", precision_text.replace("-", "\u2212")),
            ("Ordered seam pairs", f"{len(precision)} = 2 × {carrier.seam_count}"),
        )),
    )
    right.set_axis_off()
    right.set_xlim(0.0, 1.0)
    right.set_ylim(0.0, 1.0)
    for top, heading, rows in sections:
        right.text(0.02, top, heading, va="top", fontsize=14, fontweight="bold", color=INK)
        right.hlines(top - 0.065, 0.02, 0.99, color="#CBD5E1", linewidth=0.8)
        for index, (label, value) in enumerate(rows):
            y = top - 0.12 - 0.08 * index
            right.text(0.02, y, label, va="center", fontsize=12, color=INK)
            right.text(0.99, y, value, va="center", ha="right", fontsize=12,
                       fontweight="bold", color=BLUE)
    return _save(
        figure,
        "fig_ga_carrier.png",
        f"The spectral projection of the icosahedral carrier shows {carrier.port_count} numbered "
        f"ports and {carrier.seam_count} seams. The target port {target} is a dark diamond, "
        f"and its neighbours {_enumerated(neighbours)} form the proper Gaussian's Markov blanket. "
        "Their seams to the target are highlighted. Solid seams mark the front and rear seams "
        f"are dashed. The carrier has {tree_number:,} spanning trees and the agreement "
        f"sheaf has dimensions {cohomology['h0']} and {cohomology['h1']} for H0 and H1. Its "
        f"{len(group)} automorphisms include {len(rotations)} rotations and "
        f"{len(orientation_reversing)} orientation-reversing symmetries. The centred intrinsic "
        f"port variance is {marginal}. The proper Gaussian blanket size is {blanket_text} and "
        f"the precision on each of {len(precision)} ordered seam pairs is {precision_text}.",
    )


def ga_descent_figure() -> Path:
    """Compare native-law outcomes and unit-covariance Gaussian divergence levels."""
    carrier = icosahedral_carrier()
    start, end = carrier.edges[0]
    loads = tuple(2 if port == start else 0 for port in range(carrier.port_count))
    example = repair.native_repair(loads, (start, end), coin=0)
    before = (example.before[start], example.before[end])
    after = (example.after[start], example.after[end])
    census = controls.gated_descent_census()
    figure, (left, right) = plt.subplots(
        1, 2, figsize=(11.5, 4.8), gridspec_kw={"width_ratios": (0.92, 1.08)}
    )
    labels = ["Strict descent", "Move at equal energy", "No change"]
    values = [census["strict"], census["zero"], census["noop"]]
    total = sum(values)
    bars = left.barh(labels, [value / total for value in values],
                     color=[GREEN, ORANGE, GREY], height=0.56)
    for bar, value in zip(bars, values):
        left.annotate(
            f"{value:,}  ({value / total:.1%})",
            (bar.get_width(), bar.get_y() + bar.get_height() / 2),
            textcoords="offset points",
            xytext=(5, 0),
            va="center",
            fontsize=11,
            fontweight="bold",
            color=INK,
        )
    left.invert_yaxis()
    left.set_title(f"{total:,} native repair attempts", fontsize=14, fontweight="bold", pad=12)
    left.set_xlim(0, 1.24)
    left.set_xticks(np.linspace(0, 1, 5), ("0%", "25%", "50%", "75%", "100%"))
    left.tick_params(axis="y", length=0, labelsize=11)
    left.grid(False, axis="y")
    left.set_xlabel("Share of enumerated attempts", fontsize=11)
    left.text(
        0.0, -0.24, f"Seam loads {before} → {after}\nLoad energy falls by {example.energy_decrease}.",
        transform=left.transAxes, ha="left", va="top", color=INK, fontsize=10.5,
    )

    evidence = by_key()["shell-vs-free-energy"].run().values
    levels_by_shell = evidence["free_energy_levels"]
    shared_levels = evidence["levels_shared_across_shells"]
    k_values = list(levels_by_shell)
    minima = [min(levels) for levels in levels_by_shell.values()]
    maxima = [max(levels) for levels in levels_by_shell.values()]
    right.vlines(k_values, minima, maxima, color="#CBD5E1", linewidth=4, zorder=1)
    for level, shells in shared_levels.items():
        right.hlines(
            level, min(shells), max(shells), color=VERMILLION, linestyle=":", linewidth=1.6,
            zorder=2,
        )
        right.annotate(
            f"{level:g}", (max(shells), level), xytext=(9, 0), textcoords="offset points",
            va="center", fontsize=10.5, color=INK,
        )
    ordinary = [
        (k, level) for k, levels in levels_by_shell.items()
        for level in levels if level not in shared_levels
    ]
    shared = [
        (k, level) for level, shells in shared_levels.items() for k in shells
    ]
    right.scatter(
        [k for k, _ in ordinary], [level for _, level in ordinary],
        color=BLUE, s=38, zorder=4, label="Unique to one shown shell",
    )
    right.scatter(
        [k for k, _ in shared], [level for _, level in shared],
        facecolor="white", edgecolor=VERMILLION, marker="D", s=45,
        linewidth=1.65, zorder=5, label="Shared by shown shells",
    )
    right.plot(
        k_values,
        maxima,
        color=ORANGE,
        linestyle="--",
        linewidth=1.8,
        label=r"Shell maximum $4k^2$",
        zorder=3,
    )
    right.set_xticks(k_values)
    right.set_xlim(min(k_values) - 0.35, max(k_values) + 0.68)
    right.set_ylim(0, max(maxima) * 1.08)
    right.set_xlabel("IVM coordination shell $k$", fontsize=11)
    right.set_ylabel("Gaussian KL (nats)", fontsize=11)
    right.grid(False, axis="x")
    right.set_title("Gaussian costs on lattice shells", fontsize=14, fontweight="bold", pad=12)
    right.legend(loc="upper left", fontsize=10.5)
    return _save(
        figure,
        "fig_ga_descent.png",
        f"Left: shares of {total:,} native repair attempts. " + " ".join(
            f"{label}: {value:,} attempts, {value / total:.1%}."
            for label, value in zip(labels, values)
        ) + f" In one strict descent, seam loads {before} become {after}, reducing load "
        f"energy by {example.energy_decrease}. "
        f"Right: distinct Gaussian KL values for each of the first {len(k_values)} IVM shells. "
        f"Vertical lines span each shell's range. Levels shared among the displayed shells are "
        f"{_enumerated(shared_levels)} "
        "nats, shown as open diamonds connected by dotted guides across shells. "
        "A dashed curve joins the shell maxima at four k squared.",
    )


def ga_synergetics_figure() -> Path:
    """Compare the IVM solids, the carrier icosahedron and the IVM Voronoi cell."""
    figure, (left, right) = plt.subplots(
        1, 2, figsize=(11.5, 4.85), gridspec_kw={"width_ratios": (1.26, 0.74)}
    )
    ivm_volumes = geometry.derived_volume_ratios()
    ico_volume = float(geometry.icosahedron_tetravolume_exact()["tetravolume"])
    _plot_tetravolumes(left, ivm_volumes, icosahedron=ico_volume)
    left.set_title("IVM solids and the carrier", fontsize=14, fontweight="bold", pad=12)

    cell = lattice.voronoi_cell()
    density = geometry.packing_density()
    cell_volume = ivm_volumes["rhombic_dodecahedron"]
    cell_labels = ("vertices", "faces", "edges")
    cell_values = [int(density[label]) for label in cell_labels]
    faces = cell["faces"]
    edges = sorted({
        tuple(sorted((vertex, face[index - 1])))
        for face in faces for index, vertex in enumerate(face)
    })
    vertices = np.asarray(cell["vertices"], dtype=float)
    view = _clearest_view(vertices - vertices.mean(axis=0), np.array(edges))
    view /= max(np.linalg.norm(view, axis=1))
    flat = view[:, :2]
    front_edges = set()
    front_faces = []
    for face in faces:
        polygon = view[list(face)]
        a, b, c = polygon[:3]
        normal = np.cross(b - a, c - a)
        if normal[2] * (normal @ polygon.mean(axis=0)) > 0:
            right.fill(*polygon[:, :2].T, color=SKY, alpha=0.20, linewidth=0, zorder=1)
            front_edges.update(tuple(sorted((vertex, face[i - 1]))) for i, vertex in enumerate(face))
            front_faces.append((abs(normal[2]), face))
    for edge in edges:
        foreground = edge in front_edges
        right.plot(
            *flat[list(edge)].T, color=BLUE if foreground else GREY,
            linewidth=1.6 if foreground else 1, linestyle="-" if foreground else "--",
            alpha=1.0 if foreground else 0.56, zorder=2 if foreground else 1.5,
        )
    front_vertices = sorted({vertex for edge in front_edges for vertex in edge})
    rear_vertices = sorted(set(range(len(flat))) - set(front_vertices))
    right.scatter(*flat[rear_vertices].T, s=24, facecolor="white", edgecolor=GREY, alpha=0.55, zorder=2)
    right.scatter(*flat[front_vertices].T, s=33, facecolor="white", edgecolor=BLUE, zorder=3)
    _, example_face = max(front_faces, key=lambda item: item[0])
    triangle = flat[list(example_face[:3])]
    right.fill(*triangle.T, color=ORANGE, alpha=0.45, linewidth=0, zorder=1.5)
    right.plot(*triangle[[0, 2]].T, color=VERMILLION, linewidth=1.5, zorder=2.5)
    for x, y in triangle:
        right.plot((0, x), (0, y), color=VERMILLION, linestyle=":", linewidth=1.3, zorder=2.5)
    right.scatter(0, 0, s=34, marker="D", facecolor=VERMILLION, edgecolor=INK, zorder=4)
    right.set(xlim=(-1.55, 1.55), ylim=(-1.66, 1.13), aspect="equal")
    right.set_anchor("N")
    right.set_axis_off()
    right.set_title("IVM Voronoi cell", fontsize=14, fontweight="bold", pad=12)
    right.text(
        0.0, -1.07, " · ".join(f"{value} {name}" for name, value in zip(cell_labels, cell_values)),
        ha="center", va="top", fontsize=11, color=INK,
    )
    right.text(
        0.0, -1.28, f"{cell_volume} tetravolumes (exact)", ha="center", va="top",
        fontsize=12, fontweight="bold", color=INK,
    )
    right.text(
        0.0, -1.49, rf"Packing density $\pi/(3\sqrt{{2}})\approx {float(density['density']):.6f}$",
        ha="center", va="top", fontsize=11, color=INK,
    )
    return _save(
        figure,
        "fig_ga_synergetics.png",
        "Left: integer tetravolumes on a logarithmic axis: "
        + ", ".join(f"{name.replace('_', ' ')} {value}" for name, value in ivm_volumes.items())
        + f". A diamond labels the carrier icosahedron at approximately {ico_volume:.4f}. "
        + "Right: a rhombic dodecahedron represents the IVM Voronoi cell, with rear edges dashed. "
        + "An amber triangle on one face connects to the diamond at the cell centre by dotted "
        + "segments to show one tetrahedron in the volume sum. "
        + "The cell has "
        + ", ".join(f"{value} {name}" for name, value in zip(cell_labels, cell_values))
        + f", exact tetravolume {cell_volume} and packing density {float(density['density']):.6f}, "
        "equal to pi divided by three square root of two.",
    )


def pair_triangles_figure() -> Path:
    """Give each scope the same comparison layout and highlight its edge."""
    figure, axes = plt.subplots(1, 3, figsize=(11.6, 5.1))
    pairs = (OPH_FEP, FEP_SYNERGETICS, SYNERGETICS_OPH)
    pair_edges = {
        OPH_FEP: ("OPH", "FEP"),
        FEP_SYNERGETICS: ("FEP", "Synergetics"),
        SYNERGETICS_OPH: ("Synergetics", "OPH"),
    }
    vertices = {
        "OPH": (0.50, 0.72),
        "FEP": (0.18, 0.48),
        "Synergetics": (0.82, 0.48),
    }
    status_names = {
        "exact": "exact",
        "exact-on-domain": "exact on domain",
        "obstructed": "obstructed",
        "not-reproduced": "not reproduced",
    }
    grouped = by_pair()
    descriptions = []
    for axis, pair in zip(axes, pairs):
        entries = grouped[pair]
        counts = Counter(entry.kind for entry in entries)
        statuses = Counter(entry.status for entry in entries)
        axis.set(xlim=(0.0, 1.0), ylim=(0.0, 1.0))
        axis.set_axis_off()
        first, second = pair_edges[pair]
        axis.text(0.5, 0.98, f"{first} and {second}", va="top", ha="center",
                  fontsize=13, color=INK, fontweight="bold")
        for left_name, right_name in pair_edges.values():
            highlighted = pair_edges[pair] == (left_name, right_name)
            axis.plot(
                [vertices[left_name][0], vertices[right_name][0]],
                [vertices[left_name][1], vertices[right_name][1]],
                color=BLUE if highlighted else "#B8C2CC",
                linewidth=3.6 if highlighted else 1.4,
                solid_capstyle="round",
                zorder=1,
            )
        for name, (x_value, y_value) in vertices.items():
            selected = name in (first, second)
            axis.scatter(x_value, y_value, s=52, color=INK if selected else "#B8C2CC", zorder=3)
            offset = 0.065 if name == "OPH" else -0.066
            axis.text(
                x_value, y_value + offset, name, ha="center", va="center", fontsize=11, color=INK,
                fontweight="bold" if selected else "normal",
            )
        axis.text(
            0.50,
            0.555,
            str(len(entries)),
            ha="center",
            va="center",
            fontsize=17,
            fontweight="bold",
            color=INK,
            zorder=2,
        )
        # Preserve stable positions across the cards so zero counts are visible.
        for index, kind in enumerate(KINDS):
            x = 0.19 + 0.31 * index
            axis.text(x, 0.29, str(counts[kind]), ha="center", va="center", fontsize=17,
                      color=KIND_COLOURS[kind], fontweight="bold")
            axis.text(x, 0.219, kind, ha="center", va="center", color=INK, fontsize=9.5)
        axis.hlines(0.177, 0.03, 0.97, color="#CBD5E1", linewidth=0.8)
        status_phrases = [f"{statuses[key]} {label}" for key, label in status_names.items()
                          if statuses[key]]
        axis.text(
            0.50,
            0.137,
            "  ·  ".join(status_phrases[:2])
            + ("\n" + "  ·  ".join(status_phrases[2:]) if len(status_phrases) > 2 else ""),
            ha="center",
            va="top",
            fontsize=10.5,
            color=INK,
        )
        kinds_text = ", ".join(f"{kind} {counts[kind]}" for kind in KINDS)
        descriptions.append(f"{first} and {second}: {len(entries)} claims, {kinds_text}, "
                            f"with status counts of {', '.join(status_phrases)}.")
    three = grouped[THREE_WAY]
    if three:
        figure.text(
            0.50,
            0.012,
            f"{len(three)} three-way claims: " + " and ".join(entry.key for entry in three),
            ha="center",
            va="bottom",
            fontsize=11,
            color=INK,
        )
    return _save(
        figure,
        "fig_pair_triangles.png",
        "Three triangles highlight the edge for each pair and display its claim total. Counts by "
        "kind and status are listed below each triangle, with zeros shown for claim kinds. "
        + " ".join(descriptions)
        + f" The {len(three)} three-way claims are {', '.join(entry.key for entry in three)}.",
    )


def scope_matrix_figure() -> Path:
    """Align the two distributions and their shared total to one scope column."""
    grouped = by_pair()
    scopes = ("OPH and FEP", "FEP and Synergetics", "Synergetics and OPH", "Three-way")
    statuses = (*STATUS_MARKERS, "not-reproduced")
    columns = (
        "Isomorphism", "Disagreement", "Complemen-\ntarity", "Exact", "Exact on\ndomain",
        "Obstructed", "Not\nreproduced",
    )
    positions = (0.0, 1.12, 2.24, 3.61, 4.73, 5.85, 6.97)
    total_column = 8.18
    rows = [
        (
            Counter(entry.kind for entry in grouped[pair]),
            Counter(entry.status for entry in grouped[pair]),
            len(grouped[pair]),
        )
        for pair in PAIRS
    ]
    figure, axis = plt.subplots(figsize=(11.65, 3.6))
    axis.set(xlim=(-3.0, 8.85), ylim=(-0.58, 5.05))
    axis.set_axis_off()
    axis.text(-2.95, 4.03, "Scope", ha="left", va="center", fontsize=12, fontweight="bold", color=INK)
    axis.text(positions[1], 4.74, "Claim kind", ha="center", va="center",
              fontsize=14, fontweight="bold", color=INK)
    axis.text(np.mean(positions[3:]), 4.74, "Evidence status", ha="center", va="center",
              fontsize=14, fontweight="bold", color=INK)
    axis.text(total_column, 4.03, "Total", ha="center", va="center",
              fontsize=12, fontweight="bold", color=INK)
    for x, name in zip(positions, columns):
        axis.text(x, 4.03, name, ha="center", va="center", fontsize=10.5, color=INK)
    axis.hlines(3.59, -2.95, 8.7, color="#CBD5E1", linewidth=0.9)
    colourmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "registry_counts", ("#EEF4F8", "#8AC1E2")
    )
    highest = max(1, *(count for kinds, evidence, _ in rows
                       for count in (*kinds.values(), *evidence.values())))
    descriptions = []
    for y, scope, (kinds, evidence, total) in zip(reversed(range(len(rows))), scopes, rows):
        axis.text(-2.95, y, scope, ha="left", va="center", color=INK,
                  fontsize=12, fontweight="bold")
        counts = [kinds[kind] for kind in KINDS] + [evidence[status] for status in statuses]
        for x, count in zip(positions, counts):
            axis.add_patch(
                plt.Rectangle((x - 0.52, y - 0.38), 1.04, 0.76,
                              color=colourmap(count / highest), linewidth=0)
            )
            axis.text(x, y, str(count), ha="center", va="center", fontsize=14, color=INK,
                      fontweight="bold" if count else "normal")
        axis.text(total_column, y, str(total), ha="center", va="center",
                  fontsize=15, fontweight="bold", color=BLUE)
        kinds_text = ", ".join(f"{kind} {kinds[kind]}" for kind in KINDS)
        evidence_text = ", ".join(f"{status} {evidence[status]}" for status in statuses)
        descriptions.append(
            f"{scope} has {total} claims, {kinds_text} by kind and {evidence_text} by status."
        )
    return _save(
        figure,
        "fig_scope_matrix.png",
        f"One matrix groups the {len(CORRESPONDENCES)} claims by kind and evidence status, "
        "with counts and totals aligned to one row per scope. Every cell shows its count. "
        + " ".join(descriptions),
    )


def _draw_wireframe(axis, vertices, edges, centre, radius, colour) -> None:
    """Draw the projected skeleton of a solid at a fixed visual radius."""
    points = np.asarray(vertices, dtype=float)
    projected = _clearest_view(points - points.mean(axis=0), np.array(edges))[:, :2]
    span = float(np.linalg.norm(projected, axis=1).max())
    placed = projected / span * radius + np.asarray(centre)
    for a, b in edges:
        axis.plot(placed[[a, b], 0], placed[[a, b], 1], color=colour, linewidth=1.3, zorder=2)
    axis.scatter(placed[:, 0], placed[:, 1], s=14, color=colour, zorder=3)


def network_figure() -> Path:
    """The carrier read as a network, with the geometry dual pair beside it.

    Four panels, each plotted only from executed package functions: the exact
    effective-resistance ladder and its lowest spectral contribution, the
    distance shells around a port, the subgraph census and the FCC dual pair.
    """
    carrier = icosahedral_carrier()
    adjacency = [[int(value) for value in row] for row in carrier.adjacency_matrix()]
    comparison = network.carrier_resistance_geometry()
    pair_shells = comparison["shells"]
    lowest = comparison["lowest_mode_variances"]
    kirchhoff = network.kirchhoff_index(adjacency)
    distribution = network.distance_distribution(adjacency)
    intersection = network.intersection_array(adjacency)
    census = network.subgraph_census(adjacency)
    face_count = len(carrier.faces)

    cell = geometry.voronoi_cell()
    volumes = geometry.derived_volume_ratios()
    cell_tetravolume = volumes["rhombic_dodecahedron"]
    cuboctahedron = geometry.cuboctahedron_signature()
    cuboctahedron_edge_count = len(geometry.cuboctahedron_edges())
    cuboctahedron_face_counts = geometry.cuboctahedron_faces()
    cuboctahedron_tetravolume = volumes["cuboctahedron"]
    volume_ratio = cuboctahedron_tetravolume / cell_tetravolume

    figure, axes = plt.subplots(2, 2, figsize=(11.1, 7.4))
    ladder_axis, shells_axis, census_axis, dual_axis = axes.ravel()

    # --- panel 1: the resistance ladder -----------------------------------
    positions = np.arange(len(pair_shells))
    resistances = [float(entry[1]) for entry in pair_shells]
    ladder_axis.bar(positions - 0.16, resistances, width=0.28, color=BLUE,
                    label="Full variance (resistance)", zorder=2)
    ladder_axis.bar(positions + 0.16, [float(value) for value in lowest], width=0.28,
                    facecolor="white", edgecolor=VERMILLION, hatch="///",
                    label=f"Lowest modes (chord² ÷ {comparison['chord_to_lowest_variance_factor']})",
                    zorder=2)
    for position, (_, resistance, _, _) in zip(positions, pair_shells):
        ladder_axis.text(
            position - 0.16,
            float(resistance) + 0.011,
            str(resistance),
            ha="center",
            va="bottom",
            fontsize=11,
            color=INK,
        )
    ladder_axis.set_xticks(
        positions, [f"$d={step}$\n{pairs} pairs" for step, _, _, pairs in pair_shells],
        fontsize=11,
    )
    ladder_axis.set_ylabel("Gaussian contrast variance", fontsize=11)
    ladder_axis.set_ylim(0.0, max(resistances) * 1.55)
    ladder_axis.set_title(f"A  Contrast variance  ·  Kirchhoff index {kirchhoff}",
                          loc="left", fontsize=12)
    ladder_axis.legend(loc="upper left", fontsize=10)
    ladder_axis.grid(False, axis="x")

    # --- panel 2: the shells around port 0 --------------------------------
    projected_carrier, flat = _projected_carrier()
    distances = network.graph_distances(adjacency)[0]
    shell_colours = {
        0: VERMILLION,
        1: BLUE,
        2: ORANGE,
        3: GREEN,
    }
    shell_markers = ("D", "o", "s", "P")
    for port in range(projected_carrier.port_count):
        for other in projected_carrier.neighbours(port):
            if other > port:
                shells_axis.plot(
                    [flat[port, 0], flat[other, 0]],
                    [flat[port, 1], flat[other, 1]],
                    color="#9AA9B8",
                    linewidth=0.9,
                    zorder=1,
                )
    for step in range(len(distribution)):
        members = [port for port in range(projected_carrier.port_count) if distances[port] == step]
        shells_axis.scatter(
            flat[members, 0],
            flat[members, 1],
            s=188,
            color=shell_colours.get(step, GREY),
            marker=shell_markers[step],
            edgecolor=INK,
            zorder=3,
            label=f"{step}: {distribution[step]}",
        )
    for port in range(projected_carrier.port_count):
        shells_axis.text(
            flat[port, 0],
            flat[port, 1],
            str(port),
            ha="center",
            va="center",
            fontsize=9.5,
            fontweight="bold",
            color="white" if distances[port] in (0, 1) else BLACK,
            zorder=4,
        )
    shells_axis.margins(x=0.09, y=0.09)
    shells_axis.set_aspect("equal")
    shells_axis.set_axis_off()
    shells_axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 0.02),
        ncol=4,
        fontsize=10,
        title=f"Distance: ports  ·  $b={intersection[0]}$, $c={intersection[1]}$",
        title_fontsize=10,
        markerscale=0.6,
        columnspacing=0.75,
    )
    shells_axis.set_title("B  Graph-distance shells from port 0", loc="left", fontsize=12)

    # --- panel 3: the subgraph census -------------------------------------
    census_names = {
        "triangles": "triangles",
        "pentagons": "induced $5$-cycles",
        "perfect_matchings": "perfect matchings",
        "hamiltonian_cycles": "Hamiltonian cycles",
    }
    census_labels = [census_names.get(key, key.replace("_", " ")) for key in census]
    census_positions = list(range(len(census)))
    census_values = list(census.values())
    census_axis.scatter(census_values, census_positions, s=65, color=BLUE, zorder=3,
                        edgecolor="white", linewidth=0.7)
    for position, count in zip(census_positions, census_values):
        census_axis.annotate(
            f"{count:,}",
            (count, position),
            textcoords="offset points",
            xytext=(5, 0),
            ha="left",
            va="center",
            fontsize=11,
            color=INK,
        )
    census_axis.set_xscale("log")
    census_axis.axvline(face_count, color=VERMILLION, linestyle="--", linewidth=1.3)
    census_axis.text(
        face_count * 1.15,
        3.1,
        f"{face_count} faces",
        fontsize=10.5,
        color=VERMILLION,
        va="center",
    )
    census_axis.set_yticks(census_positions)
    census_axis.set_yticklabels(census_labels, fontsize=10.5)
    census_axis.set_ylim(len(census) - 0.55, -0.45)
    census_axis.set_xlim(8.0, max(census_values) * 4.0)
    census_axis.set_xlabel("Count (log scale)", fontsize=11)
    census_axis.set_title("C  Subgraph census", loc="left", fontsize=12)
    census_axis.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())

    # --- panel 4: the dual pair -------------------------------------------
    dual_axis.set_axis_off()
    dual_axis.set_title("D  FCC Voronoi cell and its dual", loc="left", fontsize=12)

    cell_edges = set()
    for face in cell["faces"]:
        for index, vertex in enumerate(face):
            following = face[(index + 1) % len(face)]
            cell_edges.add((min(vertex, following), max(vertex, following)))
    cuboctahedron_vertices = [geometry.to_xyz(site) for site in geometry.cuboctahedron_neighbors()]
    for centre, vertices, edges, colour, caption in (
        (-1.12, cell["vertices"], sorted(cell_edges), GREEN,
         f"Voronoi cell\n{cell['vertex_count']} vertices · {cell['face_count']} rhombi\n"
         f"{cell['edge_count']} edges · {cell_tetravolume} tetravolumes"),
        (1.12, cuboctahedron_vertices, geometry.cuboctahedron_edges(), BLUE,
         f"Cuboctahedron\n{cuboctahedron['vertices']} vertices\n"
         f"{len(cuboctahedron_face_counts['triangles'])} triangles + "
         f"{len(cuboctahedron_face_counts['squares'])} squares\n"
         f"{cuboctahedron_edge_count} edges · {cuboctahedron_tetravolume} tetravolumes"),
    ):
        _draw_wireframe(dual_axis, vertices, edges, (centre, 0.52), 0.69, colour)
        dual_axis.text(centre, -0.26, caption, ha="center", va="top", fontsize=10.3, color=INK)
    dual_axis.set_xlim(-2.4, 2.4)
    dual_axis.set_ylim(-1.24, 1.25)
    dual_axis.set_aspect("equal", adjustable="box")
    dual_axis.text(
        0.0,
        -1.11,
        f"Volume ratio (dual / Voronoi) = {volume_ratio}",
        ha="center",
        fontsize=10.5,
        color=INK,
    )

    return _save(
        figure,
        "fig_network.png",
        f"Four panels. A compares the full Gaussian contrast variance at the carrier's graph "
        f"distances, {_enumerated(str(entry[1]) for entry in pair_shells)}, with smaller hatched "
        f"bars for the lowest-mode contribution, approximately "
        f"{_enumerated(f'{float(value):.3f}' for value in lowest)}. The axis gives pair counts "
        "and the title gives "
        f"the Kirchhoff index {kirchhoff}. B labels the numbered carrier ports by graph distance "
        "from port zero using different marker shapes and colours. Its legend gives each shell "
        "count and the intersection array. C places the counts of triangles, induced five-cycles, "
        "matchings and Hamiltonian cycles on a log axis with a reference line at the face count. "
        "D projects the FCC Voronoi cell and its cuboctahedral dual and labels their vertices, "
        f"faces, edges and exact tetravolumes. The dual to Voronoi volume ratio is {volume_ratio}.",
    )


def vertex_shells_figure() -> Path:
    """The carrier's shells around a port against the IVM's lattice shells."""
    adjacency = [[int(value) for value in row] for row in icosahedral_carrier().adjacency_matrix()]
    carrier_shells = network.distance_distribution(adjacency)
    radii = (1, 2, 3)
    carrier_counts = carrier_shells[1:4]
    lattice_counts = tuple(lattice.shell_count(k) for k in radii)
    outward, homing = network.intersection_array(adjacency)
    figure, axis = plt.subplots(figsize=(8.8, 3.7))
    axis.hlines(radii, carrier_counts, lattice_counts, color=GREY, alpha=0.6,
                linewidth=1.5, zorder=2)
    for counts, marker, colour, alignment, offset, label in (
        (carrier_counts, "o", BLUE, "right", -8, "Carrier ports at graph distance $k$"),
        (lattice_counts, "s", GREEN, "left", 8, "IVM lattice sites on shell $k$"),
    ):
        axis.scatter(counts, radii, marker=marker, color=colour, s=85,
                     edgecolor="white", linewidth=0.7, label=label, zorder=3)
        for k, count in zip(radii, counts):
            axis.annotate(str(count), (count, k), textcoords="offset points",
                          xytext=(offset, 0), ha=alignment, va="center",
                          fontsize=12, color=INK)
    axis.set_xscale("log")
    axis.set_xticks([1, 2, 5, 10, 20, 50, 100], ["1", "2", "5", "10", "20", "50", "100"])
    axis.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    axis.set_xlim(0.7, max(lattice_counts) * 1.7)
    axis.set_yticks(radii, [f"$k={k}$" for k in radii], fontsize=11)
    axis.set_ylim(max(radii) + 0.51, min(radii) - 0.51)
    axis.set_xlabel("Number of ports or sites (log scale)")
    axis.grid(False, axis="y")
    axis.set_title(f"Carrier intersection array  ·  $b={outward}$, $c={homing}$",
                   fontsize=12, loc="left", pad=33)
    axis.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=10,
                borderaxespad=0.2)
    return _save(
        figure,
        "fig_vertex_shells.png",
        "Three rows match the graph-distance index and the lattice-shell index on a logarithmic "
        "count axis. Connected circles and squares show the two counts within each row. The "
        f"carrier's ports number {carrier_counts[0]}, {carrier_counts[1]} and {carrier_counts[2]} at "
        "graph distances one, two and three. The IVM lattice's shells hold "
        f"{lattice_counts[0]}, {lattice_counts[1]} and {lattice_counts[2]} sites at the same radii, "
        "and all three pairs differ. The "
        f"carrier is distance-regular with intersection array b = {outward}, c = {homing}, and "
        f"diameter {len(carrier_shells) - 1}.",
    )


def conjugation_figure() -> Path:
    """What Galois conjugation fixes (the counts) and what it moves (the volumes)."""
    values = as_dict()
    fixed = (
        ("spanning trees", int(values["carrier.tree_number"])),
        ("automorphisms", int(values["symmetry.group_order"])),
        ("sheaf obstructions", int(values["sheaf.h1"])),
        ("triangles", int(values["carrier.triangles"])),
    )
    exact = geometry.icosahedron_tetravolume_exact()
    volume = float(exact["tetravolume"])
    conjugate = float(exact["conjugate_tetravolume"])
    norm = exact["norm"]
    figure, (left, right) = plt.subplots(
        1, 2, figsize=(10.4, 3.95), gridspec_kw={"width_ratios": (1.2, 1.0)}
    )
    names = [name for name, _ in fixed]
    counts = [value for _, value in fixed]
    positions = np.arange(len(names))
    left.scatter(counts, positions, color=BLUE, s=85, edgecolor="white",
                 linewidth=0.7, zorder=3)
    left.set_xscale("log")
    left.set_xlim(1.0, max(counts) * 8.0)
    left.set_yticks(positions, names, fontsize=11)
    left.set_ylim(len(names) - 0.5, -0.5)
    left.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    for row, count in zip(positions, counts):
        left.annotate(
            f"{count:,}",
            (count, row),
            textcoords="offset points",
            xytext=(5, 0),
            ha="left",
            va="center",
            fontsize=10.5,
        )
    left.set_xlabel("Integer count (log scale)")
    left.set_title(r"Integer counts fixed by $\sqrt{5}\mapsto-\sqrt{5}$", loc="left")
    bars = right.barh(
        ("icosahedron", "conjugate"),
        (volume, conjugate),
        color=(BLUE, VERMILLION),
        height=0.32,
        zorder=2,
    )
    right.set_xlim(0.0, max(volume, conjugate) * 1.22)
    for bar, value in zip(bars, (volume, conjugate)):
        right.annotate(
            f"≈ {value:.5g}",
            (value, bar.get_y() + bar.get_height() / 2),
            textcoords="offset points",
            xytext=(5, 0),
            ha="left",
            va="center",
            fontsize=10.5,
        )
    right.set_ylim(1.57, -0.53)
    right.set_xlabel("Tetravolume")
    right.set_title("Tetravolume changes", loc="left")
    right.grid(False, axis="y")
    right.text(
        0.38,
        0.51,
        f"Exact norm  ·  $T\\,T'={norm}$",
        transform=right.transAxes,
        ha="center",
        va="center",
        fontsize=11,
        color=INK,
        bbox={"facecolor": "#F1F5F9", "edgecolor": SKY, "boxstyle": "round,pad=0.35"},
    )
    return _save(
        figure,
        "fig_conjugation.png",
        "Two panels. Left: points on a logarithmic axis for the counting invariants that "
        f"conjugating the square root of five fixes: the spanning-tree count {counts[0]}, the "
        f"automorphism group order {counts[1]}, the {counts[2]} sheaf obstructions and the "
        f"{counts[3]} triangles. Right: horizontal bars on a linear axis give the icosahedron's "
        f"tetravolume approximately {volume:.5g} and its conjugate approximately "
        f"{conjugate:.5g}. Their product is annotated as the exact norm {norm}.",
    )


def families_figure() -> Path:
    """The solids the phrase twelve-around-one can name, separated by their counts."""
    carrier = icosahedral_carrier()
    ico_adjacency = [[int(value) for value in row] for row in carrier.adjacency_matrix()]
    cubo_adjacency = lattice.cuboctahedron_adjacency()
    ico_degree = int(sum(ico_adjacency[0]))
    cubo_degree = int(sum(cubo_adjacency[0]))
    ico_edges = len(carrier.edges)
    cubo_edges = len(lattice.cuboctahedron_edges())
    cubo_faces = lattice.cuboctahedron_faces()
    cubo_triangles = len(cubo_faces["triangles"])
    cubo_squares = len(cubo_faces["squares"])
    ico_triangles = len(carrier.faces)
    ico_trees = invariants.spanning_tree_number(carrier)
    cubo_trees = invariants.spanning_tree_number_of_adjacency(cubo_adjacency)
    ico_volume = float(geometry.icosahedron_tetravolume_exact()["tetravolume"])
    cubo_volume = lattice.cuboctahedron_volume()
    cubo_root = math.isqrt(cubo_trees)
    figure, (left, right) = plt.subplots(
        1, 2, figsize=(10.6, 4.2), gridspec_kw={"width_ratios": (1.15, 1.0)}
    )
    names = ("Carrier", "Cuboctahedron")
    colours = (VERMILLION, BLUE)
    _, carrier_view = _carrier_view()
    cubo_vertices = [geometry.to_xyz(site) for site in lattice.cuboctahedron_neighbors()]
    columns = (2.95, 4.70)
    for x, vertices, edges, name, colour in zip(
        columns, (carrier_view, cubo_vertices), (carrier.edges, lattice.cuboctahedron_edges()),
        names, colours,
    ):
        _draw_wireframe(left, vertices, edges, (x, 3.48), 0.53, colour)
        left.text(x, 2.72, name, ha="center", color=INK, fontsize=11, fontweight="bold")
    rows = (
        ("Neighbours", ico_degree, cubo_degree),
        ("Edges", ico_edges, cubo_edges),
        ("Triangles", ico_triangles, cubo_triangles),
        ("Squares", 0, cubo_squares),
        ("Tetravolumes", f"≈ {ico_volume:.5g}", cubo_volume),
    )
    for y, (field, ico, cubo) in zip(np.linspace(2.21, 0.27, len(rows)), rows):
        left.hlines(y + 0.23, 0.08, 5.55, color=GREY, alpha=0.30, linewidth=0.7)
        left.text(0.10, y, field, color=INK, fontsize=11, va="center")
        for x, value in zip(columns, (ico, cubo)):
            left.text(x, y, str(value), color=INK, fontsize=11, va="center", ha="center")
    left.set(xlim=(0, 5.6), ylim=(0, 4.15), aspect="equal")
    left.set_axis_off()
    left.set_title("Incidence and volume", loc="left", fontsize=12)

    right.set_xscale("log")
    right.set_xlim(cubo_trees * 0.85, ico_trees * 1.8)
    right.set_ylim(1.55, -0.62)
    right.set_yticks((0, 1), names, fontsize=10.5)
    right.set_xticks((3e5, 1e6, 3e6), ("300k", "1M", "3M"))
    right.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    for value, label, row, colour, marker in (
        (ico_trees, f"{ico_trees:,}", 0, colours[0], "o"),
        (cubo_trees, f"{cubo_trees:,} = ${cubo_root}^2$", 1, colours[1], "s"),
    ):
        right.scatter(value, row, color=colour, marker=marker, s=94,
                      linewidth=0.75, edgecolor="white", zorder=3)
        right.annotate(label, (value, row), textcoords="offset points",
                       xytext=(0, 10), ha="right" if row == 0 else "left",
                       va="bottom", fontsize=11, color=INK)
    tree_ratio = Fraction(ico_trees, cubo_trees)
    right.text(0.51, 0.52, f"Carrier / cuboctahedron = {tree_ratio}",
               transform=right.transAxes, ha="center", va="center", fontsize=11,
               color=INK, bbox={"facecolor": "white", "edgecolor": "none", "pad": 3})
    right.set_title("Spanning trees of each graph", loc="left", fontsize=12)
    right.set_xlabel("Exact count (log scale)")
    return _save(
        figure,
        "fig_families.png",
        "Two panels. Left: wireframes and aligned counts compare the cuboctahedron with the carrier's "
        f"icosahedron: degree {cubo_degree} against {ico_degree}, edges {cubo_edges} against "
        f"{ico_edges}, triangles {cubo_triangles} against {ico_triangles}, and squares "
        f"{cubo_squares} against zero. The tetravolumes are {cubo_volume} and approximately "
        f"{ico_volume:.5g}. Right: the spanning trees the same "
        f"determinant counts, {ico_trees} for the carrier against {cubo_trees} for the "
        f"cuboctahedron, which is {cubo_root} squared, on a logarithmic axis. Their ratio "
        f"from carrier to cuboctahedron is {tree_ratio}.",
    )


def sandwich_figure() -> Path:
    """The surprise sandwich and the functionals that share the name free energy."""
    dictionary = by_key()["free-energy-dictionary"].run().values
    collision = by_key()["free-energy-name-collision"].run().values
    lower = float(dictionary["sandwich_lower"])
    upper = float(dictionary["sandwich_upper"])
    bayesian = float(dictionary["bayesian_surprise"])
    interface = float(collision["interface_energy_minimised_over_interior"])
    surprisal = float(collision["normalised_record_surprisal"])
    raw_potential = float(collision["unnormalised_record_potential"])
    minus = "−"

    def signed(value: float) -> str:
        return f"{value:.4g}".replace("-", minus)

    figure, (left, right) = plt.subplots(1, 2, figsize=(10.6, 3.65))
    span = upper - lower
    left.hlines(0.50, lower, upper, color=SKY, linewidth=8, zorder=1)
    left.vlines([lower, upper], 0.42, 0.58, color=BLUE, linewidth=2)
    left.scatter(bayesian, 0.50, color=VERMILLION, marker="D", s=90, edgecolor="white",
                 linewidth=1.1, zorder=3)
    left.text(
        bayesian,
        0.36,
        f"Record-change KL\n{signed(bayesian)}",
        ha="center",
        va="top",
        fontsize=10.5,
        color=VERMILLION,
    )
    for value, name in ((lower, "Lower bound"), (upper, "Upper bound")):
        left.text(
            value,
            0.64,
            f"{name}\n{signed(value)}",
            ha="center",
            va="bottom",
            fontsize=10.5,
            color=BLUE,
        )
    left.set_xlim(lower - 0.25 * span, upper + 0.25 * span)
    left.set_ylim(0.0, 1.0)
    left.set_yticks([])
    left.spines["left"].set_visible(False)
    left.grid(False)
    left.set_xlabel("Information (nats)")
    left.set_title("Bounds on record-change surprise", loc="left")

    entries = (
        ("normalised\nrecord surprisal", surprisal, GREEN),
        ("unnormalised\nrecord potential", raw_potential, BLUE),
        ("minimum\ninterface energy", interface, VERMILLION),
    )
    values = [value for _, value, _ in entries]
    right.axvline(0.0, color=BLACK, linewidth=0.8)
    for index, (_, value, color) in enumerate(entries):
        right.hlines(index, 0, value, color=color, linewidth=2)
        right.scatter(value, index, color=color, marker=("s", "D", "o")[index], s=66,
                      edgecolor="white", linewidth=0.8, zorder=3)
        right.annotate(
            signed(value),
            (value, index),
            ha="left" if value >= 0 else "right",
            va="center",
            fontsize=10,
            xytext=(7 if value >= 0 else -7, 0),
            textcoords="offset points",
        )
    right.set_xlim(min(values) * 1.51, max(values) * 1.47)
    right.set_yticks(range(len(entries)), [name.capitalize() for name, _, _ in entries])
    right.set_ylim(len(entries) - 0.45, -0.5)
    right.grid(False, axis="y")
    right.set_xlabel("Value in the carrier model")
    right.set_title("Three quantities at fixed records", loc="left")
    return _save(
        figure,
        "fig_sandwich.png",
        "Two panels read from the executed probes. Left: the sandwich interval with "
        f"its lower bound {lower:.7g} and upper bound {upper:.7g}, and the record-change KL "
        f"{bayesian:.7g} marked inside it. Right: three values at fixed records on the carrier: "
        f"minimum interface energy {interface:.7g}, unnormalised record potential "
        f"{raw_potential:.7g}, and normalised record surprisal {surprisal:.7g}.",
    )


def schedule_figure() -> Path:
    """The centered contrast between two schedules, one exact fraction per step."""
    steps_list = list(range(1, 13))
    exact_values = [controls.affine_path_contrast_exact(steps=step) for step in steps_list]
    contrasts = [float(value) for value in exact_values]
    peak = int(np.argmax(contrasts))
    figure, axis = plt.subplots(figsize=(9.4, 3.8))
    axis.plot(steps_list, contrasts, color=BLUE, marker="o", linewidth=1.65,
              markeredgecolor="white", markeredgewidth=0.8, markersize=5.5,
              linestyle=":", zorder=2)
    axis.scatter([steps_list[peak], steps_list[-1]], [contrasts[peak], contrasts[-1]],
                 marker="D", color=VERMILLION, s=73, edgecolor="white", linewidth=0.9, zorder=3)
    axis.annotate(f"Peak  {contrasts[peak]:.3g}", (steps_list[peak], contrasts[peak]),
                  textcoords="offset points", xytext=(12, 5), fontsize=10.5,
                  color=VERMILLION)
    axis.annotate(
        f"Step {steps_list[-1]} ≈ {contrasts[-1]:.5f}",
        (steps_list[-1], contrasts[-1]),
        ha="right",
        va="bottom",
        fontsize=10.5,
        color=VERMILLION,
        xytext=(steps_list[-1] - 0.25, 0.13),
        arrowprops={"arrowstyle": "-", "color": VERMILLION, "lw": 1.1},
    )
    axis.set_xticks(steps_list)
    axis.set_yscale("log")
    axis.set_yticks([0.02, 0.05, 0.10, 0.20, 0.50], ["0.02", "0.05", "0.10", "0.20", "0.50"])
    axis.minorticks_off()
    axis.set_ylim(0.015, max(contrasts) * 1.48)
    axis.grid(False, axis="x")
    axis.set_xlabel("Affine iteration step")
    axis.set_ylabel("Absolute contrast difference\n(log scale)")
    axis.set_title("Full and half-damped schedules share a fixed point", loc="left")
    return _save(
        figure,
        "fig_schedule.png",
        "The absolute difference between the full and the half-damped schedules' centered "
        "contrasts on a log scale, one point per step count from one to twelve. The schedules "
        f"share each fixed-point answer. The contrast difference peaks at step {steps_list[peak]} "
        f"and then falls to {contrasts[-1]:.7g}, exactly {exact_values[-1]} after twelve steps.",
    )


def _convex_hull(points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Andrew's monotone chain; the hull of a finite planar point set, counter-clockwise."""
    unique = sorted(set(points))
    if len(unique) < 3:
        return unique

    def cross(origin, a, b):
        return (a[0] - origin[0]) * (b[1] - origin[1]) - (a[1] - origin[1]) * (b[0] - origin[0])

    lower: list[tuple[float, float]] = []
    for point in unique:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[float, float]] = []
    for point in reversed(unique):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def shell_gauge_figure() -> Path:
    """Compare polyhedral coordination shells with spherical Gaussian KL levels."""
    bound = as_dict()
    window = 6
    shared_levels = tuple(int(level) for level in bound["shell.shared_free_energy_levels"])
    shell_range = range(1, int(bound["shell.gauge_checked_through"]) + 1)

    # --- left: the plane u3 = 0 of D3 coloured by the closed-form gauge -----
    slice_points: dict[tuple[int, int], int] = {}
    for site in geometry.sites_through_shell(window):
        halved = [int(value) // 2 for value in geometry.to_xyz(site)]
        if halved[2] == 0 and abs(halved[0]) <= window and abs(halved[1]) <= window:
            slice_points[(halved[0], halved[1])] = geometry.cuboctahedral_gauge(site)
    expected = sum(
        1
        for u1 in range(-window, window + 1)
        for u2 in range(-window, window + 1)
        if (u1 + u2) % 2 == 0
    )
    if len(slice_points) != expected:
        raise RuntimeError(f"expected {expected} lattice points in the slice, got {len(slice_points)}")
    xs = np.array([point[0] for point in slice_points])
    ys = np.array([point[1] for point in slice_points])
    gauges = np.array(list(slice_points.values()))

    figure, (left, right) = plt.subplots(
        1, 2, figsize=(11.0, 5.05), gridspec_kw={"width_ratios": (1.0, 1.18)}
    )
    cmap = plt.get_cmap("viridis", window + 1)
    scatter = left.scatter(
        xs, ys, c=gauges, cmap=cmap, vmin=-0.5, vmax=window + 0.5, s=46, edgecolor="white",
        linewidth=0.6, zorder=4,
    )
    colourbar = figure.colorbar(scatter, ax=left, ticks=range(window + 1), fraction=0.046, pad=0.02)
    colourbar.set_label("Shell index $g(u)$")
    for k in range(1, window + 1):
        ring = [point for point, value in slice_points.items() if value == k]
        hull = _convex_hull([(float(a), float(b)) for a, b in ring])
        hull.append(hull[0])
        left.plot(
            [p[0] for p in hull], [p[1] for p in hull], color=GREY, linewidth=0.9, zorder=2,
            label="shell boundary $g=k$" if k == 1 else None,
        )
    demo_k = math.isqrt(min(shared_levels) // 4)
    theta = np.linspace(0.0, 2.0 * math.pi, 400)
    circles = (
        (2 * demo_k**2, VERMILLION, "-", f"$F={2 * demo_k**2}$"),
        (4 * demo_k**2, VERMILLION, ":", f"$F={4 * demo_k**2}$"),
        (max(shared_levels), BLUE, "--", f"Shared $F={max(shared_levels)}$"),
    )
    for level, colour, style, label in circles:
        radius = math.sqrt(level / 2.0)  # F = 2 |u|^2 on the halved lattice
        left.plot(radius * np.cos(theta), radius * np.sin(theta), color=colour, linestyle=style,
                  linewidth=2.1, zorder=3, label=label)
    left.set_aspect("equal")
    left.set_xlim(-window - 0.8, window + 0.8)
    left.set_ylim(-window - 0.8, window + 0.8)
    left.set_xticks(range(-window, window + 1, 2))
    left.set_yticks(range(-window, window + 1, 2))
    left.set_xlabel("$u_1$")
    left.set_ylabel("$u_2$")
    left.set_title("A  Shells and KL contours in the $u_3=0$ slice", loc="left")
    left.grid(False)
    left.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, fontsize=10,
                handlelength=1.65, columnspacing=1.0)

    # --- right: the range of F on each shell --------------------------------
    shared: dict[Fraction, set[int]] = {}
    per_shell: dict[int, list[Fraction]] = {}
    vertex_counts = []
    for k in shell_range:
        energies = [Fraction(geometry.squared_norm(site), 2) for site in geometry.generate_shell(k)]
        levels = sorted(set(energies))
        per_shell[k] = levels
        vertex_counts.append(sum(1 for value in energies if value == 4 * k * k))
        for level in levels:
            shared.setdefault(level, set()).add(k)
    ks = np.array(list(shell_range))
    lows = np.array([float(per_shell[k][0]) for k in shell_range])
    highs = np.array([float(per_shell[k][-1]) for k in shell_range])
    right.bar(ks, highs - lows, bottom=lows, width=0.46, color=SKY, alpha=0.30, zorder=1,
              label="Min–max span")
    single_x, single_y, shared_x, shared_y = [], [], [], []
    for k in shell_range:
        for level in per_shell[k]:
            if len(shared[level]) > 1:
                shared_x.append(k)
                shared_y.append(float(level))
            else:
                single_x.append(k)
                single_y.append(float(level))
    right.scatter(single_x, single_y, s=6, color=BLUE, zorder=3, label="Unique to one shown shell")
    right.plot(shared_x, shared_y, linestyle="none", marker="_", markersize=8,
               markeredgewidth=1.15, color=VERMILLION, zorder=4, label="Shared by shown shells")
    curve = np.linspace(ks.min() - 0.3, ks.max() + 0.3, 100)
    right.plot(
        curve, 2.0 * curve**2, color=VERMILLION, linewidth=1.6, zorder=2,
        label="Lower bound $2k^2$",
    )
    right.plot(curve, 4.0 * curve**2, color=BLACK, linewidth=1.6, linestyle="--", zorder=2,
               label="Maximum $4k^2$")
    witness = max(shared_levels)
    witness_shells = sorted(shared[Fraction(witness)])
    right.hlines(witness, witness_shells[0], witness_shells[-1], color=BLUE,
                 linestyle="--", linewidth=1.65, zorder=3)
    right.annotate(f"$F={witness}$", (witness_shells[-1], witness), xytext=(9, -13),
                   textcoords="offset points", color=BLUE, fontsize=10)
    right.set_xticks(list(shell_range))
    right.set_xlim(0.45, shell_range[-1] + 0.60)
    right.set_xlabel("Coordination shell $k$")
    right.set_ylabel("Gaussian KL (nats)")
    right.set_title("B  Divergence at each coordination shell", loc="left")
    right.legend(loc="upper left", fontsize=9, ncol=2, handlelength=1.6, columnspacing=0.9)
    right.grid(False, axis="x")
    if len(set(vertex_counts)) != 1:
        raise RuntimeError(f"the maximum 4k^2 is not attained at a constant count: {vertex_counts}")
    return _save(
        figure,
        "fig_shell_gauge.png",
        "Two panels compare coordination shells with unit-covariance Gaussian divergence. "
        "Left: the lattice points of "
        f"the face-centred cubic lattice D3 in the plane u3 equal to zero with both coordinates "
        f"between minus {window} and {window}, coloured by their closed-form coordination "
        "shell; the shell boundaries are drawn as nested squares in this slice, and three "
        f"circles are overlaid for the divergence levels {2 * demo_k**2:g}, "
        f"{4 * demo_k**2:g} and {max(shared_levels):g}, so a shell is a polyhedral surface and a "
        "divergence level is a sphere that crosses several shells. Right: for shells one to "
        "twelve, a vertical bar from the smallest to the largest divergence on the shell, a mark "
        "for every distinct level with horizontal ticks marking levels shared among the "
        "displayed shells, and the "
        f"curves two k squared and four k squared; the maximum four k squared is attained at "
        f"exactly {vertex_counts[0]} sites on every shell, the twelve vertices. An upper bound "
        f"on divergence bounds the shell index. A horizontal dashed line links the shells "
        f"{_enumerated(witness_shells)} at the same divergence {witness} shown on the left.",
    )


def extended_claims_figure() -> Path:
    """Show the checked repair, record, dual, and lattice comparisons."""
    figure, axes = plt.subplots(2, 2, figsize=(11.6, 7.5))
    repair_axis, records_axis, dual_axis, shell_axis = axes.ravel()

    entry = native_dynamics.first_entry_census()
    repair = native_dynamics.repair_distance_lumping_census()
    mean_wait = int(Fraction(entry["mean_attempts"]))
    steps = range(6 * mean_wait + 1)
    transition = tuple(
        tuple(Fraction(value, repair["attempts_per_state"]) for value in row)
        for row in repair["attempt_count_matrix"]
    )
    mass = (Fraction(1), *(Fraction(0) for _ in transition[1:]))
    trajectory = []
    for _ in steps:
        trajectory.append(mass)
        mass = tuple(
            sum(probability * chance for probability, chance in zip(mass, column, strict=True))
            for column in zip(*transition, strict=True)
        )

    colours = ("#CBD5E1", BLUE, SKY, GREEN)
    repair_axis.stackplot(
        steps, *np.asarray(trajectory, dtype=float).T, colors=colours,
        labels=("Doubled", "Distance 1", "Distance 2", "Distance 3"),
    )
    stationary_at = steps[-1] + mean_wait
    bottom = 0.0
    for value, colour in zip(repair["stationary_class_masses"], colours, strict=True):
        height = float(Fraction(value))
        repair_axis.bar(stationary_at, height, width=0.58 * mean_wait, bottom=bottom, color=colour)
        bottom += height
    repair_ticks = list(range(0, steps[-1] + 1, mean_wait))
    repair_axis.set_xticks(
        [*repair_ticks, stationary_at], labels=[*(str(step) for step in repair_ticks), "Stationary"]
    )
    repair_axis.set(
        xlim=(0, stationary_at + 0.65 * mean_wait), ylim=(0, 1),
        xlabel="Native attempts", ylabel="Probability",
    )
    repair_axis.set_title("A  Repair from a doubled load", loc="left", fontweight="bold")
    repair_axis.set_yticks((0, 0.25, 0.5, 0.75, 1))
    repair_axis.grid(False, axis="x")
    repair_axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=4, fontsize=9)

    record_census = free_energy.carrier_grounded_forest_census()
    layouts = record_census["triple_layouts"]
    optimum = record_census["best_triple_determinant"]
    greedy = max(value for value, _ in record_census["greedy_triple_determinants"])
    low = min(row[2] for row in layouts)
    span = optimum - low
    records_axis.set_xlim((low - 0.11 * span) / 1e6, (optimum + 0.55 * span) / 1e6)
    records_axis.set_yticks(
        range(len(layouts)),
        labels=[f"{'–'.join(map(str, distances))}  ·  {count}" for distances, count, _, _ in layouts],
    )
    for y, (_, _, value, _) in enumerate(layouts):
        colour, marker, suffix = GREY, "o", ""
        if value == optimum:
            colour, marker, suffix = BLUE, "o", "  optimum"
        elif value == greedy:
            colour, marker, suffix = VERMILLION, "D", "  best greedy"
        records_axis.scatter(value / 1e6, y, s=65, marker=marker, color=colour, zorder=3)
        records_axis.annotate(
            f"{value:,}{suffix}", (value / 1e6, y), xytext=(8, 0),
            textcoords="offset points", va="center", fontsize=10, color=INK,
        )
    records_axis.invert_yaxis()
    records_axis.set(xlabel="Grounded determinant (millions)", ylabel="Distances · record sets")
    records_axis.set_title("B  Which three records reduce volume most?", loc="left", fontweight="bold")
    records_axis.grid(False, axis="y")

    dual_census = planar_duality.planar_duality_census()
    for y, key in enumerate(("carrier", "ivm")):
        row = dual_census[key]
        (primal,), (dual,) = row["primal_resistances"], row["dual_resistances"]
        for start, width, colour in ((0, primal, BLUE), (primal, dual, VERMILLION)):
            dual_axis.barh(y, float(width), left=float(start), height=0.49, color=colour, zorder=2)
            dual_axis.text(
                float(start + width / 2), y, str(width), va="center", ha="center",
                color="white", fontsize=13, fontweight="bold", zorder=3,
            )
    dual_axis.set_yticks((0, 1), labels=("Icosahedron and\ndodecahedron", "FCC shell and\nVoronoi cell"))
    dual_axis.set(xlim=(0, 1.02), ylim=(1.6, -0.6), xlabel="Variance across an edge")
    dual_axis.set_xticks((0, 0.25, 0.5, 0.75, 1))
    dual_axis.set_title("C  Primal and dual edge variances", loc="left", fontweight="bold")
    dual_axis.grid(False, axis="y")
    dual_axis.legend(
        handles=[plt.Rectangle((0, 0), 1, 1, color=colour, label=label)
                 for colour, label in ((BLUE, "Primal"), (VERMILLION, "Dual"))],
        loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=10,
    )

    shell = network.cuboctahedral_resistance_geometry()
    classes = shell["pair_classes"]
    full = [float(row[2]) for row in classes]
    lowest = [float(row[3]) for row in classes]
    rows = range(len(classes))
    shell_axis.hlines(list(rows), lowest, full, color="#A2AFBD", linewidth=1.4)
    shell_axis.scatter(full, list(rows), color=BLUE, s=62, label="Full variance", zorder=3)
    shell_axis.scatter(
        lowest, list(rows), color="white", edgecolor=VERMILLION, marker="D", s=64,
        linewidth=1.6, label="Lowest modes", zorder=3,
    )
    for y, (_, _, variance, projected, _) in enumerate(classes):
        for value, offset, align in ((projected, -7, "right"), (variance, 7, "left")):
            shell_axis.annotate(
                str(value), (float(value), y), xytext=(offset, 0), textcoords="offset points",
                ha=align, va="center", fontsize=10,
            )
    shell_axis.set_yticks(
        list(rows),
        labels=[f"{distance}  ·  {chord}  ({count})" for distance, chord, _, _, count in classes],
    )
    shell_axis.set(xlim=(0, max(full) + 0.13), ylim=(len(classes) - 0.52, -0.48),
                   xlabel="Unit-precision contrast variance", ylabel="Steps · chord² (pairs)")
    shell_axis.grid(False, axis="y")
    shell_axis.set_title("D  Cuboctahedral pair contrasts", loc="left", fontweight="bold")
    shell_axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=10)

    pair_details = ", ".join(
        f"graph distance {distance} and squared chord {chord}: {count} pairs "
        f"with variances {variance} and {projected}"
        for distance, chord, variance, projected, count in classes
    )
    dual_details = ", and ".join(
        f"the {name} pair contributes {row['primal_resistances'][0]} and {row['dual_resistances'][0]}"
        for name, row in dual_census.items()
    )
    return _save(
        figure, "fig_extensions.png",
        f"Panel A shows distance probabilities over {steps[-1]} native attempts from a doubled "
        "load, with the stationary distribution shown as a separate bar. First binary entry "
        "always occurs at distance one. The stationary law places positive mass on all three "
        f"binary distances, with masses {_enumerated(repair['stationary_class_masses'][1:])}. "
        f"Panel B compares the grounded determinants for all {record_census['triples_checked']} "
        f"three-record placements by sorted pair distances. The largest determinant is {optimum:,}. "
        f"Greedy extension of an optimal antipodal pair reaches {greedy:,}. "
        f"Panel C shows matched primal and dual edge variances: {dual_details}, summing to one "
        f"for each pair. Panel D compares full variance and lowest-mode variance by edge-length-two "
        f"cuboctahedral pair class: {pair_details}.",
    )


def three_way_figure() -> Path:
    """Compare the four paired tests using the exact census values."""
    figure, axes = plt.subplots(2, 2, figsize=(10.2, 6.9))
    repair_axis, dual_axis, mode_axis, plane_axis = axes.ravel()
    carrier_name, fcc_name = "Carrier", "FCC shell"
    shells = (carrier_name, fcc_name)

    repair = native_dynamics.repair_contrast_lumping_census()
    (carrier_pair,) = repair["carrier"]["distance_two_witnesses"]
    fcc_pairs = repair["cuboctahedron"]["distance_two_witnesses"]
    repair_rows = [
        (carrier_name, carrier_pair, BLUE, "o"),
        *((f"FCC · chord² {pair[1]}", pair, VERMILLION, "D") for pair in fcc_pairs),
    ]
    for y, (_, (_, _, _, _, probability), colour, marker) in enumerate(repair_rows):
        value = float(Fraction(probability))
        repair_axis.hlines(y, 0, value, color=colour, linewidth=1.7)
        repair_axis.scatter(value, y, marker=marker, color=colour, s=64, zorder=3)
        repair_axis.annotate(
            probability, (value, y), xytext=(8, 0), textcoords="offset points",
            va="center", fontsize=11, fontweight="bold",
        )
    repair_axis.set_yticks(range(len(repair_rows)), labels=[row[0] for row in repair_rows])
    repair_axis.set_xticks((0, 0.025, 0.05, 0.075, 0.10),
                           labels=("0", "2.5%", "5%", "7.5%", "10%"))
    repair_axis.set(xlim=(0, 0.107), ylim=(len(repair_rows) - 0.5, -0.5),
                    xlabel="Chance of adjacency in one attempt")
    repair_axis.set_title("A  Repair from distance two", loc="left", fontweight="bold")
    repair_axis.grid(False, axis="y")

    duality = weighted_duality.weighted_planar_duality_census()
    dual_rows = []
    for name, key in zip(shells, ("carrier", "ivm"), strict=True):
        row = duality[key]
        for policy, field in (("reciprocal", "reciprocal_dual_variance"),
                              ("copied", "copied_dual_variance")):
            dual_rows.append((f"{name} · {policy}", row["primal_variance"], row[field]))
    dual_positions = (0, 1, 2.35, 3.35)
    for y, (_, primal, dual) in zip(dual_positions, dual_rows, strict=True):
        segments = ((0, primal, BLUE, "white"), (primal, dual, ORANGE, INK))
        for left, fraction, colour, ink in segments:
            dual_axis.barh(y, float(fraction), left=float(left), height=0.58,
                           color=colour, edgecolor="white", linewidth=0.5, zorder=2)
            dual_axis.text(float(left + fraction / 2), y, str(fraction), ha="center",
                           va="center", color=ink, fontsize=11, fontweight="bold", zorder=4)
    dual_axis.axvline(1, color=INK, linewidth=1.15, linestyle="--", zorder=3)
    dual_axis.set_yticks(dual_positions, labels=[row[0] for row in dual_rows])
    dual_axis.set_xticks((0, 0.5, 1, 1.5))
    dual_axis.set(xlim=(0, 1.65), ylim=(3.85, -0.5), xlabel="Normalised variance at the changed seam")
    dual_axis.set_title("B  Weighted primal–dual balance", loc="left", fontweight="bold")
    dual_axis.grid(False, axis="y")
    dual_axis.legend(
        handles=[plt.Rectangle((0, 0), 1, 1, color=colour, label=label)
                 for colour, label in ((BLUE, "Primal"), (ORANGE, "Dual"))],
        loc="upper center", bbox_to_anchor=(0.5, -0.23), ncol=2,
    )

    modes = quadratic_modes.quadratic_mode_census()
    mode_values = [(modes[key]["variance"][0], modes[key]["variance"][2])
                   for key in ("carrier", "fcc")]
    for index, (label, colour, shift) in enumerate(
        (("Diagonal", GREEN, -0.20), ("Off-diagonal", PURPLE, 0.20)),
    ):
        values = [row[index] for row in mode_values]
        positions = [group + shift for group in range(len(shells))]
        bars = mode_axis.bar(positions, [float(value) for value in values], width=0.34,
                             color=colour, edgecolor=INK, linewidth=0.5, label=label)
        mode_axis.bar_label(bars, labels=[str(value) for value in values], padding=3,
                            fontsize=11, fontweight="bold")
    mode_axis.set_xticks(range(len(shells)), labels=shells)
    mode_axis.set_yticks((0, 0.05, 0.10, 0.15))
    mode_axis.set(xlim=(-0.55, 1.55), ylim=(0, 0.194),
                  ylabel="Variance of a unit coefficient")
    mode_axis.set_title("C  Gaussian quadratic modes", loc="left", fontweight="bold")
    mode_axis.grid(False, axis="x")
    mode_axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)

    planes = geometric_screening.central_plane_blanket_census()
    plane_values = []
    for y, key in enumerate(("icosahedron", "cuboctahedron")):
        total = planes[key]["vertex_spanned_planes"]
        screened = planes[key]["screening_planes"]
        crossed = total - screened
        plane_values.append((screened, crossed))
        for left, count, colour, label in (
            (0, screened, BLUE, "screened"),
            (screened, crossed, VERMILLION, "crossed"),
        ):
            if count:
                plane_axis.barh(y, count, left=left, color=colour, height=0.48)
                plane_axis.annotate(
                    f"{count} {label}", (left + count, y), xytext=(7, 0),
                    textcoords="offset points", va="center", fontsize=11, fontweight="bold",
                )
    plane_axis.set_yticks(range(len(shells)), labels=shells)
    plane_axis.set_xticks((0, 5, 10, 15))
    plane_axis.set(xlim=(0, 22), ylim=(1.55, -0.55), xlabel="Vertex-spanned central planes")
    plane_axis.set_title("D  Do central planes screen?", loc="left", fontweight="bold")
    plane_axis.grid(False, axis="y")

    repair_details = " and ".join(
        f"chord squared {pair[1]} gives probability {pair[4]} and variance {pair[2]}"
        for pair in fcc_pairs
    )
    dual_details = ". ".join(
        f"{name} has primal variance {row['primal_variance']}, reciprocal dual "
        f"variance {row['reciprocal_dual_variance']} and copied dual variance "
        f"{row['copied_dual_variance']}, giving sums of one and {row['copied_variance_sum']}"
        for name, row in ((carrier_name, duality["carrier"]), (fcc_name, duality["ivm"]))
    )
    mode_details = " and ".join(
        f"{name} has {diagonal} and {mixed}"
        for name, (diagonal, mixed) in zip(shells, mode_values, strict=True)
    )
    plane_details = " and ".join(
        f"{name} has {screened} screening and {crossed} crossing planes"
        for name, (screened, crossed) in zip(shells, plane_values, strict=True)
    )
    return _save(
        figure, "fig_three_way.png",
        "Panel A compares pairs at graph distance two. In one repair attempt the carrier "
        f"pair reaches adjacency with probability {carrier_pair[4]} and has contrast variance "
        f"{carrier_pair[2]}. On the FCC shell, {repair_details}, with edge length two. "
        f"Panel B increases one primal conductance to {duality['carrier']['primal_weight']} "
        f"and compares normalised variances. {dual_details}. "
        "Panel C compares the variances of unit diagonal and off-diagonal trace-free "
        f"quadratic coefficients, respectively. {mode_details}. "
        "Panel D counts vertex-spanned central planes. Screening means no seam crosses "
        f"between the strict halves. {plane_details}.",
    )


def repair_transfer_figure() -> Path:
    """Show which geometric summaries transfer repair dynamics."""
    figure, axes = plt.subplots(2, 2, figsize=(12.3, 7.2))
    cube_axis, fcc_axis, family_axis, lag_axis = axes.ravel()

    cube = contrast_quotients.distance_regular_repair_census()
    variance = tuple(map(Fraction, cube["resistances"]))
    distances = tuple(range(1, len(variance) + 1))
    cube_axis.plot(distances, list(map(float, variance)), color=BLUE, marker="o", linewidth=2.2)
    for distance, value in zip(distances, variance, strict=True):
        cube_axis.annotate(str(value), (distance, float(value)), xytext=(0, 9),
                           textcoords="offset points", ha="center", color=BLUE, fontsize=12)
    cube_axis.set_xticks(
        distances, [f"{distance}\n{size} pairs" for distance, size in
                    zip(distances, cube["class_sizes"], strict=True)],
    )
    cube_axis.set(xlim=(0.7, len(distances) + 0.3),
                  ylim=(float(min(variance)) - 0.04, float(max(variance)) + 0.06),
                  xlabel="Graph distance", ylabel="Gaussian contrast variance")
    cube_axis.set_title("A  Distance fixes contrast on the 4-cube", loc="left", fontweight="bold")

    fcc = contrast_quotients.periodic_fcc_contrast_census()
    witness_rates = tuple(Fraction(row[-1]) for row in fcc["witnesses"])
    for row_index, rate in enumerate(witness_rates):
        fcc_axis.barh(row_index, float(rate), height=0.42, color=VERMILLION)
        fcc_axis.scatter([float(rate)], [row_index], s=35, color=VERMILLION, zorder=3)
        fcc_axis.annotate(str(rate), (float(rate), row_index), xytext=(7, 0),
                          textcoords="offset points", va="center", color=VERMILLION, fontsize=12)
    fcc_axis.set_yticks(
        range(len(witness_rates)), [str(row[0]) for row in fcc["witnesses"]],
    )
    fcc_axis.set(
        xlim=(0, 1.37 * float(max(witness_rates))), ylim=(1.5, -0.5),
        xlabel=f"One-attempt probability of variance {fcc['target_variance']}",
        ylabel="Pair displacement mod 4",
    )
    fcc_axis.set_title(
        f"B  FCC pairs at variance {fcc['witnesses'][0][3]} move differently",
        loc="left", fontweight="bold",
    )

    family = tuple(contrast_quotients.bipartite_contrast_census(size) for size in (2, 3, 4))
    positions = np.arange(len(family))
    for field, label, colour, marker, offset in (
        ("stationary_next_target_given_collision", "Present only", BLUE, "o", -15),
        ("stationary_next_target_given_previous_target_and_collision",
         "Target at previous step", VERMILLION, "D", 9),
    ):
        rates = tuple(Fraction(row[field]) for row in family)
        family_axis.plot(positions, list(map(float, rates)), color=colour, marker=marker,
                         linewidth=2, label=label)
        for position, rate in zip(positions, rates, strict=True):
            family_axis.annotate(str(rate), (position, float(rate)), xytext=(0, offset),
                                 textcoords="offset points", ha="center", color=colour, fontsize=11)
    family_axis.set_xticks(
        positions, [f"$K_{{{left},{right}}}$" for left, right in (row["side_sizes"] for row in family)],
    )
    family_axis.set(xlim=(-0.3, len(family) - 0.7), ylim=(0.06, 0.20),
                    ylabel="Probability of larger-side pair")
    family_axis.set_title("C  Previous contrast changes the next step", loc="left", fontweight="bold")
    family_axis.legend(loc="upper right", frameon=False, fontsize=10.5)

    lag = root_shells.root_shell_lag_census()
    defects = tuple(lag[key]["lag_two"] - lag[key]["lag_one"]**2 for key in ("unit", "weighted"))
    for row_index, (defect, colour) in enumerate(zip(defects, (BLUE, VERMILLION), strict=True)):
        lag_axis.barh(row_index, float(defect), height=0.42, color=colour)
        lag_axis.scatter([float(defect)], [row_index], s=35, color=colour, zorder=3)
        lag_axis.annotate(str(defect), (float(defect), row_index), xytext=(7, 0),
                          textcoords="offset points", va="center", color=colour, fontsize=12)
    lag_axis.set_yticks((0, 1), ("Uniform seams", "One weighted seam"))
    lag_axis.set(xlim=(0, 1.37 * float(max(defects))), ylim=(1.5, -0.5),
                 xlabel=r"Lag-two covariance excess along the chord, $M_2-M_1^2$")
    lag_axis.set_title("D  A weighted 24-cell seam breaks closure", loc="left", fontweight="bold")

    for axis in axes.ravel():
        axis.tick_params(axis="both", labelsize=11.5)
    for axis in (fcc_axis, lag_axis):
        axis.grid(False, axis="y")
        axis.grid(True, axis="x", alpha=0.24)

    cube_details = ". ".join(
        f"Distance {distance} has variance {value} for {size} pairs"
        for distance, value, size in zip(distances, variance, cube["class_sizes"], strict=True)
    )
    fcc_details = ". ".join(
        f"Displacement {offset} has distance {distance}, {common} common neighbours, variance {value} "
        f"and target probability {rate}"
        for offset, distance, common, value, rate in fcc["witnesses"]
    )
    family_details = ". ".join(
        f"For side sizes {row['side_sizes']}, the probability is {row['stationary_next_target_given_collision']} "
        f"with the current contrast alone and {row['stationary_next_target_given_previous_target_and_collision']} "
        "with the previous target also observed"
        for row in family
    )
    return _save(
        figure, "fig_repair_transfer.png",
        f"Panel A is the 4-cube contrast ladder. {cube_details}. "
        f"Panel B shows two pairs on the period-{fcc['period']} FCC graph and the chance of reaching "
        f"variance {fcc['target_variance']} in one attempt. {fcc_details}. "
        f"Panel C shows stationary next-step probabilities on complete bipartite graphs. {family_details}. "
        f"Panel D subtracts squared lag one from measured lag two in the seam direction of the 24-cell. "
        f"Uniform seams give {defects[0]}. Weighting one seam gives measured lag two "
        f"{lag['weighted']['lag_two']} and squared lag one {lag['weighted']['lag_one']**2}, "
        f"an excess of {defects[1]}.",
    )


def centroid_memory_figure() -> Path:
    """Show common lag covariances and the conditional return witness."""
    figure, (lag_axis, return_axis) = plt.subplots(1, 2, figsize=(9.5, 3.65))
    census = centroid_dynamics.centroid_gaussian_census()
    shells = ("Carrier", "FCC shell")
    rows = (census["carrier"], census["cuboctahedron"])
    lags = np.arange(61)
    for row, name, colour, linestyle in zip(
        rows, shells, (BLUE, VERMILLION), ("-", "--"), strict=True,
    ):
        constant, radical = (float(Fraction(part)) for part in row["retention_parts"])
        retention = constant + radical * math.sqrt(5)
        lag_axis.plot(lags, retention**lags, color=colour, linestyle=linestyle,
                      linewidth=2.2, label=name)
    lag_axis.set(xlim=(0, 60), ylim=(0, 1.02), xlabel="Lag (repair attempts)",
                 ylabel="Lag covariance / variance")
    lag_axis.set_xticks((0, 15, 30, 45, 60))
    lag_axis.set_yticks((0, 0.25, 0.5, 0.75, 1))
    lag_axis.set_title("A  Shared repair and Gaussian covariance", loc="left", fontweight="bold")
    lag_axis.legend(loc="upper right")

    histories = (
        ("At zero", "next_target_given_zero", SKY, -0.19),
        ("At zero after target", "next_target_given_previous_target_and_zero", BLUE, 0.19),
    )
    for label, field, colour, shift in histories:
        values = [row[field] for row in rows]
        bars = return_axis.bar(
            [index + shift for index in range(len(shells))],
            [float(Fraction(value)) for value in values], width=0.32,
            color=colour, edgecolor=INK, linewidth=0.5, label=label,
        )
        return_axis.bar_label(bars, labels=values, padding=3, fontsize=11, fontweight="bold")
    return_axis.set_xticks(range(len(shells)), labels=shells)
    return_axis.set_yticks((0, 0.01, 0.02), labels=("0", "1%", "2%"))
    return_axis.set(xlim=(-0.55, 1.55), ylim=(0, 0.0255), ylabel="Probability in one attempt")
    return_axis.set_title("B  Return to the target centroid", loc="left", fontweight="bold")
    return_axis.grid(False, axis="x")
    return_axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)

    alt_shells = ("carrier", "FCC shell")
    retention_details = " and ".join(
        f"{row['retention']} for the {name}"
        for name, row in zip(alt_shells, rows, strict=True)
    )
    return_details = " and ".join(
        f"{row['next_target_given_zero']} and "
        f"{row['next_target_given_previous_target_and_zero']} on the {name}"
        for name, row in zip(alt_shells, rows, strict=True)
    )
    return _save(
        figure, "fig_centroid_memory.png",
        "Panel A shows common normalised lag covariance for stationary repair and the matched "
        "sampled Ornstein–Uhlenbeck process. On each shell one curve shows both processes. "
        f"The covariance is alpha to the power of the lag, with retention alpha {retention_details}. "
        "Panel B shows the probability of returning to a nonzero target on the next attempt, "
        "conditioning on a current centroid of zero alone or additionally on being at the "
        f"target in the previous attempt. The probabilities are {return_details}, respectively.",
    )


def geometry_transfer_figure() -> Path:
    """Compare topology, quadratic symmetry and two aspects of Gaussian screening."""
    figure, axes = plt.subplots(2, 2, figsize=(10.3, 6.9))
    torus_axis, root_axis, sum_axis, record_axis = axes.ravel()

    torus = surface_duality.toroidal_duality_census()
    torus_cases = (("Unit", torus["unit"]), ("Reciprocal", torus["reciprocal"]),
                   ("Copied", torus["copied"]))
    components = (
        ("Primal", "primal_diagonal", BLUE, "white"),
        ("Dual", "dual_diagonal", ORANGE, INK),
        ("Remainder", "remainder_diagonal", GREEN, INK),
    )
    positions = (0, 1.1, 2.2)
    for y, (_, case) in zip(positions, torus_cases, strict=True):
        left = Fraction(0)
        for _, field, colour, ink in components:
            value = case[field]
            start = left if value >= 0 else Fraction(0)
            midpoint = float(start + value / 2)
            torus_axis.barh(y, float(value), left=float(start), height=0.43, color=colour,
                            edgecolor="white", linewidth=0.6, zorder=2)
            if field == "remainder_diagonal" and value > 0:
                torus_axis.annotate(str(value), (midpoint, y), xytext=(0, 10),
                                    textcoords="offset points", ha="center", va="bottom",
                                    fontsize=10.5, fontweight="bold")
            else:
                torus_axis.text(midpoint, y, str(value), va="center", ha="center",
                                color=ink, fontsize=10.5, fontweight="bold", zorder=4)
            if value > 0:
                left += value
    torus_axis.set_yticks(positions, labels=[f"{name} · rank {row['remainder_rank']}"
                                            for name, row in torus_cases])
    torus_axis.axvline(0, color=INK, linewidth=1.0, zorder=3)
    torus_axis.axvline(1, color=GREY, linewidth=1.0, linestyle="--", zorder=1)
    torus_axis.set_xticks((-0.5, 0, 0.5, 1, 1.5))
    torus_axis.set(xlim=(-0.62, 1.62), ylim=(2.7, -0.63), xlabel="Normalised diagonal at the selected edge")
    torus_axis.set_title("A  The torus leaves a remainder", loc="left", fontweight="bold")
    torus_axis.grid(False, axis="y")
    torus_axis.legend(
        handles=[plt.Rectangle((0, 0), 1, 1, color=colour, label=label)
                 for label, _, colour, _ in components],
        loc="upper center", bbox_to_anchor=(0.5, -0.21), ncol=3,
    )

    quadratic = root_shells.root_shell_quadratics()
    dimensions = [case[0] for case in quadratic]
    variance_ratios = [case[5][0] / case[5][1] for case in quadratic]
    root_axis.axhline(1, color=GREY, linewidth=1.1, linestyle="--", zorder=1)
    root_axis.plot(dimensions, list(map(float, variance_ratios)), linewidth=1.6, color=BLUE, zorder=2)
    for dimension, ratio in zip(dimensions, variance_ratios, strict=True):
        isotropic = ratio == 1
        root_axis.scatter(dimension, float(ratio), marker="D" if isotropic else "o",
                          facecolor="white" if isotropic else BLUE, edgecolor=BLUE,
                          linewidth=1.5, s=58, zorder=3)
        root_axis.annotate(
            str(ratio), (dimension, float(ratio)), xytext=(0, 8), textcoords="offset points",
            ha="center", va="bottom", fontsize=11, fontweight="bold",
        )
    root_axis.set_xticks(dimensions, labels=[fr"$D_{{{value}}}$" for value in dimensions])
    root_axis.set_yticks((0.5, 0.75, 1, 1.25))
    root_axis.set(xlim=(dimensions[0] - 0.38, dimensions[-1] + 0.38), ylim=(0.5, 1.46),
                  xlabel="Root shell", ylabel="Diagonal / mixed variance")
    root_axis.set_title("B  Quadratic variance reverses order", loc="left", fontweight="bold")
    root_axis.grid(False, axis="x")

    screening = constraint_screening.fcc_sum_constraint_census()
    frequencies = (1, 2)
    plane_cases = (("Axis plane", "axis", BLUE, -0.19), ("Diagonal plane", "diagonal", SKY, 0.19))
    squared_correlations = {
        frequency: {
            key: screening[f"frequency_{frequency}_{key}"]["soft_canonical_correlation_squared"]
            for _, key, _, _ in plane_cases
        }
        for frequency in frequencies
    }
    for label, key, colour, shift in plane_cases:
        values = [squared_correlations[frequency][key] for frequency in frequencies]
        bars = sum_axis.bar([frequency + shift for frequency in frequencies], list(map(float, values)),
                            width=0.32, color=colour, edgecolor=INK, linewidth=0.5, label=label)
        sum_axis.bar_label(bars, labels=[f"{float(value):.1%}" for value in values],
                           padding=3, fontsize=10, fontweight="bold")
    shell_labels = [f"{frequency} · {screening[f'frequency_{frequency}_axis']['site_count']} sites"
                    for frequency in frequencies]
    sum_axis.set_xticks(frequencies, labels=shell_labels)
    sum_axis.set_yticks((0, 0.25, 0.5, 0.75, 1), labels=("0", "25%", "50%", "75%", "100%"))
    sum_axis.set(xlim=(0.47, 2.53), ylim=(0, 1.08), xlabel="FCC boundary frequency and size",
                 ylabel="Squared canonical correlation")
    sum_axis.set_title("C  A noisy global sum couples the sides", loc="left", fontweight="bold")
    sum_axis.grid(False, axis="x")
    sum_axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)

    record_sets = [
        (frequency, key, colour, screening[f"frequency_{frequency}_{key}"]["partition_sizes"][0],
         screening[f"frequency_{frequency}_{key}"]["site_count"])
        for frequency in frequencies for _, key, colour, _ in plane_cases
    ]
    record_positions = (0, 1, 2.35, 3.35)
    for y, (_, _, colour, recorded, total) in zip(record_positions, record_sets, strict=True):
        share = float(Fraction(recorded, total))
        record_axis.barh(y, share, height=0.48, color=colour, edgecolor=INK, linewidth=0.5)
        record_axis.annotate(f"{recorded} of {total}", (share, y), xytext=(7, 0),
                             textcoords="offset points", va="center", fontsize=11,
                             fontweight="bold")
    record_axis.set_yticks(record_positions, labels=[f"{frequency} · {key.capitalize()}"
                                                     for frequency, key, *_ in record_sets])
    record_axis.set_xticks((0, 0.15, 0.30, 0.45, 0.60), labels=("0", "15%", "30%", "45%", "60%"))
    record_axis.set(xlim=(0, 0.64), ylim=(3.85, -0.5), xlabel="Share of boundary sites on the plane",
                    ylabel="Frequency · plane")
    record_axis.set_title("D  How many sites does the plane record?", loc="left", fontweight="bold")
    record_axis.grid(False, axis="y")

    torus_details = ". ".join(
        f"{name} has diagonal entries {_enumerated(row[field] for _, field, _, _ in components)} "
        f"and remainder rank {row['remainder_rank']}"
        for name, row in torus_cases
    )
    screening_details = ". ".join(
        f"Frequency {frequency} has axis and diagonal squared correlations "
        f"{_enumerated(squared_correlations[frequency][key] for _, key, _, _ in plane_cases)}"
        for frequency in frequencies
    )
    record_details = ". ".join(
        f"At frequency {frequency}, the axis and diagonal planes record "
        f"{_enumerated(recorded for level, _, _, recorded, _ in record_sets if level == frequency)} "
        f"of {screening[f'frequency_{frequency}_axis']['site_count']} sites, respectively"
        for frequency in frequencies
    )
    return _save(
        figure, "fig_geometry_transfer.png",
        "Panel A shows primal, dual and remainder diagonals at the same torus edge. "
        f"The weighted cases set the selected primal conductance to {torus['changed_weight']}. "
        f"{torus_details}. The remainder trace is {torus['unit']['remainder_trace']} in every case. "
        "Panel B shows the ratio of diagonal to mixed quadratic coefficient variance. "
        f"For dimensions {_enumerated(dimensions)} the ratios are {_enumerated(variance_ratios)}, "
        "respectively. The modes agree at dimension four. Panel C shows how conditioning "
        "on a global sum with unit noise induces dependence across an FCC separating plane. "
        f"{screening_details}, respectively. Panel D shows the share of boundary sites on "
        f"each separating plane. {record_details}.",
    )


def jitterbug_figure() -> Path:
    """Separate the moving bars, the convex hull, and the inference model."""
    figure = plt.figure(figsize=(10.5, 6.2))
    layout = figure.add_gridspec(3, 4, height_ratios=(2.25, 0.32, 2.3))
    ico = jitterbug.ICOSAHEDRON_PARAMETER
    intermediate = Fraction(4, 5)
    landmarks = (
        (Fraction(1), "Cuboctahedron", "$t=1$"),
        (intermediate, "Intermediate hull", "$t=4/5$"),
        (ico, "Regular icosahedral hull", r"$t=(\sqrt{5}-1)/2$"),
        (Fraction(0), "Octahedral limit", "$t=0$ · paired joints"),
    )
    bars = jitterbug.jitterbug_rigid_edges()
    added = jitterbug.jitterbug_matching_edges()
    # A common oblique view separates the vertices from unrelated projected
    # edges across all four landmarks. Every panel uses bar length two.
    rotation = _rotation(math.radians(45), math.radians(35))
    for column, (parameter, title, label) in enumerate(landmarks):
        axis = figure.add_subplot(layout[0, column])
        position = float(parameter)
        points = np.asarray(jitterbug.scaled_jitterbug_vertices(parameter))
        points, inverse = np.unique(points, axis=0, return_inverse=True)
        view = points @ rotation.T
        visible_bars = sorted({tuple(sorted((inverse[a], inverse[b]))) for a, b in bars})
        visible_diagonals = [(inverse[a], inverse[b]) for a, b in added] if 0 < position < 1 else []
        segments = [
            (float(view[[a, b], 2].mean()), a, b, colour, width, linestyle)
            for edges, colour, width, linestyle in (
                (visible_bars, BLUE, 1.75 if position else 2.45, "-"),
                (visible_diagonals, VERMILLION, 1.55, (0, (3.0, 2.0))),
            )
            for a, b in edges
        ]
        for depth, a, b, colour, width, linestyle in sorted(segments):
            axis.plot(view[[a, b], 0], view[[a, b], 1], color=colour, linewidth=width,
                      linestyle=linestyle, solid_capstyle="round", alpha=0.94 if depth >= 0 else 0.48)
        axis.scatter(view[:, 0], view[:, 1], color=INK, alpha=0.95, s=15 if position else 23,
                     edgecolors="white", linewidths=0.45, zorder=4)
        axis.set(xlim=(-2.04, 2.04), ylim=(-2.04, 2.04), aspect="equal")
        axis.set_title(f"{title}\n{label}", fontsize=11.5, linespacing=1.3, fontweight="bold")
        axis.set_axis_off()

    legend_axis = figure.add_subplot(layout[1, :])
    legend_axis.set_axis_off()
    legend_axis.legend(
        handles=(
            plt.Line2D([], [], color=BLUE, linewidth=2.0, label="Persistent bars"),
            plt.Line2D([], [], color=VERMILLION, linewidth=1.8, linestyle=(0, (3.0, 2.0)),
                       label="Diagonals added by the convex hull"),
        ),
        loc="center", ncol=2, columnspacing=2.6, handlelength=2.4,
    )

    length_axis = figure.add_subplot(layout[2, :2])
    weight_axis = figure.add_subplot(layout[2, 2:])
    t_values = np.linspace(0, 1, 121)
    bar_a, bar_b = bars[0]
    diag_a, diag_b = added[0]

    def relative_diagonal(parameter: float) -> float:
        points = np.asarray(jitterbug.scaled_jitterbug_vertices(float(parameter)))
        return float(np.linalg.norm(points[diag_a] - points[diag_b])
                     / np.linalg.norm(points[bar_a] - points[bar_b]))

    length_axis.plot(t_values, [relative_diagonal(parameter) for parameter in t_values],
                     color=VERMILLION, linewidth=2.1)
    length_axis.axhline(1, color=GREY, linestyle="--", linewidth=1.05)
    length_axis.axvline(float(ico), color=GREY, linestyle="--", linewidth=1.05)
    length_axis.scatter([float(ico)], [relative_diagonal(float(ico))], color=VERMILLION,
                        edgecolors="white", linewidths=0.8, s=60, zorder=4)
    length_axis.annotate("Equal lengths", (float(ico), 1), xytext=(-8, 9),
                         textcoords="offset points", ha="right", fontsize=10)
    length_axis.set(
        xlim=(-0.025, 1.025), ylim=(0, 1.51), xlabel="Contraction parameter $t$",
        ylabel="Short diagonal / bar length",
    )
    length_axis.set_xticks((0, float(ico), 1), labels=("$0$", r"$1/\varphi$", "$1$"))
    length_axis.set_yticks((0, 0.5, 1, math.sqrt(2)), labels=("$0$", "$0.5$", "$1$", r"$\sqrt{2}$"))
    length_axis.set_title("Hull geometry", loc="left", fontweight="bold")
    length_axis.grid(False, axis="x")

    # Plot the tuned inference model only on the cuboctahedral-to-icosahedral
    # interval, where its additional positive couplings range from zero to one.
    weight_parameters = (ico, intermediate, Fraction(1))
    weights = [jitterbug_inference.jitterbug_coordinate_inference(parameter)["tuned_weight"]
               for parameter in weight_parameters]
    weight_grid = [ico + Fraction(step, 72) * (1 - ico) for step in range(73)]
    weight_line = [jitterbug_inference.jitterbug_coordinate_inference(parameter)["tuned_weight"]
                   for parameter in weight_grid]
    weight_axis.plot(list(map(float, weight_grid)), list(map(float, weight_line)),
                     color=GREEN, linewidth=2.1)
    weight_axis.scatter(list(map(float, weight_parameters)), list(map(float, weights)),
                        color=GREEN, edgecolors="white", linewidths=0.8, s=55, zorder=4)
    for parameter, value, align, offset in zip(
        weight_parameters, weights, ("left", "center", "right"),
        ((7, -9), (0, 9), (-7, 9)), strict=True,
    ):
        weight_axis.annotate(
            format_value(value), (float(parameter), float(value)), xytext=offset,
            textcoords="offset points", ha=align, va="bottom", fontsize=11, fontweight="bold",
        )
    weight_axis.set(xlim=(float(ico) - 0.026, 1.035), ylim=(-0.05, 1.09),
                    xlabel="Contraction parameter $t$", ylabel=r"Tuned conductance $w_*=1/t-t$")
    weight_axis.set_xticks(list(map(float, weight_parameters)),
                          labels=(r"$1/\varphi$", "$4/5$", "$1$"))
    weight_axis.set_yticks((0, 0.25, 0.5, 0.75, 1))
    weight_axis.set_title("Gaussian inference model", loc="left", fontweight="bold")
    weight_axis.grid(False, axis="x")

    return _save(
        figure, "fig_jitterbug.png",
        "Four views share one rotation, projection and scale. The cuboctahedron at t=1, "
        "intermediate hull at t=4/5, and regular icosahedral hull at t=1/phi each have "
        "twelve separate joints. At t=0 the joints meet in pairs at six octahedral vertices, "
        f"with two of the original {len(bars)} bars along each apparent octahedral edge. "
        f"In the intermediate and regular views, {len(added)} dashed orange lines show diagonal "
        "edges added by the convex hull of the joints. "
        "The short diagonal divided by the fixed bar length increases from zero at t=0 to "
        "one at t=1/phi and square root of two at t=1. The inference plot covers t from 1/phi "
        "to 1. The added matching conductance w equals 1/t minus t and takes the values "
        f"{_enumerated(weights)} at t=1/phi, 4/5 and 1, respectively.",
    )


def jitterbug_poster_figure() -> Path:
    """Write the manuscript image and a standalone vector poster."""
    if __package__:
        from .jitterbug_poster import draw_jitterbug_poster
    else:
        from jitterbug_poster import draw_jitterbug_poster

    figure, alt = draw_jitterbug_poster(_rotation(math.radians(45), math.radians(35)))
    FIGURES.mkdir(parents=True, exist_ok=True)
    target = FIGURES / "fig_ga_jitterbug.png"
    figure.savefig(target, dpi=300)
    with matplotlib.rc_context({"pdf.fonttype": 42}):
        figure.savefig(
            FIGURES / "OPH-FEP-Synergetics_Jitterbug_Poster.pdf",
            metadata={"Title": "Fuller's Jitterbug: From an FCC shell to the OPH carrier",
                      "Author": "Daniel Ari Friedman, Bernhard Mueller", "CreationDate": None},
        )
    plt.close(figure)
    ALT_TEXT[target.name] = alt
    return target


def write_registry() -> Path:
    """Write the accessibility registry the renderer reads for figure alt text.

    A label is the file stem with underscores turned into hyphens, so that
    `fig_ga_triangle.png` is `fig:ga-triangle` and matches the `{#fig:ga-triangle}`
    label the manuscript writes under the same figure.
    """
    FIGURES.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": "1.2",
        "generator": "scripts/make_figures.py",
        "figures": [
            {
                "label": "fig:" + name.removeprefix("fig_").removesuffix(".png").replace("_", "-"),
                "filename": name,
                "alt_text": alt,
                "alt": alt,
            }
            for name, alt in sorted(ALT_TEXT.items())
        ],
    }
    target = FIGURES / "figure_registry.json"
    target.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    return target


def main() -> int:
    """Render every figure, then the registry, and print each written path."""
    _style()
    builders = (
        carrier_figure,
        registry_figure,
        descent_figure,
        dstability_figure,
        mediator_figure,
        blanket_figure,
        orbit_figure,
        ivm_figure,
        subdivision_figure,
        cover_figure,
        ga_triangle_figure,
        ga_carrier_figure,
        ga_descent_figure,
        ga_synergetics_figure,
        jitterbug_poster_figure,
        pair_triangles_figure,
        scope_matrix_figure,
        network_figure,
        vertex_shells_figure,
        conjugation_figure,
        families_figure,
        sandwich_figure,
        schedule_figure,
        shell_gauge_figure,
        extended_claims_figure,
        three_way_figure,
        centroid_memory_figure,
        repair_transfer_figure,
        geometry_transfer_figure,
        jitterbug_figure,
    )
    unrendered = {
        name
        for name, obj in globals().items()
        if name.endswith("_figure") and callable(obj) and obj not in builders
    }
    if unrendered:
        raise RuntimeError(
            f"figure builders missing from the builders tuple: {sorted(unrendered)}"
        )
    for builder in builders:
        print(builder())
    print(write_registry())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

