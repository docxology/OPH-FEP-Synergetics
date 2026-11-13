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
    OPH_FEP,
    PAIRS,
    SYNERGETICS_OPH,
    THREE_WAY,
    by_key,
    by_pair,
)
from oph_fep_synergetics.core import (  # noqa: E402
    controls,
    free_energy,
    invariants,
    network,
    sheaf,
    symmetry,
)
from oph_fep_synergetics.core.carrier import (  # noqa: E402
    OPH_ORIENTED_BASE_FACES,
    IcosahedralCarrier,
    icosahedral_carrier,
)
from oph_fep_synergetics.geometry import lattice  # noqa: E402
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
            "savefig.dpi": 200,
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "mathtext.default": "regular",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.25,
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
    """The carrier, drawn in the spectral embedding of its own Laplacian."""
    values = as_dict()
    carrier, flat = _projected_carrier()

    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 5.0))

    for face in OPH_ORIENTED_BASE_FACES:
        polygon = flat[list(face)]
        left.fill(polygon[:, 0], polygon[:, 1], color=SKY, alpha=0.28, linewidth=0)
    for port in range(carrier.port_count):
        for other in carrier.neighbours(port):
            if other > port:
                left.plot(
                    [flat[port, 0], flat[other, 0]],
                    [flat[port, 1], flat[other, 1]],
                    color=INK,
                    linewidth=1.1,
                    zorder=2,
                )
    left.scatter(flat[:, 0], flat[:, 1], s=210, color="white", edgecolor=INK, zorder=3)
    for port in range(carrier.port_count):
        left.text(
            flat[port, 0],
            flat[port, 1],
            str(port),
            ha="center",
            va="center",
            fontsize=8.5,
            color=INK,
            zorder=4,
        )
    left.set_aspect("equal")
    left.set_axis_off()
    left.set_title(
        f"{carrier.port_count} ports, {carrier.seam_count} seams, {carrier.face_count} faces, "
        f"every port of degree {carrier.degree(0)}"
    )

    spectrum = carrier.laplacian_spectrum()
    counts = {round(float(value), 6): 0 for value in spectrum}
    for value in spectrum:
        counts[round(float(value), 6)] += 1
    marks = sorted(counts)
    right.vlines(marks, 0, [counts[m] for m in marks], color=BLUE, linewidth=2.2)
    right.scatter(marks, [counts[m] for m in marks], color=BLUE, zorder=3, s=34)
    for value in marks:
        right.annotate(
            f"{value:.4g} ×{counts[value]}",
            (value, counts[value]),
            textcoords="offset points",
            xytext=(5, 5),
            fontsize=9,
        )
    lam1_value = carrier.spectral_gap()
    normalised = lam1_value / (2.0 * carrier.seam_count)
    first = next(value for value in marks if value > 1e-9)
    right.annotate(
        f"$\\lambda_1 = 5 - \\sqrt{{5}} = {lam1_value:.6f}$\n"
        f"normalised gap $\\lambda_1/{2 * carrier.seam_count} = {normalised:.6f}$",
        xy=(first, counts[first] + 0.05),
        xytext=(0.35, 5.7),
        fontsize=9,
        arrowprops={"arrowstyle": "->", "color": VERMILLION, "lw": 1.0},
        color=VERMILLION,
    )
    right.set_xlabel("Laplacian eigenvalue")
    right.set_ylabel("multiplicity")
    right.set_ylim(0, 6.8)
    right.set_xlim(-0.5, max(marks) + 1.2)
    right.set_title("Closed-form spectrum of the carrier")

    return _save(
        figure,
        "fig_carrier.png",
        f"Two panels. Left: the {int(values['carrier.ports'])} ports of the carrier drawn in the "
        f"three-dimensional spectral embedding of its own Laplacian, with the "
        f"{int(values['carrier.seams'])} seams drawn and the {int(values['carrier.faces'])} "
        f"oriented faces shaded; every port is labelled with its index and meets "
        f"{int(values['carrier.degree'])} others, and the planar projection is chosen so that "
        "no port lies near an edge it does not belong to. "
        "Right: the Laplacian spectrum, with the zero eigenvalue once, five minus root five "
        "three times, six five times and five plus root five three times, and the normalised "
        "gap marked.",
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
    offset = flat[:, None, :] - start[None, :, :]
    fraction = np.clip(np.einsum("vej,ej->ve", offset, direction) / length_squared, 0.0, 1.0)
    nearest = start[None, :, :] + fraction[:, :, None] * direction[None, :, :]
    distance = np.linalg.norm(flat[:, None, :] - nearest, axis=2)
    own = (np.arange(flat.shape[0])[:, None] == edges[None, :, 0]) | (
        np.arange(flat.shape[0])[:, None] == edges[None, :, 1]
    )
    return float(np.min(np.where(own, np.inf, distance)))


def _projected_carrier() -> tuple[IcosahedralCarrier, np.ndarray]:
    """Return the carrier and its deterministic planar spectral projection."""
    carrier = icosahedral_carrier()
    eigenvalues, eigenvectors = np.linalg.eigh(carrier.laplacian())
    lam1 = 5.0 - math.sqrt(5.0)
    columns = [index for index, value in enumerate(eigenvalues) if abs(value - lam1) < 1e-9]
    if len(columns) != 3:
        raise RuntimeError(f"expected a three-dimensional lambda_1 eigenspace, got {len(columns)}")
    points = eigenvectors[:, columns] * math.sqrt(carrier.port_count)
    edges = np.array(
        [
            (port, other)
            for port in range(carrier.port_count)
            for other in carrier.neighbours(port)
            if other > port
        ]
    )
    best: tuple[float, np.ndarray] = (-1.0, np.zeros((carrier.port_count, 2)))
    for theta in np.linspace(0.0, math.pi, 36, endpoint=False):
        for phi in np.linspace(0.0, math.pi, 18, endpoint=False):
            flat = (points @ _rotation(theta, phi).T)[:, :2]
            score = _edge_clearance(flat, edges)
            if score > best[0]:
                best = (score, flat)
    return carrier, best[1]


def registry_figure() -> Path:
    """The whole registry as three status-marked columns."""
    kinds = ("isomorphism", "disagreement", "complementarity")
    columns = {
        kind: sorted(
            (entry for entry in CORRESPONDENCES if entry.kind == kind),
            key=lambda entry: (entry.status, entry.key),
        )
        for kind in kinds
    }
    rows = max(len(entries) for entries in columns.values())
    figure, axes = plt.subplots(1, 3, figsize=(12.0, 5.8))
    for axis, kind in zip(axes, kinds):
        entries = columns[kind]
        for row, entry in enumerate(entries):
            axis.scatter(
                [0.0],
                [row],
                marker=STATUS_MARKERS.get(entry.status, "D"),
                s=110,
                color=KIND_COLOURS[kind],
                edgecolor=INK,
                linewidth=0.9,
                zorder=3,
            )
            axis.text(0.12, row, entry.key, va="center", fontsize=9.5, family="monospace")
        axis.set_xlim(-0.2, 1.5)
        axis.set_ylim(rows - 0.4, -0.8)
        axis.set_yticks([])
        axis.set_xticks([])
        axis.grid(False)
        for spine in axis.spines.values():
            spine.set_visible(False)
        axis.set_title(f"{kind} ({len(entries)})", color=KIND_COLOURS[kind], fontsize=11)
    present = [status for status in STATUS_MARKERS if any(e.status == status for e in CORRESPONDENCES)]
    handles = [
        plt.Line2D(
            [],
            [],
            marker=STATUS_MARKERS[status],
            linestyle="",
            markerfacecolor="white",
            markeredgecolor=INK,
            markeredgewidth=1.2,
            label=status,
            markersize=9,
        )
        for status in present
    ]
    axes[1].legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.01),
        ncol=len(handles),
        fontsize=10,
        title="marker = status; colour = kind",
        title_fontsize=9,
    )
    figure.suptitle(
        f"{len(CORRESPONDENCES)} correspondences: every row is executed by "
        "oph-fep-synergetics verify",
        fontsize=11,
    )
    kind_counts = {kind: len(columns[kind]) for kind in kinds}
    return _save(
        figure,
        "fig_registry.png",
        f"The {len(CORRESPONDENCES)} correspondences arranged in three columns by kind. "
        f"{_counted(kind_counts['isomorphism'], 'isomorphism', 'isomorphisms')} on the left, "
        f"{_counted(kind_counts['disagreement'], 'disagreement', 'disagreements')} in the middle, "
        f"{_counted(kind_counts['complementarity'], 'complementarity', 'complementarities')} on "
        "the right. Each row is labelled with its registry key and carries a marker for its "
        "status: a circle for exact, a square for exact-on-domain, and a cross for obstructed. "
        "Rows are aligned across the three columns.",
    )


def descent_figure() -> Path:
    """The descent census and a mismatch-one excursion."""
    bound = as_dict()
    census = controls.gated_descent_census()
    excursion = controls.mismatch_one_excursion(steps=6)
    energies = [0.5 * sum(value * value for value in state) for state in excursion]

    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.4))
    labels = ["strict descent", "zero descent", "no-op"]
    values = [census["strict"], census["zero"], census["noop"]]
    bars = left.bar(labels, values, color=[GREEN, ORANGE, GREY])
    for bar, value in zip(bars, values):
        left.annotate(
            str(value),
            (bar.get_x() + bar.get_width() / 2, bar.get_height()),
            textcoords="offset points",
            xytext=(0, 3),
            ha="center",
            fontsize=9,
        )
    left.set_ylabel("transitions")
    left.set_title(
        f"{census['states']} states × 30 seams × 2 coins = {sum(values)} transitions"
    )
    left.set_ylim(0, max(values) * 1.15)

    matrix = np.asarray(excursion, dtype=float)
    image = right.imshow(matrix.T, aspect="auto", cmap="Blues", origin="lower")
    right.set_yticks(range(matrix.shape[1]))
    right.set_yticklabels([str(index) for index in range(matrix.shape[1])], fontsize=8)
    right.set_xlabel("step")
    right.set_ylabel("port")
    right.grid(False)
    right.set_title(f"mismatch-one excursion, energy {energies[0]:.2f} throughout")
    figure.colorbar(image, ax=right, fraction=0.03, pad=0.02, label="load")
    return _save(
        figure,
        "fig_descent.png",
        f"Two panels. Left: the shell census, a bar chart of "
        f"{int(bound['shell.strict_descent'])} strict descents, "
        f"{int(bound['shell.zero_descent'])} zero descents and {int(bound['shell.noop'])} "
        f"no-op moves out of {int(bound['shell.transitions'])} transitions. Right: a heat map "
        f"of a seven-state excursion over the {int(bound['carrier.ports'])} ports, in which "
        "the load moves between ports while the quadratic energy stays at one half throughout.",
    )


def dstability_figure() -> Path:
    """Eigenvalue loci and the sign of the third Routh entry."""
    values = as_dict()
    samples = controls.D_STABILITY_SAMPLES
    onset = float(values["d_stability.routh_onset"])
    k4_third = float(controls.routh_first_column(4.0)[2])
    figure, (left, right) = plt.subplots(1, 2, figsize=(12.0, 4.8))

    spectra = {k: controls.d_stability_eigenvalues(k) for k in samples}
    reals = np.concatenate([spectra[k].real for k in samples])
    xmax = float(reals.max()) * 3.0
    left.axvspan(0.0, xmax, color=VERMILLION, alpha=0.10, linewidth=0)
    left.text(
        xmax * 0.5,
        0.0,
        "right half-plane: unstable",
        rotation=90,
        ha="center",
        va="center",
        fontsize=9,
        color=VERMILLION,
    )
    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(math.log(min(samples)), math.log(max(samples)) * 1.12)
    handles = []
    for k in samples:
        colour = cmap(norm(math.log(k)))
        eigenvalues = spectra[k]
        stable = eigenvalues.real <= 1e-9
        left.scatter(
            eigenvalues.real[stable],
            eigenvalues.imag[stable],
            s=60,
            color=colour,
            edgecolor=INK,
            linewidth=0.7,
            zorder=3,
        )
        left.scatter(
            eigenvalues.real[~stable],
            eigenvalues.imag[~stable],
            s=90,
            marker="X",
            color=colour,
            edgecolor=INK,
            linewidth=0.9,
            zorder=3,
        )
        handles.append(
            plt.Line2D([], [], marker="o", linestyle="", color=colour, markeredgecolor=INK,
                       markersize=8, label=f"k = {k:g}")
        )
    handles.append(
        plt.Line2D([], [], marker="X", linestyle="", color="white", markeredgecolor=INK,
                   markersize=9, label="X: Re > 0")
    )
    left.axvline(0.0, color=BLACK, linewidth=1.0, linestyle="--")
    left.set_xlim(reals.min() * 1.1, xmax)
    left.set_xlabel("real part")
    left.set_ylabel("imaginary part")
    left.set_title("eigenvalues of the scaled matrix\n(right half-plane = unstable)")
    left.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        title="coupling k",
        fontsize=9,
    )

    grid = np.linspace(0.05, 5.0, 400)
    third = [controls.routh_first_column(k)[2] for k in grid]
    right.plot(grid, third, color=BLUE, linewidth=2.0)
    right.axhline(0.0, color=BLACK, linewidth=0.9)
    right.axvline(onset, color=VERMILLION, linestyle="--", linewidth=1.1)
    right.annotate(
        f"onset at $k=\\sqrt{{10}}={onset:.4f}$",
        xy=(onset, 0.0),
        xytext=(onset + 0.3, 3.6),
        fontsize=9,
        color=VERMILLION,
        arrowprops={"arrowstyle": "->", "color": VERMILLION, "lw": 1.0},
    )
    right.scatter([4.0], [k4_third], color=VERMILLION, zorder=3, s=40)
    right.annotate(
        f"$k=4$ gives ${k4_third:g}$".replace("-", "−"),
        xy=(4.0, k4_third),
        xytext=(4.15, 0.9),
        fontsize=9,
        color=VERMILLION,
        arrowprops={"arrowstyle": "->", "color": VERMILLION, "lw": 1.0},
    )
    right.set_xlabel("coupling $k$")
    right.set_ylabel("third Routh entry")
    right.set_ylim(min(third) * 1.25, 8.5)
    right.set_title(r"$2(10-k^2)/(k+2)$: negative means unstable")
    return _save(
        figure,
        "fig_dstability.png",
        "Two panels. Left: the eigenvalues of the scaled matrix in the complex plane for "
        f"{len(values['d_stability.samples'])} couplings, coloured by coupling from dark to light "
        "with a legend outside the plot, a dashed imaginary axis, and the right half-plane shaded "
        "as unstable; the complex pair that crosses into it is drawn with crosses. No lines join "
        "eigenvalues across couplings. "
        "Right: the third Routh entry two times ten minus k squared over k plus two against the "
        f"coupling, crossing zero at {format_value(values['d_stability.routh_onset'])} and equal "
        f"to {format_value(values['d_stability.routh_k4_third'])} at k equal to four.",
    )


def mediator_figure() -> Path:
    """The pairwise margin and the reduction error."""
    bound = as_dict()
    threshold = float(bound["mediator.threshold"])
    below = float(bound["mediator.epsilon_below"])
    margin_below = float(bound["mediator.margin_below"])
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.4))
    grid = np.linspace(0.0, 1.0, 200)
    margins = [controls.mediator_pairwise_margin(value) for value in grid]
    left.plot(grid, margins, color=BLUE, linewidth=2.0)
    left.axhline(0.0, color=BLACK, linewidth=0.9)
    left.axvline(threshold, color=VERMILLION, linestyle="--", linewidth=1.0)
    left.scatter([below, threshold], [margin_below, 0.0], color=VERMILLION, zorder=3)
    fraction = Fraction(threshold).limit_denominator(100)
    left.annotate(
        f"threshold $\\varepsilon={fraction}$",
        xy=(threshold, 0.0),
        xytext=(threshold + 0.1, -0.75),
        fontsize=9,
        color=VERMILLION,
        arrowprops={"arrowstyle": "->", "color": VERMILLION, "lw": 1.0},
    )
    left.annotate(
        f"$\\varepsilon={below:g}$:\nmargin ${margin_below:.1f}$".replace("-", "−"),
        xy=(below, margin_below),
        xytext=(0.02, 0.9),
        fontsize=9,
        color=VERMILLION,
        arrowprops={"arrowstyle": "->", "color": VERMILLION, "lw": 1.0},
    )
    left.set_xlabel(r"$\varepsilon$")
    left.set_ylabel("pairwise margin")
    left.set_title(r"margin $=3\varepsilon-1$")

    reduction = controls.reduction_control()
    names = ["same observation", "first after", "second after", "response distance", "worst-case floor"]
    keys = ["same_observation", "first_after", "second_after", "response_distance", "worst_case_error_lower_bound"]
    values = [float(reduction[key]) for key in keys]
    bars = right.barh(names, values, color=[BLUE, GREEN, VERMILLION, BLUE, ORANGE])
    for bar, value in zip(bars, values):
        right.annotate(
            f"{value:g}".replace("-", "−"),
            (value, bar.get_y() + bar.get_height() / 2),
            textcoords="offset points",
            xytext=(5 if value >= 0 else -5, 0),
            ha="left" if value >= 0 else "right",
            va="center",
            fontsize=9,
        )
    right.axvline(0.0, color=BLACK, linewidth=0.9)
    right.invert_yaxis()
    right.set_xlim(min(values) * 1.5, max(values) * 1.3)
    right.set_title("two histories, one observation")
    return _save(
        figure,
        "fig_mediator.png",
        "Two panels. Left: the pairwise margin three epsilon minus one against epsilon, crossing "
        f"zero at the threshold {fraction}, with the point at epsilon equal to {below:g} marked at "
        f"{margin_below:.1f}. Right: a horizontal bar chart of the reduction control, showing that "
        "the same observation scores one, the two criteria score one and minus one after the "
        "reduction, the response distance is two, and the worst-case error floor is one.",
    )


def blanket_figure() -> Path:
    """The drift that a stationary law does not fix, and free entropy production."""
    bound = as_dict()
    drift = controls.blanket_dynamics_check(0.75)
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.6))
    names = ["exterior", "blanket", "internal"]
    values = [
        float(drift["mu_drift_external"]),
        float(drift["mu_drift_blanket"]),
        float(drift["mu_drift_internal"]),
    ]
    bars = left.bar(names, values, color=[BLUE, ORANGE, VERMILLION])
    for bar, value in zip(bars, values):
        left.annotate(
            f"{value:g}".replace("-", "−"),
            (bar.get_x() + bar.get_width() / 2, value),
            textcoords="offset points",
            xytext=(0, 4 if value >= 0 else -4),
            ha="center",
            va="bottom" if value >= 0 else "top",
            fontsize=9,
        )
    left.axhline(0.0, color=BLACK, linewidth=0.9)
    span = max(values) - min(values)
    left.set_ylim(min(values) - 0.15 * span, max(values) + 0.15 * span)
    left.set_ylabel(f"drift component at $q={bound['blanket.q']:g}$")
    left.set_title(
        "stationary law factorises;\n"
        f"cross drift difference {float(drift['cross_drift_difference']):g}"
    )

    grid = np.linspace(0.0, 1.5, 200)
    production = [controls.solenoidal_check(value)["entropy_production"] for value in grid]
    right.plot(grid, production, color=BLUE, linewidth=2.0)
    marks = [float(q) for q in bound["solenoidal.q_grid"]]
    marked = [controls.solenoidal_check(value)["entropy_production"] for value in marks]
    right.scatter(marks, marked, color=VERMILLION, zorder=3, s=40)
    offsets = [((8, 6), "left"), ((-8, 8), "right"), ((-8, 6), "right")]
    for value, height, (offset, align) in zip(marks, marked, offsets):
        right.annotate(
            f"{height:.3g}",
            (value, height),
            textcoords="offset points",
            xytext=offset,
            ha=align,
            fontsize=9,
            color=VERMILLION,
        )
    right.set_xlim(-0.08, 1.65)
    right.set_xlabel("solenoidal parameter $q$")
    right.set_ylabel("entropy production")
    right.set_title(r"$2q^2$, free while every member is stationary")
    return _save(
        figure,
        "fig_blanket.png",
        f"Two panels. Left: the three drift components at q equal to "
        f"{format_value(bound['blanket.q'])}, namely {_enumerated(bound['blanket.drift'])}, "
        f"with the cross drift difference of {format_value(bound['blanket.cross_drift'])} "
        "annotated; the stationary law factorises even though the drift couples the exterior "
        "and the interior. Right: entropy production against the solenoidal parameter q, "
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

    figure = plt.figure(figsize=(11.5, 5.4))
    grid = figure.add_gridspec(2, 2, width_ratios=(1.0, 1.4))
    left = figure.add_subplot(grid[:, 0])
    top = figure.add_subplot(grid[0, 1])
    bottom = figure.add_subplot(grid[1, 1])

    left.plot(samples[:, 0], samples[:, 1], color=BLUE, linewidth=1.3)
    theta = np.linspace(0.0, 2.0 * math.pi, 400)
    for bound, style in ((r_min, "-"), (r_max, "--")):
        left.plot(
            bound * np.cos(theta),
            bound * np.sin(theta),
            color=VERMILLION,
            linestyle=style,
            linewidth=1.1,
        )
    left.scatter([0.0], [0.0], color=BLACK, s=20, zorder=3)
    angle = math.radians(-60.0)
    left.annotate(
        f"$r_{{min}}={r_min:.4f}$",
        xy=(r_min * math.cos(angle), r_min * math.sin(angle)),
        xytext=(0.0, -0.3),
        ha="center",
        fontsize=9,
        color=VERMILLION,
        arrowprops={"arrowstyle": "->", "color": VERMILLION, "lw": 1.0},
    )
    angle = math.radians(58.0)
    left.annotate(
        f"$r_{{max}}={r_max:.4f}$ (dashed)",
        xy=(r_max * math.cos(angle), r_max * math.sin(angle)),
        xytext=(0.15, 1.38),
        ha="left",
        fontsize=9,
        color=VERMILLION,
        arrowprops={"arrowstyle": "->", "color": VERMILLION, "lw": 1.0},
    )
    left.set_aspect("equal")
    left.set_xlim(-1.3, 1.9)
    left.set_ylim(-1.3, 1.6)
    left.set_xlabel("x")
    left.set_ylabel("y")
    left.set_title(f"relaxation orbit in the annulus $r_{{min}}\\leq r\\leq r_{{max}}$")

    times = np.arange(radius.size) * sample_dt
    keep = times <= window
    top.plot(times[keep], radius[keep], color=BLUE, linewidth=1.5)
    top.axhline(r_min, color=VERMILLION, linestyle="-", linewidth=1.0)
    top.axhline(r_max, color=VERMILLION, linestyle="--", linewidth=1.0)
    top.text(window * 1.01, r_max, "$r_{max}$", color=VERMILLION, fontsize=9, va="center")
    top.text(window * 1.01, r_min, "$r_{min}$", color=VERMILLION, fontsize=9, va="center")
    arrow_y = r_min - 0.05
    top.annotate(
        "",
        xy=(period, arrow_y),
        xytext=(0.0, arrow_y),
        arrowprops={"arrowstyle": "<->", "color": INK, "lw": 1.1},
    )
    top.text(0.2, arrow_y - 0.012, f"one period = {period:.3f}", ha="left", va="top",
             fontsize=9, color=INK)
    top.set_xlim(0.0, window)
    top.set_ylim(r_min - 0.17, r_max + 0.05)
    top.set_xlabel("relaxation time $t$ (continuous flow)")
    top.set_ylabel("radius")
    top.set_title("relaxation: radius over five periods")

    steps = np.arange(sync_radius.size)
    show = steps <= window
    bottom.plot(steps[show], sync_radius[show], color=ORANGE, marker="o", markersize=4.5,
                linewidth=1.0, markeredgecolor=INK, markeredgewidth=0.5)
    bottom.axhline(sync_floor, color=VERMILLION, linestyle="--", linewidth=1.0)
    bottom.text(window * 1.01, sync_floor, f"floor {sync_floor:g}", color=VERMILLION, fontsize=9,
                va="center")
    bottom.set_xlim(0.0, window)
    bottom.set_xlabel("iteration $n$ (synchronous clipped map, discrete steps)")
    bottom.set_ylabel("radius")
    bottom.set_title("synchronous clipped map: a different object, a different band")
    return _save(
        figure,
        "fig_orbit.png",
        "Three panels. Left: the relaxation orbit in the plane, a closed curve that neither "
        "escapes nor reaches the origin, with the annulus that contains it drawn as two circles "
        "labelled r minimum and r maximum. Top right: the radius of the continuous relaxation "
        "against time over five periods, oscillating between the two radii, with one period "
        "marked. Bottom right: the radius of the synchronous clipped map against its discrete "
        "iteration number, confined to a different band with its floor marked.",
    )


def ivm_figure() -> Path:
    """IVM shell sizes and tetravolume ratios."""
    shells = list(range(0, 9))
    counts = [lattice.shell_count(index) for index in shells]
    cumulative = [lattice.cumulative_count(index) for index in shells]

    figure, (left, right) = plt.subplots(1, 2, figsize=(11.5, 4.6))
    bars = left.bar(shells, counts, color=BLUE, width=0.7)
    for bar, value in zip(bars, counts):
        left.annotate(
            str(value),
            (bar.get_x() + bar.get_width() / 2, value),
            textcoords="offset points",
            xytext=(0, 3),
            ha="center",
            fontsize=9,
        )
    left.set_xticks(shells)
    left.set_xticklabels([f"{k}\n{total}" for k, total in zip(shells, cumulative)])
    left.set_xlabel("shell $k$   (second row: cumulative sites)")
    left.set_ylabel("sites in shell")
    left.set_title(r"$10k^2+2$ sites in shell $k\geq 1$")
    left.set_ylim(0, max(counts) * 1.12)

    names = list(lattice.VOLUME_RATIOS)
    volumes = [float(lattice.VOLUME_RATIOS[name]) for name in names]
    pretty = [name.replace("_", " ") for name in names]
    order = np.argsort(volumes)
    bars = right.barh([pretty[index] for index in order], [volumes[index] for index in order],
                      color=BLUE)
    right.set_xscale("log")
    for bar, index in zip(bars, order):
        right.annotate(
            f"{volumes[index]:g}",
            (volumes[index], bar.get_y() + bar.get_height() / 2),
            textcoords="offset points",
            xytext=(5, 0),
            va="center",
            fontsize=9,
        )
    ticks = sorted(set(volumes))
    right.set_xticks(ticks)
    right.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    right.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    right.set_xlabel("tetravolume, unit tetrahedron = 1 (log scale)")
    right.set_title("exact integer volumes")
    right.set_xlim(0.8, max(volumes) * 1.8)
    return _save(
        figure,
        "fig_ivm.png",
        "Two panels. Left: the site count of each IVM shell, ten k squared plus two for shells "
        "one onward, as bars for shells zero to eight, with the cumulative count printed under "
        "each shell number; the first shell holds the twelve neighbours. Right: tetravolume "
        "ratios on a logarithmic scale, with the tetrahedron at one, the cube at three, the "
        "octahedron at four, the rhombic dodecahedron at six, the cuboctahedron at twenty and "
        "the truncated octahedron at ninety-six.",
    )


def subdivision_figure() -> Path:
    """The same quadratic sequence in three places."""
    frequencies = list(range(1, 7))
    geodesic_vertices = [lattice.geodesic_vertex_count(value) for value in frequencies]
    geodesic_faces = [lattice.geodesic_face_count(value) for value in frequencies]
    goldberg_vertices = [lattice.goldberg_vertex_count(value) for value in frequencies]
    goldberg_faces = [lattice.goldberg_face_count(value) for value in frequencies]
    shell = [lattice.shell_count(value) for value in frequencies]

    figure, axis = plt.subplots(figsize=(9.0, 5.2))
    axis.plot(frequencies, geodesic_vertices, color=BLUE, linewidth=1.2, zorder=1)
    axis.plot(frequencies, geodesic_faces, color=VERMILLION, linewidth=1.2, zorder=1)
    # Coincident series are drawn as nested open markers of decreasing size so that
    # each stays visible under the larger one.
    nested = (
        (geodesic_vertices, "o", 15, BLUE, "geodesic vertices"),
        (goldberg_faces, "s", 10, GREEN, "Goldberg faces"),
        (shell, "D", 5.5, ORANGE, "IVM shell count"),
        (geodesic_faces, "^", 15, VERMILLION, "geodesic faces"),
        (goldberg_vertices, "v", 8, PURPLE, "Goldberg vertices"),
    )
    for series, marker, size, colour, label in nested:
        axis.plot(
            frequencies,
            series,
            linestyle="",
            marker=marker,
            markersize=size,
            markerfacecolor="none",
            markeredgecolor=colour,
            markeredgewidth=1.8,
            label=label,
            zorder=3,
        )
    for value in frequencies:
        axis.annotate(
            f"{lattice.geodesic_face_count(value)}",
            (value, lattice.geodesic_face_count(value)),
            textcoords="offset points",
            xytext=(0, 12),
            ha="center",
            fontsize=9,
            color=VERMILLION,
        )
        axis.annotate(
            f"{lattice.geodesic_vertex_count(value)}",
            (value, lattice.geodesic_vertex_count(value)),
            textcoords="offset points",
            xytext=(11, -15),
            ha="left",
            fontsize=9,
            color=BLUE,
        )
    axis.set_xticks(frequencies)
    axis.set_xlim(0.7, 6.6)
    axis.set_ylim(-40, max(geodesic_faces) * 1.1)
    axis.set_xlabel("frequency $f$")
    axis.set_ylabel("count")
    axis.set_title(r"$10f^2+2$ and $20f^2$: one duality, three appearances")
    axis.legend(
        loc="upper left",
        title="series that coincide are drawn nested,\nso each marker is visible",
        title_fontsize=9,
        fontsize=9,
    )
    return _save(
        figure,
        "fig_subdivision.png",
        "A line chart against subdivision frequency from one to six. Geodesic vertex count, "
        "Goldberg face count and the IVM shell count coincide on ten f squared plus two, drawn "
        "as nested open markers so that all three stay visible, while the geodesic face count "
        "and the Goldberg vertex count coincide on twenty f squared, also drawn nested; the "
        "two curves are the two sides of one duality, and the count values are annotated for "
        "the geodesic vertices and the geodesic faces.",
    )


def cover_figure() -> Path:
    """Cover art: the carrier's ports and seams on a dark field."""
    carrier, flat = _projected_carrier()
    figure = plt.figure(figsize=(7.0, 7.0))
    figure.patch.set_facecolor("#0b0f14")
    axis = figure.add_subplot(111)
    axis.set_facecolor("#0b0f14")
    axis.set_axis_off()
    for face in OPH_ORIENTED_BASE_FACES:
        axis.fill(
            [flat[port, 0] for port in face],
            [flat[port, 1] for port in face],
            color="#1d4e56",
            alpha=0.35,
            linewidth=0,
            zorder=0,
        )
    for port in range(carrier.port_count):
        for other in carrier.neighbours(port):
            if other > port:
                axis.plot(
                    [flat[port, 0], flat[other, 0]],
                    [flat[port, 1], flat[other, 1]],
                    color="#5c8f96",
                    linewidth=1.2,
                    alpha=0.9,
                    zorder=1,
                )
    axis.scatter(
        flat[:, 0],
        flat[:, 1],
        s=150,
        color="#0b0f14",
        edgecolor="#e6c27a",
        linewidth=1.3,
        zorder=3,
    )
    for port in range(carrier.port_count):
        axis.text(
            flat[port, 0],
            flat[port, 1],
            str(port),
            ha="center",
            va="center",
            fontsize=7.5,
            color="#e6c27a",
            zorder=4,
        )
    axis.set_aspect("equal")
    margin = 0.12
    axis.set_xlim(flat[:, 0].min() - margin, flat[:, 0].max() + margin)
    axis.set_ylim(flat[:, 1].min() - margin, flat[:, 1].max() + margin)
    return _save(
        figure,
        "fig_cover.png",
        "The twelve ports of the carrier as pale gold nodes on a dark field, joined by its "
        "thirty seams in teal and arranged in the planar spectral projection of the "
        "icosahedron; each port is labelled with its index and the wireframe carries no title.",
        facecolor="#0b0f14",
    )


def _kind_phrases(counts: Counter) -> list[str]:
    """The non-zero kind counts of a ledger, one phrase each."""
    names = (
        ("isomorphism", "isomorphism", "isomorphisms"),
        ("disagreement", "disagreement", "disagreements"),
        ("complementarity", "complementarity", "complementarities"),
    )
    return [_counted(counts[kind], one, many) for kind, one, many in names if counts[kind]]


def ga_triangle_figure() -> Path:
    """The three frameworks as a triangle, with each edge's ledger and its key claims."""
    figure, axis = plt.subplots(figsize=(10.0, 6.0))
    axis.set_axis_off()
    axis.set_aspect("equal")
    vertices = {
        "OPH": (0.50, 0.86),
        "FEP": (0.10, 0.12),
        "Synergetics": (0.90, 0.12),
    }
    expansions = {
        "OPH": "Observer Patch Holography",
        "FEP": "Free Energy Principle",
        "Synergetics": "quadray / IVM",
    }
    for left, right in (("OPH", "FEP"), ("FEP", "Synergetics"), ("Synergetics", "OPH")):
        axis.plot(
            [vertices[left][0], vertices[right][0]],
            [vertices[left][1], vertices[right][1]],
            color=INK,
            linewidth=2.0,
            zorder=1,
        )
    for name, (x_value, y_value) in vertices.items():
        axis.scatter([x_value], [y_value], s=100, color=INK, zorder=3)
    # Vertex labels sit outside the markers: beside the apex, and beyond the base corners.
    axis.text(0.56, 0.875, "OPH", ha="left", va="center", fontsize=13, fontweight="bold", color=INK)
    axis.text(0.56, 0.815, expansions["OPH"], ha="left", va="center", fontsize=10, color=INK)
    axis.text(0.05, 0.145, "FEP", ha="right", va="center", fontsize=13, fontweight="bold", color=INK)
    axis.text(0.05, 0.095, expansions["FEP"], ha="right", va="center", fontsize=10, color=INK)
    axis.text(0.95, 0.145, "Synergetics", ha="left", va="center", fontsize=13,
              fontweight="bold", color=INK)
    axis.text(0.95, 0.095, expansions["Synergetics"], ha="left", va="center", fontsize=10,
              color=INK)

    grouped = by_pair()
    named = [
        key
        for key in ("carrier-port-chart", "native-repair-law", "free-energy-dictionary")
        if key in by_key()
    ]
    highlighted = named if len(named) == 3 else []
    if not highlighted:
        highlighted = [
            next(entry.key for entry in grouped[pair] if entry.kind == "isomorphism")
            for pair in (OPH_FEP, FEP_SYNERGETICS, SYNERGETICS_OPH)
        ]

    line = 0.036  # data units per 10 pt text line at this scale
    placements = (
        (OPH_FEP, (0.18, 0.72), "right"),
        (FEP_SYNERGETICS, (0.50, 0.045), "center"),
        (SYNERGETICS_OPH, (0.82, 0.72), "left"),
    )
    for pair, (x_value, y_value), align in placements:
        counts = Counter(entry.kind for entry in grouped[pair])
        phrases = _kind_phrases(counts)
        axis.text(
            x_value,
            y_value,
            "\n".join(phrases),
            ha=align,
            va="top",
            fontsize=10.5,
            color=INK,
            linespacing=1.25,
            zorder=5,
        )
        keys = [key for key in highlighted if by_key()[key].pair == pair]
        if keys:
            axis.text(
                x_value,
                y_value - line * (len(phrases) + 0.6),
                "load-bearing:\n" + "\n".join(keys),
                ha=align,
                va="top",
                fontsize=9.5,
                color=BLUE,
                linespacing=1.25,
                zorder=5,
            )

    three = grouped.get(THREE_WAY, ())
    if three:
        three_counts = Counter(entry.kind for entry in three)
        three_text = (
            f"three-way: {len(three)}\n"
            + "\n".join(_kind_phrases(three_counts))
            + "\n"
            + "\n".join(entry.key for entry in three)
        )
    else:
        three_text = "three-way: no registered entry"
    axis.text(
        0.50,
        0.37,
        three_text,
        ha="center",
        va="center",
        fontsize=9.5,
        color=INK,
        linespacing=1.3,
        bbox={"facecolor": "#F1F5F9", "edgecolor": SKY, "boxstyle": "round,pad=0.35"},
        zorder=2,
    )
    axis.set_xlim(-0.45, 1.45)
    axis.set_ylim(-0.24, 0.98)
    return _save(
        figure,
        "fig_ga_triangle.png",
        "A large triangle whose vertices are OPH, Observer Patch Holography; FEP, Free Energy "
        "Principle; and Synergetics, quadray/IVM, each named beside its marker. Each edge carries "
        "its pairwise ledger of isomorphisms, disagreements and complementarities, with zero "
        "counts left out, and the three load-bearing identifications are named beside the edge "
        "they belong to. The centre gives the three-way count and its two keys, "
        "tetravolume-arithmetic and conjugation-norm.",
    )


def ga_carrier_figure() -> Path:
    """The carrier graph beside the OPH-facing and FEP-facing invariant ledger."""
    bound = as_dict()
    carrier, flat = _projected_carrier()
    figure, (left, right) = plt.subplots(
        1, 2, figsize=(11.0, 5.6), gridspec_kw={"width_ratios": (1.0, 0.9)}
    )
    for port in range(carrier.port_count):
        for other in carrier.neighbours(port):
            if other > port:
                left.plot(
                    [flat[port, 0], flat[other, 0]],
                    [flat[port, 1], flat[other, 1]],
                    color=INK,
                    linewidth=1.1,
                    zorder=2,
                )
    left.scatter(
        flat[:, 0],
        flat[:, 1],
        s=210,
        color="white",
        edgecolor=INK,
        linewidth=1.1,
        zorder=3,
    )
    for port in range(carrier.port_count):
        left.text(
            flat[port, 0],
            flat[port, 1],
            str(port),
            ha="center",
            va="center",
            fontsize=8.5,
            color=INK,
            zorder=4,
        )
    left.set_aspect("equal")
    left.set_axis_off()
    left.set_title(f"{carrier.port_count} ports, {carrier.seam_count} seams")

    cohomology = sheaf.cohomology_dimensions(sheaf.carrier_sheaf(carrier))
    group = symmetry.carrier_automorphisms(carrier)
    rotations = symmetry.rotation_subgroup(carrier)
    reflections = symmetry.reflection_subgroup(carrier)
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

    oph_rows = [
        ("spanning trees", f"{tree_number:,}"),
        (r"$\dim H^0$", f"{cohomology['h0']}"),
        (r"$\dim H^1$", f"{cohomology['h1']}"),
        ("symmetry group order", f"{len(group)}"),
        ("  rotations", f"{len(rotations)}"),
        ("  reflections", f"{len(reflections)}"),
    ]
    fep_rows = [
        (r"port variance $(L^{+})_{ii}$", f"{marginal.numerator}/{marginal.denominator}"),
        ("blanket size", blanket_text),
        ("pairwise precision", f"{precision_text}".replace("-", "\u2212")),
        ("  on pairs", f"{len(precision)}"),
    ]

    right.set_axis_off()
    right.set_xlim(0.0, 1.0)
    right.set_ylim(0.0, 1.0)
    right.text(0.0, 0.96, "OPH-facing", fontsize=12, color=BLUE, fontweight="bold", va="center")
    right.text(0.0, 0.40, "FEP-facing", fontsize=12, color=VERMILLION, fontweight="bold",
               va="center")
    for index, (label, value) in enumerate(oph_rows):
        y_value = 0.87 - 0.075 * index
        right.text(0.02, y_value, label, fontsize=10.5, color=INK, va="center")
        right.text(0.98, y_value, value, fontsize=10.5, color=BLUE, va="center", ha="right")
    for index, (label, value) in enumerate(fep_rows):
        y_value = 0.31 - 0.075 * index
        right.text(0.02, y_value, label, fontsize=10.5, color=INK, va="center")
        right.text(0.98, y_value, value, fontsize=10.5, color=VERMILLION, va="center", ha="right")
    return _save(
        figure,
        "fig_ga_carrier.png",
        f"Two panels. Left: the {int(bound['carrier.ports'])} ports and "
        f"{int(bound['carrier.seams'])} seams of the carrier in its planar "
        "spectral projection, with each port labelled by index and no two edges nearly "
        "coincident. Right: a two-part ledger of invariants, OPH-facing above and FEP-facing "
        "below, listing the spanning-tree count, the dimensions of H0 and H1, the automorphism "
        "group order with its rotation and reflection split, the port variance (a diagonal "
        "entry of the Laplacian pseudoinverse) as a fraction, the blanket size and the "
        "pairwise precision.",
    )


def ga_descent_figure() -> Path:
    """The descent census and the free-energy levels of the first five IVM shells."""
    bound = as_dict()
    census = controls.gated_descent_census()
    multiplicities = controls.path_multiplicity_values()
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.5, 5.4))
    labels = ["strict descent", "zero descent", "no-op"]
    values = [census["strict"], census["zero"], census["noop"]]
    bars = left.bar(labels, values, color=[GREEN, ORANGE, GREY])
    for bar, value in zip(bars, values):
        left.annotate(
            str(value),
            (bar.get_x() + bar.get_width() / 2, bar.get_height()),
            textcoords="offset points",
            xytext=(0, 3),
            ha="center",
            fontsize=9,
        )
    left.set_ylabel("transitions")
    left.set_title(
        f"{census['states']} states × 30 seams × 2 coins\n= {sum(values)} transitions"
    )
    left.set_ylim(0, max(values) * 1.25)
    left.text(
        0.03,
        0.97,
        "path multiplicities: " + ", ".join(str(value) for value in multiplicities),
        transform=left.transAxes,
        va="top",
        fontsize=9,
        color=INK,
    )

    k_values = list(range(1, 6))
    origin = np.zeros(3)
    identity = np.eye(3)
    shared: dict[float, list[int]] = {}
    kl_by_k: dict[int, float] = {}
    for k in k_values:
        shell = geometry.generate_shell(k)
        levels = sorted({float(0.5 * geometry.squared_norm(site)) for site in shell})
        right.scatter([k] * len(levels), levels, color=BLUE, s=30, zorder=3)
        for level in levels:
            shared.setdefault(level, []).append(k)
        target = geometry.NEIGHBOR_SQUARED_DISTANCE * k * k
        site = next(site for site in shell if geometry.squared_norm(site) == target)
        kl_by_k[k] = free_energy.gaussian_kl(origin, identity, geometry.to_xyz(site), identity)
    right.plot(
        k_values,
        [kl_by_k[k] for k in k_values],
        color=ORANGE,
        linestyle="--",
        marker="o",
        linewidth=1.4,
        label=r"Gaussian KL $=4k^2$",
        zorder=2,
    )
    for level, shells in shared.items():
        if len(shells) > 1:
            right.axhline(level, color=VERMILLION, linestyle=":", linewidth=1.1)
            right.annotate(
                f"level {level:g}: shells {', '.join(str(k) for k in shells)}",
                xy=(max(shells), level),
                xytext=(5.12, level - 1.0),
                va="top",
                fontsize=9,
                color=VERMILLION,
            )
    right.set_xticks(k_values)
    right.set_xlim(0.6, 7.6)
    right.set_xlabel("shell $k$")
    right.set_ylabel(r"free energy $F=\frac{1}{2}\|x\|^2$")
    right.set_title(
        f"levels {_enumerated(bound['shell.shared_free_energy_levels'])} "
        "are shared between shells"
    )
    right.legend(loc="upper left")
    return _save(
        figure,
        "fig_ga_descent.png",
        f"Two panels. Left: the descent census as three bars, "
        f"{int(bound['shell.strict_descent'])} strict descents, "
        f"{int(bound['shell.zero_descent'])} zero descents and {int(bound['shell.noop'])} "
        f"no-op moves, annotated with the path-multiplicity values "
        f"{_enumerated(bound['shell.multiplicities'])}. Right: the distinct free-energy levels "
        f"of the first {len(bound['ivm.shell_counts'])} IVM shells at integer shell numbers, with "
        f"levels {_enumerated(bound['shell.shared_free_energy_levels'])} marked by dotted lines "
        "where they are shared between shells, and the Gaussian KL curve annotated as four k "
        "squared.",
    )


def ga_synergetics_figure() -> Path:
    """The integral tetravolume table beside the icosahedron and its Voronoi cell."""
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.5, 5.2))
    names = list(geometry.VOLUME_RATIOS)
    volumes = [float(geometry.VOLUME_RATIOS[name]) for name in names]
    icosahedron = geometry.icosahedron_tetravolume()
    icosahedron_volume = float(icosahedron["tetravolume"])
    bar_names = [name.replace("_", " ") for name in names] + ["icosahedron"]
    bar_values = volumes + [icosahedron_volume]
    integral = [value == int(value) for value in bar_values]
    bars = left.barh(
        bar_names,
        bar_values,
        color=[BLUE if flag else VERMILLION for flag in integral],
        edgecolor=INK,
    )
    for bar, flag in zip(bars, integral):
        if not flag:
            bar.set_hatch("//")
    left.set_xscale("log")
    for bar, value, flag in zip(bars, bar_values, integral):
        left.annotate(
            f"{value:g}" if flag else f"{value:.4f}  (not an integer)",
            (value, bar.get_y() + bar.get_height() / 2),
            textcoords="offset points",
            xytext=(5, 0),
            va="center",
            fontsize=9,
            color=INK if flag else VERMILLION,
        )
    left.invert_yaxis()
    left.set_xlim(0.8, max(bar_values) * 12.0)
    left.set_xlabel("tetravolume, unit tetrahedron = 1 (log scale)")
    left.set_title("integral table solids; the icosahedron is not")
    left.legend(
        handles=[
            matplotlib.patches.Patch(facecolor=BLUE, edgecolor=INK, label="integer tetravolume"),
            matplotlib.patches.Patch(facecolor=VERMILLION, edgecolor=INK, hatch="//",
                                     label="non-integer tetravolume"),
        ],
        loc="upper right",
        fontsize=9,
    )

    cell = geometry.voronoi_cell()
    density = geometry.packing_density()
    cell_labels = ["vertices", "faces", "edges"]
    cell_values = [
        int(cell["vertex_count"]),
        int(cell["face_count"]),
        int(cell["edge_count"]),
    ]
    bars = right.bar(cell_labels, cell_values, color=[BLUE, GREEN, ORANGE])
    for bar, value in zip(bars, cell_values):
        right.annotate(
            str(value),
            (bar.get_x() + bar.get_width() / 2, bar.get_height()),
            textcoords="offset points",
            xytext=(0, 3),
            ha="center",
            fontsize=9,
        )
    right.set_ylabel("count")
    right.set_ylim(0, max(cell_values) * 1.6)
    right.set_title("Voronoi cell of the IVM")
    right.text(
        0.03,
        0.97,
        f"packing density {float(density['density']):.6f}\n"
        f"closed form        {float(density['closed_form']):.6f}",
        transform=right.transAxes,
        va="top",
        fontsize=9.5,
        color=INK,
    )
    return _save(
        figure,
        "fig_ga_synergetics.png",
        "Two panels. Left: logarithmic bars of the tetravolumes in the synergetic table, the "
        "tetrahedron, cube, octahedron, rhombic dodecahedron, cuboctahedron and truncated "
        "octahedron, all integral and drawn solid, beside the icosahedron's measured "
        "tetravolume, which is not integral and is drawn hatched and labelled as such. "
        "Right: the Voronoi cell's fourteen vertices, twelve faces and twenty-four "
        "edges, with the measured packing density annotated against its closed form.",
    )


def pair_triangles_figure() -> Path:
    """One small triangle for each pairwise scope, with its kind and status ledger."""
    figure, axes = plt.subplots(1, 3, figsize=(12.5, 4.6))
    pairs = (OPH_FEP, FEP_SYNERGETICS, SYNERGETICS_OPH)
    pair_edges = {
        OPH_FEP: ("OPH", "FEP"),
        FEP_SYNERGETICS: ("FEP", "Synergetics"),
        SYNERGETICS_OPH: ("Synergetics", "OPH"),
    }
    vertices = {
        "OPH": (0.50, 0.90),
        "FEP": (0.06, 0.16),
        "Synergetics": (0.94, 0.16),
    }
    grouped = by_pair()
    for axis, pair in zip(axes, pairs):
        entries = grouped[pair]
        counts = Counter(entry.kind for entry in entries)
        statuses = Counter(entry.status for entry in entries)
        for left_name, right_name in (
            ("OPH", "FEP"),
            ("FEP", "Synergetics"),
            ("Synergetics", "OPH"),
        ):
            highlighted = pair_edges[pair] == (left_name, right_name)
            axis.plot(
                [vertices[left_name][0], vertices[right_name][0]],
                [vertices[left_name][1], vertices[right_name][1]],
                color=ORANGE if highlighted else "#B8C2CC",
                linewidth=4.0 if highlighted else 1.2,
                zorder=1,
            )
        for name, (x_value, y_value) in vertices.items():
            axis.scatter([x_value], [y_value], s=80, color=INK, zorder=3)
            if name == "OPH":
                axis.text(x_value, y_value + 0.045, name, ha="center", va="bottom",
                          fontsize=11, fontweight="bold", color=INK)
            else:
                axis.text(x_value, y_value - 0.05, name, ha="center", va="top",
                          fontsize=11, fontweight="bold", color=INK)
        axis.text(
            0.50,
            0.40,
            "\n".join(_kind_phrases(counts)),
            ha="center",
            va="center",
            fontsize=10,
            color=INK,
            linespacing=1.4,
            zorder=2,
        )
        status_phrases = [f"{statuses[name]} {name}" for name in sorted(statuses) if statuses[name]]
        axis.text(
            0.50,
            -0.07,
            "status: " + ", ".join(status_phrases),
            ha="center",
            va="top",
            fontsize=9,
            color=INK,
            wrap=True,
        )
        axis.set_xlim(-0.12, 1.12)
        axis.set_ylim(-0.18, 1.08)
        axis.set_aspect("equal")
        axis.set_axis_off()
        first, second = pair_edges[pair]
        axis.set_title(f"{first} / {second}", fontsize=11)
    three = grouped.get(THREE_WAY, ())
    if three:
        figure.text(
            0.50,
            0.0,
            f"three-way: {len(three)}: " + ", ".join(entry.key for entry in three),
            ha="center",
            va="top",
            fontsize=10,
            color=INK,
        )
    return _save(
        figure,
        "fig_pair_triangles.png",
        "Three small triangles side by side, one for each pairwise scope: OPH and FEP, FEP and "
        "Synergetics, and Synergetics and OPH. Each triangle names all three frameworks outside "
        "its vertex markers and annotates its pair's non-zero counts by kind, with the status "
        "breakdown beneath and the pair's edge highlighted. A centre note records the three-way "
        "count and its two entries, tetravolume-arithmetic and conjugation-norm.",
    )


def _contrast(foreground: str, background: str) -> float:
    """The WCAG contrast ratio between two hex colours."""

    def luminance(hex_colour: str) -> float:
        value = hex_colour.lstrip("#")
        channels = []
        for index in (0, 2, 4):
            component = int(value[index : index + 2], 16) / 255.0
            channels.append(
                component / 12.92 if component <= 0.03928 else ((component + 0.055) / 1.055) ** 2.4
            )
        return 0.2126 * channels[0] + 0.7152 * channels[1] + 0.0722 * channels[2]

    high, low = sorted((luminance(foreground), luminance(background)), reverse=True)
    return (high + 0.05) / (low + 0.05)


def scope_matrix_figure() -> Path:
    """The registry as a dense categorical matrix in registry order."""
    pair_colours = {
        OPH_FEP: "#D9DEE3",
        FEP_SYNERGETICS: "#B4BCC5",
        SYNERGETICS_OPH: "#8F99A4",
        THREE_WAY: "#68727D",
    }
    status_colours = {
        "exact": GREEN,
        "exact-on-domain": ORANGE,
        "obstructed": VERMILLION,
        "not-reproduced": GREY,
    }

    def text_colour(hex_colour: str) -> str:
        return max(("#000000", "#FFFFFF"), key=lambda option: _contrast(option, hex_colour))

    figure, axis = plt.subplots(figsize=(10.0, 10.0))
    axis.set_axis_off()
    axis.set_xlim(0.0, 5.0)
    axis.set_ylim(0.0, len(CORRESPONDENCES) + 1.2)
    headers = ["pair", "kind", "status"]
    for column, header in enumerate(headers):
        axis.text(
            column + 2.15,
            len(CORRESPONDENCES) + 0.55,
            header,
            ha="center",
            va="center",
            fontsize=10.5,
            fontweight="bold",
            color=INK,
        )
    for row, entry in enumerate(CORRESPONDENCES):
        y_value = len(CORRESPONDENCES) - row - 0.5
        axis.text(
            0.02,
            y_value,
            entry.key,
            ha="left",
            va="center",
            fontsize=9.5,
            family="monospace",
            color=INK,
        )
        cells = (
            (pair_colours.get(entry.pair, GREY), entry.pair),
            (KIND_COLOURS.get(entry.kind, GREY), entry.kind),
            (status_colours.get(entry.status, GREY), entry.status),
        )
        for column, (colour, label) in enumerate(cells):
            axis.add_patch(
                plt.Rectangle(
                    (column + 1.70, y_value - 0.42),
                    0.90,
                    0.84,
                    facecolor=colour,
                    edgecolor="white",
                    linewidth=0.5,
                )
            )
            axis.text(
                column + 2.15,
                y_value,
                label,
                ha="center",
                va="center",
                fontsize=9.0,
                color=text_colour(colour),
            )
    axis.set_title("Registry scope matrix: pair, kind and status", fontsize=11)
    return _save(
        figure,
        "fig_scope_matrix.png",
        "A dense categorical matrix with one row per registry entry in registry order. The row "
        "label is the entry key, and three columns give its pair, in neutral greys, its kind and "
        "its status, each cell carrying its own text label. "
        f"{len(CORRESPONDENCES)} rows are shown.",
    )


def network_figure() -> Path:
    """The carrier read as a network, with the geometry dual pair beside it.

    Four panels, each plotted only from executed package functions: the exact
    effective-resistance ladder and its Kirchhoff total, the distance shells
    around a port with the intersection array, the subgraph census against the
    face count, and the Voronoi cell beside the cuboctahedron it is dual to.
    """
    carrier = icosahedral_carrier()
    adjacency = [[int(value) for value in row] for row in carrier.adjacency_matrix()]
    ladder = network.resistance_ladder(adjacency)
    kirchhoff = network.kirchhoff_index(adjacency)
    distribution = network.distance_distribution(adjacency)
    intersection = network.intersection_array(adjacency)
    census = network.subgraph_census(adjacency)
    face_count = len(carrier.faces)

    cell = geometry.voronoi_cell()
    cell_tetravolume = Fraction(geometry.voronoi_cell_tetravolume()).limit_denominator()
    cuboctahedron = geometry.cuboctahedron_signature()
    cuboctahedron_edge_count = len(geometry.cuboctahedron_edges())
    cuboctahedron_face_counts = geometry.cuboctahedron_faces()
    cuboctahedron_tetravolume = Fraction(geometry.cuboctahedron_volume())
    volume_ratio = cuboctahedron_tetravolume / cell_tetravolume

    figure, axes = plt.subplots(2, 2, figsize=(12.0, 10.0))
    ladder_axis, shells_axis, census_axis, dual_axis = axes.ravel()

    # --- panel 1: the resistance ladder -----------------------------------
    positions = list(range(len(ladder)))
    resistances = [float(entry[1]) for entry in ladder]
    ladder_axis.bar(positions, resistances, width=0.55, color=KIND_COLOURS["isomorphism"])
    for position, (step, resistance, pairs) in zip(positions, ladder):
        ladder_axis.text(
            position,
            float(resistance) + 0.004,
            f"{resistance.numerator}/{resistance.denominator}\n{_counted(pairs, 'pair', 'pairs')}",
            ha="center",
            va="bottom",
            fontsize=9,
            color=INK,
        )
    ladder_axis.set_xticks(positions)
    ladder_axis.set_xticklabels([f"d = {step}" for step, _, _ in ladder], fontsize=9)
    ladder_axis.set_ylabel("effective resistance")
    ladder_axis.set_ylim(0.0, max(resistances) * 1.45)
    ladder_axis.set_title("1. Resistance ladder: one exact value per distance", fontsize=10.5)
    ladder_axis.text(
        0.03,
        0.93,
        f"Kirchhoff index {kirchhoff}",
        transform=ladder_axis.transAxes,
        fontsize=9.5,
        color=KIND_COLOURS["disagreement"],
    )

    # --- panel 2: the shells around port 0 --------------------------------
    projected_carrier, flat = _projected_carrier()
    distances = network.graph_distances(adjacency)[0]
    shell_colours = {
        0: VERMILLION,
        1: BLUE,
        2: ORANGE,
        3: GREEN,
    }
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
            s=230,
            color=shell_colours.get(step, GREY),
            edgecolor=INK,
            zorder=3,
            label=f"distance {step}: {distribution[step]}",
        )
    for port in range(projected_carrier.port_count):
        shells_axis.text(
            flat[port, 0],
            flat[port, 1],
            str(port),
            ha="center",
            va="center",
            fontsize=8.5,
            color="white" if distances[port] in (0, 1, 3) else BLACK,
            zorder=4,
        )
    shells_axis.set_aspect("equal")
    shells_axis.set_axis_off()
    shells_axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 0.02),
        ncol=2,
        fontsize=9,
        title=f"intersection array: b = {intersection[0]}, c = {intersection[1]}",
        title_fontsize=9,
    )
    shells_axis.set_title("2. Shells around port 0, by graph distance", fontsize=10.5)

    # --- panel 3: the subgraph census -------------------------------------
    census_names = {
        "triangles": "triangles",
        "pentagons": "induced 5-cycles\n(links of the ports)",
        "perfect_matchings": "perfect matchings",
        "hamiltonian_cycles": "Hamiltonian cycles",
    }
    census_labels = [census_names.get(key, key.replace("_", " ")) for key in census]
    census_positions = list(range(len(census)))
    census_values = list(census.values())
    census_axis.barh(
        census_positions,
        census_values,
        height=0.6,
        color=[BLUE, ORANGE, GREEN, VERMILLION],
    )
    for position, count in zip(census_positions, census_values):
        census_axis.annotate(
            str(count),
            (count, position),
            textcoords="offset points",
            xytext=(5, 0),
            ha="left",
            va="center",
            fontsize=9.5,
            color=INK,
        )
    census_axis.set_xscale("log")
    census_axis.axvline(face_count, color=VERMILLION, linestyle="--", linewidth=1.3)
    census_axis.text(
        face_count * 1.06,
        -0.47,
        f"{face_count} faces",
        fontsize=9,
        color=VERMILLION,
        va="center",
    )
    census_axis.set_yticks(census_positions)
    census_axis.set_yticklabels(census_labels, fontsize=9)
    census_axis.invert_yaxis()
    census_axis.set_xlim(1.0, max(census_values) * 4.0)
    census_axis.set_xlabel("count (log scale)")
    census_axis.set_title("3. Subgraph census against the face count", fontsize=10.5)

    # --- panel 4: the dual pair -------------------------------------------
    dual_axis.set_axis_off()
    dual_axis.set_title("4. The dual pair: Voronoi cell and cuboctahedron", fontsize=10.5)

    def planar(points) -> np.ndarray:
        pts = np.asarray(points, dtype=float)
        pts = pts - pts.mean(axis=0)
        e1 = np.array([1.0, -1.0, 0.0]) / math.sqrt(2.0)
        e2 = np.array([1.0, 1.0, -2.0]) / math.sqrt(6.0)
        return np.stack([pts @ e1, pts @ e2], axis=1)

    def draw_solid(centre_x, projected, edges, colour, caption):
        span = max(float(np.linalg.norm(row)) for row in projected)
        placed = projected / span * 0.46 + np.array([centre_x, 0.32])
        for a, b in edges:
            dual_axis.plot(
                [placed[a, 0], placed[b, 0]],
                [placed[a, 1], placed[b, 1]],
                color=colour,
                linewidth=1.2,
                zorder=2,
            )
        dual_axis.scatter(placed[:, 0], placed[:, 1], s=30, color=colour, zorder=3)
        dual_axis.text(
            centre_x,
            -0.24,
            caption,
            ha="center",
            va="top",
            fontsize=9,
            color=INK,
        )

    cell_edges = set()
    for face in cell["faces"]:
        for index, vertex in enumerate(face):
            following = face[(index + 1) % len(face)]
            cell_edges.add((min(vertex, following), max(vertex, following)))
    draw_solid(
        -1.05,
        planar(cell["vertices"]),
        sorted(cell_edges),
        KIND_COLOURS["complementarity"],
        f"Voronoi cell: {cell['vertex_count']} vertices,\n"
        f"{cell['face_count']} rhombic faces, {cell['edge_count']} edges,\n"
        f"tetravolume {cell_tetravolume}",
    )
    cuboctahedron_vertices = planar(
        [geometry.to_xyz(site) for site in geometry.cuboctahedron_neighbors()]
    )
    cuboctahedron_adj = geometry.cuboctahedron_adjacency()
    cuboctahedron_edge_pairs = [
        (a, b)
        for a in range(len(cuboctahedron_adj))
        for b in range(a + 1, len(cuboctahedron_adj))
        if cuboctahedron_adj[a][b]
    ]
    draw_solid(
        1.05,
        cuboctahedron_vertices,
        cuboctahedron_edge_pairs,
        KIND_COLOURS["isomorphism"],
        f"Cuboctahedron: {cuboctahedron['vertices']} vertices,\n"
        f"{len(cuboctahedron_face_counts['triangles'])} triangles + "
        f"{len(cuboctahedron_face_counts['squares'])} squares, {cuboctahedron_edge_count} edges,\n"
        f"tetravolume {cuboctahedron_tetravolume}",
    )
    dual_axis.set_xlim(-2.4, 2.4)
    dual_axis.set_ylim(-0.78, 1.05)
    dual_axis.text(
        0.0,
        -0.70,
        f"volume ratio {cuboctahedron_tetravolume}/{cell_tetravolume} = {volume_ratio}, exact",
        ha="center",
        fontsize=9.5,
        color=INK,
    )

    return _save(
        figure,
        "fig_network.png",
        "Four panels, every number read from the executed package at draw time. Panel one is the "
        "resistance ladder, one bar for each distance class of the carrier's port graph, labelled "
        "with its exact fraction and its pair count and annotated with the Kirchhoff index. Panel "
        "two draws the twelve ports in the icosahedral spectral projection, coloured by graph "
        "distance from port zero with that shell distribution as its legend below the drawing and "
        "the intersection array in the legend title. Panel three is the subgraph census on a "
        "logarithmic axis, horizontal bars for triangles, induced five-cycles (the chordless "
        "links of the twelve ports, which are the pentagonal faces of the dual dodecahedron), "
        "perfect matchings and Hamiltonian cycles, with the icosahedron face count drawn as a "
        "reference line. Panel four compares the Voronoi cell and the cuboctahedron by their "
        "vertex, face and edge counts and their tetravolumes, ending in the exact volume ratio.",
    )


def vertex_shells_figure() -> Path:
    """The carrier's shells around a port against the IVM's lattice shells."""
    adjacency = [[int(value) for value in row] for row in icosahedral_carrier().adjacency_matrix()]
    carrier_shells = network.distance_distribution(adjacency)
    radii = (1, 2, 3)
    carrier_counts = carrier_shells[1:4]
    lattice_counts = tuple(lattice.shell_count(k) for k in radii)
    outward, homing = network.intersection_array(adjacency)
    figure, axis = plt.subplots(figsize=(9.0, 4.8))
    width = 0.38
    positions = np.array(radii, dtype=float)
    floor = 0.5  # log-axis baseline, so a count of one still has a visible bar
    carrier_bars = axis.bar(
        positions - width / 2,
        np.array(carrier_counts) - floor,
        width,
        bottom=floor,
        color=BLUE,
        label="carrier: ports at graph distance $k$",
    )
    lattice_bars = axis.bar(
        positions + width / 2,
        np.array(lattice_counts) - floor,
        width,
        bottom=floor,
        color=GREEN,
        hatch="//",
        edgecolor=INK,
        label="IVM lattice: sites on shell $k$",
    )
    axis.set_yscale("log")
    axis.set_ylim(floor, max(lattice_counts) * 3.0)
    ticks = sorted(set(carrier_counts) | set(lattice_counts))
    axis.set_yticks(ticks)
    axis.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axis.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    for bars in (carrier_bars, lattice_bars):
        for bar in bars:
            axis.annotate(
                f"{int(round(bar.get_height() + floor))}",
                (bar.get_x() + bar.get_width() / 2, bar.get_height() + floor),
                textcoords="offset points",
                xytext=(0, 3),
                ha="center",
                va="bottom",
                fontsize=9.5,
            )
    axis.set_xticks(list(radii))
    axis.set_xlim(0.45, 3.55)
    axis.set_xlabel("radius $k$ from a port")
    axis.set_ylabel("count (log scale)")
    axis.set_title(
        f"shells around a port: {carrier_counts} against {lattice_counts};\n"
        f"intersection array $b={outward}$, $c={homing}$",
        fontsize=10.5,
    )
    axis.legend(loc="upper left")
    return _save(
        figure,
        "fig_vertex_shells.png",
        "Two bar series against the integer radius from a port, on a logarithmic count axis. The "
        f"carrier's ports number {carrier_counts[0]}, {carrier_counts[1]} and {carrier_counts[2]} at "
        "graph distances one, two and three, while the IVM lattice's shells hold "
        f"{lattice_counts[0]}, {lattice_counts[1]} and {lattice_counts[2]} sites at the same radii, "
        "so the sequences disagree at the first shell and no radius repairs the reading. The "
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
    volume = float(geometry.icosahedron_tetravolume()["tetravolume"])
    norm = geometry.tetravolume_galois_norm()
    conjugate = float(norm) / volume
    figure, (left, right) = plt.subplots(
        1, 2, figsize=(11.0, 4.4), gridspec_kw={"width_ratios": (1.15, 1.0)}
    )
    names = [name for name, _ in fixed]
    counts = [value for _, value in fixed]
    bars = left.barh(names, counts, color=BLUE)
    left.set_xscale("log")
    left.set_xlim(1.0, max(counts) * 8.0)
    for bar, count in zip(bars, counts):
        left.annotate(
            f"{count:,}",
            (count, bar.get_y() + bar.get_height() / 2),
            textcoords="offset points",
            xytext=(5, 0),
            ha="left",
            va="center",
            fontsize=9.5,
        )
    left.invert_yaxis()
    left.set_xlabel("count (log scale)")
    left.set_title(r"fixed by $\sqrt{5}\to-\sqrt{5}$")
    bars = right.barh(
        ("icosahedron", "conjugate"),
        (volume, conjugate),
        color=(VERMILLION, ORANGE),
        height=0.55,
    )
    right.set_xscale("log")
    right.set_xlim(1.0, max(volume, conjugate) * 6.0)
    for bar, value in zip(bars, (volume, conjugate)):
        right.annotate(
            f"{value:.6g}",
            (value, bar.get_y() + bar.get_height() / 2),
            textcoords="offset points",
            xytext=(5, 0),
            ha="left",
            va="center",
            fontsize=9.5,
        )
    right.invert_yaxis()
    right.set_xlabel("tetravolume (log scale)")
    right.set_title("moved by the same map")
    right.text(
        0.55,
        0.5,
        f"tetravolume × conjugate = {norm}\n(the exact Galois norm)",
        transform=right.transAxes,
        ha="center",
        va="center",
        fontsize=10,
        color=INK,
        bbox={"facecolor": "#F1F5F9", "edgecolor": SKY, "boxstyle": "round,pad=0.35"},
    )
    figure.suptitle("the counting invariants are fixed, the volumes are not", fontsize=11)
    return _save(
        figure,
        "fig_conjugation.png",
        "Two panels of horizontal bars on logarithmic axes. Left: the counting invariants that "
        f"conjugating the square root of five fixes: the spanning-tree count {counts[0]}, the "
        f"automorphism group order {counts[1]}, the {counts[2]} sheaf obstructions and the "
        f"{counts[3]} triangles. Right: the icosahedron's tetravolume {volume:.6g} beside its "
        f"conjugate {conjugate:.6g}, with the product annotated as the exact norm {norm}.",
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
    ico_volume = float(geometry.icosahedron_tetravolume()["tetravolume"])
    cubo_volume = float(lattice.cuboctahedron_volume())
    cubo_root = math.isqrt(cubo_trees)
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.8))
    rows = ("degree", "edges", "triangles", "squares")
    cubo_counts = (cubo_degree, cubo_edges, cubo_triangles, cubo_squares)
    ico_counts = (ico_degree, ico_edges, ico_triangles, 0)
    positions = np.arange(len(rows))
    width = 0.38
    cubo_bars = left.barh(
        positions - width / 2,
        cubo_counts,
        width,
        color=BLUE,
        label="cuboctahedron",
    )
    ico_bars = left.barh(
        positions + width / 2,
        ico_counts,
        width,
        color=VERMILLION,
        label="carrier (icosahedron)",
    )
    left.set_yticks(positions, rows)
    left.invert_yaxis()
    left.set_xlim(0.0, max(cubo_counts + ico_counts) * 1.2)
    for bars, values in ((cubo_bars, cubo_counts), (ico_bars, ico_counts)):
        for bar, value in zip(bars, values):
            left.annotate(
                f"{value}",
                (bar.get_width(), bar.get_y() + bar.get_height() / 2),
                textcoords="offset points",
                xytext=(5, 0),
                va="center",
                ha="left",
                fontsize=9.5,
            )
    left.set_title(
        f"counts; tetravolumes {cubo_volume:.6g} against {ico_volume:.6g}",
        fontsize=10.5,
    )
    left.legend(loc="lower right")
    trees = (ico_trees, cubo_trees)
    bars = right.bar(("carrier", "cuboctahedron"), trees, color=(VERMILLION, BLUE))
    right.set_yscale("log")
    right.set_ylim(trees[1] / 4.0, trees[0] * 6.0)
    for bar, count in zip(bars, trees):
        label = f"{count:,}"
        if count == cubo_trees:
            label += f"\n$= {cubo_root}^2$"
        right.annotate(
            label,
            (bar.get_x() + bar.get_width() / 2, count),
            textcoords="offset points",
            xytext=(0, 4),
            ha="center",
            va="bottom",
            fontsize=9.5,
        )
    right.set_title("spanning trees of the seam graph, exact", fontsize=10.5)
    right.set_ylabel("spanning trees (log scale)")
    figure.suptitle("the same twelve-around-one phrase, two different solids", fontsize=11)
    return _save(
        figure,
        "fig_families.png",
        "Two panels. Left: horizontal bars comparing the cuboctahedron with the carrier's "
        f"icosahedron: degree {cubo_degree} against {ico_degree}, edges {cubo_edges} against "
        f"{ico_edges}, triangles {cubo_triangles} against {ico_triangles}, and squares "
        f"{cubo_squares} against 0, with the tetravolumes {cubo_volume:.6g} and "
        f"{ico_volume:.6g} annotated in the title. Right: the spanning trees the same "
        f"determinant counts, {ico_trees} for the carrier against {cubo_trees} for the "
        f"cuboctahedron, which is {cubo_root} squared, on a logarithmic axis. The cuboctahedron "
        "is blue and the carrier vermillion in both panels.",
    )


def sandwich_figure() -> Path:
    """The surprise sandwich and the functionals that share the name free energy."""
    dictionary = by_key()["free-energy-dictionary"].run().values
    collision = by_key()["free-energy-name-collision"].run().values
    lower = float(dictionary["sandwich_lower"])
    upper = float(dictionary["sandwich_upper"])
    bayesian = float(dictionary["bayesian_surprise"])
    numerical = float(dictionary["numerical_surprise"])
    surprisal = float(dictionary["surprisal"])
    variational = float(dictionary["variational_free_energy"])
    interface = float(collision["interface_energy_minimised_over_interior"])
    at_posterior = float(collision["variational_free_energy_at_posterior"])
    constant = float(collision["difference"])
    minus = "−"

    def signed(value: float) -> str:
        return f"{value:.4g}".replace("-", minus)

    figure, (left, right) = plt.subplots(1, 2, figsize=(11.5, 4.8))
    span = upper - lower
    left.axvspan(lower, upper, color=SKY, alpha=0.30)
    left.axvline(bayesian, color=VERMILLION, linewidth=2.2)
    left.text(
        bayesian + 0.02 * span,
        0.40,
        f"Bayesian surprise\n{signed(bayesian)}",
        ha="left",
        va="center",
        fontsize=9.5,
        color=VERMILLION,
    )
    for value, name, align, offset in (
        (lower, "lower bound", "right", -0.02 * span),
        (upper, "upper bound", "left", 0.02 * span),
    ):
        left.axvline(value, color=BLUE, linewidth=1.2, linestyle="--")
        left.text(
            value + offset,
            0.80,
            f"{name}\n{signed(value)}",
            ha=align,
            va="center",
            fontsize=9.5,
            color=BLUE,
        )
    left.set_xlim(lower - 0.5 * span, upper + 0.5 * span)
    left.set_ylim(0.0, 1.0)
    left.set_yticks([])
    left.set_xlabel("surprisal units")
    left.set_title("the sandwich brackets the Bayesian value", fontsize=10.5)

    # The two bars that the model constant relates sit on adjacent top rows.
    entries = (
        ("numerical", numerical, BLUE),
        ("Bayesian", bayesian, BLUE),
        ("surprisal", surprisal, GREEN),
        ("variational", variational, GREEN),
        ("VFE at posterior", at_posterior, VERMILLION),
        ("interface", interface, VERMILLION),
    )
    names = [name for name, _, _ in entries]
    values = [value for _, value, _ in entries]
    bars = right.barh(names, values, color=[colour for _, _, colour in entries])
    right.axvline(0.0, color=BLACK, linewidth=0.8)
    for bar, value in zip(bars, values):
        right.annotate(
            signed(value),
            (value, bar.get_y() + bar.get_height() / 2),
            ha="left" if value >= 0 else "right",
            va="center",
            fontsize=9,
            xytext=(5 if value >= 0 else -5, 0),
            textcoords="offset points",
        )
    top_row = len(entries) - 1
    bracket_y = top_row + 0.75
    for x_value, row in ((interface, top_row), (at_posterior, top_row - 1)):
        right.plot([x_value, x_value], [row + 0.4, bracket_y], color=VERMILLION, linewidth=1.0,
                   linestyle=":")
    right.annotate(
        "",
        xy=(at_posterior, bracket_y),
        xytext=(interface, bracket_y),
        arrowprops={"arrowstyle": "<->", "color": VERMILLION, "lw": 1.3},
    )
    right.text(
        (interface + at_posterior) / 2,
        bracket_y + 0.12,
        f"$\\Delta$ = model constant = {signed(constant)}",
        ha="center",
        va="bottom",
        fontsize=9.5,
        color=VERMILLION,
    )
    right.set_xlim(min(values) * 1.5, max(values) * 1.3)
    right.set_ylim(-0.6, bracket_y + 0.8)
    right.set_title("values called free energy", fontsize=10.5)
    right.tick_params(axis="y", labelsize=9.5)
    figure.suptitle("what the probes report about the values called free energy", fontsize=11)
    return _save(
        figure,
        "fig_sandwich.png",
        "Two panels read from the executed probes. Left: the sandwich interval with "
        f"its lower bound {lower:.7g} and upper bound {upper:.7g}, and the Bayesian surprise "
        f"{bayesian:.7g} marked inside it. Right: horizontal bars for the functionals the "
        f"literature calls free energy in the same linear-Gaussian model: the numerical "
        f"surprise {numerical:.7g}, the Bayesian surprise {bayesian:.7g}, the surprisal "
        f"{surprisal:.7g}, the variational free energy {variational:.7g}, the variational free "
        f"energy at the posterior {at_posterior:.7g} and the interface energy minimised over "
        f"the interior {interface:.7g}; a bracket above the last two bars marks that the "
        f"interior and posterior values differ by the model constant {constant:.7g}.",
    )


def schedule_figure() -> Path:
    """The centered contrast between two schedules, one exact fraction per step."""
    steps_list = list(range(1, 13))
    contrasts = [
        abs(float(controls.affine_path_contrast_exact(steps=step)))
        for step in steps_list
    ]
    exact = controls.affine_path_contrast_exact(steps=12)
    figure, axis = plt.subplots(figsize=(8.8, 4.4))
    bars = axis.bar([str(step) for step in steps_list], contrasts, color=VERMILLION)
    for bar, value in zip(bars, contrasts):
        if value > 0:
            axis.annotate(
                f"{value:.4g}",
                (bar.get_x() + bar.get_width() / 2, value),
                textcoords="offset points",
                xytext=(0, 3),
                ha="center",
                va="bottom",
                fontsize=8,
            )
    axis.annotate(
        f"after twelve steps: {exact}",
        (0.97, 0.88),
        ha="right",
        va="top",
        fontsize=9,
        color=INK,
        xycoords="axes fraction",
    )
    axis.set_xlabel("steps of the affine iteration")
    axis.set_ylabel("absolute centered-contrast difference")
    axis.set_title("full and half-damped schedules agree on answers, not on contrasts", fontsize=10.5)
    return _save(
        figure,
        "fig_schedule.png",
        "The absolute difference between the full and the half-damped schedules' centered "
        "contrasts, one exact value per step count from one to twelve. The schedules share "
        f"every answer while the contrast difference decays to {contrasts[-1]:.7g}, which is "
        f"exactly {exact} after twelve steps.",
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
    """A coordination shell is a polyhedral surface; a free-energy level is a sphere."""
    bound = as_dict()
    window = 6
    shared_levels = [float(level) for level in bound["shell.shared_free_energy_levels"]]
    top_shell = geometry.MAX_SHELL if hasattr(geometry, "MAX_SHELL") else 12
    shell_range = range(1, 13)

    # --- left: the plane u3 = 0 of D3 coloured by the closed-form gauge -----
    slice_points: dict[tuple[int, int], int] = {}
    for site in geometry.sites_through_shell(window):
        halved = [int(round(value)) // 2 for value in geometry.to_xyz(site)]
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
        1, 2, figsize=(12.0, 5.8), gridspec_kw={"width_ratios": (1.0, 1.05)}
    )
    cmap = plt.get_cmap("viridis", window + 1)
    scatter = left.scatter(
        xs, ys, c=gauges, cmap=cmap, vmin=-0.5, vmax=window + 0.5, s=70, edgecolor=INK,
        linewidth=0.5, zorder=4,
    )
    colourbar = figure.colorbar(scatter, ax=left, ticks=range(window + 1), fraction=0.046, pad=0.02)
    colourbar.set_label("coordination shell $g(u)$")
    for k in range(1, window + 1):
        ring = [point for point, value in slice_points.items() if value == k]
        hull = _convex_hull([(float(a), float(b)) for a, b in ring])
        hull.append(hull[0])
        left.plot(
            [p[0] for p in hull], [p[1] for p in hull], color=INK, linewidth=1.1, zorder=2,
            label="shell boundary $g=k$" if k == 1 else None,
        )
    demo_k = math.isqrt(int(round(min(shared_levels) / 4.0)))
    theta = np.linspace(0.0, 2.0 * math.pi, 400)
    circles = (
        (2.0 * demo_k**2, VERMILLION, "-", f"$F={2 * demo_k**2:g}=2k^2$, $k={demo_k}$"),
        (4.0 * demo_k**2, ORANGE, "-", f"$F={4 * demo_k**2:g}=4k^2$, $k={demo_k}$"),
        (max(shared_levels), PURPLE, "--", f"$F={max(shared_levels):g}$, shared level"),
    )
    for level, colour, style, label in circles:
        radius = math.sqrt(level / 2.0)  # F = 2 |u|^2 on the halved lattice
        left.plot(radius * np.cos(theta), radius * np.sin(theta), color=colour, linestyle=style,
                  linewidth=2.0, zorder=3, label=f"free-energy level {label}")
    left.set_aspect("equal")
    left.set_xlim(-window - 0.8, window + 0.8)
    left.set_ylim(-window - 0.8, window + 0.8)
    left.set_xticks(range(-window, window + 1, 2))
    left.set_yticks(range(-window, window + 1, 2))
    left.set_xlabel("$u_1$")
    left.set_ylabel("$u_2$")
    left.set_title("plane $u_3=0$ of the lattice $D_3$:\nshells are squares here, levels are circles")
    left.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, fontsize=8.5)

    # --- right: the range of F on each shell --------------------------------
    shared: dict[float, set[int]] = {}
    per_shell: dict[int, list[float]] = {}
    vertex_counts = []
    for k in shell_range:
        energies = [0.5 * geometry.squared_norm(site) for site in geometry.generate_shell(k)]
        levels = sorted(set(energies))
        per_shell[k] = levels
        vertex_counts.append(sum(1 for value in energies if value == 4.0 * k * k))
        for level in levels:
            shared.setdefault(level, set()).add(k)
    ks = np.array(list(shell_range))
    lows = np.array([per_shell[k][0] for k in shell_range])
    highs = np.array([per_shell[k][-1] for k in shell_range])
    right.bar(ks, highs - lows, bottom=lows, width=0.55, color=SKY, alpha=0.55, zorder=1,
              label="range of $F$ on shell")
    single_x, single_y, shared_x, shared_y = [], [], [], []
    for k in shell_range:
        for level in per_shell[k]:
            if len(shared[level]) > 1:
                shared_x.append(k)
                shared_y.append(level)
            else:
                single_x.append(k)
                single_y.append(level)
    right.scatter(single_x, single_y, s=14, color=BLUE, zorder=3, label="level on one shell")
    right.scatter(shared_x, shared_y, s=34, color=ORANGE, edgecolor=INK, linewidth=0.6, zorder=4,
                  label="level shared by shells")
    curve = np.linspace(ks.min() - 0.3, ks.max() + 0.3, 100)
    right.plot(curve, 2.0 * curve**2, color=VERMILLION, linewidth=1.6, zorder=2, label="$2k^2$")
    right.plot(curve, 4.0 * curve**2, color=BLACK, linewidth=1.6, linestyle="--", zorder=2,
               label="$4k^2$")
    right.set_xticks(list(shell_range))
    right.set_xlim(0.4, 12.6)
    right.set_xlabel("coordination shell $k$")
    right.set_ylabel(r"free energy $F=\frac{1}{2}\|x\|^2$")
    right.set_title("free-energy range on each shell: neighbouring shells overlap")
    right.legend(loc="upper left", fontsize=8.5, ncol=2)
    if len(set(vertex_counts)) != 1:
        raise RuntimeError(f"the maximum 4k^2 is not attained at a constant count: {vertex_counts}")
    right.text(
        0.02,
        0.50,
        f"$F\\leq 4k^2$, equality only at\nthe {vertex_counts[0]} vertices;\n"
        r"$F\geq 2k^2$, so $F\leq f$ forces" "\n"
        r"$k\leq\sqrt{f/2}$: free energy bounds" "\n"
        "the shell, but neither is a\nfunction of the other",
        transform=right.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        color=INK,
        bbox={"facecolor": "white", "edgecolor": SKY, "boxstyle": "round,pad=0.35"},
        zorder=5,
    )
    del top_shell
    return _save(
        figure,
        "fig_shell_gauge.png",
        "Two panels on the coordination shell and the free energy. Left: the lattice points of "
        f"the face-centred cubic lattice D3 in the plane u3 equal to zero with both coordinates "
        f"between minus {window} and {window}, coloured by their closed-form coordination "
        "shell; the shell boundaries are drawn as nested squares in this slice, and three "
        f"circles are overlaid for the free-energy levels {2 * demo_k**2:g}, "
        f"{4 * demo_k**2:g} and {max(shared_levels):g}, so a shell is a polyhedral surface and a "
        "free-energy level is a sphere that crosses several shells. Right: for shells one to "
        "twelve, a vertical bar from the smallest to the largest free energy on the shell, a dot "
        "for every distinct level with the levels shared between shells highlighted, and the "
        f"curves two k squared and four k squared; the maximum four k squared is attained at "
        f"exactly {vertex_counts[0]} sites on every shell, the twelve vertices. A free energy "
        "bounded above bounds the shell, but shell and level are not functions of one another.",
    )


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
        pair_triangles_figure,
        scope_matrix_figure,
        network_figure,
        vertex_shells_figure,
        conjugation_figure,
        families_figure,
        sandwich_figure,
        schedule_figure,
        shell_gauge_figure,
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

