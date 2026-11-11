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

KIND_COLOURS = {
    "isomorphism": "#0072b2",
    "disagreement": "#d55e00",
    "complementarity": "#009e73",
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
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
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
    carrier = icosahedral_carrier()
    eigenvalues, eigenvectors = np.linalg.eigh(carrier.laplacian())
    lam1 = 5.0 - math.sqrt(5.0)
    columns = [index for index, value in enumerate(eigenvalues) if abs(value - lam1) < 1e-9]
    if len(columns) != 3:
        raise RuntimeError(f"expected a three-dimensional lambda_1 eigenspace, got {len(columns)}")
    points = eigenvectors[:, columns] * math.sqrt(carrier.port_count)

    # The eigenspace is the icosahedron up to an orthogonal transform, so choose a
    # rotation whose planar projection separates the ports.  The search is a fixed
    # grid, which keeps the drawing deterministic.
    best: tuple[float, np.ndarray] = (-1.0, np.zeros((carrier.port_count, 2)))
    for theta in np.linspace(0.0, math.pi, 24, endpoint=False):
        for phi in np.linspace(0.0, math.pi, 12, endpoint=False):
            rotation = _rotation(theta, phi)
            projected = points @ rotation.T
            flat = projected[:, :2]
            gaps = [
                float(np.linalg.norm(flat[i] - flat[j]))
                for i in range(carrier.port_count)
                for j in range(i + 1, carrier.port_count)
            ]
            score = min(gaps)
            if score > best[0]:
                best = (score, flat)
    flat = best[1]

    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 5.0))

    for face in OPH_ORIENTED_BASE_FACES:
        polygon = flat[list(face)]
        left.fill(polygon[:, 0], polygon[:, 1], color="#c9d6e4", alpha=0.55, linewidth=0)
    for port in range(carrier.port_count):
        for other in carrier.neighbours(port):
            if other > port:
                left.plot(
                    [flat[port, 0], flat[other, 0]],
                    [flat[port, 1], flat[other, 1]],
                    color="#243b53",
                    linewidth=1.1,
                    zorder=2,
                )
    left.scatter(flat[:, 0], flat[:, 1], s=190, color="white", edgecolor="#102a43", zorder=3)
    for port in range(carrier.port_count):
        left.text(
            flat[port, 0],
            flat[port, 1],
            str(port),
            ha="center",
            va="center",
            fontsize=7.5,
            color="#102a43",
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
    right.vlines(marks, 0, [counts[m] for m in marks], color="#1f4e79", linewidth=2.2)
    right.scatter(marks, [counts[m] for m in marks], color="#1f4e79", zorder=3, s=28)
    for value in marks:
        right.annotate(
            f"{value:.4g} x{counts[value]}",
            (value, counts[value]),
            textcoords="offset points",
            xytext=(4, 4),
            fontsize=7.5,
        )
    lam1_value = carrier.spectral_gap()
    normalised = lam1_value / (2.0 * carrier.seam_count)
    right.annotate(
        f"$\\lambda_1 = 5 - \\sqrt{{5}} = {lam1_value:.6f}$\n"
        f"normalised gap $\\lambda_1/60 = {normalised:.6f}$",
        xy=(lam1, 3.05),
        xytext=(0.35, 5.9),
        fontsize=8,
        arrowprops={"arrowstyle": "->", "color": "#a4262c", "lw": 1.0},
        color="#a4262c",
    )
    right.set_xlabel("Laplacian eigenvalue")
    right.set_ylabel("multiplicity")
    right.set_ylim(0, 6.8)
    right.set_title("Closed-form spectrum of the carrier")

    return _save(
        figure,
        "fig_carrier.png",
        f"Two panels. Left: the {int(values['carrier.ports'])} ports of the carrier drawn in the "
        f"three-dimensional spectral embedding of its own Laplacian, with the "
        f"{int(values['carrier.seams'])} seams drawn and the {int(values['carrier.faces'])} "
        f"oriented faces shaded; every port is labelled with its index and meets "
        f"{int(values['carrier.degree'])} others. "
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


def _projected_carrier() -> tuple[IcosahedralCarrier, np.ndarray]:
    """Return the carrier and its deterministic planar spectral projection."""
    carrier = icosahedral_carrier()
    eigenvalues, eigenvectors = np.linalg.eigh(carrier.laplacian())
    lam1 = 5.0 - math.sqrt(5.0)
    columns = [index for index, value in enumerate(eigenvalues) if abs(value - lam1) < 1e-9]
    if len(columns) != 3:
        raise RuntimeError(f"expected a three-dimensional lambda_1 eigenspace, got {len(columns)}")
    points = eigenvectors[:, columns] * math.sqrt(carrier.port_count)
    best: tuple[float, np.ndarray] = (-1.0, np.zeros((carrier.port_count, 2)))
    for theta in np.linspace(0.0, math.pi, 24, endpoint=False):
        for phi in np.linspace(0.0, math.pi, 12, endpoint=False):
            rotation = _rotation(theta, phi)
            projected = points @ rotation.T
            flat = projected[:, :2]
            gaps = [
                float(np.linalg.norm(flat[i] - flat[j]))
                for i in range(carrier.port_count)
                for j in range(i + 1, carrier.port_count)
            ]
            score = min(gaps)
            if score > best[0]:
                best = (score, flat)
    return carrier, best[1]


def registry_figure() -> Path:
    """The whole registry as three status-marked columns."""
    figure, axes = plt.subplots(1, 3, figsize=(12.0, 5.2))
    for axis, kind in zip(axes, ("isomorphism", "disagreement", "complementarity")):
        entries = [entry for entry in CORRESPONDENCES if entry.kind == kind]
        entries = sorted(entries, key=lambda entry: (entry.status, entry.key))
        positions = list(range(len(entries)))
        for row, entry in zip(positions, entries):
            axis.scatter(
                [0.0],
                [row],
                marker=STATUS_MARKERS.get(entry.status, "o"),
                s=70,
                color=KIND_COLOURS[kind],
            )
            axis.text(0.14, row, entry.key, va="center", fontsize=9.4, family="monospace")
        axis.set_xlim(-0.25, 1.5)
        axis.set_ylim(-0.8, max(len(entries), 1) - 0.2)
        axis.invert_yaxis()
        axis.set_yticks([])
        axis.set_xticks([])
        axis.grid(False)
        for spine in axis.spines.values():
            spine.set_visible(False)
        axis.set_title(f"{kind} ({len(entries)})", color=KIND_COLOURS[kind], fontsize=10)
    handles = [
        plt.Line2D([], [], marker=marker, linestyle="", color="#4a4a4a", label=status, markersize=7)
        for status, marker in STATUS_MARKERS.items()
    ]
    figure.legend(handles=handles, loc="lower center", ncol=3, fontsize=8)
    figure.suptitle(
        f"{len(CORRESPONDENCES)} correspondences: every row is executed by "
        "oph-fep-synergetics verify",
        fontsize=10,
    )
    kind_counts = {
        kind: sum(1 for entry in CORRESPONDENCES if entry.kind == kind)
        for kind in ("isomorphism", "disagreement", "complementarity")
    }
    return _save(
        figure,
        "fig_registry.png",
        f"The {len(CORRESPONDENCES)} correspondences arranged in three columns by kind. "
        f"{_counted(kind_counts['isomorphism'], 'isomorphism', 'isomorphisms')} on the left, "
        f"{_counted(kind_counts['disagreement'], 'disagreement', 'disagreements')} in the middle, "
        f"{_counted(kind_counts['complementarity'], 'complementarity', 'complementarities')} on "
        "the right. Each row is labelled with its registry key and carries a marker for its "
        "status: a circle for exact, a square for exact-on-domain, and a cross for obstructed.",
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
    bars = left.bar(labels, values, color=["#2f6b3a", "#c98a00", "#7a7a7a"])
    for bar, value in zip(bars, values):
        left.annotate(
            str(value),
            (bar.get_x() + bar.get_width() / 2, bar.get_height()),
            textcoords="offset points",
            xytext=(0, 3),
            ha="center",
            fontsize=8,
        )
    left.set_ylabel("transitions")
    left.set_title(
        f"{census['states']} states x 30 seams x 2 coins = {sum(values)} transitions"
    )
    left.set_ylim(0, max(values) * 1.15)

    matrix = np.asarray(excursion, dtype=float)
    image = right.imshow(matrix.T, aspect="auto", cmap="Blues", origin="lower")
    right.set_yticks(range(matrix.shape[1]))
    right.set_yticklabels([str(index) for index in range(matrix.shape[1])], fontsize=7)
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
    onset = math.sqrt(10.0)
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.6))
    for k in samples:
        eigenvalues = controls.d_stability_eigenvalues(k)
        unstable = bool(np.any(eigenvalues.real > 1e-9))
        colour = "#a4262c" if unstable else "#1f4e79"
        left.scatter(eigenvalues.real, eigenvalues.imag, s=42, color=colour, label=f"k = {k:g}")
        left.plot(eigenvalues.real, eigenvalues.imag, linestyle=":", color=colour, linewidth=0.8)
    left.axvline(0.0, color="black", linewidth=0.9)
    left.set_xlabel("real part")
    left.set_ylabel("imaginary part")
    left.set_title("eigenvalues of the scaled matrix")
    left.legend(fontsize=7, loc="lower right")

    grid = np.linspace(0.05, 5.0, 400)
    third = [controls.routh_first_column(k)[2] for k in grid]
    right.plot(grid, third, color="#1f4e79", linewidth=1.8)
    right.axhline(0.0, color="black", linewidth=0.9)
    right.axvline(onset, color="#a4262c", linestyle="--", linewidth=1.0)
    right.annotate(
        f"onset at k = sqrt(10) = {onset:.4f}",
        xy=(onset, 0.0),
        xytext=(onset + 0.35, 2.6),
        fontsize=8,
        color="#a4262c",
        arrowprops={"arrowstyle": "->", "color": "#a4262c", "lw": 1.0},
    )
    right.scatter([4.0], [controls.routh_first_column(4.0)[2]], color="#a4262c", zorder=3)
    right.annotate("k = 4 gives -2", xy=(4.0, -2.0), xytext=(2.4, -5.2), fontsize=8)
    right.set_xlabel("coupling k")
    right.set_ylabel("third Routh entry")
    right.set_title("2 (10 - k^2) / (k + 2)")
    return _save(
        figure,
        "fig_dstability.png",
        "Two panels. Left: the eigenvalue loci of the scaled matrix in the complex plane for "
        f"{len(values['d_stability.samples'])} couplings, with a dashed imaginary axis; the pair "
        "that crosses into the right half plane is drawn in red and the stable cases in blue. "
        "Right: the third Routh entry two times ten minus k squared over k plus two against the "
        f"coupling, crossing zero at {format_value(values['d_stability.routh_onset'])} and equal "
        f"to {format_value(values['d_stability.routh_k4_third'])} at k equal to four.",
    )


def mediator_figure() -> Path:
    """The pairwise margin and the reduction error."""
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.4))
    grid = np.linspace(0.0, 1.0, 200)
    margins = [controls.mediator_pairwise_margin(value) for value in grid]
    left.plot(grid, margins, color="#1f4e79", linewidth=1.8)
    left.axhline(0.0, color="black", linewidth=0.9)
    left.axvline(1.0 / 3.0, color="#a4262c", linestyle="--", linewidth=1.0)
    left.scatter([0.2, 1.0 / 3.0], [controls.mediator_pairwise_margin(0.2), 0.0], color="#a4262c", zorder=3)
    left.annotate("epsilon = 1/3", xy=(1.0 / 3.0, 0.0), xytext=(0.42, 0.55), fontsize=8, color="#a4262c")
    left.set_xlabel("epsilon")
    left.set_ylabel("pairwise margin")
    left.set_title("3 epsilon - 1")

    reduction = controls.reduction_control()
    names = ["same observation", "first after", "second after", "response distance", "worst-case floor"]
    keys = ["same_observation", "first_after", "second_after", "response_distance", "worst_case_error_lower_bound"]
    values = [float(reduction[key]) for key in keys]
    bars = right.barh(names, values, color=["#1f4e79", "#2f6b3a", "#a4262c", "#1f4e79", "#a4262c"])
    for bar, value in zip(bars, values):
        right.annotate(
            f"{value:g}",
            (value, bar.get_y() + bar.get_height() / 2),
            textcoords="offset points",
            xytext=(4 if value >= 0 else -18, 0),
            va="center",
            fontsize=8,
        )
    right.axvline(0.0, color="black", linewidth=0.9)
    right.set_xlim(min(values) * 1.4, max(values) * 1.25)
    right.set_title("two histories, one observation")
    return _save(
        figure,
        "fig_mediator.png",
        "Two panels. Left: the pairwise margin three epsilon minus one against epsilon, crossing "
        "zero at the threshold one third, with the point at epsilon equal to 0.2 marked at minus "
        "0.4. Right: a horizontal bar chart of the reduction control, showing that the same "
        "observation scores one, the two criteria score one and minus one after the reduction, "
        "the response distance is two, and the worst-case error floor is one.",
    )


def blanket_figure() -> Path:
    """The drift that a stationary law does not fix, and free entropy production."""
    bound = as_dict()
    drift = controls.blanket_dynamics_check(0.75)
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.4))
    names = ["exterior", "blanket", "internal"]
    values = [
        float(drift["mu_drift_external"]),
        float(drift["mu_drift_blanket"]),
        float(drift["mu_drift_internal"]),
    ]
    bars = left.bar(names, values, color=["#1f4e79", "#c98a00", "#a4262c"])
    for bar, value in zip(bars, values):
        left.annotate(
            f"{value:g}",
            (bar.get_x() + bar.get_width() / 2, value),
            textcoords="offset points",
            xytext=(0, 3 if value >= 0 else -11),
            ha="center",
            fontsize=8,
        )
    left.axhline(0.0, color="black", linewidth=0.9)
    left.set_ylabel("drift component at q = 0.75")
    left.set_title(
        f"stationary law factorises, cross drift {float(drift['cross_drift_difference']):g}"
    )

    grid = np.linspace(0.0, 1.5, 200)
    production = [controls.solenoidal_check(value)["entropy_production"] for value in grid]
    right.plot(grid, production, color="#1f4e79", linewidth=1.8)
    marks = [0.0, 0.5, 1.5]
    right.scatter(
        marks,
        [controls.solenoidal_check(value)["entropy_production"] for value in marks],
        color="#a4262c",
        zorder=3,
    )
    for value in marks:
        right.annotate(
            f"{controls.solenoidal_check(value)['entropy_production']:.3g}",
            (value, controls.solenoidal_check(value)["entropy_production"]),
            textcoords="offset points",
            xytext=(4, -10),
            fontsize=8,
            color="#a4262c",
        )
    right.set_xlabel("solenoidal parameter q")
    right.set_ylabel("entropy production")
    right.set_title("2 q^2, free while every member is stationary")
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
    sync = controls.two_node_sync_orbit(steps=500)

    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.8))
    left.plot(samples[:, 0], samples[:, 1], color="#1f4e79", linewidth=1.2)
    theta = np.linspace(0.0, 2.0 * math.pi, 400)
    for bound, style in ((measured["min_radius"], "-"), (measured["max_radius"], "--")):
        left.plot(
            bound * np.cos(theta),
            bound * np.sin(theta),
            color="#a4262c",
            linestyle=style,
            linewidth=0.9,
        )
    left.scatter([0.0], [0.0], color="black", s=18, zorder=3)
    left.set_aspect("equal")
    left.set_xlabel("x")
    left.set_ylabel("y")
    left.set_title(
        f"orbit confined to [{measured['min_radius']:.4f}, {measured['max_radius']:.4f}]"
    )

    times = np.arange(radius.size) * sample_dt
    right.plot(times, radius, color="#1f4e79", linewidth=1.2, label="relaxation")
    sync_radius = np.linalg.norm(sync, axis=1)
    right.plot(
        np.arange(sync_radius.size) * 1.0,
        sync_radius,
        color="#c98a00",
        linewidth=1.0,
        label="synchronous clipped map",
    )
    right.axhline(measured["min_radius"], color="#a4262c", linestyle="--", linewidth=0.9)
    right.annotate(
        f"period {measured['period']:.3f}",
        xy=(measured["period"], measured["min_radius"]),
        xytext=(measured["period"] + 1.2, measured["min_radius"] + 0.12),
        fontsize=8,
        color="#a4262c",
        arrowprops={"arrowstyle": "->", "color": "#a4262c", "lw": 1.0},
    )
    right.set_xlabel("time")
    right.set_ylabel("radius")
    right.set_title("radius against time")
    right.legend(fontsize=7, loc="lower right")
    return _save(
        figure,
        "fig_orbit.png",
        "Two panels. Left: the relaxation orbit in the plane, a closed curve that neither escapes "
        "nor reaches the origin, with the annulus that contains it drawn at the measured inner "
        "and outer radii. Right: the radius against time for the relaxation orbit and for the "
        "synchronous clipped map, which is confined to a different band; the relaxation period is "
        "marked.",
    )


def ivm_figure() -> Path:
    """IVM shell sizes and tetravolume ratios."""
    shells = list(range(0, 9))
    counts = [lattice.shell_count(index) for index in shells]
    cumulative = [lattice.cumulative_count(index) for index in shells]

    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.4))
    bars = left.bar([str(index) for index in shells], counts, color="#1f4e79", label="shell")
    for bar, value in zip(bars, counts):
        left.annotate(
            str(value),
            (bar.get_x() + bar.get_width() / 2, value),
            textcoords="offset points",
            xytext=(0, 3),
            ha="center",
            fontsize=7.5,
        )
    twin = left.twinx()
    twin.plot(
        [str(index) for index in shells],
        cumulative,
        color="#c98a00",
        marker="o",
        markersize=4,
        linewidth=1.4,
        label="cumulative",
    )
    twin.set_ylabel("cumulative sites", color="#c98a00")
    twin.grid(False)
    left.set_xlabel("shell k")
    left.set_ylabel("sites in shell")
    left.set_title("10k^2 + 2 sites per shell")
    left.set_ylim(0, max(counts) * 1.15)

    names = list(lattice.VOLUME_RATIOS)
    volumes = [float(lattice.VOLUME_RATIOS[name]) for name in names]
    pretty = [name.replace("_", " ") for name in names]
    order = np.argsort(volumes)
    bars = right.barh(
        [pretty[index] for index in order],
        [volumes[index] for index in order],
        color=["#2f6b3a" if volumes[index] < 10 else "#a4262c" for index in order],
    )
    right.set_xscale("log")
    for bar, index in zip(bars, order):
        right.annotate(
            f"{volumes[index]:g}",
            (volumes[index], bar.get_y() + bar.get_height() / 2),
            textcoords="offset points",
            xytext=(4, 0),
            va="center",
            fontsize=8,
        )
    right.set_xlabel("tetravolume, unit tetrahedron = 1 (log scale)")
    right.set_title("exact rational volumes")
    right.set_xlim(0.8, max(volumes) * 1.8)
    return _save(
        figure,
        "fig_ivm.png",
        "Two panels. Left: the site count of each IVM shell, ten k squared plus two, as bars for "
        "shells zero to eight, with the cumulative count overlaid on a second axis; the first "
        "shell holds the twelve neighbours. Right: tetravolume ratios on a logarithmic scale, "
        "with the tetrahedron at one, the cube at three, the octahedron at four, the rhombic "
        "dodecahedron at six, the cuboctahedron at twenty and the truncated octahedron at "
        "ninety-six.",
    )


def subdivision_figure() -> Path:
    """The same quadratic sequence in three places."""
    frequencies = list(range(1, 7))
    geodesic_vertices = [lattice.geodesic_vertex_count(value) for value in frequencies]
    geodesic_faces = [lattice.geodesic_face_count(value) for value in frequencies]
    goldberg_vertices = [lattice.goldberg_vertex_count(value) for value in frequencies]
    goldberg_faces = [lattice.goldberg_face_count(value) for value in frequencies]
    shell = [lattice.shell_count(value) for value in frequencies]

    figure, axis = plt.subplots(figsize=(8.2, 4.6))
    axis.plot(frequencies, geodesic_vertices, marker="o", color="#1f4e79", label="geodesic vertices")
    axis.plot(
        frequencies,
        goldberg_faces,
        marker="x",
        linestyle="--",
        color="#2f6b3a",
        label="Goldberg faces",
    )
    axis.plot(
        frequencies,
        shell,
        marker="s",
        linestyle=":",
        color="#c98a00",
        label="IVM shell count",
    )
    axis.plot(
        frequencies,
        geodesic_faces,
        marker="^",
        color="#a4262c",
        label="geodesic faces = Goldberg vertices",
    )
    for value in frequencies:
        axis.annotate(
            f"{lattice.geodesic_vertex_count(value)}",
            (value, lattice.geodesic_vertex_count(value)),
            textcoords="offset points",
            xytext=(0, 6),
            ha="center",
            fontsize=7,
            color="#1f4e79",
        )
    axis.set_xlabel("frequency f")
    axis.set_ylabel("count")
    axis.set_title("10f^2 + 2 and 20f^2: one duality, three appearances")
    axis.legend(fontsize=7.5, loc="upper left")
    return _save(
        figure,
        "fig_subdivision.png",
        "A line chart against subdivision frequency from one to six. Geodesic vertex count, "
        "Goldberg face count and the IVM shell count coincide on ten f squared plus two, while "
        "the geodesic face count and the Goldberg vertex count coincide on twenty f squared; the "
        "two curves are the two sides of one duality, and the count values are annotated for the "
        "geodesic vertices.",
    )


def cover_figure() -> Path:
    """Cover art: the carrier's ports and seams on a dark field."""
    carrier, flat = _projected_carrier()
    figure = plt.figure(figsize=(7.0, 7.0))
    figure.patch.set_facecolor("#0b0f14")
    axis = figure.add_subplot(111)
    axis.set_facecolor("#0b0f14")
    axis.set_axis_off()
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


def ga_triangle_figure() -> Path:
    """The three frameworks as a triangle, with each edge's ledger and its key claims."""
    figure, axis = plt.subplots(figsize=(10.0, 7.0))
    axis.set_axis_off()
    axis.set_aspect("equal")
    vertices = {
        "OPH": (0.50, 0.94),
        "FEP": (0.08, 0.10),
        "Synergetics": (0.92, 0.10),
    }
    edges = (
        ("OPH", "FEP", OPH_FEP, (0.22, 0.55), "center"),
        ("FEP", "Synergetics", FEP_SYNERGETICS, (0.50, 0.035), "center"),
        ("Synergetics", "OPH", SYNERGETICS_OPH, (0.78, 0.55), "center"),
    )
    for left, right, _pair, _position, _align in edges:
        axis.plot(
            [vertices[left][0], vertices[right][0]],
            [vertices[left][1], vertices[right][1]],
            color="#243b53",
            linewidth=1.8,
            zorder=1,
        )
    expansions = {
        "OPH": "Observer Patch Holography",
        "FEP": "Free Energy Principle",
        "Synergetics": "quadray/IVM",
    }
    for name, (x_value, y_value) in vertices.items():
        axis.scatter([x_value], [y_value], s=110, color="#102a43", zorder=3)
        axis.text(
            x_value,
            y_value,
            name,
            ha="center",
            va="center",
            fontsize=10.5,
            color="white",
            zorder=4,
        )
        if name == "OPH":
            expansion_y, vertical = y_value + 0.065, "bottom"
        else:
            expansion_y, vertical = y_value - 0.065, "top"
        axis.text(
            x_value,
            expansion_y,
            expansions[name],
            ha="center",
            va=vertical,
            fontsize=7.2,
            color="#243b53",
            zorder=4,
        )

    grouped = by_pair()
    named = [
        key
        for key in ("carrier-port-chart", "native-repair-law", "free-energy-dictionary")
        if key in by_key()
    ]
    highlighted = named if len(named) == 3 else []
    if not highlighted:
        highlighted = [
            next(
                entry.key
                for entry in grouped[pair]
                if entry.kind == "isomorphism"
            )
            for _left, _right, pair, _position, _align in edges
        ]

    for left, right, pair, position, align in edges:
        entries = grouped[pair]
        counts = Counter(entry.kind for entry in entries)
        ledger = ", ".join(
            (
                _counted(counts["isomorphism"], "isomorphism", "isomorphisms"),
                _counted(counts["disagreement"], "disagreement", "disagreements"),
                _counted(counts["complementarity"], "complementarity", "complementarities"),
            )
        )
        keys = [key for key in highlighted if by_key()[key].pair == pair]
        axis.text(
            position[0],
            position[1],
            ledger,
            ha=align,
            va="center",
            fontsize=8.2,
            color="#243b53",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.86, "pad": 1.2},
            zorder=5,
        )
        if keys:
            axis.text(
                position[0],
                position[1] - 0.045,
                "load-bearing: " + ", ".join(keys),
                ha=align,
                va="top",
                fontsize=7.4,
                color="#c98a00",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.86, "pad": 1.0},
                zorder=6,
            )

    three = grouped.get(THREE_WAY, ())
    if three:
        three_counts = Counter(entry.kind for entry in three)
        three_text = (
            f"three-way: {len(three)} "
            f"({_counted(three_counts['isomorphism'], 'isomorphism', 'isomorphisms')}, "
            f"{_counted(three_counts['disagreement'], 'disagreement', 'disagreements')}, "
            f"{_counted(three_counts['complementarity'], 'complementarity', 'complementarities')})\n"
            + ", ".join(entry.key for entry in three)
        )
    else:
        three_text = "three-way: no registered entry"
    axis.text(
        0.50,
        0.45,
        three_text,
        ha="center",
        va="center",
        fontsize=8.6,
        color="#102a43",
        bbox={"facecolor": "#f2f5f8", "edgecolor": "#c9d6e4", "pad": 2.0},
        zorder=2,
    )
    axis.set_xlim(-0.10, 1.10)
    axis.set_ylim(-0.12, 1.10)
    return _save(
        figure,
        "fig_ga_triangle.png",
        "A large triangle whose vertices are OPH, Observer Patch Holography; FEP, Free Energy "
        "Principle; and Synergetics, quadray/IVM. Each edge carries its pairwise ledger of "
        "isomorphisms, disagreements and complementarities, and the three load-bearing "
        "identifications are named beside the edge they belong to. The centre gives the three-way "
        "count and its two keys, tetravolume-arithmetic and conjugation-norm.",
    )


def ga_carrier_figure() -> Path:
    """The carrier graph beside the OPH-facing and FEP-facing invariant ledger."""
    bound = as_dict()
    carrier, flat = _projected_carrier()
    figure, (left, right) = plt.subplots(1, 2, figsize=(10.0, 7.0))
    for port in range(carrier.port_count):
        for other in carrier.neighbours(port):
            if other > port:
                left.plot(
                    [flat[port, 0], flat[other, 0]],
                    [flat[port, 1], flat[other, 1]],
                    color="#243b53",
                    linewidth=1.1,
                    zorder=2,
                )
    left.scatter(
        flat[:, 0],
        flat[:, 1],
        s=130,
        color="white",
        edgecolor="#102a43",
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
            fontsize=7.2,
            color="#102a43",
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
        ("spanning-tree number", f"{tree_number}"),
        ("dim H0", f"{cohomology['h0']}"),
        ("dim H1", f"{cohomology['h1']}"),
        (
            "group order",
            f"{len(group)} = {len(rotations)} rot + {len(reflections)} refl",
        ),
    ]
    fep_rows = [
        ("marginal variance", f"{marginal.numerator}/{marginal.denominator}"),
        ("blanket size", blanket_text),
        ("pairwise precision", f"{precision_text} on {len(precision)} pairs"),
    ]

    right.set_axis_off()
    right.set_xlim(0.0, 1.0)
    right.set_ylim(0.0, 1.0)
    right.text(0.02, 0.95, "OPH-facing", fontsize=10.5, color="#1f4e79", fontweight="bold")
    right.text(0.54, 0.95, "FEP-facing", fontsize=10.5, color="#a4262c", fontweight="bold")
    for index, (label, value) in enumerate(oph_rows):
        y_value = 0.86 - 0.11 * index
        right.text(0.02, y_value, label, fontsize=8.4, color="#243b53", va="center")
        right.text(0.44, y_value, value, fontsize=8.4, color="#1f4e79", va="center", ha="right")
    for index, (label, value) in enumerate(fep_rows):
        y_value = 0.86 - 0.11 * index
        right.text(0.54, y_value, label, fontsize=8.4, color="#243b53", va="center")
        right.text(0.98, y_value, value, fontsize=8.4, color="#a4262c", va="center", ha="right")
    return _save(
        figure,
        "fig_ga_carrier.png",
        f"Two panels. Left: the {int(bound['carrier.ports'])} ports and "
        f"{int(bound['carrier.seams'])} seams of the carrier in its planar "
        "spectral projection, with each port labelled by index. Right: a two-column ledger of "
        "invariants, OPH-facing on the left and FEP-facing on the right, listing the spanning-tree "
        "number, the dimensions of H0 and H1, the automorphism group order with its rotation and "
        "reflection split, the marginal variance as a fraction, the blanket size and the pairwise "
        "precision.",
    )


def ga_descent_figure() -> Path:
    """The descent census and the free-energy levels of the first five IVM shells."""
    bound = as_dict()
    census = controls.gated_descent_census()
    multiplicities = controls.path_multiplicity_values()
    figure, (left, right) = plt.subplots(1, 2, figsize=(10.0, 7.0))
    labels = ["strict descent", "zero descent", "no-op"]
    values = [census["strict"], census["zero"], census["noop"]]
    bars = left.bar(labels, values, color=["#2f6b3a", "#c98a00", "#7a7a7a"])
    for bar, value in zip(bars, values):
        left.annotate(
            str(value),
            (bar.get_x() + bar.get_width() / 2, bar.get_height()),
            textcoords="offset points",
            xytext=(0, 3),
            ha="center",
            fontsize=8,
        )
    left.set_ylabel("transitions")
    left.set_title(
        f"{census['states']} states x 30 seams x 2 coins = {sum(values)} transitions"
    )
    left.set_ylim(0, max(values) * 1.18)
    left.text(
        0.02,
        0.95,
        "path multiplicities: " + ", ".join(str(value) for value in multiplicities),
        transform=left.transAxes,
        fontsize=8,
        color="#243b53",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.2},
    )

    k_values = list(range(1, 6))
    origin = np.zeros(3)
    identity = np.eye(3)
    shared: dict[float, list[int]] = {}
    kl_by_k: dict[int, float] = {}
    for k in k_values:
        shell = geometry.generate_shell(k)
        levels = sorted({float(0.5 * geometry.squared_norm(site)) for site in shell})
        right.scatter([k] * len(levels), levels, color="#1f4e79", s=26, zorder=3)
        for level in levels:
            shared.setdefault(level, []).append(k)
        target = geometry.NEIGHBOR_SQUARED_DISTANCE * k * k
        site = next(site for site in shell if geometry.squared_norm(site) == target)
        kl_by_k[k] = free_energy.gaussian_kl(origin, identity, geometry.to_xyz(site), identity)
    right.plot(
        k_values,
        [kl_by_k[k] for k in k_values],
        color="#c98a00",
        linestyle="--",
        marker="o",
        linewidth=1.2,
        label="Gaussian KL 4k^2",
    )
    for level, shells in shared.items():
        if len(shells) > 1:
            right.axhline(level, color="#a4262c", linestyle=":", linewidth=0.9)
            right.annotate(
                f"level {level:g} on shells {', '.join(str(k) for k in shells)}",
                xy=(max(shells), level),
                xytext=(max(shells) + 0.12, level + 0.6),
                fontsize=7.8,
                color="#a4262c",
            )
    right.set_xlabel("shell k")
    right.set_ylabel("0.5 times squared norm")
    right.set_title("levels 36 and 52 are shared between shells")
    right.legend(fontsize=7.5, loc="upper left")
    return _save(
        figure,
        "fig_ga_descent.png",
        f"Two panels. Left: the descent census as three bars, "
        f"{int(bound['shell.strict_descent'])} strict descents, "
        f"{int(bound['shell.zero_descent'])} zero descents and {int(bound['shell.noop'])} "
        f"no-op moves, annotated with the path-multiplicity values "
        f"{_enumerated(bound['shell.multiplicities'])}. Right: the distinct free-energy levels "
        f"of the first {len(bound['ivm.shell_counts'])} IVM shells, with levels "
        f"{_enumerated(bound['shell.shared_free_energy_levels'])} visibly shared between shells "
        "and the Gaussian KL curve annotated as four k squared.",
    )


def ga_synergetics_figure() -> Path:
    """The integral tetravolume table beside the icosahedron and its Voronoi cell."""
    figure, (left, right) = plt.subplots(1, 2, figsize=(10.0, 7.0))
    names = list(geometry.VOLUME_RATIOS)
    volumes = [float(geometry.VOLUME_RATIOS[name]) for name in names]
    icosahedron = geometry.icosahedron_tetravolume()
    icosahedron_volume = float(icosahedron["tetravolume"])
    icosahedron_integral = bool(icosahedron["integral"])
    bar_names = [name.replace("_", " ") for name in names] + ["icosahedron"]
    bar_values = volumes + [icosahedron_volume]
    colours = ["#2f6b3a" if value == int(value) else "#a4262c" for value in bar_values]
    bars = left.barh(bar_names, bar_values, color=colours)
    left.set_xscale("log")
    for bar, value in zip(bars, bar_values):
        left.annotate(
            f"{value:g}",
            (value, bar.get_y() + bar.get_height() / 2),
            textcoords="offset points",
            xytext=(4, 0),
            va="center",
            fontsize=8,
        )
    left.set_xlabel("tetravolume, unit tetrahedron = 1 (log scale)")
    left.set_title("integral table solids, non-integral icosahedron")
    left.text(
        0.02,
        0.04,
        f"icosahedron integral: {icosahedron_integral}",
        transform=left.transAxes,
        fontsize=8,
        color="#a4262c" if not icosahedron_integral else "#2f6b3a",
    )

    cell = geometry.voronoi_cell()
    density = geometry.packing_density()
    cell_labels = ["vertices", "faces", "edges"]
    cell_values = [
        int(cell["vertex_count"]),
        int(cell["face_count"]),
        int(cell["edge_count"]),
    ]
    bars = right.bar(cell_labels, cell_values, color=["#1f4e79", "#2f6b3a", "#c98a00"])
    for bar, value in zip(bars, cell_values):
        right.annotate(
            str(value),
            (bar.get_x() + bar.get_width() / 2, bar.get_height()),
            textcoords="offset points",
            xytext=(0, 3),
            ha="center",
            fontsize=8,
        )
    right.set_ylabel("count")
    right.set_title("Voronoi cell of the IVM")
    right.text(
        0.02,
        0.95,
        f"packing density {float(density['density']):.6f} vs closed form "
        f"{float(density['closed_form']):.6f}",
        transform=right.transAxes,
        fontsize=8,
        color="#243b53",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.2},
    )
    return _save(
        figure,
        "fig_ga_synergetics.png",
        "Two panels. Left: logarithmic bars of the tetravolumes in the synergetic table, the "
        "tetrahedron, cube, octahedron, rhombic dodecahedron, cuboctahedron and truncated "
        "octahedron, all integral, beside the icosahedron's measured tetravolume, which is not "
        "integral. Right: the Voronoi cell's fourteen vertices, twelve faces and twenty-four "
        "edges, with the measured packing density annotated against its closed form.",
    )


def pair_triangles_figure() -> Path:
    """One small triangle for each pairwise scope, with its kind and status ledger."""
    figure, axes = plt.subplots(1, 3, figsize=(11.0, 4.5))
    pairs = (OPH_FEP, FEP_SYNERGETICS, SYNERGETICS_OPH)
    pair_edges = {
        OPH_FEP: ("OPH", "FEP"),
        FEP_SYNERGETICS: ("FEP", "Synergetics"),
        SYNERGETICS_OPH: ("Synergetics", "OPH"),
    }
    vertices = {
        "OPH": (0.50, 0.92),
        "FEP": (0.08, 0.12),
        "Synergetics": (0.92, 0.12),
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
                color="#c98a00" if highlighted else "#9aa9b8",
                linewidth=3.0 if highlighted else 1.0,
                zorder=1,
            )
        for name, (x_value, y_value) in vertices.items():
            axis.scatter([x_value], [y_value], s=55, color="#102a43", zorder=3)
            axis.text(
                x_value,
                y_value,
                name,
                ha="center",
                va="center",
                fontsize=7.2,
                color="white",
                zorder=4,
            )
        lines = [
            f"isomorphisms {counts['isomorphism']}",
            f"disagreements {counts['disagreement']}",
            f"complementarities {counts['complementarity']}",
            f"exact {statuses['exact']}",
            f"exact-on-domain {statuses['exact-on-domain']}",
            f"obstructed {statuses['obstructed']}",
            f"not-reproduced {statuses['not-reproduced']}",
        ]
        axis.text(
            0.50,
            0.48,
            "\n".join(lines),
            ha="center",
            va="center",
            fontsize=7.0,
            color="#243b53",
            zorder=2,
        )
        axis.set_xlim(0.0, 1.0)
        axis.set_ylim(0.0, 1.0)
        axis.set_aspect("equal")
        axis.set_axis_off()
        axis.set_title(pair.replace("-", " / "), fontsize=9.5)
    three = grouped.get(THREE_WAY, ())
    if three:
        figure.text(
            0.50,
            0.015,
            "three-way: " + str(len(three)) + " " + ", ".join(entry.key for entry in three),
            ha="center",
            fontsize=7.8,
            color="#243b53",
        )
    return _save(
        figure,
        "fig_pair_triangles.png",
        "Three small triangles side by side, one for each pairwise scope: OPH and FEP, FEP and "
        "Synergetics, and Synergetics and OPH. Each triangle names all three frameworks and "
        "annotates its pair's counts by kind and its status breakdown, with the pair's edge "
        "highlighted. A centre note records the three-way count and its two entries, "
        "tetravolume-arithmetic and conjugation-norm.",
    )


def scope_matrix_figure() -> Path:
    """The registry as a dense categorical matrix in registry order."""
    pair_colours = {
        OPH_FEP: "#1f4e79",
        FEP_SYNERGETICS: "#c98a00",
        SYNERGETICS_OPH: "#2f6b3a",
        THREE_WAY: "#7a7a7a",
    }
    status_colours = {
        "exact": "#2f6b3a",
        "exact-on-domain": "#c98a00",
        "obstructed": "#a4262c",
        "not-reproduced": "#7a7a7a",
    }
    probe_colours = {"own": "#1f4e79", "shared": "#7a7a7a"}

    def text_colour(hex_colour: str) -> str:
        value = hex_colour.lstrip("#")
        red, green, blue = (int(value[index : index + 2], 16) for index in (0, 2, 4))
        luminance = 0.299 * red + 0.587 * green + 0.114 * blue
        return "black" if luminance > 150 else "white"

    figure, axis = plt.subplots(figsize=(11.0, 10.0))
    axis.set_axis_off()
    axis.set_xlim(0.0, 5.70)
    axis.set_ylim(0.0, len(CORRESPONDENCES) + 1.2)
    probe_names = [entry.probe.__name__ for entry in CORRESPONDENCES]
    unique_probes = len(set(probe_names)) == len(probe_names)
    headers = ["pair", "kind", "status", "own probe"]
    for column, header in enumerate(headers):
        axis.text(
            column + 2.15,
            len(CORRESPONDENCES) + 0.55,
            header,
            ha="center",
            va="center",
            fontsize=9.5,
            fontweight="bold",
            color="#102a43",
        )
    for row, entry in enumerate(CORRESPONDENCES):
        y_value = len(CORRESPONDENCES) - row - 0.5
        axis.text(
            0.02,
            y_value,
            entry.key,
            ha="left",
            va="center",
            fontsize=9.0,
            family="monospace",
            color="#102a43",
        )
        probe_label = "own" if callable(entry.probe) and unique_probes else "shared"
        cells = (
            (pair_colours.get(entry.pair, "#7a7a7a"), entry.pair),
            (KIND_COLOURS.get(entry.kind, "#7a7a7a"), entry.kind),
            (status_colours.get(entry.status, "#7a7a7a"), entry.status),
            (probe_colours[probe_label], probe_label),
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
    axis.set_title("Registry scope matrix: pair, kind, status and probe ownership", fontsize=10)
    return _save(
        figure,
        "fig_scope_matrix.png",
        "A dense categorical matrix with one row per registry entry in registry order. The row "
        "label is the entry key, and four coloured columns give its pair, kind, status and "
        "whether the probe is the entry's own executable probe. "
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

    figure, axes = plt.subplots(2, 2, figsize=(11.5, 9.0))
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
            fontsize=7.6,
            color="#102a43",
        )
    ladder_axis.set_xticks(positions)
    ladder_axis.set_xticklabels([f"d = {step}" for step, _, _ in ladder], fontsize=8)
    ladder_axis.set_ylabel("effective resistance")
    ladder_axis.set_ylim(0.0, max(resistances) * 1.45)
    ladder_axis.set_title("1. The resistance ladder: one exact resistance per distance class")
    ladder_axis.text(
        0.03,
        0.93,
        f"Kirchhoff index {kirchhoff}",
        transform=ladder_axis.transAxes,
        fontsize=8.2,
        color=KIND_COLOURS["disagreement"],
    )

    # --- panel 2: the shells around port 0 --------------------------------
    projected_carrier, flat = _projected_carrier()
    distances = network.graph_distances(adjacency)[0]
    shell_colours = {
        0: KIND_COLOURS["disagreement"],
        1: KIND_COLOURS["isomorphism"],
        2: "#c98a00",
        3: KIND_COLOURS["complementarity"],
    }
    for port in range(projected_carrier.port_count):
        for other in projected_carrier.neighbours(port):
            if other > port:
                shells_axis.plot(
                    [flat[port, 0], flat[other, 0]],
                    [flat[port, 1], flat[other, 1]],
                    color="#9aa9b8",
                    linewidth=0.9,
                    zorder=1,
                )
    for step in range(len(distribution)):
        members = [port for port in range(projected_carrier.port_count) if distances[port] == step]
        shells_axis.scatter(
            flat[members, 0],
            flat[members, 1],
            s=190,
            color=shell_colours.get(step, "#7a7a7a"),
            edgecolor="#102a43",
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
            fontsize=7.2,
            color="white",
            zorder=4,
        )
    shells_axis.set_aspect("equal")
    shells_axis.set_axis_off()
    shells_axis.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=7.4,
        title=f"intersection array\nb = {intersection[0]}, c = {intersection[1]}",
        title_fontsize=7.4,
    )
    shells_axis.set_title("2. The shells around port 0, coloured by graph distance")

    # --- panel 3: the subgraph census -------------------------------------
    census_positions = list(range(len(census)))
    census_values = list(census.values())
    census_axis.bar(
        census_positions,
        census_values,
        width=0.6,
        color=[
            KIND_COLOURS["isomorphism"],
            "#c98a00",
            KIND_COLOURS["complementarity"],
            KIND_COLOURS["disagreement"],
        ],
    )
    for position, count in zip(census_positions, census_values):
        census_axis.text(
            position,
            count,
            str(count),
            ha="center",
            va="bottom",
            fontsize=8.0,
            color="#102a43",
        )
    census_axis.axhline(face_count, color=KIND_COLOURS["disagreement"], linestyle="--", linewidth=1.2)
    census_axis.annotate(
        f"faces {face_count}",
        xy=(census_positions[-1] + 0.3, face_count),
        xytext=(census_positions[-1] + 0.42, face_count + max(census_values) * 0.03),
        fontsize=7.8,
        color=KIND_COLOURS["disagreement"],
        ha="left",
    )
    census_axis.set_xticks(census_positions)
    census_axis.set_xticklabels(
        ["tri.", "pent.", "match.", "Ham."], fontsize=8.4
    )
    census_axis.set_ylabel("count")
    census_axis.set_ylim(0.0, max(census_values) * 1.15)
    census_axis.set_title("3. The subgraph census, against the face count")

    # --- panel 4: the dual pair -------------------------------------------
    dual_axis.set_axis_off()
    dual_axis.set_title("4. The dual pair: the Voronoi cell and the cuboctahedron")

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
                linewidth=1.1,
                zorder=2,
            )
        dual_axis.scatter(placed[:, 0], placed[:, 1], s=26, color=colour, zorder=3)
        dual_axis.text(
            centre_x,
            -0.24,
            caption,
            ha="center",
            va="top",
            fontsize=7.8,
            color="#102a43",
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
        f"Voronoi cell: {cell['vertex_count']} vertices, {cell['face_count']} rhombic\n"
        f"faces, {cell['edge_count']} edges, tetravolume {cell_tetravolume}",
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
        f"Cuboctahedron: {cuboctahedron['vertices']} vertices, "
        f"{len(cuboctahedron_face_counts['triangles'])} triangles +\n"
        f"{len(cuboctahedron_face_counts['squares'])} squares, {cuboctahedron_edge_count} edges, "
        f"tetravolume {cuboctahedron_tetravolume}",
    )
    dual_axis.set_xlim(-2.1, 2.1)
    dual_axis.set_ylim(-0.78, 1.05)
    dual_axis.text(
        0.0,
        -0.70,
        f"volume ratio {cuboctahedron_tetravolume}/{cell_tetravolume} = {volume_ratio}, exact",
        ha="center",
        fontsize=8.6,
        color="#102a43",
    )

    return _save(
        figure,
        "fig_network.png",
        "Four panels, every number read from the executed package at draw time. Panel one is the "
        "resistance ladder, one bar for each distance class of the carrier's port graph, labelled "
        "with its exact fraction and its pair count and annotated with the Kirchhoff index. Panel "
        "two draws the twelve ports in the icosahedral spectral projection, coloured by graph "
        "distance from port zero with that shell distribution as its legend and the intersection "
        "array annotated. Panel three is the subgraph census, bars for triangles, pentagons, "
        "perfect matchings and Hamiltonian cycles with the icosahedron face count drawn as a "
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
    figure, axis = plt.subplots(figsize=(8.4, 4.4))
    width = 0.38
    positions = np.arange(len(radii))
    carrier_bars = axis.bar(
        positions - width / 2,
        carrier_counts,
        width,
        color=KIND_COLOURS["isomorphism"],
        label="carrier: ports at graph distance k",
    )
    lattice_bars = axis.bar(
        positions + width / 2,
        lattice_counts,
        width,
        color=KIND_COLOURS["complementarity"],
        label="IVM lattice: sites on shell k",
    )
    axis.set_yscale("log")
    axis.set_ylim(0.8, max(lattice_counts) * 1.6)
    for bars in (carrier_bars, lattice_bars):
        for bar in bars:
            axis.annotate(
                f"{int(bar.get_height())}",
                (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                ha="center",
                va="bottom",
                fontsize=8,
            )
    axis.set_xticks(positions, [f"radius {k}" for k in radii])
    axis.set_ylabel("count (log scale)")
    axis.set_title(
        f"shells around a port: {carrier_counts} against {lattice_counts}; "
        f"b = {outward}, c = {homing}"
    )
    axis.legend(fontsize=7.5, loc="upper left")
    return _save(
        figure,
        "fig_vertex_shells.png",
        "Two bar series against the radius from a port, on a logarithmic count axis. The "
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
        ("spanning-tree number", int(values["carrier.tree_number"])),
        ("automorphisms", int(values["symmetry.group_order"])),
        ("sheaf obstructions", int(values["sheaf.h1"])),
        ("triangles", int(values["carrier.triangles"])),
    )
    volume = float(geometry.icosahedron_tetravolume()["tetravolume"])
    norm = geometry.tetravolume_galois_norm()
    conjugate = float(norm) / volume
    figure, (left, right) = plt.subplots(1, 2, figsize=(10.0, 4.4))
    names = [name for name, _ in fixed]
    counts = [value for _, value in fixed]
    bars = left.bar(names, counts, color="#1f4e79")
    left.set_yscale("log")
    left.set_ylim(1.0, max(counts) * 4.0)
    for bar, count in zip(bars, counts):
        left.annotate(
            f"{count}",
            (bar.get_x() + bar.get_width() / 2, count),
            ha="center",
            va="bottom",
            fontsize=8,
        )
    left.set_title("fixed by sqrt(5) -> -sqrt(5)")
    left.tick_params(axis="x", labelrotation=12)
    bars = right.bar(("icosahedron", "conjugate"), (volume, conjugate), color=("#a4262c", "#c98a00"))
    right.set_yscale("log")
    right.set_ylim(1.0, max(volume, conjugate) * 4.0)
    for bar, value in zip(bars, (volume, conjugate)):
        right.annotate(
            f"{value:.6g}",
            (bar.get_x() + bar.get_width() / 2, value),
            ha="center",
            va="bottom",
            fontsize=8,
        )
    right.set_title(f"moved; the product is the exact norm {norm}")
    right.set_ylabel("tetravolume (log scale)")
    figure.suptitle("the counting invariants are fixed, the volumes are not", fontsize=10)
    return _save(
        figure,
        "fig_conjugation.png",
        "Two panels on logarithmic axes. Left: the counting invariants that conjugating the "
        f"square root of five fixes: the spanning-tree number {counts[0]}, the automorphism "
        f"group order {counts[1]}, the {counts[2]} sheaf obstructions and the {counts[3]} triangles. "
        f"Right: the icosahedron's tetravolume {volume:.6g} beside its conjugate {conjugate:.6g}; "
        f"the two multiply to the exact norm {norm}.",
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
    cubo_triangles = sum(1 for face in cubo_faces.values() if len(face) == 3)
    cubo_squares = sum(1 for face in cubo_faces.values() if len(face) == 4)
    ico_triangles = len(carrier.faces)
    ico_trees = invariants.spanning_tree_number(carrier)
    cubo_trees = invariants.spanning_tree_number_of_adjacency(cubo_adjacency)
    ico_volume = float(geometry.icosahedron_tetravolume()["tetravolume"])
    cubo_volume = float(lattice.cuboctahedron_volume())
    figure, (left, right) = plt.subplots(1, 2, figsize=(10.0, 4.6))
    rows = ("degree", "edges", "triangles", "squares")
    cubo_counts = (cubo_degree, cubo_edges, cubo_triangles, cubo_squares)
    ico_counts = (ico_degree, ico_edges, ico_triangles, 0)
    positions = np.arange(len(rows))
    width = 0.38
    cubo_bars = left.barh(
        positions - width / 2,
        cubo_counts,
        width,
        color=KIND_COLOURS["isomorphism"],
        label="cuboctahedron",
    )
    ico_bars = left.barh(
        positions + width / 2,
        ico_counts,
        width,
        color=KIND_COLOURS["disagreement"],
        label="carrier (icosahedron)",
    )
    left.set_yticks(positions, rows)
    left.invert_yaxis()
    left.set_xlim(0.0, max(cubo_counts + ico_counts) * 1.35)
    for bars, values in ((cubo_bars, cubo_counts), (ico_bars, ico_counts)):
        for bar, value in zip(bars, values):
            left.annotate(
                f"{value}" if value else "none",
                (bar.get_width(), bar.get_y() + bar.get_height() / 2),
                va="center",
                ha="left",
                fontsize=7.5,
            )
    left.set_title(
        f"counts; tetravolumes {cubo_volume:.6g} against {ico_volume:.6g}",
        fontsize=10,
    )
    left.legend(fontsize=7.5, loc="lower right")
    trees = (ico_trees, cubo_trees)
    bars = right.bar(("carrier", "cuboctahedron"), trees, color=("#1f4e79", "#c98a00"))
    right.set_yscale("log")
    right.set_ylim(trees[1] / 4.0, trees[0] * 4.0)
    for bar, count in zip(bars, trees):
        right.annotate(
            f"{count:,}",
            (bar.get_x() + bar.get_width() / 2, count),
            ha="center",
            va="bottom",
            fontsize=8,
        )
    right.annotate(
        f"{math.isqrt(cubo_trees)}^2",
        (1.0, cubo_trees),
        ha="center",
        va="top",
        fontsize=7.5,
        color="#102a43",
    )
    right.set_title("spanning trees of the seam graph, exact", fontsize=10)
    right.set_ylabel("spanning trees (log scale)")
    figure.suptitle("the same twelve-around-one phrase, two different solids", fontsize=10)
    return _save(
        figure,
        "fig_families.png",
        "Two panels. Left: horizontal bars comparing the cuboctahedron with the carrier's "
        f"icosahedron: degree {cubo_degree} against {ico_degree}, edges {cubo_edges} against "
        f"{ico_edges}, triangles {cubo_triangles} against {ico_triangles}, and squares "
        f"{cubo_squares} against none, with the tetravolumes {cubo_volume:.6g} and "
        f"{ico_volume:.6g} annotated in the title. Right: the spanning trees the same "
        f"determinant counts, {ico_trees} for the carrier against {cubo_trees} for the "
        f"cuboctahedron, which is {math.isqrt(cubo_trees)} squared, on a logarithmic axis.",
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
    figure, (left, right) = plt.subplots(1, 2, figsize=(10.0, 4.2))
    left.axvspan(lower, upper, color="#dfe9f3", label="the sandwich")
    left.axvline(bayesian, color="#a4262c", linewidth=2.0)
    left.annotate(
        f"Bayesian surprise {bayesian:.7g}",
        (bayesian, 0.5),
        ha="center",
        va="center",
        fontsize=8,
        color="#a4262c",
    )
    left.set_xlim(lower - 0.08, upper + 0.08)
    for value, name in ((lower, "lower bound"), (upper, "upper bound")):
        left.axvline(value, color="#1f4e79", linewidth=1.0, linestyle="--")
        left.annotate(
            f"{name} {value:.7g}",
            (value, 0.9),
            ha="left" if value < 0.3 else "right",
            va="center",
            fontsize=7.5,
            color="#1f4e79",
        )
    left.set_ylim(0.0, 1.0)
    left.set_yticks([])
    left.set_xlabel("surprisal units")
    left.set_title("the sandwich brackets the Bayesian value", fontsize=10)
    names = ("numerical", "Bayesian", "interface", "surprisal", "variational", "VFE at posterior")
    values = (numerical, bayesian, interface, surprisal, variational, at_posterior)
    colours = ["#1f4e79", "#1f4e79", "#a4262c", "#2f6b3a", "#2f6b3a", "#a4262c"]
    bars = right.barh(names, values, color=colours)
    right.axvline(0.0, color="black", linewidth=0.8)
    for bar, value in zip(bars, values):
        right.annotate(
            f"{value:.7g}",
            (value, bar.get_y() + bar.get_height() / 2),
            ha="left" if value >= 0 else "right",
            va="center",
            fontsize=7.2,
            xytext=(4 if value >= 0 else -4, 0),
            textcoords="offset points",
        )
    right.set_xlim(min(values) * 1.55, max(values) * 1.25)
    right.annotate(
        "",
        xy=(at_posterior, 5.0),
        xytext=(interface, 2.0),
        arrowprops={"arrowstyle": "<->", "color": "#a4262c", "lw": 1.1},
    )
    right.annotate(
        f"$\\Delta$ = model constant\n{constant:.7g}",
        xy=(-1.02, 3.5),
        ha="center",
        va="center",
        fontsize=7.6,
        color="#a4262c",
    )
    right.set_title("five readings, one name", fontsize=10)
    right.tick_params(axis="y", labelsize=8)
    figure.suptitle("what the probes report about the values called free energy", fontsize=10)
    return _save(
        figure,
        "fig_sandwich.png",
        "Two panels read from the executed probes. Left: the sandwich interval with "
        f"its lower bound {lower:.7g} and upper bound {upper:.7g}, and the Bayesian surprise "
        f"{bayesian:.7g} marked inside it. Right: horizontal bars for the functionals the "
        f"literature calls free energy in the same linear-Gaussian model: the numerical "
        f"surprise {numerical:.7g}, the Bayesian surprise {bayesian:.7g}, the interface energy "
        f"minimised over the interior {interface:.7g}, the surprisal {surprisal:.7g}, the "
        f"variational free energy {variational:.7g}, and the variational free energy at the "
        f"posterior {at_posterior:.7g}; the interior and posterior values differ by the model "
        f"constant {constant:.7g}.",
    )


def schedule_figure() -> Path:
    """The centered contrast between two schedules, one exact fraction per step."""
    steps_list = list(range(1, 13))
    contrasts = [
        abs(float(controls.affine_path_contrast_exact(steps=step)))
        for step in steps_list
    ]
    exact = controls.affine_path_contrast_exact(steps=12)
    figure, axis = plt.subplots(figsize=(8.4, 4.2))
    bars = axis.bar([str(step) for step in steps_list], contrasts, color="#a4262c")
    for bar, value in zip(bars, contrasts):
        if value > 0:
            axis.annotate(
                f"{value:.4g}",
                (bar.get_x() + bar.get_width() / 2, value),
                ha="center",
                va="bottom",
                fontsize=7,
            )
    axis.annotate(
        f"after twelve steps: {exact}",
        (0.97, 0.88),
        ha="right",
        va="top",
        fontsize=8,
        color="#102a43",
        xycoords="axes fraction",
    )
    axis.set_xlabel("steps of the affine iteration")
    axis.set_ylabel("absolute centered-contrast difference")
    axis.set_title("full and half-damped schedules agree on answers, not on contrasts", fontsize=10)
    return _save(
        figure,
        "fig_schedule.png",
        "The absolute difference between the full and the half-damped schedules' centered "
        "contrasts, one exact value per step count from one to twelve. The schedules share "
        f"every answer while the contrast difference decays to {contrasts[-1]:.7g}, which is "
        f"exactly {exact} after twelve steps.",
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
