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
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from oph_fep_synergetics.analysis.registry import CORRESPONDENCES  # noqa: E402
from oph_fep_synergetics.core import controls  # noqa: E402
from oph_fep_synergetics.core.carrier import (  # noqa: E402
    OPH_ORIENTED_BASE_FACES,
    icosahedral_carrier,
)
from oph_fep_synergetics.geometry import lattice  # noqa: E402

FIGURES: Path = ROOT / "output" / "figures"

KIND_COLOURS = {
    "isomorphism": "#1f4e79",
    "disagreement": "#a4262c",
    "complementarity": "#2f6b3a",
}
STATUS_MARKERS = {"exact": "o", "exact-on-domain": "s", "obstructed": "X"}
ALT_TEXT: dict[str, str] = {}


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


def _save(figure: plt.Figure, name: str, alt: str) -> Path:
    FIGURES.mkdir(parents=True, exist_ok=True)
    target = FIGURES / name
    figure.tight_layout()
    figure.savefig(target, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    ALT_TEXT[name] = alt
    return target


def carrier_figure() -> Path:
    """The carrier, drawn in the spectral embedding of its own Laplacian."""
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
    gap = carrier.spectral_gap()
    right.annotate(
        f"normalised gap {gap:.6f}",
        xy=(lam1, 3),
        xytext=(lam1 + 1.2, 5.4),
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
        "Two panels. Left: the twelve ports of the carrier drawn in the three-dimensional "
        "spectral embedding of its own Laplacian, with the thirty seams drawn and the twenty "
        "oriented faces shaded; every port is labelled with its index and meets five others. "
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
            axis.text(0.14, row, entry.key, va="center", fontsize=7.4, family="monospace")
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
        "`oph-fep-synergetics verify`",
        fontsize=10,
    )
    return _save(
        figure,
        "fig_registry.png",
        "The twenty-three correspondences arranged in three columns by kind. Eight isomorphisms "
        "on the left, eleven disagreements in the middle, four complementarities on the right. "
        "Each row is labelled with its registry key and carries a marker for its status: a "
        "circle for exact, a square for exact-on-domain, and a cross for obstructed.",
    )


def descent_figure() -> Path:
    """The descent census and a mismatch-one excursion."""
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
        "Two panels. Left: the shell census, a bar chart of one hundred and twenty strict "
        "descents, six hundred zero descents and three thousand nine hundred and sixty no-op "
        "moves out of four thousand six hundred and eighty transitions. Right: a heat map of a "
        "seven-state excursion over the twelve ports, in which the load moves between ports "
        "while the quadratic energy stays at one half throughout.",
    )


def dstability_figure() -> Path:
    """Eigenvalue loci and the sign of the third Routh entry."""
    samples = (0.5, 1.0, 2.0, 4.0, 7.0)
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
        "Two panels. Left: the eigenvalue loci of the scaled matrix in the complex plane for five "
        "couplings, with a dashed imaginary axis; the pair that crosses into the right half plane "
        "is drawn in red and the stable cases in blue. Right: the third Routh entry two times ten "
        "minus k squared over k plus two against the coupling, crossing zero at the square root "
        "of ten and equal to minus two at k equal to four.",
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
        "Two panels. Left: the three drift components at q equal to 0.75, namely 1.5, minus 0.25 "
        "and minus 2, with the cross drift difference of 1.5 annotated; the stationary law "
        "factorises even though the drift couples the exterior and the interior. Right: entropy "
        "production against the solenoidal parameter q, rising as two q squared, with the marked "
        "values 0, 0.5 and 4.5 at q equal to 0, 0.5 and 1.5.",
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


def write_registry() -> Path:
    """Write the accessibility registry the renderer reads for figure alt text."""
    FIGURES.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": "1.2",
        "generator": "scripts/make_figures.py",
        "figures": [
            {"label": f"fig:{name.removeprefix('fig_').removesuffix('.png')}", "filename": name, "alt": alt}
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
    )
    for builder in builders:
        print(builder())
    print(write_registry())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
