"""Source-bound values quoted by the manuscript.

Every entry carries the value the prose may quote, the exact form where the
source gives one, and the source that establishes it.  ``MANUSCRIPT_VALUES`` is
the single place a number enters the prose, so a change in the code or a re-pin
of a source shows up here rather than inside a sentence.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class BoundValue:
    """One value the manuscript may quote, with its source."""

    key: str
    value: object
    source: str
    exact: str | None = None
    note: str | None = None


def _value(key: str) -> BoundValue:
    for entry in MANUSCRIPT_VALUES:
        if entry.key == key:
            return entry
    raise KeyError(key)


def as_dict() -> dict[str, object]:
    """Return the values keyed by identifier, for JSON projection."""
    return {entry.key: entry.value for entry in MANUSCRIPT_VALUES}


def format_value(value: object) -> str:
    """Format a bound value so that it reads naturally inside a sentence.

    Integers stay integers, floats keep six significant digits, sequences become
    parenthesised tuples, and mappings become ``name value`` pairs.  The point is
    that the manuscript never formats a number itself: it injects a token and the
    token carries this formatting.
    """
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.6g}"
    if isinstance(value, dict):
        return ", ".join(
            f"{str(key).replace('_', ' ')} {format_value(item)}" for key, item in value.items()
        )
    if isinstance(value, (tuple, list)):
        return "(" + ", ".join(format_value(item) for item in value) + ")"
    return str(value)


def token_for(key: str) -> str:
    """Return the ``TOKEN`` name a bound value is injected under."""
    return key.upper().replace(".", "_")


def as_tokens() -> dict[str, str]:
    """Return the manuscript token mapping, keyed by ``TOKEN`` name."""
    return {token_for(entry.key): format_value(entry.value) for entry in MANUSCRIPT_VALUES}


MANUSCRIPT_VALUES: tuple[BoundValue, ...] = (
    BoundValue("carrier.ports", 12, "OPH carrier (book/44-appendix-a.md:45)", "twelve"),
    BoundValue("carrier.seams", 30, "OPH carrier (book/44-appendix-a.md:45)", "thirty"),
    BoundValue("carrier.faces", 20, "OPH carrier (book/44-appendix-a.md:45)", "twenty"),
    BoundValue(
        "carrier.degree",
        5,
        "OPH carrier: 'Every port meets five triangles' (book/44-appendix-a.md:45)",
        "five",
    ),
    BoundValue(
        "carrier.lambda1",
        2.763932022500210,
        "OPH carrier Laplacian spectrum",
        "5 - sqrt(5)",
    ),
    BoundValue(
        "carrier.shell_gap",
        0.04606553370833684,
        "'Agreement and Surprise' Section 6.2: the spectral gap on the shell and the full carrier",
        "(5 - sqrt(5)) / 60",
    ),
    BoundValue(
        "shell.total_two",
        78,
        "OPH model/observer/native.py (states) and the paper Section 6.2",
        "C(13, 11)",
    ),
    BoundValue("shell.binary", 66, "OPH model/observer/native.py (states)", None),
    BoundValue(
        "shell.transitions",
        4680,
        "78 configurations x 30 seams x 2 coins",
        "78 * 30 * 2",
    ),
    BoundValue(
        "shell.strict_descent",
        120,
        "gated_descent_census over the 78-state shell",
        None,
        note="broken-seam repairs that strictly lower the quadratic energy",
    ),
    BoundValue(
        "shell.zero_descent",
        600,
        "gated_descent_census over the 78-state shell",
        None,
        note="mismatch-one moves that leave the quadratic energy unchanged; the executable content of PR-01",
    ),
    BoundValue(
        "shell.noop",
        3960,
        "gated_descent_census over the 78-state shell",
        None,
    ),
    BoundValue(
        "shell.multiplicities",
        (1, 2, 50, 52),
        "path_multiplicity_census; the paper's 'm in {1, 2, 50, 52}'",
        "{1, 2, 50, 52}",
        note="native transition multiplicities; the conditional surprisal of a path is a sum of log(60/m) terms",
    ),
    BoundValue(
        "schedule.affine_contrast",
        0.023040849975,
        "'Agreement and Surprise' Section 4.4",
        "921633999 / 40000000000",
        note="the centered-contrast difference between full and half-damped iteration after twelve steps",
    ),
    BoundValue(
        "d_stability.routh_k4",
        (1.0, 6.0, -2.0, 24.0),
        "'Agreement and Surprise' Section 4.1",
        "(1, 6, -2, 24)",
        note="the Routh first column at k = 4; two right-half-plane roots",
    ),
    BoundValue(
        "d_stability.lyapunov_determinant",
        4.2592592592592595,
        "'Agreement and Surprise' Proposition 10",
        "115 / 27",
    ),
    BoundValue(
        "two_node.relaxation_min_radius",
        0.84,
        "'Agreement and Surprise' Proposition 15: 'at distance at least 0.84'",
        None,
        note="measured 0.8401 for the relaxation orbit; the source reports the bound, this package the measurement",
    ),
    BoundValue(
        "two_node.relaxation_period",
        7.58,
        "'Agreement and Surprise' Proposition 15: 'period 7.58 in units of the relaxation time'",
        None,
        note="measured 7.575 at dt = 1e-3",
    ),
    BoundValue(
        "two_node.sync_min_radius",
        1.0,
        "two_node_orbit_min_radius",
        None,
        note="the synchronous clipped map is a different object from the relaxation orbit",
    ),
    BoundValue(
        "ivm.neighbor_moves",
        12,
        "QuadMath src/quadmath/core/quadray.py: IVM_NEIGHBOR_MOVES",
        "the distinct permutations of (2, 1, 1, 0)",
    ),
    BoundValue(
        "ivm.neighbor_squared_distance",
        8,
        "QuadMath src/quadmath/lattice/omni_numbering.py",
        None,
        note="center-to-neighbor distance 2*sqrt(2) in the scale-one embedding",
    ),
    BoundValue(
        "ivm.shell_counts",
        (1, 12, 42, 92, 162),
        "QuadMath shell_count(k) = 10k^2 + 2",
        "10k^2 + 2",
    ),
    BoundValue(
        "ivm.cumulative_counts",
        (1, 13, 55, 147, 309),
        "QuadMath cumulative_count",
        None,
    ),
    BoundValue(
        "ivm.volume_ratios",
        {"tetrahedron": 1, "cube": 3, "octahedron": 4, "rhombic_dodecahedron": 6, "cuboctahedron": 20, "truncated_octahedron": 96},
        "QuadMath docs/manuscript/03_quadray_methods.md",
        None,
        note="unit IVM tetrahedron = 1",
    ),
    BoundValue(
        "mediator.threshold",
        0.3333333333333333,
        "'Agreement and Surprise' Proposition 8",
        "1/3",
        note="a pairwise representation of the reduced boundary exists exactly when epsilon >= 1/3",
    ),
    BoundValue(
        "directed_gradient.direct",
        (1.0, 0.0),
        "'Agreement and Surprise' Appendix B",
        "-J^{-1} c",
    ),
    BoundValue(
        "directed_gradient.adjoint",
        (1.0, 2.0),
        "'Agreement and Surprise' Appendix B",
        "-J^{-T} c",
    ),
    BoundValue(
        "fep_formal.topics",
        168,
        "FEP_Formal config/catalogue_metadata.yaml, roster fep-001 .. fep-168",
        "168",
    ),
    BoundValue(
        "paper.five_surprises",
        ("numerical", "bayesian", "surprisal", "critic", "meaningful"),
        "'Agreement and Surprise' Table 7",
        None,
        note="the first three coincide up to curvature in the linear Gaussian corner; the last two need a task",
    ),
    BoundValue(
        "d_stability.routh_onset",
        3.1622776601683795,
        "'Agreement and Surprise' Section 4.1 characteristic polynomial, derived here",
        "sqrt(10)",
        note="the third Routh entry is 2(10 - k^2)/(k + 2), so for k > 0 stability holds exactly when k < sqrt(10)",
    ),
    BoundValue(
        "two_node.relaxation_min_radius_measured",
        0.8400707833642703,
        "two_node_relaxation_orbit, measured in this package",
        None,
        note="the source states the bound 0.84; this is the measurement the prose reports beside it",
    ),
    BoundValue(
        "two_node.relaxation_max_radius",
        1.0533912226419835,
        "two_node_relaxation_orbit, measured in this package",
        None,
        note="the outer radius of the annulus that contains the relaxation orbit",
    ),
    BoundValue(
        "two_node.relaxation_period_measured",
        7.575,
        "two_node_relaxation_orbit, measured in this package at dt = 1e-3",
        None,
        note="the source states 7.58; this is the measurement the prose reports beside it",
    ),
    BoundValue("registry.entries", 23, "analysis/registry.py CORRESPONDENCES", "23"),
    BoundValue("registry.isomorphisms", 8, "analysis/registry.py CORRESPONDENCES", "8"),
    BoundValue("registry.disagreements", 11, "analysis/registry.py CORRESPONDENCES", "11"),
    BoundValue("registry.complementarities", 4, "analysis/registry.py CORRESPONDENCES", "4"),
)

MANUSCRIPT_TOKENS: dict[str, str] = as_tokens()
