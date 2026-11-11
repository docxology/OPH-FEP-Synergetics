"""Source-bound values quoted by the manuscript.

Every entry carries the value the prose may quote, the exact form where the
source gives one, and the source that establishes it.  ``MANUSCRIPT_VALUES`` is
the single place a number enters the prose, so a change in the code or a re-pin
of a source shows up here rather than inside a sentence.

Where a value can be *derived* from the package it is derived rather than
restated: the registry counts are read from the registry, the tetravolume from
the polyhedron it measures, the packing density from the Voronoi cell.  A
transcribed copy is a number that can silently stop matching its source.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

from ..analysis import registry as registry_module
from ..core import cohomology_dimensions, controls as control_module, icosahedral_carrier
from ..core import invariants as invariant_module
from ..core import network as network_module
from ..core import symmetry as symmetry_module
from ..core.carrier import OPH_NEIGHBOURS
from ..core.sheaf import carrier_sheaf
from ..geometry import (
    VOLUME_RATIOS,
    angle_defects,
    cuboctahedron_adjacency,
    cuboctahedron_edges,
    cuboctahedron_faces,
    cuboctahedron_volume,
    generate_shell,
    icosahedron_tetravolume,
    neighbor_moves,
    packing_density,
    shell_count,
    squared_norm,
    tetravolume_galois_norm,
    tetravolume_ratio_squared_exact,
    voronoi_cell,
)


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
        return f"{value:.7g}"
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


def _carrier() -> object:
    return icosahedral_carrier()


def _registry_counts() -> dict[str, int]:
    """Count the registry by kind and by pairwise scope, from the registry itself."""
    counts = {"entries": len(registry_module.CORRESPONDENCES)}
    for kind in (
        registry_module.ISOMORPHISM,
        registry_module.DISAGREEMENT,
        registry_module.COMPLEMENTARITY,
    ):
        counts[kind] = sum(1 for entry in registry_module.CORRESPONDENCES if entry.kind == kind)
    for pair in registry_module.PAIRS:
        counts[pair] = sum(1 for entry in registry_module.CORRESPONDENCES if entry.pair == pair)
    return counts


def _shell_free_energy_levels() -> tuple[int, ...]:
    """Count the distinct values of ``||q||^2 / 2`` on each of the first five shells."""
    return tuple(len({squared_norm(site) for site in generate_shell(k)}) for k in range(1, 6))


def _shared_free_energy_levels() -> tuple[float, ...]:
    """Return the free-energy levels that occur on more than one of the first five shells."""
    shells: dict[float, set[int]] = {}
    for k in range(1, 6):
        for site in generate_shell(k):
            shells.setdefault(0.5 * float(squared_norm(site)), set()).add(k)
    return tuple(sorted(level for level, carriers in shells.items() if len(carriers) > 1))


def _network_adjacency() -> list[list[int]]:
    """Return the carrier's adjacency as exact integers for the network readings."""
    return [[int(value) for value in row] for row in icosahedral_carrier().adjacency_matrix()]


def _resistance_readings() -> dict[str, object]:
    """Read the resistance ladder, the pair counts and the Kirchhoff index."""
    adjacency = _network_adjacency()
    ladder = network_module.resistance_ladder(adjacency)
    return {
        "adjacent": str(ladder[0][1]),
        "distance_two": str(ladder[1][1]),
        "antipodal": str(ladder[2][1]),
        "pair_counts": tuple(pairs for _, _, pairs in ladder),
        "kirchhoff_index": int(network_module.kirchhoff_index(adjacency)),
    }


def _shell_readings() -> dict[str, object]:
    """Read the carrier's shells around a port against the lattice's own shells."""
    adjacency = _network_adjacency()
    shells = network_module.distance_distribution(adjacency)
    onward, back = network_module.intersection_array(adjacency)
    return {
        "carrier": shells[1:],
        "lattice": tuple(shell_count(step) for step in range(1, len(shells))),
        "diameter": len(shells) - 1,
        "intersection": "(" + ", ".join(str(value) for value in onward) + "; "
        + ", ".join(str(value) for value in back) + ")",
    }


def _defect_readings() -> dict[str, object]:
    """Read the angle defect of the carrier's triangulation, in units of pi."""
    carrier = icosahedral_carrier()
    defects = angle_defects(carrier.faces, carrier.port_count)
    return {
        "per_port": str(defects["per_vertex"][0]),
        "total": str(defects["total"]),
        "euler_characteristic_times_two": defects["euler_characteristic_times_two"],
    }


def _cuboctahedron_readings() -> dict[str, object]:
    """Read the counts that separate the cuboctahedron from the carrier."""
    adjacency = [list(row) for row in cuboctahedron_adjacency()]
    faces = cuboctahedron_faces()
    return {
        "spanning_trees": invariant_module.spanning_tree_number_of_adjacency(adjacency),
        "triangles": len(faces["triangles"]),
        "squares": len(faces["squares"]),
        "edges": len(cuboctahedron_edges()),
        "matchings": network_module.subgraph_census(adjacency)["perfect_matchings"],
    }


def _volume_ratio() -> str:
    """Return the dual pair's tetravolume ratio, exactly.

    The cell is the rhombic dodecahedron, whose tetravolume the declared table
    gives as six, so the ratio against the cuboctahedron's twenty is 10/3 and is
    read from the table rather than transcribed.
    """
    return str(Fraction(cuboctahedron_volume()) / VOLUME_RATIOS["rhombic_dodecahedron"])


def _conjugate_tetravolume() -> float:
    """Return the tetravolume of the solid whose squared volume is the conjugate.

    Conjugating ``sqrt 5`` in ``175 + 75 sqrt 5`` gives a positive square, and its
    root is the volume of the conjugate solid. The two volumes multiply to the
    exact norm 50, which is what makes the conjugation a usable check rather than
    a curiosity.
    """
    first, second = tetravolume_ratio_squared_exact()
    return math.sqrt(float(first) - float(second) * math.sqrt(5.0))


_COUNTS = _registry_counts()
_COHOMOLOGY = cohomology_dimensions(carrier_sheaf(icosahedral_carrier()))
_EQUIVARIANCE = symmetry_module.law_equivariance()
_ICOSAHEDRON = icosahedron_tetravolume()
_DENSITY = packing_density()
_CELL = voronoi_cell()
_SHELL_LEVELS = _shell_free_energy_levels()
_SHARED_LEVELS = _shared_free_energy_levels()
#: The network readings each solve exact rational systems, so they are read once.
_RESISTANCE = _resistance_readings()
_SHELLS = _shell_readings()
_DEFECTS = _defect_readings()
_CUBOCTAHEDRON = _cuboctahedron_readings()
_CENSUS = network_module.subgraph_census(_network_adjacency())


def _probe_values(key: str) -> dict[str, object]:
    """Run one registered probe; probe-derived values quote what verify reports."""
    result = registry_module.by_key()[key].run()
    if not result.passed:
        raise AssertionError(f"probe {key} failed while binding manuscript values")
    return dict(result.values)


def _chart_values() -> tuple[BoundValue, ...]:
    """The derived port chart, one bound value per port, for the adjacency table."""
    return tuple(
        BoundValue(
            f"carrier.chart_{index}",
            ", ".join(str(port) for port in row),
            "core/carrier.py OPH_NEIGHBOURS: the port chart derived from OPH_ORIENTED_BASE_FACES",
            None,
            note=f"the five neighbours of port {index}",
        )
        for index, row in enumerate(OPH_NEIGHBOURS)
    )


#: Probe outputs the prose quotes directly. Running the probes here keeps the
#: bound values and the verified values from drifting apart: a change in a
#: probe changes the prose at the next render or fails the suite.
_VE = _probe_values("variational-equivalence")
_FE = _probe_values("free-energy-dictionary")
_COLLISION = _probe_values("free-energy-name-collision")
_CERT = _probe_values("certificate-species")
_SPECTRUM = _probe_values("spectral-gap")
_PR01 = _probe_values("pr01-obstruction")
_BLANKET = _probe_values("stationary-blanket")
_SOLENODIAL = _probe_values("solenoidal-freedom")
_DSTAB = _probe_values("d-stability")
_MEDIATOR = _probe_values("mediator-obstruction")
_REDUCTION = _probe_values("reduction-criterion")
_FOURCYCLE = _probe_values("four-cycle-local-minimum")
_EP = _probe_values("equilibrium-propagation")
_REPLICATOR = _probe_values("replicator-potential")
_GMRF = _probe_values("gmrf-independence")
_EIGEN = _probe_values("eigenspace-irreducibility")


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
    BoundValue(
        "registry.entries",
        _COUNTS["entries"],
        "analysis/registry.py CORRESPONDENCES",
        str(_COUNTS["entries"]),
    ),
    BoundValue(
        "registry.isomorphisms",
        _COUNTS[registry_module.ISOMORPHISM],
        "analysis/registry.py CORRESPONDENCES",
        str(_COUNTS[registry_module.ISOMORPHISM]),
    ),
    BoundValue(
        "registry.disagreements",
        _COUNTS[registry_module.DISAGREEMENT],
        "analysis/registry.py CORRESPONDENCES",
        str(_COUNTS[registry_module.DISAGREEMENT]),
    ),
    BoundValue(
        "registry.complementarities",
        _COUNTS[registry_module.COMPLEMENTARITY],
        "analysis/registry.py CORRESPONDENCES",
        str(_COUNTS[registry_module.COMPLEMENTARITY]),
    ),
    BoundValue(
        "registry.oph_fep",
        _COUNTS[registry_module.OPH_FEP],
        "analysis/registry.py CORRESPONDENCES, pairwise scope oph-fep",
        None,
    ),
    BoundValue(
        "registry.fep_synergetics",
        _COUNTS[registry_module.FEP_SYNERGETICS],
        "analysis/registry.py CORRESPONDENCES, pairwise scope fep-synergetics",
        None,
    ),
    BoundValue(
        "registry.synergetics_oph",
        _COUNTS[registry_module.SYNERGETICS_OPH],
        "analysis/registry.py CORRESPONDENCES, pairwise scope synergetics-oph",
        None,
    ),
    BoundValue(
        "registry.three_way",
        _COUNTS[registry_module.THREE_WAY],
        "analysis/registry.py CORRESPONDENCES, pairwise scope three-way",
        None,
    ),
    BoundValue(
        "carrier.tree_number",
        invariant_module.spanning_tree_number(),
        "core/invariants.py, Kirchhoff: the exact reduced-Laplacian determinant",
        "5184000",
    ),
    BoundValue(
        "gmrf.marginal_variance",
        float(invariant_module.marginal_variance()),
        "core/invariants.py, the reciprocal spectrum summed exactly in the rationals",
        "7/36",
    ),
    BoundValue(
        "gmrf.blanket_size",
        invariant_module.blanket_sizes()[0],
        "core/invariants.py: the Markov blanket of a port is its neighbours",
        "five",
    ),
    BoundValue(
        "gmrf.pairwise_precision",
        min(invariant_module.pairwise_precision().values()),
        "core/invariants.py: every seam's pairwise precision entry",
        "-1",
    ),
    BoundValue(
        "sheaf.coboundary_rank",
        _COHOMOLOGY["coboundary_rank"],
        "core/sheaf.py cohomology_dimensions, the rank of the coboundary",
        None,
    ),
    BoundValue(
        "sheaf.h0",
        _COHOMOLOGY["h0"],
        "core/sheaf.py: the global sections of the carrier's scalar sheaf",
        "one",
    ),
    BoundValue(
        "sheaf.h1",
        _COHOMOLOGY["h1"],
        "core/sheaf.py: edge assignments modulo coboundaries",
        None,
    ),
    BoundValue(
        "sheaf.euler_characteristic",
        _COHOMOLOGY["euler_characteristic"],
        "core/sheaf.py: twelve nodes minus thirty edges",
        "12 - 30",
    ),
    BoundValue(
        "symmetry.group_order",
        len(symmetry_module.carrier_automorphisms()),
        "core/symmetry.py carrier_automorphisms, by exhaustive backtracking",
        "120",
    ),
    BoundValue(
        "symmetry.rotations",
        len(symmetry_module.rotation_subgroup()),
        "core/symmetry.py rotation_subgroup, from OPH's oriented faces",
        "60",
    ),
    BoundValue(
        "symmetry.reflections",
        len(symmetry_module.reflection_subgroup()),
        "core/symmetry.py reflection_subgroup",
        "60",
    ),
    BoundValue(
        "symmetry.configurations",
        len(icosahedral_carrier().shell_configurations(total=2)),
        "OPH shell configurations at total load two",
        "C(13, 11)",
    ),
    BoundValue(
        "symmetry.orbits",
        _EQUIVARIANCE["orbits"],
        "core/symmetry.py law_equivariance: orbits of configurations under the group",
        None,
    ),
    BoundValue(
        "symmetry.outcome_multisets",
        _EQUIVARIANCE["distinct_outcome_multisets"],
        "core/symmetry.py law_equivariance: distinct descent-outcome multisets",
        None,
    ),
    BoundValue(
        "icosahedron.tetravolume",
        _ICOSAHEDRON["tetravolume"],
        "geometry/solids.py icosahedron_tetravolume, measured from the oriented surface",
        "(5/2)(3 sqrt 2 + sqrt 10)",
    ),
    BoundValue(
        "icosahedron.tetravolume_squared",
        _ICOSAHEDRON["tetravolume_squared"],
        "geometry/solids.py, the measured square against the exact one",
        "175 + 75 sqrt 5",
    ),
    BoundValue(
        "icosahedron.faces_reoriented",
        _ICOSAHEDRON["faces_reoriented"],
        "geometry/solids.py orient_faces_outward: faces enumerated without an orientation",
        None,
    ),
    BoundValue(
        "icosahedron.degree",
        icosahedral_carrier().degree(0),
        "OPH carrier: every port meets five triangles",
        "five",
    ),
    BoundValue(
        "cuboctahedron.degree",
        max(sum(row) for row in cuboctahedron_adjacency()),
        "geometry/lattice.py cuboctahedron_adjacency: every vertex is four-valent",
        "four",
    ),
    BoundValue(
        "cuboctahedron.tetravolume",
        float(VOLUME_RATIOS["cuboctahedron"]),
        "QuadMath polyhedra table read through geometry.VOLUME_RATIOS",
        "20",
    ),
    BoundValue(
        "tetravolume.integral_solids",
        tuple(int(value) for value in VOLUME_RATIOS.values()),
        "geometry.VOLUME_RATIOS: the tetravolumes of the IVM's integral solids",
        "1, 3, 4, 6, 20, 96",
    ),
    BoundValue(
        "voronoi.vertices",
        _CELL["vertex_count"],
        "geometry/lattice.py voronoi_cell, the dual of the twelve neighbours",
        None,
    ),
    BoundValue(
        "voronoi.faces",
        _CELL["face_count"],
        "geometry/lattice.py voronoi_cell",
        None,
    ),
    BoundValue(
        "voronoi.edges",
        _CELL["edge_count"],
        "geometry/lattice.py voronoi_cell",
        None,
    ),
    BoundValue(
        "packing.density",
        _DENSITY["density"],
        "geometry/lattice.py packing_density, measured from the cell and the ball",
        "pi / (3 sqrt 2)",
    ),
    BoundValue(
        "packing.cell_tetravolume",
        _DENSITY["cell_tetravolume"],
        "geometry/lattice.py: the cell's volume in tetrahedra of the prime vector",
        "6",
    ),
    BoundValue(
        "packing.ball_tetravolume",
        _DENSITY["ball_tetravolume"],
        "geometry/lattice.py: the inscribed ball's volume in the same tetrahedra",
        "pi sqrt 2",
    ),
    BoundValue(
        "shell.free_energy_level_counts",
        _SHELL_LEVELS,
        "geometry.generate_shell: distinct values of half the squared norm per shell",
        None,
    ),
    BoundValue(
        "shell.shared_free_energy_levels",
        _SHARED_LEVELS,
        "geometry.generate_shell: free-energy levels carried by more than one shell",
        None,
    ),
    BoundValue(
        "kl.one_step",
        0.5 * float(squared_norm(neighbor_moves()[0])),
        "core/free_energy.py gaussian_kl at the IVM's first neighbour",
        "half of 8",
    ),
    BoundValue(
        "kl.step_squared_distance",
        float(squared_norm(neighbor_moves()[0])),
        "geometry.NEIGHBOR_SQUARED_DISTANCE, the IVM's squared neighbour distance",
        "8",
    ),
    BoundValue(
        "resistance.adjacent",
        _RESISTANCE["adjacent"],
        "core/network.py resistance_ladder: exact rational solve of L x = e_i - e_j on a seam",
        "11/30",
    ),
    BoundValue(
        "resistance.distance_two",
        _RESISTANCE["distance_two"],
        "core/network.py resistance_ladder: ports at graph distance two",
        "7/15",
    ),
    BoundValue(
        "resistance.antipodal",
        _RESISTANCE["antipodal"],
        "core/network.py resistance_ladder: ports at graph distance three, across the carrier",
        "1/2",
    ),
    BoundValue(
        "resistance.pair_counts",
        _RESISTANCE["pair_counts"],
        "core/network.py resistance_ladder: how many port pairs each distance class holds",
        None,
    ),
    BoundValue(
        "resistance.kirchhoff_index",
        _RESISTANCE["kirchhoff_index"],
        "core/network.py kirchhoff_index: the 66 resistances summed, against 12 times 7/3",
        "28",
    ),
    BoundValue(
        "carrier.shell_sizes",
        _SHELLS["carrier"],
        "core/network.py distance_distribution: shells of ports around one port",
        None,
    ),
    BoundValue(
        "lattice.shell_sizes",
        _SHELLS["lattice"],
        "geometry.shell_count: the IVM's own shell sizes for the same radii",
        None,
    ),
    BoundValue(
        "carrier.diameter",
        _SHELLS["diameter"],
        "core/network.py graph_distances: the carrier's graph diameter",
        "three",
    ),
    BoundValue(
        "carrier.intersection_array",
        _SHELLS["intersection"],
        "core/network.py intersection_array: neighbours outward and back, shell by shell",
        None,
    ),
    BoundValue(
        "carrier.triangles",
        _CENSUS["triangles"],
        "core/network.py subgraph_census: three-cycles, exactly the twenty faces",
        "twenty",
    ),
    BoundValue(
        "carrier.pentagons",
        _CENSUS["pentagons"],
        "core/network.py subgraph_census: five-cycles, the dual dodecahedron's vertices",
        "twelve",
    ),
    BoundValue(
        "carrier.perfect_matchings",
        _CENSUS["perfect_matchings"],
        "core/network.py subgraph_census: perfect matchings of the custody graph",
        "125",
    ),
    BoundValue(
        "carrier.hamiltonian_cycles",
        _CENSUS["hamiltonian_cycles"],
        "core/network.py subgraph_census: labelled Hamiltonian cycles",
        None,
    ),
    BoundValue(
        "carrier.defect_per_port",
        _DEFECTS["per_port"],
        "geometry/solids.py angle_defects: the defect at a port, in units of pi",
        "1/3 pi",
    ),
    BoundValue(
        "carrier.total_defect",
        _DEFECTS["total"],
        "geometry/solids.py angle_defects: the total defect against 2 pi times chi",
        "4 pi",
    ),
    BoundValue(
        "cuboctahedron.spanning_trees",
        _CUBOCTAHEDRON["spanning_trees"],
        "core/invariants.py spanning_tree_number_of_adjacency: the cuboctahedron's own count",
        "576 squared",
    ),
    BoundValue(
        "cuboctahedron.spanning_trees_root",
        math.isqrt(int(_CUBOCTAHEDRON["spanning_trees"])),
        "core/invariants.py: the square root of the cuboctahedron's spanning-tree count",
        None,
    ),
    BoundValue(
        "cuboctahedron.triangles",
        _CUBOCTAHEDRON["triangles"],
        "geometry/lattice.py cuboctahedron_faces: the eight triangles",
        "eight",
    ),
    BoundValue(
        "cuboctahedron.squares",
        _CUBOCTAHEDRON["squares"],
        "geometry/lattice.py cuboctahedron_faces: the six squares",
        "six",
    ),
    BoundValue(
        "cuboctahedron.matchings",
        _CUBOCTAHEDRON["matchings"],
        "core/network.py subgraph_census on the cuboctahedron",
        None,
    ),
    BoundValue(
        "dual.volume_ratio",
        _volume_ratio(),
        "geometry.VOLUME_RATIOS: the cuboctahedron's twenty against the cell's six",
        "10/3",
    ),
    BoundValue(
        "tetravolume.galois_norm",
        str(tetravolume_galois_norm()),
        "geometry/solids.py tetravolume_galois_norm: T times its conjugate, exactly",
        "50",
    ),
    BoundValue(
        "icosahedron.tetravolume_conjugate",
        _conjugate_tetravolume(),
        "geometry/solids.py: the conjugate solid's tetravolume, from the conjugated square",
        None,
    ),
    BoundValue(
        "gmrf.reciprocal_spectrum_sum",
        str(invariant_module.reciprocal_spectrum_sum(icosahedral_carrier())),
        "core/invariants.py reciprocal_spectrum_sum: the exact sum of reciprocal eigenvalues",
        "7/3",
    ),
    BoundValue(
        "cuboctahedron.vertices",
        len(cuboctahedron_adjacency()),
        "geometry/lattice.py cuboctahedron_adjacency: the twelve neighbouring sites",
        "twelve",
    ),
    BoundValue(
        "cuboctahedron.edges",
        _CUBOCTAHEDRON["edges"],
        "geometry/lattice.py cuboctahedron_edges: the cube's line graph's edges",
        "twenty-four",
    ),
    BoundValue(
        "cuboctahedron.faces",
        int(_CUBOCTAHEDRON["triangles"]) + int(_CUBOCTAHEDRON["squares"]),
        "geometry/lattice.py cuboctahedron_faces: eight triangles and six squares",
        "fourteen",
    ),
    # --- probe-derived values: the prose quotes what the probes report ---
    BoundValue(
        "spectral.lambda_max",
        float(_SPECTRUM["spectrum"][-1]),
        "probe:spectral-gap:spectrum:max",
        "5 + sqrt(5)",
    ),
    BoundValue(
        "certificate.second_eigenvalue",
        float(_CERT["second_eigenvalue"]),
        "probe:certificate-species:second_eigenvalue",
        "1 - (5 - sqrt(5)) / 60",
    ),
    BoundValue(
        "certificate.infinity_norm",
        float(_CERT["induced_infinity_norm"]),
        "probe:certificate-species:induced_infinity_norm",
        "1",
    ),
    BoundValue(
        "ve.iterations",
        int(_VE["iterations"]),
        "probe:variational-equivalence:iterations",
        None,
        note="the perceive fixed-point iteration's steps to the harmonic extension",
    ),
    BoundValue(
        "ve.residual",
        float(_VE["residual"]),
        "probe:variational-equivalence:residual",
        None,
    ),
    BoundValue(
        "ve.excess",
        float(_VE["excess_free_energy_at_posterior"]),
        "probe:variational-equivalence:excess_free_energy_at_posterior",
        "0",
    ),
    BoundValue(
        "fe.numerical",
        float(_FE["numerical_surprise"]),
        "probe:free-energy-dictionary:numerical_surprise",
        None,
    ),
    BoundValue(
        "fe.bayesian",
        float(_FE["bayesian_surprise"]),
        "probe:free-energy-dictionary:bayesian_surprise",
        None,
    ),
    BoundValue(
        "fe.surprisal",
        float(_FE["surprisal"]),
        "probe:free-energy-dictionary:surprisal",
        None,
    ),
    BoundValue(
        "fe.variational",
        float(_FE["variational_free_energy"]),
        "probe:free-energy-dictionary:variational_free_energy",
        None,
    ),
    BoundValue(
        "fe.sandwich_lower",
        float(_FE["sandwich_lower"]),
        "probe:free-energy-dictionary:sandwich_lower",
        None,
    ),
    BoundValue(
        "fe.sandwich_upper",
        float(_FE["sandwich_upper"]),
        "probe:free-energy-dictionary:sandwich_upper",
        None,
    ),
    BoundValue(
        "fe.interface",
        float(_COLLISION["interface_energy_minimised_over_interior"]),
        "probe:free-energy-name-collision:interface_energy_minimised_over_interior",
        None,
    ),
    BoundValue(
        "fe.variational_at_posterior",
        float(_COLLISION["variational_free_energy_at_posterior"]),
        "probe:free-energy-name-collision:variational_free_energy_at_posterior",
        None,
    ),
    BoundValue(
        "fe.difference",
        float(_COLLISION["difference"]),
        "probe:free-energy-name-collision:difference",
        None,
    ),
    BoundValue(
        "fe.model_constant",
        float(_COLLISION["expected_model_constant"]),
        "probe:free-energy-name-collision:expected_model_constant",
        None,
    ),
    BoundValue(
        "pr01.coin_zero_current",
        float(_PR01["coin_zero_current"]),
        "probe:pr01-obstruction:coin_zero_current",
        None,
    ),
    BoundValue(
        "pr01.coin_one_current",
        float(_PR01["coin_one_current"]),
        "probe:pr01-obstruction:coin_one_current",
        None,
    ),
    BoundValue(
        "pr01.coin_zero_decrease",
        float(_PR01["coin_zero_decrease"]),
        "probe:pr01-obstruction:coin_zero_decrease",
        "0",
    ),
    BoundValue(
        "pr01.coin_one_decrease",
        float(_PR01["coin_one_decrease"]),
        "probe:pr01-obstruction:coin_one_decrease",
        "0",
    ),
    BoundValue(
        "pr01.walk_energy",
        float(_PR01["excursion_energy"]),
        "probe:pr01-obstruction:excursion_energy",
        "1/2",
    ),
    BoundValue(
        "blanket.q",
        float(_BLANKET["q"]),
        "probe:stationary-blanket:q",
        None,
    ),
    BoundValue(
        "blanket.drift",
        (
            float(_BLANKET["mu_drift_external"]),
            float(_BLANKET["mu_drift_blanket"]),
            float(_BLANKET["mu_drift_internal"]),
        ),
        "probe:stationary-blanket:mu_drift:row",
        None,
    ),
    BoundValue(
        "gmrf.marginal_variance_exact",
        "7/36",
        "core/invariants.py marginal_variance: the exact rational marginal variance of a port",
        "7/36",
    ),
    BoundValue(
        "blanket.cross_drift",
        float(_BLANKET["cross_drift_difference"]),
        "probe:stationary-blanket:cross_drift_difference",
        None,
    ),
    BoundValue(
        "solenoidal.q_grid",
        tuple(float(q) for q in _SOLENODIAL["q_grid"]),
        "probe:solenoidal-freedom:q_grid",
        None,
    ),
    BoundValue(
        "solenoidal.entropy_production",
        tuple(float(value) for value in _SOLENODIAL["entropy_production"]),
        "probe:solenoidal-freedom:entropy_production",
        None,
    ),
    BoundValue(
        "d_stability.real_stable",
        float(min(_DSTAB["eigenvalue_real_parts"])),
        "probe:d-stability:eigenvalue_real_parts:min",
        None,
    ),
    BoundValue(
        "d_stability.real_unstable",
        float(max(_DSTAB["eigenvalue_real_parts"])),
        "probe:d-stability:eigenvalue_real_parts:max",
        None,
    ),
    BoundValue(
        "d_stability.imag_unstable",
        float(abs(max(_DSTAB["eigenvalue_imag_parts"], key=abs))),
        "probe:d-stability:eigenvalue_imag_parts:absmax",
        None,
    ),
    BoundValue(
        "d_stability.lyapunov_residual",
        float(_DSTAB["lyapunov_residual"]),
        "probe:d-stability:lyapunov_residual",
        None,
    ),
    BoundValue(
        "d_stability.routh_k4_third",
        float(_DSTAB["routh_first_column_at_k4"][2]),
        "probe:d-stability:routh_first_column_at_k4:2",
        "-2",
    ),
    BoundValue(
        "d_stability.samples",
        tuple(float(value) for value in control_module.D_STABILITY_SAMPLES),
        "core/controls.py D_STABILITY_SAMPLES: the eigenvalue-locus sample grid",
        None,
    ),
    BoundValue(
        "mediator.epsilon_below",
        float(_MEDIATOR["below_threshold_epsilon"]),
        "probe:mediator-obstruction:below_threshold_epsilon",
        None,
    ),
    BoundValue(
        "mediator.margin_below",
        float(_MEDIATOR["margin_below_threshold"]),
        "probe:mediator-obstruction:margin_below_threshold",
        None,
    ),
    BoundValue(
        "mediator.margin_at",
        float(_MEDIATOR["margin_at_threshold"]),
        "probe:mediator-obstruction:margin_at_threshold",
        "0",
    ),
    BoundValue(
        "reduction.same_observation",
        float(_REDUCTION["same_observation"]),
        "probe:reduction-criterion:same_observation",
        "1",
    ),
    BoundValue(
        "reduction.first_after",
        float(_REDUCTION["first_after"]),
        "probe:reduction-criterion:first_after",
        None,
    ),
    BoundValue(
        "reduction.second_after",
        float(_REDUCTION["second_after"]),
        "probe:reduction-criterion:second_after",
        None,
    ),
    BoundValue(
        "reduction.response_distance",
        float(_REDUCTION["response_distance"]),
        "probe:reduction-criterion:response_distance",
        "2",
    ),
    BoundValue(
        "reduction.worst_case",
        float(_REDUCTION["worst_case_error_lower_bound"]),
        "probe:reduction-criterion:worst_case_error_lower_bound",
        "1",
    ),
    BoundValue(
        "four_cycle.violations",
        int(_FOURCYCLE["trap_violations"]),
        "probe:four-cycle-local-minimum:trap_violations",
        "2",
    ),
    BoundValue(
        "ep.closed_form",
        str(_EP["closed_form_contrast"]),
        "probe:equilibrium-propagation:closed_form_contrast",
        None,
        note="kept at full printed precision: the closed-form and relaxed contrasts differ in the eighth digit",
    ),
    BoundValue(
        "ep.relaxed",
        str(_EP["relaxed_contrast"]),
        "probe:equilibrium-propagation:relaxed_contrast",
        None,
        note="kept at full printed precision for the same reason",
    ),
    BoundValue(
        "ep.exact_gradient",
        float(_EP["exact_gradient"]),
        "probe:equilibrium-propagation:exact_gradient",
        "-1.25",
    ),
    BoundValue(
        "ep.closed_form_error",
        float(_EP["closed_form_error"]),
        "probe:equilibrium-propagation:closed_form_error",
        None,
    ),
    BoundValue(
        "ep.gradient_distance",
        float(_EP["distance_to_gradient"]),
        "probe:equilibrium-propagation:distance_to_gradient",
        None,
    ),
    BoundValue(
        "replicator.rate",
        float(_REPLICATOR["potential_rate_direct"]),
        "probe:replicator-potential:potential_rate_direct",
        None,
    ),
    BoundValue(
        "replicator.field_norm",
        float(_REPLICATOR["field_norm"]),
        "probe:replicator-potential:field_norm",
        None,
    ),
    BoundValue(
        "gmrf.ordered_pairs",
        int(_GMRF["pairwise_precision_entries"]),
        "probe:gmrf-independence:pairwise_precision_entries",
        None,
    ),
    BoundValue(
        "symmetry.eigenspace_dimensions",
        tuple(sorted(int(value) for value in _EIGEN["dimensions"])[1:]),
        "probe:eigenspace-irreducibility:dimensions:sortedtail",
        None,
        note="the non-trivial eigenspace dimensions beside the trivial one",
    ),
    BoundValue(
        "not_reproduced.port_zero_occupation",
        "35/264",
        "'Agreement and Surprise' (private draft): recorded not-reproduced; the record/selector protocol is not specified in its text",
        None,
        note="not-reproduced",
    ),
    BoundValue(
        "not_reproduced.move_probabilities",
        "13/15 and 1/60",
        "'Agreement and Surprise' (private draft): recorded not-reproduced; the probability space and the conditioning event are not determined by the text",
        None,
        note="not-reproduced",
    ),
    *_chart_values(),
)

MANUSCRIPT_TOKENS: dict[str, str] = as_tokens()
