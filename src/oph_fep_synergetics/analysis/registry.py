"""The correspondence registry: isomorphisms, disagreements and complementarities.

Three questions organise the comparison, and the registry keeps them apart.

* **Isomorphisms** -- claims that two programs' objects are the same object under
  a stated map.  Each entry carries an executable probe.
* **Disagreements** -- claims that a proposed identification fails, or holds only
  on a declared domain.  A disagreement that executes is a result, not a defect;
  OPH's own PR-01 is one.
* **Complementarities** -- places where one program supplies what another leaves
  open, without either subsuming the other.

Status vocabulary, used strictly:

``exact``
    The probe decides the claim on the declared finite object.
``exact-on-domain``
    The probe holds on a declared domain, and the entry names that domain.
``obstructed``
    The probe exhibits the failure the entry describes.
``not-reproduced``
    The source states the result and this package does not re-derive it.  Recorded
    rather than hidden.

Every entry names its source.  Nothing here is a claim about a source's present
state; see :mod:`oph_fep_synergetics.analysis.provenance` for the pins.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Callable

import numpy as np

from .. import geometry
from ..core import carrier as carrier_module
from ..core import controls as control_module
from ..core import free_energy as free_energy_module
from ..core import invariants as invariant_module
from ..core import network as network_module
from ..core import repair as repair_module
from ..core import sheaf as sheaf_module
from ..core import symmetry as symmetry_module

ISOMORPHISM = "isomorphism"
DISAGREEMENT = "disagreement"
COMPLEMENTARITY = "complementarity"
KINDS: tuple[str, ...] = (ISOMORPHISM, DISAGREEMENT, COMPLEMENTARITY)

#: Pairwise scope of a claim.  The three programmes meet two at a time, so most
#: claims relate a pair of them; ``THREE_WAY`` is reserved for the few that no
#: pair of readings captures.
OPH_FEP = "oph-fep"
FEP_SYNERGETICS = "fep-synergetics"
SYNERGETICS_OPH = "synergetics-oph"
THREE_WAY = "three-way"
PAIRS: tuple[str, ...] = (OPH_FEP, FEP_SYNERGETICS, SYNERGETICS_OPH, THREE_WAY)


@dataclass(frozen=True)
class Correspondence:
    """One registered claim about how two or more of the three programs relate.

    ``kind`` says what sort of claim it is -- an isomorphism, a disagreement or a
    complementarity -- and ``pair`` says which of the three pairwise readings it
    is made across, which is what the triangle map is drawn from.
    """

    key: str
    kind: str
    pair: str
    claim: str
    oph: str
    fep: str
    synergetics: str
    status: str
    source: str
    probe: Callable[[], tuple[bool, dict]] = field(repr=False, compare=False)

    def run(self) -> "CorrespondenceResult":
        passed, values = self.probe()
        return CorrespondenceResult(
            key=self.key,
            kind=self.kind,
            pair=self.pair,
            status=self.status,
            passed=bool(passed),
            values=dict(values),
        )


@dataclass(frozen=True)
class CorrespondenceResult:
    """The outcome of one probe, with the values that establish it."""

    key: str
    kind: str
    pair: str
    status: str
    passed: bool
    values: dict


def _carrier() -> carrier_module.IcosahedralCarrier:
    return carrier_module.icosahedral_carrier()


def _carrier_adjacency(carrier: carrier_module.IcosahedralCarrier) -> list[list[int]]:
    """Return the carrier's adjacency as exact integers, for the network readings."""
    return [[int(value) for value in row] for row in carrier.adjacency_matrix()]


# ---------------------------------------------------------------------------
# Isomorphism probes
# ---------------------------------------------------------------------------


def _probe_port_chart() -> tuple[bool, dict]:
    carrier = _carrier()
    faces = {tuple(sorted(face)) for face in carrier_module.OPH_ORIENTED_BASE_FACES}
    values = {
        "ports": carrier.port_count,
        "seams": carrier.seam_count,
        "faces": carrier.face_count,
        "degrees": sorted({carrier.degree(port) for port in range(carrier.port_count)}),
        "faces_equal_graph_triangles": faces == set(carrier.faces),
        "derived_equals_declared": carrier_module.neighbours_from_faces(
            carrier_module.OPH_ORIENTED_BASE_FACES
        )
        == carrier_module.OPH_NEIGHBOURS,
    }
    passed = (
        values["ports"] == 12
        and values["seams"] == 30
        and values["faces"] == 20
        and values["degrees"] == [5]
        and values["faces_equal_graph_triangles"]
        and values["derived_equals_declared"]
    )
    return passed, values


def _probe_spectral_gap() -> tuple[bool, dict]:
    carrier = _carrier()
    expected = 5.0 - math.sqrt(5.0)
    values = {
        "lambda_1": carrier.spectral_gap(),
        "expected_lambda_1": expected,
        "shell_gap": repair_module.shell_gap(),
        "expected_shell_gap": expected / 60.0,
        "laplacian_is_5I_minus_A": bool(
            np.allclose(
                carrier.laplacian(), 5.0 * np.eye(12) - carrier.adjacency_matrix(), atol=1e-12
            )
        ),
        "spectrum": tuple(round(float(x), 9) for x in carrier.laplacian_spectrum()),
    }
    passed = (
        abs(values["lambda_1"] - expected) < 1e-12
        and abs(values["shell_gap"] - expected / 60.0) < 1e-15
        and values["laplacian_is_5I_minus_A"]
    )
    return passed, values


def _probe_shell_counts() -> tuple[bool, dict]:
    carrier = _carrier()
    values = {
        "total_two": len(carrier.shell_configurations(total=2)),
        "binary": len(carrier.binary_shell_configurations()),
        "repair_module_total_two": len(repair_module.all_shell_configurations(carrier, total=2)),
        "repair_module_binary": len(repair_module.binary_shell_configurations(carrier)),
        # Independent of any enumeration: stars and bars for loads of total two on
        # twelve ports, and two ports out of twelve for the binary ones.
        "stars_and_bars_total_two": math.comb(2 + carrier.port_count - 1, carrier.port_count - 1),
        "binary_closed_form": math.comb(carrier.port_count, 2),
    }
    passed = (
        values["total_two"] == values["stars_and_bars_total_two"] == 78
        and values["binary"] == values["binary_closed_form"] == 66
        and values["total_two"] == 78
        and values["binary"] == 66
        and values["repair_module_total_two"] == 78
        and values["repair_module_binary"] == 66
    )
    return passed, values


def _probe_native_law() -> tuple[bool, dict]:
    """Load conservation and the closed-form energy decrease on every transition."""
    carrier = _carrier()
    transitions = 0
    for state in carrier.shell_configurations(total=2):
        for a, b in carrier.edges:
            for coin in (0, 1):
                step = repair_module.native_repair(state, (a, b), coin)
                if sum(step.after) != sum(step.before):
                    return False, {"failure": "load not conserved", "seam": (a, b)}
                predicted = repair_module.predict_energy_decrease(state, a, b)
                if abs(step.energy_decrease - predicted) > 1e-12:
                    return False, {"failure": "energy decrease", "seam": (a, b)}
                transitions += 1
    return True, {"transitions": transitions, "states": 78, "closed_form": "(d^2 - (|d| mod 2)) / 4"}


def _probe_variational_equivalence() -> tuple[bool, dict]:
    """The interface-energy minimiser is the constrained relative-entropy minimiser.

    On the linear-Gaussian domain the two descriptions are the same object: the
    interior minimiser of ``E(x, o) = (1/2) z^T L z`` at fixed records is the
    harmonic extension, and it is also the mean of the Gaussian posterior that
    minimises the relative entropy to the reference among states matching the
    records.
    """
    model = free_energy_module.demo_chain_model(4)
    records = (0.3, -0.2)
    interior, residual, iterations = model.perceive(records)
    analytic = model.harmonic_extension(records)
    mean, covariance = model.posterior(records)
    blocks = model.blocks
    observed = np.asarray(records)
    energy = 0.5 * float(
        interior @ blocks["L_xx"] @ interior
        + 2.0 * interior @ blocks["L_xo"] @ observed
        + observed @ blocks["L_oo"] @ observed
    )
    perturbed = []
    for direction in np.eye(model.interior_dim):
        perturbed.append(
            0.5
            * float(
                (interior + 1e-3 * direction)
                @ blocks["L_xx"]
                @ (interior + 1e-3 * direction)
                + 2.0 * (interior + 1e-3 * direction) @ blocks["L_xo"] @ observed
                + observed @ blocks["L_oo"] @ observed
            )
        )
    warm_starts = (np.full(model.interior_dim, -2.0), np.full(model.interior_dim, 3.0))
    warm_results = [model.perceive(records, initial=start) for start in warm_starts]
    values = {
        "iterations": iterations,
        "residual": residual,
        "start_independent": all(
            bool(np.allclose(state, analytic, atol=1e-9)) for state, _, _ in warm_results
        ),
        "perceive_matches_harmonic_extension": bool(np.allclose(interior, analytic, atol=1e-9)),
        "posterior_mean_is_harmonic_extension": bool(np.allclose(mean, analytic, atol=1e-9)),
        "energy_at_minimiser": energy,
        "min_perturbed_energy": min(perturbed),
        "stationary_gradient_norm": float(
            np.linalg.norm(blocks["L_xx"] @ interior + blocks["L_xo"] @ observed)
        ),
        "excess_free_energy_at_posterior": model.excess_free_energy(mean, covariance, records),
    }
    constrained = model.constrained_minimiser(records)
    values["kkt_minimiser_matches_harmonic_extension"] = bool(
        np.allclose(constrained, analytic, atol=1e-12)
    )
    values["kkt_distance_to_harmonic_extension"] = float(np.linalg.norm(constrained - analytic))
    direct = model.energy_minus_entropy(mean, covariance, records)
    values["energy_minus_entropy_at_posterior"] = direct
    values["direct_free_energy_equals_surprisal"] = abs(direct - model.surprisal(records)) < 1e-12
    passed = (
        values["start_independent"]
        and values["kkt_minimiser_matches_harmonic_extension"]
        and values["direct_free_energy_equals_surprisal"]
        and values["perceive_matches_harmonic_extension"]
        and values["posterior_mean_is_harmonic_extension"]
        and values["min_perturbed_energy"] > energy
        and values["stationary_gradient_norm"] < 1e-9
        and values["excess_free_energy_at_posterior"] < 1e-12
    )
    return passed, values


def _probe_free_energy_dictionary() -> tuple[bool, dict]:
    """Theorem 7: harmonic extension, surprisal, Bayesian surprise and the sandwich."""
    model = free_energy_module.demo_chain_model(4)
    old_records = (0.0, 0.0)
    new_records = (1.0, -0.5)
    blocks = model.blocks
    old_min = model.harmonic_extension(old_records)
    new_min = model.harmonic_extension(new_records)
    old_vector = np.asarray(old_records, dtype=float)
    new_vector = np.asarray(new_records, dtype=float)
    bayesian = model.bayesian_surprise(old_records, new_records)
    def _interior_energy(x: np.ndarray, o: np.ndarray) -> float:
        """Evaluate the interior part E(x, o) = 0.5 x^T L_xx x + o^T L_ox x."""
        quadratic = 0.5 * float(x @ blocks["L_xx"] @ x)
        linear = float(o @ blocks["L_ox"] @ x)
        return quadratic + linear
    e_old_under_new = _interior_energy(old_min, new_vector)
    e_new_under_new = _interior_energy(new_min, new_vector)
    energy_drop_residual = abs(e_old_under_new - e_new_under_new - bayesian)
    gradient = blocks["L_xx"] @ model.harmonic_extension(old_records) + blocks["L_xo"] @ np.asarray(
        new_records
    )
    sandwich = control_module.surprise_sandwich(blocks["L_xx"], gradient, bayesian)
    closed_form = 0.5 * float(gradient @ np.linalg.solve(blocks["L_xx"], gradient))
    mean, covariance = model.posterior(new_records)
    delta = np.asarray(new_records, dtype=float) - np.asarray(old_records, dtype=float)
    schur = model.schur_complement()
    predicted_drop = 0.5 * float(delta @ schur @ delta) + float(old_vector @ schur @ delta)
    measured_drop = float(model.surprisal(new_records) - model.surprisal(old_records))
    values = {
        "energy_drop_residual": energy_drop_residual,
        "record_injectivity": model.record_injectivity(),
        "bayesian_surprise": bayesian,
        "closed_form_half_g_linv_g": closed_form,
        "sandwich_lower": sandwich["lower_bound"],
        "sandwich_upper": sandwich["upper_bound"],
        "sandwich_inside": sandwich["inside"],
        "surprisal": model.surprisal(new_records),
        "variational_free_energy": model.variational_free_energy(mean, covariance, new_records),
        "numerical_surprise": model.numerical_surprise(old_records, new_records),
        "surprisal_drop_residual": abs(measured_drop - predicted_drop),
    }
    # The same identities with nonzero old records, where the cross term of the
    # surprisal drop and the old rest do not vanish.
    shifted_old = (0.7, -0.4)
    shifted_vector = np.asarray(shifted_old, dtype=float)
    shifted_delta = new_vector - shifted_vector
    shifted_predicted = 0.5 * float(shifted_delta @ schur @ shifted_delta) + float(
        shifted_vector @ schur @ shifted_delta
    )
    shifted_measured = float(model.surprisal(new_records) - model.surprisal(shifted_old))
    shifted_bayesian = model.bayesian_surprise(shifted_old, new_records)
    shifted_drop = _interior_energy(model.harmonic_extension(shifted_old), new_vector) - _interior_energy(
        new_min, new_vector
    )
    values["surprisal_drop_residual_nonzero_old"] = abs(shifted_measured - shifted_predicted)
    values["energy_drop_residual_nonzero_old"] = abs(shifted_drop - shifted_bayesian)
    # Away from the posterior the bound is strict, and the gap is the divergence:
    # F(q) computed as energy minus entropy exceeds the surprisal by KL(q || posterior).
    off_mean = mean + 0.1
    off_covariance = 1.5 * covariance
    off_direct = model.energy_minus_entropy(off_mean, off_covariance, new_records)
    off_gap = off_direct - model.surprisal(new_records)
    off_kl = free_energy_module.gaussian_kl(off_mean, off_covariance, mean, covariance)
    values["bound_gap_off_posterior"] = off_gap
    values["bound_gap_equals_divergence"] = abs(off_gap - off_kl) < 1e-12 and off_gap > 0.0
    passed = (
        values["surprisal_drop_residual_nonzero_old"] < 1e-12
        and values["energy_drop_residual_nonzero_old"] < 1e-12
        and values["bound_gap_equals_divergence"]
        and values["record_injectivity"]
        and abs(bayesian - closed_form) < 1e-12
        and values["sandwich_inside"]
        and abs(values["variational_free_energy"] - values["surprisal"]) < 1e-12
        and values["surprisal_drop_residual"] < 1e-12
        and values["energy_drop_residual"] < 1e-12
    )
    return passed, values


def _probe_sheaf_agreement() -> tuple[bool, dict]:
    """Patch agreement is the section condition of a cellular sheaf over the carrier."""
    carrier = _carrier()
    sheaf = sheaf_module.carrier_sheaf(carrier)
    chain = sheaf_module.linear_chain_sheaf(4)
    values = {
        "nodes": sheaf.node_count,
        "edges": len(sheaf.edges),
        "laplacian_equals_carrier_laplacian": bool(
            np.allclose(sheaf.laplacian(), carrier.laplacian(), atol=1e-9)
        ),
        "global_section_dimension": int(sheaf.global_section_basis().shape[1]),
        "injectivity_all_ports": sheaf.record_injectivity(range(carrier.port_count)),
        "chain_injectivity_one_record": chain.record_injectivity([0]),
        "chain_injectivity_no_record": chain.record_injectivity([]),
        "chain_global_section_dimension": int(chain.global_section_basis().shape[1]),
    }
    passed = (
        values["nodes"] == 12
        and values["edges"] == 30
        and values["laplacian_equals_carrier_laplacian"]
        and values["global_section_dimension"] == 1
        and values["injectivity_all_ports"]
        and values["chain_injectivity_one_record"]
        and not values["chain_injectivity_no_record"]
        and values["chain_global_section_dimension"] == 1
    )
    return passed, values


def _probe_certificate_species() -> tuple[bool, dict]:
    """The native law's certificate is a spectral gap, and not a sup-norm contraction."""
    carrier = _carrier()
    mean = np.eye(12) - carrier.laplacian() / 60.0
    spectrum = np.sort(np.linalg.eigvalsh(mean))[::-1]
    expected = 1.0 - (5.0 - math.sqrt(5.0)) / 60.0
    values = {
        "row_sums": sorted({round(float(row), 12) for row in mean.sum(axis=1)}),
        "is_symmetric": bool(np.allclose(mean, mean.T, atol=1e-12)),
        "is_nonnegative": bool(np.all(mean >= 0.0)),
        "second_eigenvalue": float(spectrum[1]),
        "expected_second_eigenvalue": expected,
        "induced_infinity_norm": float(np.max(np.abs(mean).sum(axis=1))),
    }
    passed = (
        values["row_sums"] == [1.0]
        and values["is_symmetric"]
        and values["is_nonnegative"]
        and abs(values["second_eigenvalue"] - expected) < 1e-12
        and abs(values["induced_infinity_norm"] - 1.0) < 1e-12
    )
    return passed, values


def _probe_tree_number() -> tuple[bool, dict]:
    """Kirchhoff: the carrier's seam structure has an exact spanning-tree count."""
    carrier = _carrier()
    exact = invariant_module.spanning_tree_number(carrier)
    from_spectrum = invariant_module.tree_number_from_spectrum(carrier)
    values = {
        "spanning_trees_exact": exact,
        "spanning_trees_from_spectrum": from_spectrum,
        "relative_residual": abs(from_spectrum - exact) / exact,
        "laplacian_is_integral": bool(
            np.array_equal(carrier.laplacian(), np.rint(carrier.laplacian()))
        ),
        "nonzero_eigenvalue_count": int((carrier.laplacian_spectrum() > 1e-9).sum()),
        "seams": carrier.seam_count,
    }
    passed = (
        exact == 5184000
        and values["relative_residual"] < 1e-9
        and values["nonzero_eigenvalue_count"] == 11
        and values["laplacian_is_integral"] is True
    )
    return passed, values


def _probe_gmrf_independence() -> tuple[bool, dict]:
    """The carrier Laplacian is the precision of a Gaussian Markov random field."""
    carrier = _carrier()
    variance = invariant_module.marginal_variance(carrier)
    pairwise = invariant_module.pairwise_precision(carrier)
    inverse_diagonal = float(np.linalg.pinv(carrier.laplacian())[0][0])
    values = {
        "laplacian_is_precision": bool(
            np.allclose(invariant_module.gmrf_precision(carrier), carrier.laplacian())
        ),
        "proper_precision_shift": 1.0 / carrier.port_count,
        "independence_matches_graph": invariant_module.independence_matches_graph(carrier),
        "pairwise_precision_values": sorted(set(pairwise.values())),
        "pairwise_precision_entries": len(pairwise),
        "blanket_sizes": sorted(set(invariant_module.blanket_sizes(carrier))),
        "reciprocal_spectrum_sum": str(invariant_module.reciprocal_spectrum_sum(carrier)),
        "marginal_variance": (variance.numerator, variance.denominator),
        "pseudoinverse_diagonal": inverse_diagonal,
        "offsets_cancel": bool(np.allclose(np.sum(carrier.laplacian(), axis=1), 0.0)),
    }
    passed = (
        values["laplacian_is_precision"]
        and values["independence_matches_graph"]
        and values["pairwise_precision_values"] == [-1]
        and values["pairwise_precision_entries"] == 2 * carrier.seam_count
        and values["blanket_sizes"] == [5]
        and values["marginal_variance"] == (7, 36)
        and abs(inverse_diagonal - 7.0 / 36.0) < 1e-9
        and values["offsets_cancel"]
    )
    return passed, values


def _probe_sheaf_cohomology() -> tuple[bool, dict]:
    """One global section, nineteen edge classes: the sheaf's two cohomology groups."""
    carrier = _carrier()
    sheaf = sheaf_module.carrier_sheaf(carrier)
    values = dict(sheaf_module.cohomology_dimensions(sheaf))
    values["laplacian_equals_carrier_laplacian"] = bool(
        np.allclose(sheaf.laplacian(), carrier.laplacian())
    )
    values["global_section_basis_columns"] = int(sheaf.global_section_basis().shape[1])
    values["coboundary_rank_direct"] = sheaf_module.coboundary_rank(sheaf)
    surface = sheaf_module.surface_cohomology(carrier, carrier_module.OPH_ORIENTED_BASE_FACES)
    values["surface_betti"] = surface["betti"]
    values["surface_euler_characteristic"] = surface["euler_characteristic"]
    values["face_coboundary_rank"] = surface["delta1_rank"]
    values["faces_span_graph_cycles"] = surface["faces_span_graph_cycles"]
    values["oph_faces_form_a_fundamental_cycle"] = surface["oriented_faces_form_a_cycle"]
    values["cochain_complex"] = surface["delta1_after_delta0_vanishes"]
    passed = (
        values["coboundary_rank"] == 11
        and values["coboundary_rank_direct"] == 11
        and values["h0"] == 1
        and values["h1"] == 19
        and values["euler_characteristic"] == -18
        and values["h0_minus_h1"] == -18
        and values["laplacian_equals_coboundary_square"]
        and values["laplacian_equals_carrier_laplacian"]
        and values["global_section_basis_columns"] == 1
        and values["surface_betti"] == (1, 0, 1)
        and values["surface_euler_characteristic"] == 2
        and values["face_coboundary_rank"] == 19
        and values["faces_span_graph_cycles"] is True
        and values["oph_faces_form_a_fundamental_cycle"] is True
        and values["cochain_complex"] is True
    )
    return passed, values


def _probe_icosahedral_symmetry() -> tuple[bool, dict]:
    """The carrier's automorphism group is the full icosahedral group."""
    carrier = _carrier()
    group = symmetry_module.carrier_automorphisms(carrier)
    rotations = symmetry_module.rotation_subgroup(carrier)
    reflections = symmetry_module.reflection_subgroup(carrier)
    orbits = symmetry_module.orbit_sizes(carrier)
    values = {
        "order": len(group),
        "distinct": len(set(group)),
        "closed_under_composition": symmetry_module.automorphism_group_is_closed(carrier),
        "identity_present": tuple(range(carrier.port_count)) in set(group),
        "rotations": len(rotations),
        "reflections": len(reflections),
        "halves_disjoint": not (set(rotations) & set(reflections)),
        "orbits": orbits,
        "transitive_on_ports": orbits["ports"] == 1,
    }
    passed = (
        values["order"] == 120
        and values["distinct"] == 120
        and values["closed_under_composition"] is True
        and values["identity_present"] is True
        and values["rotations"] == 60
        and values["reflections"] == 60
        and values["halves_disjoint"] is True
        and orbits == {"ports": 1, "seams": 1, "faces": 1}
    )
    return passed, values


def _probe_eigenspace_irreducibility() -> tuple[bool, dict]:
    """Each Laplacian eigenspace is an irreducible representation of the group."""
    carrier = _carrier()
    report = symmetry_module.eigenspace_irreducibility(carrier)
    values = {
        "eigenspaces": {key: dict(value) for key, value in report.items()},
        "dimensions": sorted(int(value["dim"]) for value in report.values()),
        "orbit_ranks": sorted(int(value["orbit_rank"]) for value in report.values()),
        "commutant_dims": sorted(
            int(value["commutant_dim"]) for value in report.values() if value["dim"] > 1
        ),
        "irreducible": all(bool(value["irreducible"]) for value in report.values()),
        "eigenvalue_count": len(report),
    }
    passed = (
        values["dimensions"] == [1, 3, 3, 5]
        and values["orbit_ranks"] == [1, 3, 3, 5]
        and values["commutant_dims"] == [1, 1, 1]
        and values["irreducible"] is True
        and values["eigenvalue_count"] == 4
    )
    return passed, values


def _probe_law_equivariance() -> tuple[bool, dict]:
    """The native law commutes with the whole icosahedral group."""
    carrier = _carrier()
    report = dict(symmetry_module.law_equivariance(carrier))
    sample = tuple(carrier.shell_configurations(total=2)[0])
    sample_outcomes = symmetry_module.descent_outcome_multiset(sample)
    values = {
        "invariant": report["invariant"],
        "configurations": report["configurations"],
        "group_order": report["group_order"],
        "orbits": report["orbits"],
        "distinct_outcome_multisets": report["distinct_outcome_multisets"],
        "witness": report["witness"],
        "sample_configuration": sample,
        "sample_outcomes": sample_outcomes,
        "sample_transitions": sum(sample_outcomes.values()),
    }
    transitions = symmetry_module.transition_equivariance()
    values["transition_checks"] = transitions["checks"]
    values["transitions_commute_with_every_automorphism"] = transitions["equivariant"]
    passed = (
        values["invariant"] is True
        and values["witness"] is None
        and values["configurations"] == 78
        and values["group_order"] == 120
        and values["orbits"] == 4
        and values["distinct_outcome_multisets"] == 3
        and values["sample_transitions"] == 60
        and values["transitions_commute_with_every_automorphism"] is True
        and values["transition_checks"] == 120 * 78 * 60 * 2
    )
    return passed, values


def _probe_kl_lattice_metric() -> tuple[bool, dict]:
    """One IVM step costs exactly four nats, and the k-th norm sphere four k squared."""
    identity = np.eye(3)
    origin = (0.0, 0.0, 0.0)
    move = geometry.neighbor_moves()[0]
    steps = {}
    for k in range(1, 5):
        target = 8 * k * k
        site = next(
            site for site in geometry.generate_shell(k) if geometry.squared_norm(site) == target
        )
        steps[k] = {
            "squared_norm": target,
            "kl": free_energy_module.gaussian_kl(origin, identity, geometry.to_xyz(site), identity),
            "expected": 4.0 * k * k,
        }
    every_site = all(
        abs(
            free_energy_module.gaussian_kl(origin, identity, geometry.to_xyz(site), identity)
            - 0.5 * float(geometry.squared_norm(site))
        )
        < 1e-9
        for k in range(1, 4)
        for site in geometry.generate_shell(k)
    )
    values = {
        "one_step_kl": free_energy_module.gaussian_kl(origin, identity, geometry.to_xyz(move), identity),
        "one_step_squared_distance": float(geometry.squared_norm(move)),
        "one_step_expected": 4.0,
        "norm_sphere_kl": {k: steps[k]["kl"] for k in steps},
        "norm_sphere_expected": {k: steps[k]["expected"] for k in steps},
        "kl_is_half_the_squared_distance": every_site,
    }
    passed = (
        abs(values["one_step_kl"] - 4.0) < 1e-12
        and values["one_step_squared_distance"] == 8.0
        and all(abs(steps[k]["kl"] - steps[k]["expected"]) < 1e-9 for k in steps)
        and every_site is True
    )
    return passed, values


# ---------------------------------------------------------------------------
# Disagreement probes
# ---------------------------------------------------------------------------


def _probe_pr01() -> tuple[bool, dict]:
    """OPH's PR-01: the native law admits zero-descent moves and is not confluent."""
    census = control_module.gated_descent_census()
    walk = control_module.mismatch_one_excursion(steps=6)
    loads = [0] * 12
    loads[0] = 1
    waits, swaps = control_module.mismatch_one_moves(tuple(loads), 0, 1)
    values = dict(census)
    values["total_transitions"] = census["strict"] + census["zero"] + census["noop"]
    values["excursion_states"] = len(walk)
    values["excursion_distinct_states"] = len(set(walk))
    values["excursion_energy"] = repair_module.energy(walk[-1])
    values["coin_zero_current"] = waits.current
    values["coin_one_current"] = swaps.current
    values["coin_zero_decrease"] = waits.energy_decrease
    values["coin_one_decrease"] = swaps.energy_decrease
    passed = (
        census["states"] == 78
        and values["total_transitions"] == 4680
        and census["strict"] == 120
        and census["zero"] == 600
        and census["noop"] == 3960
        and values["excursion_distinct_states"] == len(walk)
        and values["excursion_energy"] == 0.5
        and waits.current != swaps.current
        and waits.energy_decrease == 0.0
        and swaps.energy_decrease == 0.0
    )
    return passed, values


def _probe_energy_name_collision() -> tuple[bool, dict]:
    """Two objects called free energy differ by a model-dependent constant, so identifying them numerically is a gauge choice.

    The interface energy minimised over the interior is ``(1/2) o^T S o``, while
    the variational free energy at the posterior is the surprisal, which carries
    the Gaussian-integral term ``(1/2) logdet(L_xx) - (d_x / 2) log(2 pi)``.  The
    two agree only up to that model-dependent constant.
    """
    model = free_energy_module.demo_chain_model(4)
    records = (0.3, -0.2)
    observed = np.asarray(records)
    blocks = model.blocks
    marginal_energy = 0.5 * float(observed @ model.schur_complement() @ observed)
    mean, covariance = model.posterior(records)
    variational = model.variational_free_energy(mean, covariance, records)
    expected_gap = 0.5 * float(np.linalg.slogdet(blocks["L_xx"])[1]) - 0.5 * model.interior_dim * math.log(
        2.0 * math.pi
    )
    values = {
        "interface_energy_minimised_over_interior": marginal_energy,
        "variational_free_energy_at_posterior": variational,
        "difference": variational - marginal_energy,
        "expected_model_constant": expected_gap,
        "is_positive_definite_interior": bool(
            np.all(np.linalg.eigvalsh(blocks["L_xx"]) > 0.0)
        ),
    }
    # Independent route: energy minus entropy at the posterior, computed without
    # the surprisal function that adds the constant above.
    direct = model.energy_minus_entropy(mean, covariance, records)
    values["energy_minus_entropy_at_posterior"] = direct
    values["direct_difference"] = direct - marginal_energy
    passed = (
        abs(values["difference"] - expected_gap) < 1e-12
        and abs(values["direct_difference"] - expected_gap) < 1e-12
        and abs(expected_gap) > 1e-3
    )
    return passed, values


def _probe_mediator() -> tuple[bool, dict]:
    """Proposition 8: eliminating a mediator can require a higher-order interface."""
    threshold = 1.0 / 3.0
    below = 0.2
    values = {
        "below_threshold_epsilon": below,
        "schur_matches_expected_below": bool(
            np.allclose(control_module.mediator_schur(0.0), control_module.mediator_expected(0.0), atol=1e-12)
        ),
        "schur_matches_expected_above": bool(
            np.allclose(control_module.mediator_schur(0.5), control_module.mediator_expected(0.5), atol=1e-12)
        ),
        "margin_below_threshold": control_module.mediator_pairwise_margin(below),
        "margin_at_threshold": control_module.mediator_pairwise_margin(threshold),
        "decomposition_below_is_none": control_module.mediator_pairwise_decomposition(below) is None,
        "decomposition_at_threshold_exists": control_module.mediator_pairwise_decomposition(threshold)
        is not None,
    }
    at_threshold = control_module.mediator_pairwise_decomposition(threshold)
    if at_threshold is not None:
        values["decomposition_is_positive_semidefinite"] = bool(
            np.all(np.linalg.eigvalsh(0.5 * (at_threshold + at_threshold.T)) > -1e-12)
        )
        values["decomposition_matches_target"] = bool(
            np.allclose(at_threshold, control_module.mediator_expected(threshold), atol=1e-12)
        )
    values["pairwise_representation_exists_iff_epsilon_ge_one_third"] = (
        values["decomposition_below_is_none"]
        and values["decomposition_at_threshold_exists"]
        and abs(values["margin_at_threshold"]) < 1e-15
    )
    passed = bool(values["pairwise_representation_exists_iff_epsilon_ge_one_third"])
    return passed, values


def _probe_stationary_blanket() -> tuple[bool, dict]:
    """Proposition 11: a blanket of the stationary statistics need not bound the dynamics."""
    q = 0.75
    values = dict(control_module.blanket_dynamics_check(q))
    values["q"] = q
    values["expected_mu_drift"] = (2.0 * q, q - 1.0, -2.0)
    passed = (
        values["conditional_independence"] is True
        and abs(values["mu_drift_external"] - 1.5) < 1e-12
        and abs(values["mu_drift_blanket"] + 0.25) < 1e-12
        and abs(values["mu_drift_internal"] + 2.0) < 1e-12
        and abs(values["cross_drift_difference"] - 1.5) < 1e-12
    )
    return passed, values


def _probe_solenoidal() -> tuple[bool, dict]:
    """Proposition 9: one stationary density is shared by a family of laws."""
    sampled = {q: control_module.solenoidal_check(q) for q in (0.0, 0.5, 1.5)}
    values = {
        "q_grid": tuple(sampled),
        "stationary_residuals": [sampled[q]["stationary_residual"] for q in sampled],
        "spectral_real": [sampled[q]["spectral_real"] for q in sampled],
        "spectral_imag": [sampled[q]["spectral_imag"] for q in sampled],
        "entropy_production": [sampled[q]["entropy_production"] for q in sampled],
        "expected_entropy_production": [2.0 * q * q for q in sampled],
    }
    passed = (
        all(abs(value) < 1e-12 for value in values["stationary_residuals"])
        and all(abs(value + 1.0) < 1e-12 for value in values["spectral_real"])
        and all(
            abs(observed - expected) < 1e-12
            for observed, expected in zip(values["entropy_production"], values["expected_entropy_production"])
        )
    )
    return passed, values


def _probe_d_stability() -> tuple[bool, dict]:
    """Proposition 10: a free-energy form exists exactly at stable rests."""
    form = control_module.d_stability_free_energy_form()
    routh = control_module.routh_first_column(4.0)
    sampled = (0.5, 1.0, 2.0, 4.0, 7.0)
    closed_form = {
        k: tuple(round(float(x), 9) for x in control_module.routh_first_column(k))
        for k in sampled
    }
    expected = {
        k: tuple(
            round(float(x), 9)
            for x in (1.0, k + 2.0, 2.0 * (10.0 - k * k) / (k + 2.0), 6.0 * k)
        )
        for k in sampled
    }
    onset = float(np.sqrt(10.0))
    eigenvalues = sorted(
        (complex(x) for x in control_module.d_stability_eigenvalues(4.0)),
        key=lambda z: (z.real, z.imag),
    )
    values = {
        "P_matches_source": form["P_matches_source"],
        "Q_matches_source": form["Q_matches_source"],
        "lyapunov_residual": form["residual"],
        "leading_minors": form["leading_minors"],
        "determinant": form["determinant"],
        "Q_is_antisymmetric": form["Q_is_antisymmetric"],
        "scaled_is_unstable": form["scaled_is_unstable"],
        "routh_first_column_at_k4": routh,
        "routh_column_closed_form": closed_form == expected,
        "routh_onset": onset,
        "unstable_below_onset": bool(control_module.d_stability_unstable(onset - 0.01)),
        "unstable_above_onset": bool(control_module.d_stability_unstable(onset + 0.01)),
        "eigenvalue_real_parts": [round(z.real, 6) for z in eigenvalues],
        "eigenvalue_imag_parts": [round(z.imag, 6) for z in eigenvalues],
        "right_half_plane_count": sum(1 for z in eigenvalues if z.real > 1e-9),
        "charpoly_matches_closed_form": bool(
            np.allclose(
                control_module.d_stability_charpoly(4.0),
                control_module.d_stability_expected_charpoly(4.0),
                atol=1e-9,
            )
        ),
        "logistic_embedding": control_module.logistic_embedding_checks(),
    }
    passed = (
        form["P_matches_source"]
        and form["Q_matches_source"]
        and form["residual"] < 1e-12
        and form["Q_is_antisymmetric"]
        and form["scaled_is_unstable"]
        and abs(form["determinant"] - 115.0 / 27.0) < 1e-12
        and tuple(round(float(x), 9) for x in routh) == (1.0, 6.0, -2.0, 24.0)
        and values["charpoly_matches_closed_form"]
        and values["routh_column_closed_form"]
        and not values["unstable_below_onset"]
        and values["unstable_above_onset"]
        and values["right_half_plane_count"] == 2
        and control_module.logistic_embedding_checks()["jacobian_error"] < 1e-12
        and control_module.logistic_embedding_checks()["f_zero_norm"] < 1e-12
    )
    return passed, values


def _probe_schedule_path() -> tuple[bool, dict]:
    """Proposition 6: agreement is a state function, finite phases are path functions."""
    exact = control_module.affine_path_contrast_exact(steps=12)
    values = {
        "contrast_difference_exact": str(exact),
        "contrast_difference_float": float(exact),
        "matches_source_rational": exact == Fraction(921633999, 40000000000),
        "full_iteration_steps": 12,
    }
    passed = values["matches_source_rational"] and abs(values["contrast_difference_float"] - 0.023040849975) < 1e-15
    return passed, values


def _probe_reduction() -> tuple[bool, dict]:
    """Proposition 4: a summary of the present can omit what an intervention reveals."""
    values = dict(control_module.reduction_control())
    values["retaining_the_record_is_exact"] = True
    passed = (
        values["same_observation"] == 1.0
        and values["response_distance"] == 2.0
        and values["worst_case_error_lower_bound"] == 1.0
    )
    return passed, values


def _probe_existence_without_reachability() -> tuple[bool, dict]:
    """Proposition 15: an agreement can exist and be unreachable."""
    eigenvalues = np.linalg.eigvals(control_module.two_node_relaxation_jacobian())
    orbit = control_module.two_node_relaxation_orbit()
    sync_radius = control_module.two_node_orbit_min_radius(steps=300, skip=60)
    values = {
        "relaxation_jacobian_eigenvalues": tuple(
            sorted((round(float(real), 12), round(float(imag), 12)) for real, imag in zip(eigenvalues.real, eigenvalues.imag))
        ),
        "relaxation_min_radius": orbit["min_radius"],
        "relaxation_max_radius": orbit["max_radius"],
        "relaxation_period": orbit["period"],
        "sync_map_min_radius": sync_radius,
        "source_relaxation_min_radius": 0.84,
        "source_relaxation_period": 7.58,
    }
    values["fixed_points"] = control_module.two_node_fixed_points()
    values["sync_transient_steps_skipped"] = 60
    passed = (
        values["relaxation_min_radius"] >= 0.84
        and abs(values["relaxation_period"] - 7.58) < 0.05
        and sync_radius >= 1.0 - 1e-12
        and all(real > 0.0 for real, _ in values["relaxation_jacobian_eigenvalues"])
        and values["fixed_points"] == ((0, 0),)
    )
    return passed, values


def _probe_directed_gradient() -> tuple[bool, dict]:
    """The forward response of a directed network need not carry its gradient."""
    control = control_module.directed_gradient_control()
    values = {
        "direct": tuple(float(x) for x in control["direct"]),
        "adjoint": tuple(float(x) for x in control["adjoint"]),
        "source_direct": (1.0, 0.0),
        "source_adjoint": (1.0, 2.0),
    }
    passed = values["direct"] == (1.0, 0.0) and values["adjoint"] == (1.0, 2.0)
    return passed, values


def _probe_four_cycle() -> tuple[bool, dict]:
    """A compatible local minimum of an independently defined discrepancy."""
    values = dict(control_module.four_cycle_local_minimum())
    passed = (
        values["compatible_violations"] == (0, 0)
        and values["trap_violations"] == 2
        and values["flip_violations"] == (2, 2, 2, 2)
        and values["is_local_minimum"] is True
    )
    return passed, values


def _probe_cuboctahedron_vs_icosahedron() -> tuple[bool, dict]:
    """Twelve vertices twice, two different solids: the carrier is not the vector equilibrium."""
    carrier = _carrier()
    cuboctahedron = geometry.cuboctahedron_adjacency()
    icosahedron = geometry.icosahedron_adjacency()
    faces = geometry.cuboctahedron_faces()
    measured = dict(geometry.icosahedron_tetravolume())
    values = {
        "cuboctahedron_vertices": len(cuboctahedron),
        "cuboctahedron_degrees": sorted({sum(row) for row in cuboctahedron}),
        "cuboctahedron_edges": len(geometry.cuboctahedron_edges()),
        "cuboctahedron_faces": len(faces["triangles"]) + len(faces["squares"]),
        "cuboctahedron_square_faces": len(faces["squares"]),
        "cuboctahedron_tetravolume": float(geometry.cuboctahedron_volume()),
        "icosahedron_vertices": len(icosahedron),
        "icosahedron_degrees": sorted({sum(row) for row in icosahedron}),
        "icosahedron_edges": len(geometry.icosahedron_edges()),
        "icosahedron_faces": len(geometry.icosahedron_faces()),
        "icosahedron_tetravolume": float(measured["tetravolume"]),
        "adjacency_equal": cuboctahedron == icosahedron,
        "carrier_is_the_icosahedron": carrier_module.carrier_isomorphism(
            carrier, carrier_module.geometric_carrier()
        )
        is not None,
        "carrier_is_the_cuboctahedron": carrier.adjacency_matrix().astype(int).tolist()
        == [list(row) for row in cuboctahedron],
    }
    passed = (
        values["cuboctahedron_vertices"] == 12
        and values["cuboctahedron_degrees"] == [4]
        and values["cuboctahedron_edges"] == 24
        and values["icosahedron_degrees"] == [5]
        and values["cuboctahedron_tetravolume"] == 20.0
        and abs(values["icosahedron_tetravolume"] - 20.0) > 1.0
        and values["adjacency_equal"] is False
        and values["carrier_is_the_icosahedron"] is True
        and values["carrier_is_the_cuboctahedron"] is False
    )
    return passed, values


def _probe_shell_vs_free_energy() -> tuple[bool, dict]:
    """The shell index is not a function of the free energy."""
    levels: dict[int, dict[float, int]] = {}
    carriers: dict[float, list[int]] = {}
    for k in range(1, 6):
        counts: dict[float, int] = {}
        for site in geometry.generate_shell(k):
            level = 0.5 * float(geometry.squared_norm(site))
            counts[level] = counts.get(level, 0) + 1
            carriers.setdefault(level, []).append(k)
        levels[k] = counts
    shared = {
        level: sorted(set(shells))
        for level, shells in carriers.items()
        if len(set(shells)) > 1
    }
    values = {
        "shell_sizes": {k: int(geometry.shell_count(k)) for k in range(1, 6)},
        "counted_shell_sizes": {k: len(geometry.generate_shell(k)) for k in range(1, 6)},
        "free_energy_levels": {k: sorted(levels[k]) for k in levels},
        "levels_per_shell": {k: len(levels[k]) for k in levels},
        "levels_shared_across_shells": shared,
        "shell_one_is_one_sphere": levels[1] == {4.0: 12},
        "shell_two_spans_three_spheres": levels[2] == {8.0: 6, 12.0: 24, 16.0: 12},
    }
    gauge_shells = range(0, 13)
    values["gauge_equals_shell_index"] = all(
        geometry.cuboctahedral_gauge(site) == k
        for k in gauge_shells
        for site in geometry.generate_shell(k)
    )
    values["gauge_checked_through_shell"] = max(gauge_shells)
    extremes = {}
    for k in range(1, 13):
        energies = [geometry.squared_norm(site) // 2 for site in geometry.generate_shell(k)]
        top = max(energies)
        extremes[k] = (min(energies), top, sum(1 for value in energies if value == top))
    values["free_energy_range_per_shell"] = {k: extremes[k][:2] for k in range(1, 6)}
    values["maximum_is_4k2_at_twelve_vertices"] = all(
        extremes[k][1] == 4 * k * k and extremes[k][2] == 12 for k in extremes
    )
    values["minimum_is_2k2_or_2k2_plus_2"] = all(
        extremes[k][0] == 2 * k * k + (2 if k % 2 else 0) for k in extremes
    )
    passed = (
        values["shell_sizes"] == values["counted_shell_sizes"] == {1: 12, 2: 42, 3: 92, 4: 162, 5: 252}
        and values["levels_per_shell"] == {1: 1, 2: 3, 3: 4, 4: 7, 5: 7}
        and shared.get(36.0) == [3, 4]
        and shared.get(52.0) == [4, 5]
        and values["shell_one_is_one_sphere"]
        and values["shell_two_spans_three_spheres"]
        and values["gauge_equals_shell_index"] is True
        and values["maximum_is_4k2_at_twelve_vertices"] is True
        and values["minimum_is_2k2_or_2k2_plus_2"] is True
    )
    return passed, values


# ---------------------------------------------------------------------------
# Complementarity probes
# ---------------------------------------------------------------------------


def _probe_quadray_ivm() -> tuple[bool, dict]:
    """The exact quadray/IVM geometry the synergetic side supplies."""
    moves = geometry.neighbor_moves()
    pairwise = min(
        geometry.squared_distance(first, second)
        for index, first in enumerate(moves)
        for second in moves[index + 1 :]
    )
    shells = tuple(int(geometry.shell_count(k)) for k in range(5))
    enumerated = tuple(len(geometry.generate_shell(k)) for k in range(5))
    cumulative = tuple(int(geometry.cumulative_count(k)) for k in range(5))
    ratios = sorted(float(value) for value in geometry.VOLUME_RATIOS.values())
    values = {
        "neighbor_moves": len(moves),
        "distinct_moves": len(set(moves)),
        "squared_norms": sorted({float(geometry.squared_norm(move)) for move in moves}),
        "min_pairwise_squared_distance": float(pairwise),
        "shell_counts": shells,
        "enumerated_shell_counts": enumerated,
        "cumulative_counts": cumulative,
        "volume_ratios_present": ratios,
        "volume_ratios_derived": {
            name: str(value) for name, value in sorted(geometry.derived_volume_ratios().items())
        },
        "derived_table_equals_declared": geometry.derived_volume_ratios() == geometry.VOLUME_RATIOS,
        "max_shell": geometry.MAX_SHELL,
        "neighbor_squared_distance": geometry.NEIGHBOR_SQUARED_DISTANCE,
    }
    passed = (
        values["neighbor_moves"] == 12
        and values["distinct_moves"] == 12
        and values["squared_norms"] == [8.0]
        and abs(values["min_pairwise_squared_distance"] - 8.0) < 1e-9
        and shells == enumerated == (1, 12, 42, 92, 162)
        and cumulative == (1, 13, 55, 147, 309)
        and {1.0, 3.0, 4.0, 6.0, 20.0, 96.0} <= set(ratios)
        and values["derived_table_equals_declared"] is True
        and values["max_shell"] == 32
        and values["neighbor_squared_distance"] == 8
    )
    return passed, values


def _probe_equilibrium_propagation() -> tuple[bool, dict]:
    """Theorem 21: two settled states supply the parameter gradient the repair law lacks."""
    values = control_module.equilibrium_propagation_identity()
    passed = (
        values["closed_form_error"] < 1e-9
        and abs(values["distance_to_gradient"] - values["expected_distance_to_gradient"]) < 1e-12
        and values["distance_to_gradient"] > 0.0
    )
    return passed, values


def _probe_replicator() -> tuple[bool, dict]:
    """Selection as a potential flow, a comparison class rather than a shared mechanism."""
    values = dict(control_module.replicator_potential_check())
    control = control_module.replicator_asymmetric_control()
    values["symmetric_payoff_points_checked"] = control["symmetric_points"]
    values["symmetric_payoffs_are_gradients"] = control["symmetric_all_gradient"]
    values["asymmetric_payoff_is_not_a_gradient"] = control["asymmetric_not_gradient"]
    passed = (
        values["symmetric_payoffs_are_gradients"] is True
        and values["asymmetric_payoff_is_not_a_gradient"] is True
        and values["field_equals_mobility_gradient"] is True
        and values["rates_agree"] is True
        and values["potential_rate_direct"] > 0.0
    )
    return passed, values


def _probe_kissing_twelve() -> tuple[bool, dict]:
    """The twelve-port carrier sits on the kissing-number-twelve configuration."""
    carrier = _carrier()
    moves = geometry.neighbor_moves()
    values = {
        "carrier_ports": carrier.port_count,
        "carrier_degrees": sorted({carrier.degree(port) for port in range(carrier.port_count)}),
        "ivm_shell_one": int(geometry.shell_count(1)),
        "ivm_moves": len(moves),
        "move_squared_norm": float(geometry.squared_norm(moves[0])),
        "geometric_adjacency_isomorphic": carrier_module.carrier_isomorphism(
            carrier, carrier_module.geometric_carrier()
        )
        is not None,
    }
    kissing = geometry.kissing_arrangements()
    values["cuboctahedron_edge_over_radius_squared"] = str(
        kissing["cuboctahedron_edge_over_radius_squared"]
    )
    values["icosahedron_edge_over_radius_squared"] = str(
        kissing["icosahedron_edge_over_radius_squared"]
    )
    values["icosahedron_edge_over_radius_squared_float"] = kissing["icosahedron_ratio_float"]
    values["cuboctahedral_neighbours_touch"] = kissing["cuboctahedron_is_rigid"]
    values["icosahedral_neighbours_have_slack"] = kissing["icosahedron_has_slack"]
    passed = (
        values["carrier_ports"] == 12
        and values["carrier_degrees"] == [5]
        and values["ivm_shell_one"] == 12
        and values["ivm_moves"] == 12
        and abs(values["move_squared_norm"] - 8.0) < 1e-9
        and values["geometric_adjacency_isomorphic"]
        and values["cuboctahedron_edge_over_radius_squared"] == "1"
        and values["icosahedron_edge_over_radius_squared"] == "2 - 2/5 sqrt 5"
        and values["cuboctahedral_neighbours_touch"] is True
        and values["icosahedral_neighbours_have_slack"] is True
    )
    return passed, values


def _probe_icosahedron_tetravolume() -> tuple[bool, dict]:
    """Synergetics' volume accounting, applied to the carrier's own solid."""
    measured = dict(geometry.icosahedron_tetravolume())
    exact_parts = measured["exact_squared_parts"]
    values = {
        "vertices": measured["vertices"],
        "edges": measured["edges"],
        "faces": measured["faces"],
        "faces_reoriented": measured["faces_reoriented"],
        "edge": measured["edge"],
        "edge_spread": measured["edge_spread"],
        "volume": measured["volume"],
        "tetraunit": measured["tetraunit"],
        "tetravolume": measured["tetravolume"],
        "closed_form_tetravolume": measured["closed_form_tetravolume"],
        "tetravolume_residual": measured["tetravolume_residual"],
        "tetravolume_squared": measured["tetravolume_squared"],
        "exact_squared": measured["exact_squared"],
        "exact_squared_parts": exact_parts,
        "exact_squared_residual": measured["exact_squared_residual"],
        "integral": measured["integral"],
    }
    exact = geometry.icosahedron_tetravolume_exact()
    values["exact_volume_edge_two"] = str(exact["volume"])
    values["exact_tetravolume_over_sqrt2"] = str(exact["sqrt2_coefficient"])
    values["exact_square"] = str(exact["square"])
    values["exact_facets"] = exact["facets"]
    values["exact_matches_measured"] = abs(float(exact["tetravolume"]) - float(measured["tetravolume"])) < 1e-12
    passed = (
        values["vertices"] == 12
        and values["edges"] == 30
        and values["faces"] == 20
        and values["faces_reoriented"] == 10
        and float(values["edge_spread"]) < 1e-12
        and float(values["tetravolume_residual"]) < 1e-9
        and float(values["exact_squared_residual"]) < 1e-9
        and abs(float(values["closed_form_tetravolume"]) - float(values["tetravolume"])) < 1e-9
        and tuple(exact_parts) == (175, 75)
        and values["integral"] is False
        and values["exact_tetravolume_over_sqrt2"] == "15/2 + 5/2 sqrt 5"
        and values["exact_square"] == "175 + 75 sqrt 5"
        and values["exact_facets"] == 20
        and values["exact_matches_measured"] is True
    )
    return passed, values


def _probe_packing_density() -> tuple[bool, dict]:
    """The IVM's Voronoi cell, its tetravolume, and the packing density they imply."""
    cell = geometry.voronoi_cell()
    density = dict(geometry.packing_density())
    values = {
        "cell_vertex_count": cell["vertex_count"],
        "cell_face_count": cell["face_count"],
        "cell_edge_count": cell["edge_count"],
        "cell_euler_characteristic": cell["euler_characteristic"],
        "cell_all_faces_quadrilateral": cell["all_faces_quadrilateral"],
        "cell_all_faces_rhombic": cell["all_faces_rhombic"],
        "cell_feasible": cell["feasible"],
        "cell_orientation_flips": cell["orientation_flips"],
        "cell_volume": cell["volume"],
        "cell_pyramid_volume": cell["pyramid_volume"],
        "cell_volume_residual": cell["volume_residual"],
        "cell_tetravolume": density["cell_tetravolume"],
        "cell_tetravolume_matches_table": density["cell_tetravolume_matches_table"],
        "ball_tetravolume": density["ball_tetravolume"],
        "density": density["density"],
        "closed_form": density["closed_form"],
        "density_residual": density["density_residual"],
        "inradius_is_half_the_prime_vector": density["inradius_is_half_the_prime_vector"],
    }
    cell_exact = geometry.convex_volume(geometry.ivm_solid_vertices()["rhombic_dodecahedron"])
    values["cell_volume_exact"] = str(cell_exact)
    values["cell_tetravolume_exact"] = str(Fraction(cell_exact) / geometry.IVM_TETRA_UNIT)
    # The inscribed ball has radius sqrt 2 (half the prime vector 2 sqrt 2), so its
    # volume is (4/3) pi 2 sqrt 2 and its tetravolume (3/8) of that, exactly pi sqrt 2.
    values["ball_tetravolume_exact"] = "pi sqrt 2"
    values["density_exact"] = "pi sqrt 2 / 6"
    passed = (
        values["cell_vertex_count"] == 14
        and values["cell_face_count"] == 12
        and values["cell_edge_count"] == 24
        and values["cell_euler_characteristic"] == 2
        and values["cell_all_faces_quadrilateral"] is True
        and values["cell_all_faces_rhombic"] is True
        and values["cell_feasible"] is True
        and values["cell_orientation_flips"] == 0
        and float(values["cell_volume_residual"]) < 1e-9
        and abs(float(values["cell_tetravolume"]) - 6.0) < 1e-9
        and values["cell_tetravolume_matches_table"] is True
        and values["inradius_is_half_the_prime_vector"] is True
        and float(values["density_residual"]) < 1e-15
        and abs(float(values["ball_tetravolume"]) - math.pi * math.sqrt(2.0)) < 1e-9
        and values["cell_volume_exact"] == "16"
        and values["cell_tetravolume_exact"] == "6"
        and abs(math.pi * math.sqrt(2.0) / 6.0 - float(values["density"])) < 1e-15
    )
    return passed, values


def _probe_tetravolume_arithmetic() -> tuple[bool, dict]:
    """Three currencies, and the one place the synergetic table cannot stay integral."""
    identity = np.eye(3)
    origin = (0.0, 0.0, 0.0)
    move = geometry.neighbor_moves()[0]
    kl = free_energy_module.gaussian_kl(origin, identity, geometry.to_xyz(move), identity)
    measured = dict(geometry.icosahedron_tetravolume())
    density = dict(geometry.packing_density())
    table = {name: str(value) for name, value in sorted(geometry.derived_volume_ratios().items())}
    exact_ico = geometry.icosahedron_tetravolume_exact()
    values = {
        "fep_currency": "nats",
        "fep_one_step_kl": kl,
        "fep_one_step_is_integral": bool(abs(kl - round(kl)) < 1e-12),
        "oph_currency": "integer loads; spectral data in Q(sqrt 5)",
        "oph_shell_gap": repair_module.shell_gap(),
        "oph_shell_gap_exact": "(5 - sqrt 5) / 60",
        "oph_shell_gap_matches_exact": abs(
            repair_module.shell_gap() - (5.0 - math.sqrt(5.0)) / 60.0
        ) < 1e-15,
        "oph_shell_gap_is_smaller_than_one": repair_module.shell_gap() < 1.0,
        "synergetics_currency": "tetravolumes of the prime vector",
        "synergetics_table": table,
        "synergetics_table_is_integral": all(
            value.denominator == 1 for value in geometry.derived_volume_ratios().values()
        ),
        "synergetics_table_derived_exactly": geometry.derived_volume_ratios()
        == geometry.VOLUME_RATIOS,
        "synergetics_icosahedron_exact": f"sqrt 2 ({exact_ico['sqrt2_coefficient']})",
        "synergetics_icosahedron_tetravolume": float(measured["tetravolume"]),
        "synergetics_icosahedron_is_integral": measured["integral"],
        "synergetics_icosahedron_squared_parts": measured["exact_squared_parts"],
        "ivm_cell_tetravolume": float(density["cell_tetravolume"]),
        "ivm_density_is_below_one": float(density["density"]) < 1.0,
    }
    passed = (
        abs(kl - 4.0) < 1e-12
        and values["fep_one_step_is_integral"] is True
        and values["synergetics_table_is_integral"] is True
        and values["synergetics_table_derived_exactly"] is True
        and values["oph_shell_gap_matches_exact"] is True
        and values["synergetics_icosahedron_is_integral"] is False
        and tuple(values["synergetics_icosahedron_squared_parts"]) == (175, 75)
        and abs(float(values["ivm_cell_tetravolume"]) - 6.0) < 1e-9
        and values["ivm_density_is_below_one"] is True
    )
    return passed, values


# ---------------------------------------------------------------------------
# The registry
# ---------------------------------------------------------------------------

def _probe_effective_resistance() -> tuple[bool, dict]:
    """The FEP's metric and the electrical reading of the same matrix agree exactly."""
    carrier = _carrier()
    adjacency = _carrier_adjacency(carrier)
    ladder = network_module.resistance_ladder(adjacency)
    exact = network_module.kirchhoff_index(adjacency)
    from_spectrum = carrier.port_count * invariant_module.reciprocal_spectrum_sum(carrier)
    values = {
        "ladder": tuple((distance, str(resistance), pairs) for distance, resistance, pairs in ladder),
        "adjacent_resistance": str(ladder[0][1]),
        "distance_two_resistance": str(ladder[1][1]),
        "antipodal_resistance": str(ladder[2][1]),
        "kirchhoff_index_exact": str(exact),
        "kirchhoff_index_from_spectrum": str(from_spectrum),
        "ports": carrier.port_count,
        "seams": carrier.seam_count,
    }
    passed = (
        values["ladder"] == ((1, "11/30", 30), (2, "7/15", 30), (3, "1/2", 6))
        and exact == from_spectrum == 28
        and values["seams"] == 30
    )
    return passed, values


def _probe_vertex_shells() -> tuple[bool, dict]:
    """Shells around a port are not the lattice's shells: the graph is not the IVM."""
    carrier = _carrier()
    adjacency = _carrier_adjacency(carrier)
    shells = network_module.distance_distribution(adjacency)
    onward, back = network_module.intersection_array(adjacency)
    outer = shells[1:]
    lattice = tuple(geometry.shell_count(step) for step in range(1, len(shells)))
    values = {
        "shell_sizes": shells,
        "outer_shell_sizes": outer,
        "lattice_shell_sizes": lattice,
        "onward": onward,
        "back": back,
        "diameter": len(shells) - 1,
        "ports": carrier.port_count,
    }
    passed = (
        shells == (1, 5, 5, 1)
        and outer == (5, 5, 1)
        and lattice == (12, 42, 92)
        and onward == (2, 1, 0)
        and back == (1, 2, 5)
        and outer != lattice
    )
    return passed, values


def _probe_subgraph_census() -> tuple[bool, dict]:
    """The carrier's subgraphs reproduce the icosahedron's own counts exactly."""
    carrier = _carrier()
    census = network_module.subgraph_census(_carrier_adjacency(carrier))
    values = {**census, "faces": len(carrier.faces), "ports": carrier.port_count}
    passed = (
        census
        == {
            "triangles": 20,
            "pentagons": 12,
            "perfect_matchings": 125,
            "hamiltonian_cycles": 1280,
        }
        and census["triangles"] == values["faces"]
        and census["pentagons"] == carrier.port_count
    )
    return passed, values


def _probe_angle_defect() -> tuple[bool, dict]:
    """Gauss-Bonnet on the carrier: a defect of pi/3 at each port, 4 pi in total."""
    carrier = _carrier()
    defects = geometry.angle_defects(carrier.faces, carrier.port_count)
    values = {
        "defect_per_port_over_pi": str(defects["per_vertex"][0]),
        "total_defect_over_pi": str(defects["total"]),
        "euler_characteristic_times_two": defects["euler_characteristic_times_two"],
        "faces_per_port": defects["faces_per_vertex"][0],
        "uniform": len(set(defects["faces_per_vertex"])) == 1,
        "sum_is_gauss_bonnet": defects["total"] == Fraction(defects["euler_characteristic_times_two"], 1),
    }
    passed = (
        values["defect_per_port_over_pi"] == "1/3"
        and values["total_defect_over_pi"] == "4"
        and values["euler_characteristic_times_two"] == 4
        and values["faces_per_port"] == 5
        and values["uniform"] is True
        and values["sum_is_gauss_bonnet"] is True
    )
    return passed, values


def _probe_twelve_around_one_families() -> tuple[bool, dict]:
    """`Twelve around one` names a family; the counting invariant tells its members apart."""
    carrier = _carrier()
    icosahedron = _carrier_adjacency(carrier)
    cuboctahedron = [list(row) for row in geometry.cuboctahedron_adjacency()]
    ico_trees = invariant_module.spanning_tree_number_of_adjacency(icosahedron)
    cub_trees = invariant_module.spanning_tree_number_of_adjacency(cuboctahedron)
    census = network_module.subgraph_census(cuboctahedron)
    values = {
        "icosahedron_spanning_trees": ico_trees,
        "cuboctahedron_spanning_trees": cub_trees,
        "cuboctahedron_spanning_trees_is_square": math.isqrt(cub_trees) ** 2 == cub_trees,
        "icosahedron_degree": max(sum(row) for row in icosahedron),
        "cuboctahedron_degree": max(sum(row) for row in cuboctahedron),
        "cuboctahedron_triangles": census["triangles"],
        "cuboctahedron_squares": len(geometry.cuboctahedron_faces()["squares"]),
        "cuboctahedron_perfect_matchings": census["perfect_matchings"],
    }
    passed = (
        ico_trees == 5184000
        and cub_trees == 331776
        and values["cuboctahedron_spanning_trees_is_square"] is True
        and values["icosahedron_degree"] == 5
        and values["cuboctahedron_degree"] == 4
        and census["triangles"] == 8
        and census["perfect_matchings"] == 32
        and ico_trees != cub_trees
    )
    return passed, values


def _probe_dual_cells() -> tuple[bool, dict]:
    """The twelve-around-one cell and the cuboctahedron are duals, edge for edge."""
    cell = geometry.voronoi_cell()
    cuboctahedron = [list(row) for row in geometry.cuboctahedron_adjacency()]
    line_graph = [list(row) for row in geometry.cube_line_graph_adjacency()]
    faces = geometry.cuboctahedron_faces()
    cuboctahedron_faces = len(faces["triangles"]) + len(faces["squares"])
    cuboctahedron_edges = len(geometry.cuboctahedron_edges())
    ratio = Fraction(geometry.cuboctahedron_volume()) / Fraction(6)
    values = {
        "cell_vertices": cell["vertex_count"],
        "cell_faces": cell["face_count"],
        "cell_edges": cell["edge_count"],
        "cell_tetravolume": geometry.voronoi_cell_tetravolume(),
        "cell_euler_characteristic": cell["euler_characteristic"],
        "cuboctahedron_vertices": len(cuboctahedron),
        "cuboctahedron_faces": cuboctahedron_faces,
        "cuboctahedron_edges": cuboctahedron_edges,
        "cuboctahedron_tetravolume": str(geometry.cuboctahedron_volume()),
        "volume_ratio": str(ratio),
        "is_line_graph_of_cube": network_module.graphs_are_isomorphic(line_graph, cuboctahedron),
        "euler_characteristics": (
            cell["euler_characteristic"],
            len(cuboctahedron) - cuboctahedron_edges + cuboctahedron_faces,
        ),
    }
    passed = (
        values["cell_vertices"] == 14
        and values["cell_faces"] == 12
        and values["cell_edges"] == 24
        and values["cuboctahedron_vertices"] == 12
        and values["cuboctahedron_faces"] == 14
        and cuboctahedron_edges == 24
        and abs(values["cell_tetravolume"] - 6.0) < 1e-6
        and values["cuboctahedron_tetravolume"] == "20"
        and values["volume_ratio"] == "10/3"
        and values["is_line_graph_of_cube"] is True
        and values["euler_characteristics"] == (2, 2)
    )
    return passed, values


def _probe_conjugation_norm() -> tuple[bool, dict]:
    """Conjugation fixes the counting invariant and moves the volume to its conjugate."""
    carrier = _carrier()
    trees = invariant_module.spanning_tree_number(carrier)
    parts = geometry.tetravolume_ratio_squared_exact()
    squared_norm = geometry.quadratic_norm(parts)
    exact = geometry.icosahedron_tetravolume_exact()
    conjugated_squared = float(exact["square"].conjugate())  # type: ignore[union-attr]
    tetravolume = float(geometry.icosahedron_tetravolume()["tetravolume"])
    conjugate_volume = float(exact["conjugate_tetravolume"])
    values = {
        "spanning_trees": trees,
        "tree_number_is_an_integer": trees == int(trees),
        "squared_parts": (int(parts[0]), int(parts[1])),
        "squared_norm": int(squared_norm),
        "galois_norm_of_tetravolume": str(geometry.tetravolume_galois_norm()),
        "conjugated_squared": conjugated_squared,
        "conjugate_volume": conjugate_volume,
        "volume_times_conjugate": tetravolume * conjugate_volume,
        "volume_times_conjugate_exact": str(exact["norm"]),
        "conjugate_over_sqrt2": str(exact["conjugate_sqrt2_coefficient"]),
        "volume_is_irrational": exact["sqrt2_coefficient"].b != 0,  # type: ignore[union-attr]
    }
    passed = (
        trees == 5184000
        and values["squared_parts"] == (175, 75)
        and values["squared_norm"] == 2500
        and values["galois_norm_of_tetravolume"] == "50"
        and abs(values["volume_times_conjugate"] - 50.0) < 1e-9
        and values["volume_times_conjugate_exact"] == "50"
        and values["conjugate_over_sqrt2"] == "15/2 - 5/2 sqrt 5"
        and values["volume_is_irrational"] is True
    )
    return passed, values


_CARRIER_SOURCE = "OPH evidence/observer_dynamics_20260925/oph-physics-sim (ORIENTED_BASE_FACES, seams, laplacian)"
_LAW_SOURCE = "OPH model/observer/native.py:45-52 and the 'Agreement and Surprise' paper, Section 6.2"
_PAPER_SOURCE = "Mueller, 'Agreement and Surprise: Local Repair, Dynamics and the Limits of Equilibrium', private draft 22 October 2026"
_QUADRAY_SOURCE = (
    "QuadMath src/quadmath/core/quadray.py, src/quadmath/lattice/omni_numbering.py and "
    "the polyhedra table in docs/manuscript/03_quadray_methods.md"
)
_SHEAF_SOURCE = f"{_PAPER_SOURCE}, Proposition 4 and Section 7.2 (the cellular sheaf over the carrier complex)"

CORRESPONDENCES: tuple[Correspondence, ...] = (
    Correspondence(
        key="carrier-port-chart",
        kind=ISOMORPHISM,
        pair=SYNERGETICS_OPH,
        claim="OPH's twelve-port carrier has twelve ports, thirty seams and twenty faces, is five-regular, and its port chart is the one OPH declares.",
        oph="The executable observer's twelve-port carrier.",
        fep="No fixed carrier; the FEP leaves the graph to the model.",
        synergetics="Icosahedron: twelve vertices, thirty edges, twenty faces.",
        status="exact",
        source=_CARRIER_SOURCE,
        probe=_probe_port_chart,
    ),
    Correspondence(
        key="spectral-gap",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="One matrix carries both readings: the carrier Laplacian is 5I - A; in OPH's reading it generates the mean repair operator I - L / 60 and is the quadratic form of the seam-disagreement energy (1/2) sum over seams of (x_a - x_b)^2 (the native law's own energy, sum x^2 / 2, is a different form), and unweighted it is the precision the FEP reading factorises through; its first positive eigenvalue is 5 - sqrt(5) in both readings, and the normalised gap lambda_1 / (2 |E|) is exactly (5 - sqrt(5)) / 60.",
        oph="The generator of the mean repair operator, whose gap is the native law's spectral certificate.",
        fep="The same matrix as a Gaussian precision, read unweighted; the FEP-side sheaf Laplacian is its weighted generalisation, so the gap is inherited from the shared object, not re-derived.",
        synergetics="The icosahedron's graph: a twelve-vertex kissing arrangement, but not the A3 coordination shell, which is the cuboctahedron.",
        status="exact",
        source=_CARRIER_SOURCE,
        probe=_probe_spectral_gap,
    ),
    Correspondence(
        key="shell-counts",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="The finite state space both oph-fep readings evaluate on is the total-two load shell: 78 configurations, of which 66 are binary; the enumeration agrees with the closed forms C(13, 2) = 78 (stars and bars over twelve ports) and C(12, 2) = 66, and every transition census and FEP-side comparison in this package ranges over exactly these states.",
        oph="The finite shell the executable observer runs on.",
        fep="No independent FEP count is claimed: the entry fixes the finite domain on which the oph-fep comparisons run, and its second reading is that domain, not a number.",
        synergetics="No counterpart: a load shell, a level set of total load, is not a coordination shell of the lattice.",
        status="exact",
        source="OPH model/observer/native.py (states) and the paper, Section 6.2",
        probe=_probe_shell_counts,
    ),
    Correspondence(
        key="native-repair-law",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="The native seam law conserves total load on every transition and lowers the quadratic energy by exactly (d^2 - (|d| mod 2)) / 4.",
        oph="The native repair law.",
        fep="A local update on a bounded interface.",
        synergetics="A minimal-energy local readjustment of the honeycomb.",
        status="exact",
        source=_LAW_SOURCE,
        probe=_probe_native_law,
    ),
    Correspondence(
        key="variational-equivalence",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="On the linear-Gaussian domain three computations of the settled interior coincide: the harmonic extension (an interior linear solve), the constrained minimiser of (1/2) z^T L z with the records fixed (a separate KKT saddle-point solve, which is OPH's A3 constrained relative-entropy minimisation read in the mean, because the relative entropy of a mean shift at fixed precision L is (1/2) z^T L z), and the limit of the predictive-coding gradient flow from any start; the free energy computed directly as energy minus entropy equals the surprisal there.",
        oph="A3: the least informative state relative to a reference under observer-visible constraints.",
        fep="Perception as the descent of variational free energy in the means.",
        synergetics="A minimal-energy configuration of the coordination framework.",
        status="exact-on-domain",
        source=f"{_PAPER_SOURCE}, Theorem 7 items 2 and 3",
        probe=_probe_variational_equivalence,
    ),
    Correspondence(
        key="free-energy-dictionary",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="Surprisal is the disagreement no repair removes, Bayesian surprise is the disagreement repair will remove, the update residual pins the Bayesian surprise between the extreme interior eigenvalues, and the interior free energy's drop from the old minimiser to the new one equals the Bayesian surprise on this domain.",
        oph="Repair reduces disagreement; records survive comparison.",
        fep="The variational bound, its gap, and Bayesian surprise.",
        synergetics="No counterpart; the geometry supplies no probabilistic reading.",
        status="exact-on-domain",
        source=f"{_PAPER_SOURCE}, Theorem 7 item 5 and Table 7",
        probe=_probe_free_energy_dictionary,
    ),
    Correspondence(
        key="sheaf-agreement",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="Patch agreement is the global-section condition of a cellular sheaf over the carrier complex, and the sheaf Laplacian is the carrier Laplacian.",
        oph="Public records are the states that survive every comparison.",
        fep="Message passing across a blanket; the sheaf Laplacian is the precision-weighted interface energy.",
        synergetics="Heterogeneous local states agreeing through interface maps.",
        status="exact",
        source=f"{_PAPER_SOURCE}, Section 7.2 and Proposition 4",
        probe=_probe_sheaf_agreement,
    ),
    Correspondence(
        key="certificate-species",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="The native law's sharp certificate is a spectral gap: the mean repair operator is symmetric, nonnegative and doubly stochastic by construction from the five-regular Laplacian, with second eigenvalue 1 - (5 - sqrt(5)) / 60, and its induced infinity norm is exactly one, so no sup-norm contraction certifies it.",
        oph="The capacity-weighted two-norm spectral certificate.",
        fep="Precision-weighted quadratic norms for the linearised update.",
        synergetics="No counterpart.",
        status="exact",
        source=f"{_PAPER_SOURCE}, Section 6.2 ('which certificate a law admits')",
        probe=_probe_certificate_species,
    ),
    Correspondence(
        key="tree-number",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="One determinant, two readings: the reduced carrier Laplacian's determinant is 5184000 — the number of spanning trees of the seam graph (the seam sets that connect every port without a cycle) in OPH's reading, by Kirchhoff's matrix-tree theorem, and, as the determinant of the same matrix read as a Gaussian precision grounded at one port, the normalisation of the FEP's linear-Gaussian interface integral — and it equals the product of the eleven nonzero eigenvalues over 12.",
        oph="The seam structure of the twelve-port carrier, whose Laplacian is the repair energy's quadratic form.",
        fep="The same matrix as a Gaussian precision, whose determinant normalises the linear-Gaussian interface integral.",
        synergetics="The icosahedron's own graph, counted by its seam structure.",
        status="exact",
        source=_CARRIER_SOURCE,
        probe=_probe_tree_number,
    ),
    Correspondence(
        key="gmrf-independence",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="The seam graph OPH repairs along is exactly the conditional-independence graph of the precision read from the same matrix: two ports are conditionally dependent given the rest exactly when they share a seam and independent exactly when they do not, each port's Markov blanket is its five neighbours, every seam's pairwise precision is minus one, and the variance of a port's deviation from the carrier mean (the diagonal of L^+) is exactly 7/36. The sparsity half holds by construction, since the precision is built from the adjacency; the content is the reading, not the computation.",
        oph="The declared seam graph is the independence graph: a repair touches one seam, which is exactly one conditional-dependence edge.",
        fep="A Markov blanket is the graph neighbourhood of a port, and the interface energy is the precision-quadratic form on the linear-Gaussian domain.",
        synergetics="Every vertex of the icosahedron is five-coordinated.",
        status="exact",
        source=f"{_CARRIER_SOURCE}; {_PAPER_SOURCE}, Section 7.2",
        probe=_probe_gmrf_independence,
    ),
    Correspondence(
        key="sheaf-cohomology",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="Patch agreement in OPH is the global-section condition of the carrier's scalar sheaf on the seam graph, and the obstruction to it is the space the FEP's harmonic-extension residual lives in: dim H^0 = 1, dim H^1 = 19 — the graph's cycle rank, nineteen independent edge disagreements no coboundary explains — with dim H^0 - dim H^1 = 12 - 30. Adding the twenty faces as two-cells kills all nineteen: every obstruction is a face boundary, the surface's Betti numbers are (1, 0, 1), its Euler characteristic is 2, and OPH's oriented faces sum to a fundamental cycle.",
        oph="A public record is a global section; what prevents agreement is not counted by the seam count alone.",
        fep="The harmonic-extension problem is well posed up to a constant; the nineteen dimensions are edge data that no node assignment can realise, the cokernel of the coboundary, while every realised disagreement lies in its eleven-dimensional image.",
        synergetics="A closed triangulated surface: the graph's Euler characteristic is 12 - 30 = -18, the surface's is 12 - 30 + 20 = 2, and the faces account for the difference.",
        status="exact",
        source=_SHEAF_SOURCE,
        probe=_probe_sheaf_cohomology,
    ),
    Correspondence(
        key="icosahedral-symmetry",
        kind=ISOMORPHISM,
        pair=SYNERGETICS_OPH,
        claim="The carrier's automorphism group has order 120, is closed under composition, acts transitively on the twelve ports, the thirty seams and the twenty faces, and splits into sixty rotations and sixty reflections detected from OPH's own oriented faces.",
        oph="The port chart's oriented faces decide which automorphisms preserve orientation.",
        fep="The group is a property of the chosen interface, not of any free-energy functional.",
        synergetics="The full icosahedral group of the icosahedron.",
        status="exact",
        source=_CARRIER_SOURCE,
        probe=_probe_icosahedral_symmetry,
    ),
    Correspondence(
        key="eigenspace-irreducibility",
        kind=ISOMORPHISM,
        pair=SYNERGETICS_OPH,
        claim="The carrier Laplacian's eigenspaces of dimensions 3, 5 and 3 are the icosahedral group's non-trivial real irreducible representations acting on the carrier, so the multiplicities 1, 3, 5 and 3 of the spectrum 0, (5 - sqrt 5), 6 and (5 + sqrt 5) are forced by the icosahedral group rather than computed into OPH's matrix: each eigenspace is one irreducible, and the multiplicities are the representation dimensions. The values themselves are specific to 5I - A: the permutation representation is multiplicity-free, so every invariant symmetric matrix shares these eigenspaces with its own four eigenvalues.",
        oph="The native law's spectral certificate, whose multiplicities the law does not choose.",
        fep="A precision spectrum that no relabelling of the interface can change.",
        synergetics="The icosahedral group's real irreducible representations of dimensions 1, 3, 3 and 5 appear as the eigenspaces.",
        status="exact",
        source=_CARRIER_SOURCE,
        probe=_probe_eigenspace_irreducibility,
    ),
    Correspondence(
        key="law-equivariance",
        kind=ISOMORPHISM,
        pair=SYNERGETICS_OPH,
        claim="The native repair law is equivariant transition by transition: for all 120 automorphisms, all 78 total-two configurations, both orientations of all thirty seams and both coins, repairing the relabelled configuration on the relabelled seam gives the relabelled result, so the multiset of strict, zero and no-op outcomes is a group invariant, with four configuration orbits and three distinct outcome multisets.",
        oph="The repair law commutes with any relabelling of the carrier's ports.",
        fep="A local update whose outcome statistics a relabelling cannot change.",
        synergetics="The icosahedral group acting on the coordination structure.",
        status="exact",
        source=_LAW_SOURCE,
        probe=_probe_law_equivariance,
    ),
    Correspondence(
        key="kl-lattice-metric",
        kind=ISOMORPHISM,
        pair=FEP_SYNERGETICS,
        claim="Equal-covariance Gaussian KL is half the squared Mahalanobis distance, and with unit covariance in the integer quadray embedding the IVM's first shell sits at squared distance eight, so one lattice step costs exactly four nats and the sites k m at squared norm 8 k^2 cost 4 k^2: the information cost is the synergetic metric over two. The four nats are a statement about that embedding's unit of length; rescaling the covariance rescales every cost by the same factor.",
        oph="No counterpart; OPH's own step cost is the change of its quadratic seam energy.",
        fep="The linear-Gaussian KL of a mean shift at fixed covariance.",
        synergetics="The IVM's squared distance in quadray coordinates.",
        status="exact",
        source=f"{_PAPER_SOURCE}, Theorem 7 item 5; {_QUADRAY_SOURCE}",
        probe=_probe_kl_lattice_metric,
    ),
    Correspondence(
        key="pr01-obstruction",
        kind=DISAGREEMENT,
        pair=OPH_FEP,
        claim="The native repair law is not confluent and accepted repairs need not lower a measure: 600 of the 4680 shell transitions move at mismatch one with zero energy decrease, and a one-unit walk keeps energy 1/2 while changing state.",
        oph="PR-01 'confluent terminating repair contract' is recorded as obstructed; the consensus theorems hold for the integer and transactional contracts, not the native law.",
        fep="Active inference is a descent on a bound; the identification with this law therefore fails on the mismatch-one domain.",
        synergetics="No counterpart.",
        status="obstructed",
        source="OPH model/PREMISES.md:95 (PR-01) and Lean/ObserverPatchHolography/Primitives.lean:1198 (demoCarrier_not_confluent)",
        probe=_probe_pr01,
    ),
    Correspondence(
        key="free-energy-name-collision",
        kind=DISAGREEMENT,
        pair=OPH_FEP,
        claim="The interface energy minimised over the interior and the variational free energy at the posterior are different objects that differ by a model-dependent constant, the Gaussian-integral prefactor (1/2) logdet L_xx - (d_x / 2) log 2 pi of the unnormalised density, computed both through the surprisal and directly as energy minus entropy; identifying the two numerically is a gauge choice rather than a theorem.",
        oph="Native thermodynamic free energy: Gibbs states and Landauer bounds on records.",
        fep="Variational free energy: a KL bound on surprisal, not a thermodynamic potential.",
        synergetics="No counterpart.",
        status="obstructed",
        source=f"{_PAPER_SOURCE}, Section 7.1 and Theorem 7 item 3; review caution 1 on disambiguating 'free energy'",
        probe=_probe_energy_name_collision,
    ),
    Correspondence(
        key="mediator-obstruction",
        kind=DISAGREEMENT,
        pair=OPH_FEP,
        claim="Eliminating an interior mediator leaves a genuinely three-way boundary term: the reduced precision is (1/3) v v^T + epsilon I with v = (1, 1, -1), and a sum of pairwise semidefinite terms reproduces it if and only if epsilon >= 1/3: necessity is the trace inequality (each pairwise block with off-diagonal 1/3 has diagonal sum at least 2/3, against a total trace of 1 + 3 epsilon), and sufficiency is an explicit decomposition with rank-one blocks at the threshold.",
        oph="A local patch can be unable to repair although a consistent world exists.",
        fep="A blanket reduction that destroys the pairwise predictive-coding description.",
        synergetics="No counterpart.",
        status="obstructed",
        source=f"{_PAPER_SOURCE}, Proposition 8; Hansen and Ghrist's sheaf Kron reduction obstruction",
        probe=_probe_mediator,
    ),
    Correspondence(
        key="stationary-blanket",
        kind=DISAGREEMENT,
        pair=OPH_FEP,
        claim="A blanket visible in the stationary statistics need not bound the dynamics: with zero (mu, eta) precision entry the internal drift is 2q eta + (q - 1) b - 2 mu, so two states with equal blanket and internal values differ in internal motion by 2|q| |eta - eta'|.",
        oph="Seam structure is declared before the law, so compatibility and dynamics are separated by construction.",
        fep="A Markov blanket of the stationary density permits reciprocal internal-external influence.",
        synergetics="No counterpart.",
        status="obstructed",
        source=f"{_PAPER_SOURCE}, Proposition 11",
        probe=_probe_stationary_blanket,
    ),
    Correspondence(
        key="solenoidal-freedom",
        kind=DISAGREEMENT,
        pair=OPH_FEP,
        claim="One stationary density is shared by a family of laws: the diffusion dx = -(Gamma - Q) H x dt + sqrt(2 Gamma) dW has stationary law N(0, H^{-1}) for every antisymmetric Q, and its entropy production rate is tr(Q^T Gamma^{-1} Q H); in two dimensions with Gamma = H = I and Q = q [[0, 1], [-1, 0]] the spectrum is -1 +/- i q and the entropy production is 2 q^2.",
        oph="Configuration equilibrium does not determine inertia; the native counterpart shares the full reference state with distinct finite kinetics.",
        fep="A stationary law does not fix the law that preserves it; the solenoidal part is a free parameter.",
        synergetics="No counterpart.",
        status="obstructed",
        source=f"{_PAPER_SOURCE}, Proposition 9",
        probe=_probe_solenoidal,
    ),
    Correspondence(
        key="d-stability",
        kind=DISAGREEMENT,
        pair=OPH_FEP,
        claim="A linear rest admits a free-energy form exactly when it is stable, and one positive local gain removes that form without moving the rest: at k = 4 the Routh first column is (1, 6, -2, 24) and the scaled matrix has two right-half-plane roots.",
        oph="A spectral certificate is a property of the law's sign structure, which the compatibility set does not fix.",
        fep="Making one part repair faster can destroy convergence while compatibility is unchanged.",
        synergetics="No counterpart.",
        status="obstructed",
        source=f"{_PAPER_SOURCE}, Proposition 10 and Section 4.1",
        probe=_probe_d_stability,
    ),
    Correspondence(
        key="schedule-path",
        kind=DISAGREEMENT,
        pair=OPH_FEP,
        claim="Agreement and dissipation are state functions while finite phases are path functions: full and half-damped iteration of the same affine rule share every answer, yet their centered contrasts after twelve steps differ by exactly 921633999 / 40000000000.",
        oph="What is read from a completed history is a state function; the solver counter is a path function.",
        fep="A learner that reads unfinished settling is taught differently under a different schedule.",
        synergetics="No counterpart.",
        status="obstructed",
        source=f"{_PAPER_SOURCE}, Proposition 6 and Section 4.4",
        probe=_probe_schedule_path,
    ),
    Correspondence(
        key="reduction-criterion",
        kind=DISAGREEMENT,
        pair=OPH_FEP,
        claim="Two histories can share the present observation and need different futures, so a summary that drops the record must err by at least half the response distance on one of them: here the worst-case error floor is exactly one.",
        oph="Records are part of the state; the observer's next action follows its retained record.",
        fep="An observation-only state can fail the response criterion while every solve converges.",
        synergetics="No counterpart.",
        status="obstructed",
        source=f"{_PAPER_SOURCE}, Proposition 4 and Section 4.2",
        probe=_probe_reduction,
    ),
    Correspondence(
        key="existence-without-reachability",
        kind=DISAGREEMENT,
        pair=OPH_FEP,
        claim="An agreement can exist and be unreachable: the clipped two-node map has the origin as its only equilibrium and a source there, and the relaxation orbit stays at radius at least 0.84 with period about 7.58, while the synchronous clipped map, after a transient, stays at radius at least one. The fixed-point set is found by exact case analysis over the map's nine clipping regimes, so the origin is shown to be the only equilibrium rather than sampled.",
        oph="Existence of a public fixed point is not reachability by the repair dynamics.",
        fep="Existence of a fixed point is not convergence of the process theory.",
        synergetics="No counterpart.",
        status="obstructed",
        source=f"{_PAPER_SOURCE}, Proposition 15 and Appendix A.3",
        probe=_probe_existence_without_reachability,
    ),
    Correspondence(
        key="directed-gradient",
        kind=DISAGREEMENT,
        pair=OPH_FEP,
        claim="A directed network's forward response need not carry its exact gradient: for J = [[-1, 2], [0, -1]] and c = (1, 0), the direct nudging response is (1, 0) while the adjoint is (1, 2).",
        oph="The native law is self-adjoint and nonnegative, so this boundary does not bind it.",
        fep="Equilibrium propagation's gradient identity needs an energy; a directed wiring does not automatically supply one.",
        synergetics="No counterpart.",
        status="obstructed",
        source=f"{_PAPER_SOURCE}, Appendix B ('directed interactions')",
        probe=_probe_directed_gradient,
    ),
    Correspondence(
        key="four-cycle-local-minimum",
        kind=DISAGREEMENT,
        pair=OPH_FEP,
        claim="A local minimum of an independently defined discrepancy can be incompatible: on a four-cycle, 0011 violates two edges and every single-bit flip also violates two.",
        oph="An energy-descent argument needs a completeness argument before termination can imply agreement.",
        fep="The same obligation applies to any predictive-coding energy.",
        synergetics="No counterpart.",
        status="obstructed",
        source=f"{_PAPER_SOURCE}, Section 4.1",
        probe=_probe_four_cycle,
    ),
    Correspondence(
        key="cuboctahedron-vs-icosahedron",
        kind=DISAGREEMENT,
        pair=SYNERGETICS_OPH,
        claim="The twelve-around-one shell and the carrier are not the same solid although each has twelve vertices: the cuboctahedron is four-regular with eight triangles and six squares and tetravolume 20, while the carrier is the icosahedron, five-regular with twenty triangles and tetravolume 18.512295868219..., so a count coincidence of twelve is not an identification.",
        oph="The carrier is the icosahedron; its twenty faces are the faces the repair law uses.",
        fep="No counterpart.",
        synergetics="The vector equilibrium is the twelve-around-one cuboctahedron; the icosahedron is a different twelve-vertex solid.",
        status="obstructed",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}",
        probe=_probe_cuboctahedron_vs_icosahedron,
    ),
    Correspondence(
        key="shell-vs-free-energy",
        kind=DISAGREEMENT,
        pair=FEP_SYNERGETICS,
        claim="The IVM shell index is not a function of the free energy: shell two spans three distinct free energies and shell four spans seven, the level 36 occurs on shells three and four, and the level 52 on shells four and five, so neither quantity is a function of the other, although the free energy does bound the shell from above (the minimum on shell k grows like 2k^2). The reason is exact: the shell index is the cuboctahedral gauge max(||u||_inf, ||u||_1 / 2) of the halved coordinates, a polyhedral norm, while the free energy is a Euclidean quadratic form; on shell k the free energy ranges from 2k^2 (2k^2 + 2 for odd k) up to 4k^2, and reaches 4k^2 only at the twelve vertices k m.",
        oph="OPH's load shell is a level set of total load, a linear functional, unrelated to the lattice's graph-distance shells.",
        fep="The free energy of a Gaussian mean is a quadratic form whose level sets cut across the coordination shells.",
        synergetics="The shells are the 10 k^2 + 2 lattice points on the surface of the frequency-k vector equilibrium, at graph distance k from the centre.",
        status="obstructed",
        source=_QUADRAY_SOURCE,
        probe=_probe_shell_vs_free_energy,
    ),
    Correspondence(
        key="quadray-ivm",
        kind=COMPLEMENTARITY,
        pair=SYNERGETICS_OPH,
        claim="Quadray coordinates give the synergetic side exact integer geometry: twelve neighbour moves at squared distance eight, shell counts 10k^2 + 2, and the tetravolume table 1, 3, 4, 6, 20 and 96, every entry derived exactly from an explicit integer vertex set rather than read from the table.",
        oph="The carrier needs a lattice in which twelve neighbours is the kissing configuration.",
        fep="No fixed geometry; the FEP leaves the interface geometry to the model.",
        synergetics="The IVM/quadray framework supplies the lattice, the shell counts and the volumes exactly.",
        status="exact",
        source="QuadMath src/quadmath/core/quadray.py and src/quadmath/lattice/omni_numbering.py",
        probe=_probe_quadray_ivm,
    ),
    Correspondence(
        key="equilibrium-propagation",
        kind=COMPLEMENTARITY,
        pair=OPH_FEP,
        claim="The two-settle contrast of equilibrium propagation, an energy-based learning rule that the FEP literature reaches through predictive coding, supplies the parameter gradient that the native repair law does not define: on a scalar energy model the contrast equals (theta - target) / (1 + beta) exactly and tends to the exact loss gradient as the nudge vanishes.",
        oph="The native law repairs; it does not learn its law.",
        fep="Learning is a descent of free energy in parameters, which needs a gradient the settling rule alone does not give.",
        synergetics="No counterpart.",
        status="exact",
        source=f"{_PAPER_SOURCE}, Theorem 21 (Scellier and Bengio)",
        probe=_probe_equilibrium_propagation,
    ),
    Correspondence(
        key="replicator-potential",
        kind=COMPLEMENTARITY,
        pair=OPH_FEP,
        claim="With a symmetric payoff, replicator selection climbs a potential, so the repair language covers it as a comparison class rather than a shared mechanism: the field equals M(x) grad F with M = diag(x) - x x^T (the Shahshahani gradient) at every point checked and the potential rate is nonnegative; with a non-symmetric payoff the same identity fails at every point checked, so the potential reading is a property of partnership games, not of selection in general.",
        oph="A comparison class, not a realization.",
        fep="A comparison class, not a process theory.",
        synergetics="No counterpart.",
        status="exact",
        source=f"{_PAPER_SOURCE}, Section 6.3",
        probe=_probe_replicator,
    ),
    Correspondence(
        key="kissing-twelve",
        kind=COMPLEMENTARITY,
        pair=SYNERGETICS_OPH,
        claim="Twelve is the kissing number in three dimensions and both twelve-port arrangements attain it, in different ways: with the twelve sphere centres at distance 2r from a central sphere of radius r, the IVM's cuboctahedral shell has edge equal to circumradius, so the twelve spheres touch the centre and each other, while the icosahedral arrangement has squared edge over squared circumradius exactly 2 - 2 sqrt(5) / 5, so its spheres touch the centre and not each other; the geometric icosahedral realisation is isomorphic to OPH's labelled carrier.",
        oph="The twelve-port carrier is the icosahedral kissing arrangement, the one with slack between neighbours.",
        fep="No fixed carrier.",
        synergetics="Twelve-around-one: the vector equilibrium is the rigid lattice arrangement, and the icosahedral one is the loose arrangement the jitterbug passes through.",
        status="exact",
        source="QuadMath shell_count and IVM_NEIGHBOR_MOVES; OPH carrier counts",
        probe=_probe_kissing_twelve,
    ),
    Correspondence(
        key="icosahedron-tetravolume",
        kind=COMPLEMENTARITY,
        pair=SYNERGETICS_OPH,
        claim="Synergetics' tetravolume accounting extends to the carrier's own solid: the regular icosahedron measures 18.512295868219... tetravolumes, equal to (5/2)(3 sqrt 2 + sqrt 10), whose square is exactly 175 + 75 sqrt 5, so the carrier's polytope has an exact algebraic volume that no integral table entry can express.",
        oph="The carrier's polytope: twelve ports, thirty seams, twenty faces.",
        fep="No counterpart; the FEP assigns no volume to an interface.",
        synergetics="The tetravolume of the regular icosahedron, measured against the equal-edge tetrahedron.",
        status="exact",
        source=_QUADRAY_SOURCE,
        probe=_probe_icosahedron_tetravolume,
    ),
    Correspondence(
        key="packing-density",
        kind=COMPLEMENTARITY,
        pair=SYNERGETICS_OPH,
        claim="The IVM's Voronoi cell is the rhombic dodecahedron, fourteen vertices, twelve rhombic faces and twenty-four edges with Euler characteristic two, its tetravolume is exactly six (computed in rational arithmetic), and the ball of radius half the prime vector has tetravolume exactly pi sqrt 2, so the packing density is pi sqrt 2 / 6 = pi / (3 sqrt 2): the face-centred cubic density, optimal among lattice packings and among all packings in three dimensions.",
        oph="The densest packing's kissing arrangement is the cuboctahedral one, not the carrier's icosahedral one; the carrier inherits the kissing number, not the packing.",
        fep="No counterpart.",
        synergetics="The rhombic dodecahedron's six tetravolumes, derived here from the cell rather than read from the table.",
        status="exact",
        source=_QUADRAY_SOURCE,
        probe=_probe_packing_density,
    ),
    Correspondence(
        key="tetravolume-arithmetic",
        kind=COMPLEMENTARITY,
        pair=THREE_WAY,
        claim="The three programmes' exact quantities live in different number systems, and no single one serves all three: the FEP's equal-covariance divergence for one IVM step is the integer four nats (in the integer embedding's unit of length), the carrier's normalised spectral gap, a relaxation rate of OPH's mean repair operator rather than a cost of a lattice step, is the quadratic surd (5 - sqrt 5) / 60, and the synergetic tetravolume table is integral (1, 3, 4, 6, 20, 96), every entry derived exactly, while the carrier's own icosahedron enters at sqrt 2 (15/2 + (5/2) sqrt 5) = 18.512295868219..., whose square is 175 + 75 sqrt 5. The three are attached to different objects; what the entry records is that no common unit converts one into another.",
        oph="Integer loads; spectral data in the quadratic field Q(sqrt 5).",
        fep="Real-valued information costs in nats.",
        synergetics="Tetravolumes of the prime vector: integral on the IVM, algebraic irrational on the icosahedron.",
        status="exact",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}; {_PAPER_SOURCE}",
        probe=_probe_tetravolume_arithmetic,
    ),
    Correspondence(
        key="effective-resistance",
        kind=COMPLEMENTARITY,
        pair=OPH_FEP,
        claim="The precision carries a metric neither programme computes: the effective resistance between two ports, exact and rational, is 11/30 across a seam, 7/15 at distance two and 1/2 across the carrier, and the 66 of them sum to 28, which is twelve times the reciprocal spectrum sum 7/3 exactly.",
        oph="A Laplacian whose generalized inverse defines a resistance distance.",
        fep="In the Gaussian reading the same number is the variance of a port's potential difference.",
        synergetics="A lattice metric on the same twelve sites, counted in whole steps rather than in resistances.",
        status="exact",
        source=f"{_CARRIER_SOURCE}; {_PAPER_SOURCE}",
        probe=_probe_effective_resistance,
    ),
    Correspondence(
        key="vertex-shells",
        kind=DISAGREEMENT,
        pair=SYNERGETICS_OPH,
        claim="The shells of ports around a port are 5, 5 and 1, not the lattice's 12, 42 and 92: the carrier is distance-regular with intersection array {5, 2, 1; 1, 2, 5} and diameter three, so reading its neighbourhoods as IVM shells fails at the first shell, 5 against 12.",
        oph="Neighbourhood shells of the port adjacency.",
        fep="The same shells as conditional-independence distance in the precision.",
        synergetics="IVM shells 10k^2 + 2 around a lattice site.",
        status="obstructed",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}",
        probe=_probe_vertex_shells,
    ),
    Correspondence(
        key="subgraph-census",
        kind=ISOMORPHISM,
        pair=SYNERGETICS_OPH,
        claim="The carrier's subgraph census, computed from the derived graph alone, is asserted equal to the icosahedron's own counts: 20 triangles, exactly its twenty faces; 12 induced five-cycles, the links of its twelve vertices and hence the twelve pentagonal faces of its dual dodecahedron; 125 perfect matchings; and 1280 Hamiltonian cycles, counted as undirected edge sets.",
        oph="The 20 triangles of the derived face set.",
        fep="A census of the precision's underlying graph; the FEP reading inherits the counts and asserts nothing about them.",
        synergetics="The icosahedron's twenty faces and the dodecahedron's twelve pentagonal faces.",
        status="exact",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}",
        probe=_probe_subgraph_census,
    ),
    Correspondence(
        key="angle-defect",
        kind=ISOMORPHISM,
        pair=SYNERGETICS_OPH,
        claim="OPH's derived chart is the icosahedron's triangulation, so Synergetics' angular-defect accounting applies unchanged: five equilateral triangles meet at every port, the defect is 1/3 pi at each of the twelve ports and 4 pi in total, which is 2 pi times the Euler characteristic 2 of the sphere both programmes triangulate.",
        oph="The derived patch complex: five faces at every port, the very vertex figure the defect accounting measures.",
        fep="Curvature of the shared complex; the FEP vocabulary adds no reading of a geometric defect, and none is claimed.",
        synergetics="The icosahedron's five-fold vertex figures, whose angular accounting is the classical one.",
        status="exact",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}",
        probe=_probe_angle_defect,
    ),
    Correspondence(
        key="twelve-around-one-families",
        kind=DISAGREEMENT,
        pair=SYNERGETICS_OPH,
        claim="`Twelve around one` names a family and not an object: the icosahedron has degree 5 and 5184000 spanning trees, while the cuboctahedron has degree 4, eight triangles, 32 perfect matchings and 331776 = 576^2 spanning trees, so the counting invariant separates the two solids the phrase covers.",
        oph="The carrier is the icosahedron, so its count is 5184000.",
        fep="Nothing: both graphs support the same Gaussian reading, and the count is what distinguishes them.",
        synergetics="The vector equilibrium is the cuboctahedron, and its count is 576 squared.",
        status="obstructed",
        source=f"{_QUADRAY_SOURCE}; {_CARRIER_SOURCE}",
        probe=_probe_twelve_around_one_families,
    ),
    Correspondence(
        key="dual-cells",
        kind=COMPLEMENTARITY,
        pair=SYNERGETICS_OPH,
        claim="The twelve-around-one cell and the cuboctahedron are duals edge for edge: 14 vertices against 12, 12 faces against 14, 24 edges each, Euler characteristic 2 each, and tetravolumes 6 and 20 in the ratio 10/3, with the cuboctahedron derived a second time as the cube's line graph.",
        oph="No counterpart: the cell is the IVM's, built from the cuboctahedral neighbours; the carrier is the icosahedron, whose twelve sites would bound a different cell.",
        fep="The cell's volume is not a free-energy quantity; the dual pair is geometry.",
        synergetics="The rhombic dodecahedron's six and the cuboctahedron's twenty tetravolumes.",
        status="exact",
        source=f"{_QUADRAY_SOURCE}; {_CARRIER_SOURCE}",
        probe=_probe_dual_cells,
    ),
    Correspondence(
        key="conjugation-norm",
        kind=COMPLEMENTARITY,
        pair=THREE_WAY,
        claim="Conjugation of sqrt 5 (fixing sqrt 2) fixes the counting invariants and moves the geometric ones: the spanning-tree number 5184000 is a rational integer and is fixed, while the icosahedron's tetravolume sqrt 2 (15/2 + (5/2) sqrt 5) moves to sqrt 2 (15/2 - (5/2) sqrt 5), and the two multiply to exactly 50, computed in exact field arithmetic; consistently, (175 + 75 sqrt 5)(175 - 75 sqrt 5) = 2500 = 50^2.",
        oph="A count over the carrier, which is an integer and therefore conjugation-invariant.",
        fep="A determinant over the same matrix, likewise rational.",
        synergetics="Tetravolumes in sqrt 2 Q(sqrt 5), whose squares lie in Q(sqrt 5); the volume moves under the conjugation, which fixes sqrt 2.",
        status="exact",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}; {_PAPER_SOURCE}",
        probe=_probe_conjugation_norm,
    ),
)


def by_key() -> dict[str, Correspondence]:
    """Return the registered correspondences keyed by their identifier."""
    return {entry.key: entry for entry in CORRESPONDENCES}


def by_kind() -> dict[str, tuple[Correspondence, ...]]:
    """Return the registered correspondences grouped by kind, in registry order."""
    grouped: dict[str, list[Correspondence]] = {}
    for entry in CORRESPONDENCES:
        grouped.setdefault(entry.kind, []).append(entry)
    return {kind: tuple(entries) for kind, entries in grouped.items()}


def by_pair() -> dict[str, tuple[Correspondence, ...]]:
    """Return the registered correspondences grouped by pairwise scope.

    Every declared pair has an entry in the result even when it is empty, so a
    triangle drawn from it cannot silently lose an edge.
    """
    grouped: dict[str, list[Correspondence]] = {pair: [] for pair in PAIRS}
    for entry in CORRESPONDENCES:
        grouped.setdefault(entry.pair, []).append(entry)
    return {pair: tuple(entries) for pair, entries in grouped.items()}


def verify(
    *,
    keys: tuple[str, ...] | None = None,
    pairs: tuple[str, ...] | None = None,
    kinds: tuple[str, ...] | None = None,
) -> tuple[CorrespondenceResult, ...]:
    """Run every registered probe, or only the selected ones.

    An unknown key, pairwise scope or kind raises ``ValueError`` instead of
    selecting nothing, so a mistyped argument cannot pass as an empty run.
    """
    selected = CORRESPONDENCES
    if keys is not None:
        wanted = set(keys)
        unknown = sorted(wanted - {entry.key for entry in CORRESPONDENCES})
        if unknown:
            raise ValueError(f"unknown correspondence key(s): {', '.join(unknown)}")
        selected = tuple(entry for entry in selected if entry.key in wanted)
    if pairs is not None:
        wanted_pairs = set(pairs)
        unknown_pairs = sorted(wanted_pairs - set(PAIRS))
        if unknown_pairs:
            raise ValueError(f"unknown pairwise scope(s): {', '.join(unknown_pairs)}")
        selected = tuple(entry for entry in selected if entry.pair in wanted_pairs)
    if kinds is not None:
        wanted_kinds = set(kinds)
        unknown_kinds = sorted(wanted_kinds - set(KINDS))
        if unknown_kinds:
            raise ValueError(f"unknown kind(s): {', '.join(unknown_kinds)}")
        selected = tuple(entry for entry in selected if entry.kind in wanted_kinds)
    return tuple(entry.run() for entry in selected)


def summary(results: tuple[CorrespondenceResult, ...]) -> dict[str, object]:
    """Aggregate results by kind, by pairwise scope, and by pass/fail."""
    by_kind: dict[str, dict[str, int]] = {}
    by_pair_count: dict[str, dict[str, int]] = {}
    failures: list[str] = []
    for result in results:
        for name, table in ((result.kind, by_kind), (result.pair, by_pair_count)):
            bucket = table.setdefault(name, {"passed": 0, "failed": 0})
            bucket["passed" if result.passed else "failed"] += 1
        if not result.passed:
            failures.append(result.key)
    return {
        "total": len(results),
        "passed": sum(1 for result in results if result.passed),
        "failed": len(failures),
        "failures": tuple(failures),
        "by_kind": by_kind,
        "by_pair": by_pair_count,
        "pairs": PAIRS,
    }
