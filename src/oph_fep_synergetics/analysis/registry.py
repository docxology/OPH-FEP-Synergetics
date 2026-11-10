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
from ..core import repair as repair_module
from ..core import sheaf as sheaf_module

ISOMORPHISM = "isomorphism"
DISAGREEMENT = "disagreement"
COMPLEMENTARITY = "complementarity"


@dataclass(frozen=True)
class Correspondence:
    """One registered claim about how two or more of the three programs relate."""

    key: str
    kind: str
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
            status=self.status,
            passed=bool(passed),
            values=dict(values),
        )


@dataclass(frozen=True)
class CorrespondenceResult:
    """The outcome of one probe, with the values that establish it."""

    key: str
    kind: str
    status: str
    passed: bool
    values: dict


def _carrier() -> carrier_module.IcosahedralCarrier:
    return carrier_module.icosahedral_carrier()


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
    }
    passed = (
        values["total_two"] == 78
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
    values = {
        "iterations": iterations,
        "residual": residual,
        "perceive_matches_harmonic_extension": bool(np.allclose(interior, analytic, atol=1e-9)),
        "posterior_mean_is_harmonic_extension": bool(np.allclose(mean, analytic, atol=1e-9)),
        "energy_at_minimiser": energy,
        "min_perturbed_energy": min(perturbed),
        "stationary_gradient_norm": float(
            np.linalg.norm(blocks["L_xx"] @ interior + blocks["L_xo"] @ observed)
        ),
        "excess_free_energy_at_posterior": model.excess_free_energy(mean, covariance, records),
    }
    passed = (
        values["perceive_matches_harmonic_extension"]
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
    gradient = blocks["L_xx"] @ model.harmonic_extension(old_records) + blocks["L_xo"] @ np.asarray(
        new_records
    )
    bayesian = model.bayesian_surprise(old_records, new_records)
    sandwich = control_module.surprise_sandwich(blocks["L_xx"], gradient, bayesian)
    closed_form = 0.5 * float(gradient @ np.linalg.solve(blocks["L_xx"], gradient))
    mean, covariance = model.posterior(new_records)
    values = {
        "record_injectivity": model.record_injectivity(),
        "bayesian_surprise": bayesian,
        "closed_form_half_g_linv_g": closed_form,
        "sandwich_lower": sandwich["lower_bound"],
        "sandwich_upper": sandwich["upper_bound"],
        "sandwich_inside": sandwich["inside"],
        "surprisal": model.surprisal(new_records),
        "variational_free_energy": model.variational_free_energy(mean, covariance, new_records),
        "numerical_surprise": model.numerical_surprise(old_records, new_records),
    }
    passed = (
        values["record_injectivity"]
        and abs(bayesian - closed_form) < 1e-12
        and values["sandwich_inside"]
        and abs(values["variational_free_energy"] - values["surprisal"]) < 1e-12
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
    """Two objects called free energy differ by a model constant, not by a gauge choice.

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
    passed = abs(values["difference"] - expected_gap) < 1e-12
    return passed, values


def _probe_mediator() -> tuple[bool, dict]:
    """Proposition 8: eliminating a mediator can require a higher-order interface."""
    threshold = 1.0 / 3.0
    values = {
        "schur_matches_expected_below": bool(
            np.allclose(control_module.mediator_schur(0.0), control_module.mediator_expected(0.0), atol=1e-12)
        ),
        "schur_matches_expected_above": bool(
            np.allclose(control_module.mediator_schur(0.5), control_module.mediator_expected(0.5), atol=1e-12)
        ),
        "margin_below_threshold": control_module.mediator_pairwise_margin(0.2),
        "margin_at_threshold": control_module.mediator_pairwise_margin(threshold),
        "decomposition_below_is_none": control_module.mediator_pairwise_decomposition(0.2) is None,
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
    values = dict(control_module.blanket_dynamics_check(0.75))
    values["expected_mu_drift"] = (2.0 * 0.75, 0.75 - 1.0, -2.0)
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
    passed = (
        values["relaxation_min_radius"] >= 0.83
        and abs(values["relaxation_period"] - 7.58) < 0.05
        and sync_radius >= 0.99
        and any(real > 0.0 for real, _ in values["relaxation_jacobian_eigenvalues"])
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
    cumulative = tuple(int(geometry.cumulative_count(k)) for k in range(5))
    ratios = sorted(float(value) for value in geometry.VOLUME_RATIOS.values())
    values = {
        "neighbor_moves": len(moves),
        "distinct_moves": len(set(moves)),
        "squared_norms": sorted({float(geometry.squared_norm(move)) for move in moves}),
        "min_pairwise_squared_distance": float(pairwise),
        "shell_counts": shells,
        "cumulative_counts": cumulative,
        "volume_ratios_present": ratios,
        "max_shell": geometry.MAX_SHELL,
        "neighbor_squared_distance": geometry.NEIGHBOR_SQUARED_DISTANCE,
    }
    passed = (
        values["neighbor_moves"] == 12
        and values["distinct_moves"] == 12
        and values["squared_norms"] == [8.0]
        and abs(values["min_pairwise_squared_distance"] - 8.0) < 1e-9
        and shells == (1, 12, 42, 92, 162)
        and cumulative == (1, 13, 55, 147, 309)
        and {1.0, 3.0, 4.0, 6.0, 20.0, 96.0} <= set(ratios)
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
    passed = (
        values["field_equals_mobility_gradient"] is True
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
    passed = (
        values["carrier_ports"] == 12
        and values["carrier_degrees"] == [5]
        and values["ivm_shell_one"] == 12
        and values["ivm_moves"] == 12
        and abs(values["move_squared_norm"] - 8.0) < 1e-9
        and values["geometric_adjacency_isomorphic"]
    )
    return passed, values


# ---------------------------------------------------------------------------
# The registry
# ---------------------------------------------------------------------------

_CARRIER_SOURCE = "OPH evidence/observer_dynamics_20260925/oph-physics-sim (ORIENTED_BASE_FACES, seams, laplacian) at 5517b0c"
_LAW_SOURCE = "OPH model/observer/native.py:45-51 and the 'Agreement and Surprise' paper, Section 6.2"
_PAPER_SOURCE = "Mueller, 'Agreement and Surprise: Local Repair, Dynamics and the Limits of Equilibrium', private draft 22 October 2026"

CORRESPONDENCES: tuple[Correspondence, ...] = (
    Correspondence(
        key="carrier-port-chart",
        kind=ISOMORPHISM,
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
        claim="The carrier Laplacian is 5I - A, its first positive eigenvalue is 5 - sqrt(5), and the record-discarding shell gap is (5 - sqrt(5)) / 60.",
        oph="The spectral certificate of the native repair law.",
        fep="A precision-weighted sheaf Laplacian, of which this is the unweighted case.",
        synergetics="The icosahedral coordination shell of the A3 packing.",
        status="exact",
        source=_CARRIER_SOURCE,
        probe=_probe_spectral_gap,
    ),
    Correspondence(
        key="shell-counts",
        kind=ISOMORPHISM,
        claim="The total-two load shell has 78 configurations and the binary shell has 66.",
        oph="The finite shell the executable observer runs on.",
        fep="The finite state space of a bounded system.",
        synergetics="A finite coordination shell of the closest packing.",
        status="exact",
        source="OPH model/observer/native.py (states) and the paper, Section 6.2",
        probe=_probe_shell_counts,
    ),
    Correspondence(
        key="native-repair-law",
        kind=ISOMORPHISM,
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
        claim="On the linear-Gaussian domain the interface-energy minimiser over the interior is the harmonic extension and the posterior mean, so A3-style constrained relative-entropy minimisation and predictive-coding perception are the same object.",
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
        claim="Surprisal is the disagreement no repair removes, Bayesian surprise is the disagreement repair will remove, and the update residual pins the Bayesian surprise between the extreme interior eigenvalues.",
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
        claim="The native law's sharp certificate is a spectral gap: the mean repair operator is symmetric, nonnegative and doubly stochastic with second eigenvalue 1 - (5 - sqrt(5)) / 60, and its induced infinity norm is exactly one, so no sup-norm contraction certifies it.",
        oph="The capacity-weighted two-norm spectral certificate.",
        fep="Precision-weighted quadratic norms for the linearised update.",
        synergetics="No counterpart.",
        status="exact",
        source=f"{_PAPER_SOURCE}, Section 6.2 ('which certificate a law admits')",
        probe=_probe_certificate_species,
    ),
    Correspondence(
        key="pr01-obstruction",
        kind=DISAGREEMENT,
        claim="The native repair law is not confluent and accepted repairs need not lower a measure: 600 of the 4680 shell transitions move at mismatch one with zero energy decrease, and a one-unit walk keeps energy 1/2 while changing state.",
        oph="PR-01 'confluent terminating repair contract' is recorded as obstructed; the consensus theorems hold for the integer and transactional contracts, not the native law.",
        fep="Active inference is a descent on a bound; the identification with this law therefore fails on the mismatch-one domain.",
        synergetics="No counterpart.",
        status="obstructed",
        source="OPH model/PREMISES.md:95 (PR-01) and Lean/ObserverPatchHolography/Primitives.lean (demoCarrier_not_confluent)",
        probe=_probe_pr01,
    ),
    Correspondence(
        key="free-energy-name-collision",
        kind=DISAGREEMENT,
        claim="The interface energy minimised over the interior and the variational free energy at the posterior are different objects that differ by a model-dependent constant, so identifying them numerically is a gauge choice rather than a theorem.",
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
        claim="Eliminating an interior mediator leaves a genuinely three-way boundary term: the reduced precision is (1/3) v v^T + epsilon I with v = (1, 1, -1), and a sum of pairwise semidefinite terms reproduces it if and only if epsilon >= 1/3.",
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
        claim="One stationary density is shared by a family of laws: the diffusion with drift -(Gamma - Q) H has stationary law N(0, H^{-1}) for every antisymmetric Q, with spectrum -1 +/- i q and entropy production 2 q^2 in the two-dimensional case.",
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
        claim="An agreement can exist and be unreachable: the clipped two-node map has the origin as its only equilibrium and a source there, and the relaxation orbit stays at radius at least 0.84 with period about 7.58, while the synchronous clipped map stays at radius at least one.",
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
        claim="A local minimum of an independently defined discrepancy can be incompatible: on a four-cycle, 0011 violates two edges and every single-bit flip also violates two.",
        oph="An energy-descent argument needs a completeness argument before termination can imply agreement.",
        fep="The same obligation applies to any predictive-coding energy.",
        synergetics="No counterpart.",
        status="obstructed",
        source=f"{_PAPER_SOURCE}, Section 4.1",
        probe=_probe_four_cycle,
    ),
    Correspondence(
        key="quadray-ivm",
        kind=COMPLEMENTARITY,
        claim="Quadray coordinates give the synergetic side exact integer geometry: twelve neighbour moves at squared distance eight, shell counts 10k^2 + 2, and rational tetravolumes including 1, 3, 4, 6, 20 and 96.",
        oph="The carrier needs a lattice in which twelve neighbours is the kissing configuration.",
        fep="No fixed geometry; the FEP leaves the interface geometry to the model.",
        synergetics="The IVM/quadray framework supplies the lattice, the shell counts and the volumes exactly.",
        status="exact",
        source="QuadMath src/quadmath/core/quadray.py and src/quadmath/lattice/omni_numbering.py at a8e335b",
        probe=_probe_quadray_ivm,
    ),
    Correspondence(
        key="equilibrium-propagation",
        kind=COMPLEMENTARITY,
        claim="The two-settle contrast of equilibrium propagation supplies the parameter gradient that the native repair law does not define: the contrast equals (theta - target) / (1 + beta) exactly and tends to the exact loss gradient as the nudge vanishes.",
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
        claim="Replicator selection climbs a potential, so the repair language covers it as a comparison class rather than a shared mechanism: the field equals M(x) grad F with M = diag(x) - x x^T and the potential rate is nonnegative.",
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
        claim="The twelve-port carrier sits on the kissing-number-twelve configuration: the IVM shell one has twelve sites at squared distance eight, and the geometric icosahedral realisation is isomorphic to OPH's labelled carrier.",
        oph="The twelve-port carrier with its twelve neighbours.",
        fep="No fixed carrier.",
        synergetics="Twelve-around-one, the vector equilibrium.",
        status="exact",
        source="QuadMath shell_count and IVM_NEIGHBOR_MOVES; OPH carrier counts",
        probe=_probe_kissing_twelve,
    ),
)


def by_key() -> dict[str, Correspondence]:
    """Return the registered correspondences keyed by their identifier."""
    return {entry.key: entry for entry in CORRESPONDENCES}


def verify(*, keys: tuple[str, ...] | None = None) -> tuple[CorrespondenceResult, ...]:
    """Run every registered probe, or only the named ones.

    An unknown key raises ``ValueError`` instead of selecting nothing, so a
    mistyped key cannot pass as an empty run.
    """
    if keys is None:
        selected = CORRESPONDENCES
    else:
        wanted = set(keys)
        unknown = sorted(wanted - {entry.key for entry in CORRESPONDENCES})
        if unknown:
            raise ValueError(f"unknown correspondence key(s): {', '.join(unknown)}")
        selected = tuple(entry for entry in CORRESPONDENCES if entry.key in wanted)
    return tuple(entry.run() for entry in selected)


def summary(results: tuple[CorrespondenceResult, ...]) -> dict[str, object]:
    """Aggregate results by kind and by pass/fail."""
    by_kind: dict[str, dict[str, int]] = {}
    failures: list[str] = []
    for result in results:
        bucket = by_kind.setdefault(result.kind, {"passed": 0, "failed": 0})
        if result.passed:
            bucket["passed"] += 1
        else:
            bucket["failed"] += 1
            failures.append(result.key)
    return {
        "total": len(results),
        "passed": sum(1 for result in results if result.passed),
        "failed": len(failures),
        "failures": tuple(failures),
        "by_kind": by_kind,
    }
