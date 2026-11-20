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
from itertools import combinations, product
from typing import Callable

import numpy as np

from .. import geometry
from ..core import carrier as carrier_module
from ..core import centroid_dynamics
from ..core import controls as control_module
from ..core import free_energy as free_energy_module
from ..core import geometric_screening
from ..core import invariants as invariant_module
from ..core import native_dynamics
from ..core import network as network_module
from ..core import planar_duality
from ..core import quadratic_modes
from ..core import repair as repair_module
from ..core import sheaf as sheaf_module
from ..core import symmetry as symmetry_module
from ..core import weighted_duality
from ..geometry.lattice import ivm_quadray_cosets

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
                if step.energy_decrease != predicted:
                    return False, {"failure": "energy decrease", "seam": (a, b)}
                transitions += 1
    return True, {"transitions": transitions, "states": 78, "closed_form": "(d^2 - (|d| mod 2)) / 4"}


_RECORD_PRIOR_PRECISION = 0.5


def _carrier_gaussian_model() -> free_energy_module.GaussianSheafModel:
    return free_energy_module.GaussianSheafModel(
        sheaf_module.carrier_sheaf(_carrier()), (0, 3),
        record_prior_precision=_RECORD_PRIOR_PRECISION,
    )


def _chain_gaussian_model() -> free_energy_module.GaussianSheafModel:
    return free_energy_module.demo_chain_model(4, record_prior_precision=_RECORD_PRIOR_PRECISION)


def _probe_variational_equivalence() -> tuple[bool, dict]:
    passed, values = _check_variational_equivalence(_carrier_gaussian_model())
    chain_passed, chain = _check_variational_equivalence(_chain_gaussian_model())
    # Integrate the non-Gaussian posterior of a single oriented seam, clamped
    # at zero, using midpoints. Each nonzero displacement has positive energy,
    # so the unique minimiser is zero. The half-Gaussian integrals give the
    # mean by a separate route.
    dx = 0.001
    midpoints = np.arange(-12.0 + dx / 2, 12.0, dx)
    weights = np.exp(-free_energy_module.asymmetric_seam_energy(midpoints))
    mean = float((midpoints @ weights) / weights.sum())
    expected_mean = (0.5 - 1.0) / (math.sqrt(math.pi / 2) * (1.0 + 1.0 / math.sqrt(2)))
    candidates = np.append(midpoints, 0.0)
    minimiser = float(candidates[np.argmin(free_energy_module.asymmetric_seam_energy(candidates))])
    diverges = abs(mean - expected_mean) < 1e-6 and abs(mean - minimiser) > 0.1
    values.update({
        "chain_control_passes": bool(chain_passed),
        "chain_control_interior_dimension": chain["interior_dimension"],
        "nonquadratic_energy_minimiser": minimiser,
        "nonquadratic_posterior_mean": mean,
        "nonquadratic_expected_posterior_mean": expected_mean,
        "nonquadratic_mean_differs_from_minimiser": diverges,
    })
    return passed and chain_passed and diverges, values


def _check_variational_equivalence(
    model: free_energy_module.GaussianSheafModel,
) -> tuple[bool, dict]:
    """The interface-energy minimiser is the constrained relative-entropy minimiser.

    On the linear-Gaussian domain the two descriptions are the same object: the
    interior minimiser of ``E(x, o) = (1/2) z^T L z`` at fixed records is the
    harmonic extension, and it is also the mean of the Gaussian posterior that
    minimises the relative entropy to the reference among states matching the
    records.
    """
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
        "record_ports": model.record_nodes,
        "interior_dimension": model.interior_dim,
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
        "excess_free_energy_at_posterior_residual": model.excess_free_energy(mean, covariance, records),
    }
    constrained = model.constrained_minimiser(records)
    values["kkt_minimiser_matches_harmonic_extension"] = bool(
        np.allclose(constrained, analytic, atol=1e-12)
    )
    values["kkt_harmonic_extension_residual"] = float(np.linalg.norm(constrained - analytic))
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
        and abs(values["excess_free_energy_at_posterior_residual"]) < 1e-12
    )
    return passed, values


def _probe_free_energy_dictionary() -> tuple[bool, dict]:
    """Theorem 7: harmonic extension, surprisal, Bayesian surprise and the sandwich."""
    passed, values = _check_free_energy_dictionary(_carrier_gaussian_model())
    chain_passed, _ = _check_free_energy_dictionary(_chain_gaussian_model())
    values["chain_control_passes"] = bool(chain_passed)
    # Ground one endpoint of the asymmetric seam and read the other. A
    # symmetric record prior and the log normaliser cancel in the difference
    # between records of opposite sign and equal magnitude.
    records = np.array([-1.0, 1.0])
    energies = free_energy_module.asymmetric_seam_energy(records)
    gap = float(energies[1] - energies[0])
    values["nonquadratic_record_surprisal_difference"] = gap
    return passed and chain_passed and gap == 0.5, values


def _check_free_energy_dictionary(
    model: free_energy_module.GaussianSheafModel,
) -> tuple[bool, dict]:
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
    old_mean, old_covariance = model.posterior(old_records)
    posterior_kl = free_energy_module.gaussian_kl(old_mean, old_covariance, mean, covariance)
    delta = np.asarray(new_records, dtype=float) - np.asarray(old_records, dtype=float)
    precision = model.schur_complement() + _RECORD_PRIOR_PRECISION * np.eye(model.record_dim)
    predicted_drop = 0.5 * float(delta @ precision @ delta) + float(old_vector @ precision @ delta)
    measured_drop = float(model.record_surprisal(new_records) - model.record_surprisal(old_records))
    values = {
        "record_ports": model.record_nodes,
        "interior_dimension": model.interior_dim,
        "energy_drop_residual": energy_drop_residual,
        "record_injectivity": model.record_injectivity(),
        "bayesian_surprise": bayesian,
        "closed_form_half_g_linv_g": closed_form,
        "kl_between_posteriors": posterior_kl,
        "sandwich_lower": sandwich["lower_bound"],
        "sandwich_upper": sandwich["upper_bound"],
        "sandwich_inside": sandwich["inside"],
        "surprisal": model.record_surprisal(new_records),
        "variational_free_energy": model.variational_free_energy(mean, covariance, new_records),
        "numerical_surprise": model.numerical_surprise(old_records, new_records),
        "surprisal_drop_residual": abs(measured_drop - predicted_drop),
    }
    # The same identities with nonzero old records, where the cross term of the
    # surprisal drop and the old rest do not vanish.
    shifted_old = (0.7, -0.4)
    shifted_vector = np.asarray(shifted_old, dtype=float)
    shifted_delta = new_vector - shifted_vector
    shifted_predicted = 0.5 * float(shifted_delta @ precision @ shifted_delta) + float(
        shifted_vector @ precision @ shifted_delta
    )
    shifted_measured = model.record_surprisal(new_records) - model.record_surprisal(shifted_old)
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
        and abs(bayesian - posterior_kl) < 1e-12
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
        "injectivity_all_ports": all(
            sheaf.record_injectivity((port,)) for port in range(carrier.port_count)
        ),
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
    # The matrix above is built from the Laplacian.  The native law itself,
    # averaged over every seam and both coins, must give the same operator:
    # x - L x / 60 in exact rationals on each of the 78 configurations.
    laplacian = [[int(round(entry)) for entry in row] for row in carrier.laplacian()]
    states = carrier.shell_configurations(total=2)
    values["mean_native_step_configurations"] = len(states)
    values["mean_native_step_is_generator_step"] = all(
        repair_module.mean_native_step(state, carrier.edges)
        == tuple(
            state[port]
            - Fraction(sum(laplacian[port][other] * state[other] for other in range(12)), 60)
            for port in range(12)
        )
        for state in states
    )
    passed = (
        values["mean_native_step_is_generator_step"]
        and values["mean_native_step_configurations"] == 78
        and values["row_sums"] == [1.0]
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
    # Pair the conjugate eigenvalues 5 ± sqrt(5) before forming the product.
    from_closed_spectrum = Fraction((5**2 - 5)**3 * 6**5, carrier.port_count)
    grounded = carrier.laplacian()[1:, 1:]
    cholesky = np.linalg.cholesky(grounded)
    log_prefactor = grounded.shape[0] * math.log(2 * math.pi) / 2
    log_integral = log_prefactor - float(np.log(np.diag(cholesky)).sum())
    log_integral_from_count = log_prefactor - math.log(exact) / 2
    values = {
        "spanning_trees_exact": exact,
        "closed_form_spectral_product_matches_count": from_closed_spectrum == exact,
        "spanning_trees_from_spectrum": from_spectrum,
        "relative_residual": abs(from_spectrum - exact) / exact,
        "grounded_gaussian_log_integrals_agree": math.isclose(
            log_integral, log_integral_from_count, rel_tol=0.0, abs_tol=1e-12
        ),
        "laplacian_is_integral": bool(
            np.array_equal(carrier.laplacian(), np.rint(carrier.laplacian()))
        ),
        "nonzero_eigenvalue_count": int((carrier.laplacian_spectrum() > 1e-9).sum()),
        "seams": carrier.seam_count,
    }
    passed = (
        exact == 5184000
        and values["closed_form_spectral_product_matches_count"]
        and values["relative_residual"] < 1e-9
        and values["grounded_gaussian_log_integrals_agree"]
        and values["nonzero_eigenvalue_count"] == 11
        and values["laplacian_is_integral"] is True
    )
    return passed, values


def _probe_gmrf_independence() -> tuple[bool, dict]:
    """Check the proper Gaussian's graph and the centered intrinsic variance."""
    carrier = _carrier()
    laplacian = carrier.laplacian()
    proper = invariant_module.proper_precision(carrier)
    shift = 1.0 / carrier.port_count
    variance = invariant_module.marginal_variance(carrier)
    pairwise = invariant_module.pairwise_precision(carrier)
    inverse_diagonal = float(np.linalg.pinv(laplacian)[0][0])
    covariance = np.linalg.inv(proper)
    ports = tuple(range(carrier.port_count))
    conditioned_edges = set()
    for pair in combinations(ports, 2):
        others = tuple(port for port in ports if port not in pair)
        cross = covariance[np.ix_(pair, others)]
        conditioned = covariance[np.ix_(pair, pair)] - cross @ np.linalg.solve(
            covariance[np.ix_(others, others)], cross.T
        )
        if abs(conditioned[0, 1]) > 1e-12:
            conditioned_edges.add(pair)
    values = {
        "laplacian_is_precision": bool(
            np.allclose(invariant_module.gmrf_precision(carrier), laplacian)
        ),
        "proper_precision_shift": shift,
        "proper_precision_shift_is_diagonal": bool(
            np.allclose(proper - laplacian, shift * np.eye(carrier.port_count), rtol=0.0, atol=1e-12)
        ),
        "proper_precision_is_positive_definite": bool(np.linalg.eigvalsh(proper)[0] > 0.0),
        "independence_matches_graph": invariant_module.independence_matches_graph(carrier),
        "conditional_covariance_matches_seams": conditioned_edges == set(carrier.edges),
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
        and values["proper_precision_shift_is_diagonal"]
        and values["proper_precision_is_positive_definite"]
        and values["independence_matches_graph"]
        and values["conditional_covariance_matches_seams"]
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
    """The ungated law is confluent but nonterminating; the gate changes both."""
    census = control_module.gated_descent_census()
    walk = control_module.mismatch_one_excursion(steps=6)
    loads = [0] * 12
    loads[0] = 1
    waits, swaps = control_module.mismatch_one_moves(tuple(loads), 0, 1)
    values = dict(census)
    values["total_transitions"] = census["strict"] + census["zero"] + census["noop"]
    values["excursion_states"] = len(walk)
    values["excursion_distinct_states"] = len(set(walk))
    values["excursion_energy"] = float(repair_module.energy(walk[-1]))
    values["coin_zero_current"] = waits.current
    values["coin_one_current"] = swaps.current
    values["coin_zero_decrease"] = float(waits.energy_decrease)
    values["coin_one_decrease"] = float(swaps.energy_decrease)
    # The one-attempt law out of (1, 1, 0, ..., 0): a uniform seam and a fair coin,
    # sixty equally likely draws.  Route one enumerates the draws; route two counts
    # the seams that meet exactly one of the two adjacent occupied ports.
    carrier = _carrier()
    pair = [0] * 12
    pair[0] = pair[1] = 1
    pair_state = tuple(pair)
    kernel = control_module.native_transition_probabilities(pair_state)
    hops = {target: mass for target, mass in kernel.items() if target != pair_state}
    draws = 2 * carrier.seam_count
    one_unit_seams = (carrier.degree(0) - 1) + (carrier.degree(1) - 1)
    values["pair_ports_share_a_seam"] = 1 in carrier.neighbours(0)
    values["stay_probability"] = str(kernel[pair_state])
    values["hop_probability"] = str(max(hops.values()))
    values["hop_targets"] = len(hops)
    values["stay_probability_closed_form"] = str(1 - Fraction(one_unit_seams, draws))
    values["hops_share_one_probability"] = len(set(hops.values())) == 1
    values["pair_draws_with_energy_decrease"] = sum(
        1
        for a, b in carrier.edges
        for coin in (0, 1)
        if repair_module.native_repair(pair_state, (a, b), coin).energy_decrease != 0
    )
    rewrite = native_dynamics.native_rewrite_census()
    ensemble = native_dynamics.binary_ensemble_census()
    values.update(rewrite)
    values.update(ensemble)
    passed = (
        values["pair_ports_share_a_seam"]
        and kernel[pair_state] == Fraction(13, 15) == 1 - Fraction(one_unit_seams, draws)
        and set(hops.values()) == {Fraction(1, 60)}
        and len(hops) == one_unit_seams == 8
        and sum(kernel.values()) == 1
        and values["pair_draws_with_energy_decrease"] == 0
        and census["states"] == 78
        and values["total_transitions"] == 4680
        and census["strict"] == 120
        and census["zero"] == 600
        and census["noop"] == 3960
        and values["excursion_distinct_states"] == len(walk)
        and values["excursion_energy"] == 0.5
        and waits.current != swaps.current
        and waits.energy_decrease == 0
        and swaps.energy_decrease == 0
        and rewrite["ungated_normal_forms"] == 0
        and rewrite["ungated_has_cycle"]
        and rewrite["binary_is_closed_connected_class"]
        and rewrite["all_reach_binary_class"]
        and rewrite["ungated_confluent"]
        and rewrite["gated_strict_descent"]
        and rewrite["gated_normal_forms"] == 66
        and rewrite["gated_normal_forms_are_binary"]
        and rewrite["gated_normal_forms_per_doubled_state"] == 5
        and rewrite["gated_normal_form_count_is_uniform"]
        and ensemble["binary_chain_states"] == 66
        and ensemble["binary_chain_draws"] == draws
        and ensemble["binary_chain_stationary_dimension"] == 1
        and all(value for value in ensemble.values() if isinstance(value, bool))
    )
    return passed, values


def _probe_energy_name_collision() -> tuple[bool, dict]:
    """The raw Gaussian integral and the probability density have different offsets."""
    model = _carrier_gaussian_model()
    records = (0.3, -0.2)
    observed = np.asarray(records)
    blocks = model.blocks
    marginal_energy = model.minimum_interface_energy(records)
    raw_potential = model.unnormalised_record_potential(records)
    surprisal = model.record_surprisal(records)
    mean, covariance = model.posterior(records)
    variational = model.variational_free_energy(mean, covariance, records)
    expected_gap = 0.5 * float(np.linalg.slogdet(blocks["L_xx"])[1]) - 0.5 * model.interior_dim * math.log(
        2.0 * math.pi
    )
    direct = model.energy_minus_entropy(mean, covariance, records)
    prior_energy = 0.5 * _RECORD_PRIOR_PRECISION * float(observed @ observed)

    # Two series chains have the same two-record Schur complement, despite
    # having different interior dimensions and determinants. Unit conductance
    # over three seams and conductance 2/3 over two seams both give 1/3.
    chain = _chain_gaussian_model()
    restriction = math.sqrt(2.0 / 3.0) * np.ones((1, 1))
    weighted_sheaf = sheaf_module.CellularSheaf(
        stalk_dims=(1, 1, 1),
        edges=tuple(
            sheaf_module.EdgeRestriction(i, i + 1, restriction, restriction)
            for i in range(2)
        ),
    )
    weighted_chain = free_energy_module.GaussianSheafModel(
        weighted_sheaf, (0, 2), record_prior_precision=_RECORD_PRIOR_PRECISION,
    )
    chain_schur_agrees = bool(np.allclose(
        chain.schur_complement(), weighted_chain.schur_complement(), atol=1e-12, rtol=0.0
    ))
    chain_raw_difference = (
        chain.unnormalised_record_potential(records)
        - weighted_chain.unnormalised_record_potential(records)
    )
    chain_surprisal_difference = (
        chain.record_surprisal(records) - weighted_chain.record_surprisal(records)
    )
    old_records = (0, 0)
    contrast_records = (1, -1)
    controls = (chain, weighted_chain)
    record_change_divergences = tuple(
        control.bayesian_surprise(old_records, contrast_records) for control in controls
    )
    posterior_divergences = tuple(
        free_energy_module.gaussian_kl(
            *control.posterior(contrast_records), *control.posterior(old_records)
        )
        for control in controls
    )
    # Linear interpolation on the unit-weight path gives the two opposed
    # interior shifts. With zero boundary shifts the squared edge differences
    # are in the ratio 1:4:1. The midpoint of the weighted path stays fixed.
    first_interior_shift = Fraction(contrast_records[0], len(chain.interior_nodes) + 1)
    exact_record_change_divergences = (3 * first_interior_shift**2, Fraction(0))
    values = {
        "record_ports": model.record_nodes,
        "interior_dimension": model.interior_dim,
        "interface_energy_minimised_over_interior": marginal_energy,
        "unnormalised_record_potential": raw_potential,
        "normalised_record_surprisal": surprisal,
        "variational_free_energy_at_posterior": variational,
        "difference": raw_potential - marginal_energy,
        "expected_model_constant": expected_gap,
        "record_prior_precision": _RECORD_PRIOR_PRECISION,
        "record_prior_energy": prior_energy,
        "normalised_record_offset": surprisal - marginal_energy - prior_energy,
        "is_positive_definite_interior": bool(
            np.all(np.linalg.eigvalsh(blocks["L_xx"]) > 0.0)
        ),
        "energy_minus_entropy_at_posterior": direct,
        "chain_control_schur_agrees": chain_schur_agrees,
        "chain_control_raw_potential_difference": chain_raw_difference,
        "chain_control_normalised_surprisal_residual": chain_surprisal_difference,
        "chain_control_contrast_records": contrast_records,
        "chain_control_record_change_divergences": record_change_divergences,
        "chain_control_posterior_divergences": posterior_divergences,
        "chain_control_exact_record_change_divergences": tuple(
            str(value) for value in exact_record_change_divergences
        ),
    }
    passed = (
        values["is_positive_definite_interior"]
        and abs(values["difference"] - expected_gap) < 1e-12
        and abs(direct - surprisal) < 1e-12
        and abs(variational - surprisal) < 1e-12
        and abs(expected_gap) > 1e-3
        and chain_schur_agrees
        and abs(chain_raw_difference) > 0.1
        and abs(chain_surprisal_difference) < 1e-12
        and np.allclose(
            record_change_divergences, posterior_divergences, atol=1e-12, rtol=0.0
        )
        and np.allclose(
            record_change_divergences,
            [float(value) for value in exact_record_change_divergences],
            atol=1e-12,
            rtol=0.0,
        )
    )
    return passed, values


def _probe_two_record_resistance() -> tuple[bool, dict]:
    """Compare Schur marginalisation with a rational grounded-network solve."""
    carrier = _carrier()
    sheaf = sheaf_module.carrier_sheaf(carrier)
    adjacency = _carrier_adjacency(carrier)
    distances = network_module.graph_distances(adjacency)
    difference_form = np.array([[1.0, -1.0], [-1.0, 1.0]])
    resistance_classes = network_module.resistance_ladder(adjacency)
    ladder = {distance: resistance for distance, resistance, _ in resistance_classes}
    pair_counts = tuple(count for _, _, count in resistance_classes)
    max_error = 0.0
    matches_ladder = ladder == {1: Fraction(11, 30), 2: Fraction(7, 15), 3: Fraction(1, 2)}
    for source in range(carrier.port_count):
        for target in range(source + 1, carrier.port_count):
            distance = distances[source][target]
            model = free_energy_module.GaussianSheafModel(sheaf, (source, target))
            expected = difference_form / float(ladder[distance])
            error = np.max(np.abs(model.schur_complement() - expected))
            max_error = max(max_error, float(error))
    values = {
        "resistances_by_distance": tuple(str(value) for value in ladder.values()),
        "pairs_by_distance": pair_counts,
        "max_schur_residual": max_error,
        "exact_resistance_ladder_matches": matches_ladder,
    }
    return matches_ladder and pair_counts == (30, 30, 6) and max_error < 1e-12, values


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
    # The same test on the native law: OPH's two-event observer chooses its second
    # seam from the first record.  The enumeration over all paths and a case
    # analysis that enumerates nothing must agree.
    feedback = control_module.record_feedback_occupation()
    values["record_feedback_paths"] = feedback["paths"]
    values["record_feedback_occupation"] = str(feedback["occupation"])
    values["record_feedback_case_analysis"] = str(feedback["case_analysis"])
    values["record_feedback_second_seam_mass"] = str(feedback["second_seam_mass"])
    values["selector_disconnected_occupation"] = str(feedback["selector_disconnected"])
    values["selector_independent_occupation"] = str(feedback["selector_independent"])
    passed = (
        feedback["paths"] == 264
        and feedback["occupation"] == feedback["case_analysis"] == Fraction(35, 264)
        and feedback["second_seam_mass"] == Fraction(5, 66)
        and feedback["selector_disconnected"] == feedback["selector_independent"] == Fraction(1, 6)
        and values["same_observation"] == 1.0
        and values["response_distance"] == 2.0
        and values["worst_case_error_lower_bound"] == 1.0
    )
    return passed, values


def _probe_existence_without_reachability() -> tuple[bool, dict]:
    """Proposition 15: an agreement can exist and be unreachable."""
    eigenvalues = np.linalg.eigvals(control_module.two_node_relaxation_jacobian())
    orbit = control_module.two_node_relaxation_orbit()
    sync_radius = control_module.two_node_orbit_min_radius(steps=300, skip=60)
    boundary_states = set(product((-1, 0, 1), repeat=2)) - {(0, 0)}
    boundary_orbit = control_module.two_node_sync_orbit(
        start=(1.0, 0.0), steps=len(boundary_states)
    )
    boundary_cycle_verified = (
        {tuple(state) for state in boundary_orbit[:-1]} == boundary_states
        and np.array_equal(boundary_orbit[0], boundary_orbit[-1])
    )
    values = {
        "relaxation_jacobian_eigenvalues": tuple(
            sorted((round(float(real), 12), round(float(imag), 12)) for real, imag in zip(eigenvalues.real, eigenvalues.imag))
        ),
        "relaxation_min_radius": orbit["min_radius"],
        "relaxation_max_radius": orbit["max_radius"],
        "relaxation_period": orbit["period"],
        "sync_map_min_radius": sync_radius,
        "sync_boundary_cycle_verified": boundary_cycle_verified,
        "source_relaxation_min_radius": 0.84,
        "source_relaxation_period": 7.58,
    }
    values["fixed_points"] = control_module.two_node_fixed_points()
    values["sync_transient_steps_skipped"] = 60
    passed = (
        values["relaxation_min_radius"] >= 0.84
        and abs(values["relaxation_period"] - 7.58) < 0.05
        and sync_radius >= 1.0 - 1e-12
        and values["sync_boundary_cycle_verified"]
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
    derived_volumes = geometry.derived_volume_ratios()
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
            name: str(value) for name, value in sorted(derived_volumes.items())
        },
        "derived_table_equals_declared": derived_volumes == geometry.VOLUME_RATIOS,
        "max_shell": geometry.MAX_SHELL,
        "neighbor_squared_distance": geometry.NEIGHBOR_SQUARED_DISTANCE,
        "quadray_cosets": ivm_quadray_cosets(),
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
        and values["quadray_cosets"] == {
            "minimum_squared_distances": (0, 3, 4, 3),
            "nearest_counts": (1, 4, 6, 4),
            "nearest_pair_squared_distance_counts": (
                (), ((8, 6),), ((8, 12), (16, 3)), ((8, 6),),
            ),
        }
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
    values["cuboctahedral_neighbours_touch"] = kissing["cuboctahedron_neighbours_touch"]
    values["cuboctahedral_outer_contact_degrees"] = sorted(
        {sum(row) for row in geometry.cuboctahedron_adjacency()}
    )
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
        and values["cuboctahedral_outer_contact_degrees"] == [4]
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
    exact_vertices = geometry.ivm_solid_vertices()["rhombic_dodecahedron"]
    cell_exact = geometry.convex_volume(exact_vertices)
    values["cell_vertices_match_exact_construction"] = set(cell["vertices"]) == set(exact_vertices)
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
        and values["cell_vertices_match_exact_construction"] is True
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
    derived_volumes = geometry.derived_volume_ratios()
    table = {name: str(value) for name, value in sorted(derived_volumes.items())}
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
            value.denominator == 1 for value in derived_volumes.values()
        ),
        "synergetics_table_derived_exactly": derived_volumes == geometry.VOLUME_RATIOS,
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
        and onward == (5, 2, 1)
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



def _probe_repair_is_conditioning() -> tuple[bool, dict]:
    """On the binary shell each native repair is exact conditioning, and a variational step.

    Route one is the integer rule under its fair coin.  Route two is the uniform
    binary reference conditioned on the datum the repair leaves visible, built
    without the rule.  The free energy over a fibre is then evaluated at the
    repair's law, at biased coins and at the two deterministic rules, and
    compared with the surprisal of the datum plus the divergence from the
    posterior, which is a separate computation.
    """
    census = control_module.repair_conditioning_census()
    entropy = control_module.repair_entropy_identity()
    reference = census["binary_states"]
    one_unit_fibre = 2
    surprisal_one_unit = math.log(reference) - math.log(one_unit_fibre)
    surprisal_settled = math.log(reference)
    at_repair = control_module.fibre_free_energy((0.5, 0.5), reference)
    at_settled = control_module.fibre_free_energy((1.0,), reference)
    biased_ok = True
    for weight in (0.6, 0.75, 0.9, 1.0):
        law = (weight, 1.0 - weight)
        excess = control_module.fibre_free_energy(law, reference) - surprisal_one_unit
        divergence = sum(q * math.log(q / 0.5) for q in law if q > 0.0)
        biased_ok = biased_ok and abs(excess - divergence) < 1e-12 and excess > 0.0
    deterministic = control_module.fibre_free_energy((1.0, 0.0), reference)
    values = dict(census)
    values["kernel_equals_conditional"] = census["mismatches"] == 0
    values["one_unit_cases"] = census["fibre_sizes"].get(2, 0)
    values["settled_cases"] = census["fibre_sizes"].get(1, 0)
    values["repair_is_idempotent"] = census["non_idempotent"] == 0
    values["repair_is_reversible"] = census["asymmetric"] == 0
    values["surprisal_one_unit"] = f"log {reference // one_unit_fibre}"
    values["surprisal_settled"] = f"log {reference}"
    values["free_energy_at_repair_is_surprisal"] = (
        abs(at_repair - surprisal_one_unit) < 1e-12 and abs(at_settled - surprisal_settled) < 1e-12
    )
    values["biased_coin_excess_is_divergence"] = biased_ok
    values["deterministic_excess"] = f"log {one_unit_fibre}"
    values["deterministic_rule_pays_log_two"] = (
        abs(deterministic - surprisal_one_unit - math.log(one_unit_fibre)) < 1e-12
    )
    values["entropy_identity_seams"] = entropy["seams"]
    values["entropy_identity_holds"] = entropy["max_residual"] < 1e-12
    values["relative_entropy_never_rises"] = entropy["smallest_drop"] >= 0.0
    passed = (
        census["binary_states"] == 66
        and census["cases"] == 3960
        and values["kernel_equals_conditional"]
        and census["fibre_sizes"] == {1: 2760, 2: 1200}
        and values["repair_is_idempotent"]
        and values["repair_is_reversible"]
        and values["free_energy_at_repair_is_surprisal"]
        and values["biased_coin_excess_is_divergence"]
        and values["deterministic_rule_pays_log_two"]
        and values["entropy_identity_seams"] == 30
        and values["entropy_identity_holds"]
        and values["relative_entropy_never_rises"]
        and census["off_shell_cases"] == 360
        and census["off_shell_failures"] == 60
        and census["off_shell_binary_reference_holds"] == 60
        and census["off_shell_null_datum"] == 300
    )
    return passed, values


_CARRIER_SOURCE = "OPH evidence/observer_dynamics_20260925/oph-physics-sim (ORIENTED_BASE_FACES, seams, laplacian)"
_LAW_SOURCE = "OPH model/observer/native.py:45-52 and the 'Agreement and Surprise' paper, Section 6.2"
_PAPER_SOURCE = "Mueller, 'Agreement and Surprise: Local Repair, Dynamics and the Limits of Equilibrium', preprint of 22 October 2026 (PhilPapers MUEAAS-2)"
_QUADRAY_SOURCE = (
    "QuadMath src/quadmath/core/quadray.py, src/quadmath/lattice/omni_numbering.py and "
    "the polyhedra table in docs/manuscript/03_quadray_methods.md"
)
_SHEAF_SOURCE = f"{_PAPER_SOURCE}, Section 7.2 and Theorem 7 (the cellular sheaf over the carrier complex)"

def _probe_binary_conditional_variance() -> tuple[bool, dict]:
    """Compare exact Dirichlet and conditional-variance matrix coefficients."""
    values = native_dynamics.binary_dirichlet_census()
    return values == {
        "binary_states": 66,
        "seams_checked": 30,
        "two_state_fibres": 300,
        "coefficient_matrices_agree": True,
        "adjacent_indicator_energy": "1/495",
    }, values


def _probe_carrier_spectral_chords() -> tuple[bool, dict]:
    """Compare exact chord shells with full and projected Gaussian variance."""
    result = network_module.carrier_resistance_geometry()
    values = {
        **result,
        "shells": tuple(
            (step, str(resistance), str(chord), count)
            for step, resistance, chord, count in result["shells"]
        ),
        "lowest_eigenvalue": str(result["lowest_eigenvalue"]),
        "coordinate_norm_squared": str(result["coordinate_norm_squared"]),
        "chord_to_lowest_variance_factor": str(result["chord_to_lowest_variance_factor"]),
        "lowest_mode_variances": tuple(str(value) for value in result["lowest_mode_variances"]),
    }
    passed = (
        result["shells"] == (
            (1, Fraction(11, 30), 4, 30),
            (2, Fraction(7, 15), geometry.Surd5(6, 2), 30),
            (3, Fraction(1, 2), geometry.Surd5(10, 2), 6),
        )
        and values["ranks"] == (11, 3)
        and result["lowest_eigenvalue"] == geometry.Surd5(5, -1)
        and values["lowest_multiplicity"] == 3
        and result["coordinate_norm_squared"] == geometry.Surd5(10, 2)
        and result["chord_to_lowest_variance_factor"] == 40
        and result["lowest_mode_variances"] == (
            Fraction(1, 10),
            geometry.Surd5(Fraction(3, 20), Fraction(1, 20)),
            geometry.Surd5(Fraction(1, 4), Fraction(1, 20)),
        )
        and values["full_variance_proportional"] is False
    )
    return passed, values


def _probe_repair_first_entry() -> tuple[bool, dict]:
    values = native_dynamics.first_entry_census()
    return values == {
        "doubled_states": 12,
        "binary_states": 66,
        "attempts_per_state": 60,
        "entry_blocks_match": True,
        "survival_per_attempt": "5/6",
        "hazard_per_attempt": "1/6",
        "mean_attempts": "6",
        "variance_attempts": "30",
        "entry_pairs_per_doubled_state": 5,
        "uniform_doubled_entry_pairs": 30,
        "uniform_doubled_entry_pair_mass": "1/30",
        "stationary_pair_mass": "1/66",
        "entry_to_stationary_support_ratio": "11/5",
        "entry_to_stationary_total_variation": "6/11",
        "entry_one_port_occupancy": "1/6",
        "stationary_one_port_occupancy": "1/6",
        "entry_one_port_marginals_match": True,
    }, values


def _probe_repair_distance_lumping() -> tuple[bool, dict]:
    values = native_dynamics.repair_distance_lumping_census()
    return values == {
        "class_sizes": (12, 30, 30, 6),
        "attempts_per_state": 60,
        "attempt_count_matrix": ((50, 10, 0, 0), (0, 56, 4, 0), (0, 4, 54, 2), (0, 0, 10, 50)),
        "all_rows_match": True,
        "class_uniform_lift": True,
        "stationary_class_masses": ("0", "5/11", "5/11", "1/11"),
        "stationarity_verified": True,
        "binary_eigenvalue_sum": "5/3",
        "binary_eigenvalue_product": "311/450",
        "binary_eigenvalue_discriminant": "1/75",
    }, values


def _probe_grounded_forest_precision() -> tuple[bool, dict]:
    result = free_energy_module.carrier_grounded_forest_census()
    values = {key: result[key] for key in (
        "pairs_checked", "triples_checked", "pair_layouts", "triple_layouts", "matches_contrast_covariance",
    )}
    passed = (
        values["matches_contrast_covariance"] is True
        and values["pairs_checked"] == 66
        and values["triples_checked"] == 220
        and tuple(row[:3] for row in values["pair_layouts"]) == (
            ((1,), 30, 1900800), ((2,), 30, 2419200), ((3,), 6, 2592000),
        )
        and tuple(row[:3] for row in values["triple_layouts"]) == (
            ((1, 1, 1), 20, 522720), ((1, 1, 2), 60, 604800),
            ((1, 2, 2), 60, 712800), ((1, 2, 3), 60, 743040),
            ((2, 2, 2), 20, 846720),
        )
    )
    return passed, values


def _probe_greedy_record_placement() -> tuple[bool, dict]:
    result = free_energy_module.carrier_grounded_forest_census()
    values = {key: result[key] for key in (
        "best_pair_example", "best_pair_determinant", "best_triple_example",
        "best_triple_determinant", "worst_triple_example", "worst_triple_determinant",
        "greedy_triple_example", "greedy_triple_determinants", "optimum_to_greedy_ratio",
    )}
    values["greedy_triple_determinant"] = result["greedy_triple_determinants"][0][0]
    return result["matches_contrast_covariance"] is True and values == {
        "best_pair_example": (0, 3),
        "best_pair_determinant": 2592000,
        "best_triple_example": (0, 2, 8),
        "best_triple_determinant": 846720,
        "worst_triple_example": (0, 1, 5),
        "worst_triple_determinant": 522720,
        "greedy_triple_example": (0, 1, 2),
        "greedy_triple_determinant": 743040,
        "greedy_triple_determinants": ((743040, 60),),
        "optimum_to_greedy_ratio": "49/43",
    }, values


def _probe_planar_dual_covariance() -> tuple[bool, dict]:
    result = planar_duality.planar_duality_census()
    values = {
        "families": {
            name: {
                **row,
                "primal_resistances": tuple(map(str, row["primal_resistances"])),
                "dual_resistances": tuple(map(str, row["dual_resistances"])),
            }
            for name, row in result.items()
        },
        "matrix_entries": sum(row["matrix_entries_checked"] for row in result.values()),
    }
    for name in ("carrier", "ivm"):
        row = values["families"][name]
        values[f"{name}_resistances"] = row["primal_resistances"] + row["dual_resistances"]
    return values["families"] == {
        "carrier": {
            "edge_count": 30, "cut_ranks": (11, 19), "tree_counts": (5184000, 5184000),
            "primal_resistances": ("11/30",), "dual_resistances": ("19/30",),
            "matrix_entries_checked": 900, "covariances_sum_to_identity": True,
        },
        "ivm": {
            "edge_count": 24, "cut_ranks": (11, 13), "tree_counts": (331776, 331776),
            "primal_resistances": ("11/24",), "dual_resistances": ("13/24",),
            "matrix_entries_checked": 576, "covariances_sum_to_identity": True,
        },
    }, values


def _probe_cuboctahedral_spectral_chords() -> tuple[bool, dict]:
    result = network_module.cuboctahedral_resistance_geometry()
    values = {
        "spectrum": result["spectrum"],
        "ranks": result["ranks"],
        "chord_to_lowest_variance_factor": int(result["chord_to_lowest_variance_factor"]),
        "lowest_mode_variances": tuple(str(row[3]) for row in result["pair_classes"]),
        "kirchhoff_index": int(result["kirchhoff_index"]),
    }
    passed = (
        result["chord_to_lowest_variance_factor"] == 32
        and result["kirchhoff_index"] == 37
        and values["spectrum"] == ((0, 1), (2, 3), (4, 3), (6, 5))
        and values["ranks"] == (11, 3)
        and values["lowest_mode_variances"] == ("1/8", "1/4", "3/8", "1/2")
    )
    return passed, values


def _probe_cuboctahedral_distance_split() -> tuple[bool, dict]:
    result = network_module.cuboctahedral_resistance_geometry()
    classes = tuple(
        (distance, int(chord), str(resistance), count)
        for distance, chord, resistance, _, count in result["pair_classes"]
    )
    distance_two = tuple(row for row in classes if row[0] == 2)
    values = {
        "pair_classes": classes,
        "distance_two_chords": tuple(row[1] for row in distance_two),
        "distance_two_resistances": tuple(row[2] for row in distance_two),
    }
    return classes == ((1, 4, "11/24", 24), (2, 8, "7/12", 12),
                       (2, 12, "5/8", 24), (3, 16, "2/3", 6)), values


def _probe_repair_contrast_quotient() -> tuple[bool, dict]:
    values = native_dynamics.repair_contrast_lumping_census()
    expected = {
        "carrier": (60, ((56, 4, 0), (4, 54, 2), (0, 10, 50)), ("5/11", "5/11", "1/11")),
        "cuboctahedron": (
            48, ((44, 2, 2, 0), (4, 40, 4, 0), (2, 2, 42, 2), (0, 0, 8, 40)),
            ("4/11", "2/11", "4/11", "1/11"),
        ),
    }
    checks = ("contrast_lumpable", "chord_partition_matches", "all_rows_match_token_moves",
              "class_uniform_lift", "stationarity_verified")
    return all(
        (row["attempts_per_state"], row["attempt_count_matrix"], row["stationary_class_masses"])
        == expected[name] and all(row[key] for key in checks)
        for name, row in values.items()
    ), values


def _probe_repair_distance_counterexample() -> tuple[bool, dict]:
    cases = native_dynamics.repair_contrast_lumping_census()
    carrier, fcc = cases["carrier"], cases["cuboctahedron"]
    witness_rates = tuple(Fraction(row[-1]) for row in fcc["distance_two_witnesses"])
    sizes = tuple(row[3] for row in fcc["classes"] if row[0] == 2)
    weights = tuple(Fraction(size, sum(sizes)) for size in sizes)
    stationary_rate = sum(weight * rate for weight, rate in zip(weights, witness_rates, strict=True))
    history_rate = sum(
        weight * rate**2 for weight, rate in zip(weights, witness_rates, strict=True)
    ) / stationary_rate
    values = {
        "carrier_distance_lumpable": carrier["distance_lumpable"],
        "fcc_distance_lumpable": fcc["distance_lumpable"],
        "fcc_distance_two_witnesses": fcc["distance_two_witnesses"],
        "stationary_adjacency_probability": stationary_rate,
        "adjacency_probability_after_previous_adjacency": history_rate,
    }
    return (
        carrier["distance_lumpable"] and not fcc["distance_lumpable"]
        and fcc["distance_two_adjacency_counts_match"]
        and tuple(row[1:] for row in fcc["distance_two_witnesses"])
        == (("8", "7/12", 2, "1/12"), ("12", "5/8", 1, "1/24"))
        and (stationary_rate, history_rate) == (Fraction(1, 18), Fraction(1, 16))
    ), values


def _probe_reciprocal_dual_weights() -> tuple[bool, dict]:
    values = weighted_duality.weighted_planar_duality_census()
    expected = {
        "carrier": (900, Fraction(44, 63), Fraction(19, 63), (10886400, 2721600)),
        "ivm": (576, Fraction(44, 57), Fraction(13, 57), (787968, 196992)),
    }
    return all(
        (row["entries_checked"], row["primal_variance"], row["reciprocal_dual_variance"],
         row["weighted_tree_counts"]) == expected[name]
        and row["reciprocal_defect_entries"] == row["rescaled_reciprocal_defect_entries"] == 0
        and row["rank_one_checks_agree"]
        for name, row in values.items()
    ), values


def _probe_copied_dual_weight_obstruction() -> tuple[bool, dict]:
    results = weighted_duality.weighted_planar_duality_census()
    values = {name: {key: row[key] for key in (
        "edge", "primal_weight", "copied_variance_sum", "copied_defect_entries",
        "copied_cut_pairing_nonzeros", "reciprocal_raw_variance_sum",
    )} for name, row in results.items()}
    return (
        (values["carrier"]["copied_variance_sum"], values["ivm"]["copied_variance_sum"])
        == (Fraction(2872, 1827), Fraction(1912, 1197))
        and (values["carrier"]["copied_defect_entries"], values["ivm"]["copied_defect_entries"])
        == (676, 484)
        and all(row["primal_weight"] == 4 and row["copied_cut_pairing_nonzeros"] == 4
                for row in values.values())
    ), values


def _probe_quadratic_even_field() -> tuple[bool, dict]:
    values = quadratic_modes.quadratic_mode_census()
    expected = {
        "carrier": ((Fraction(8, 5),) * 5, (6,) * 5, (Fraction(48, 5),) * 5),
        "fcc": ((1, 1, 2, 2, 2), (6, 6, 4, 4, 4), (6, 6, 8, 8, 8)),
    }
    checks = ("linear_frame_isotropic", "quadratic_modes_centered", "frame_orthogonal",
              "precision_orthogonal", "edge_and_laplacian_precision_agree",
              "laplacian_preserves_modes")
    return all(
        (row["frame"], row["laplacian_eigenvalues"], row["precision"]) == expected[name]
        and all(row[key] for key in checks)
        for name, row in values.items()
    ), values


def _probe_quadratic_isotropy_obstruction() -> tuple[bool, dict]:
    cases = quadratic_modes.quadratic_mode_census()
    values = {
        name: {
            "isotropic_linear_frame": row["linear_frame_isotropic"],
            "diagonal_and_mixed_variances": (row["variance"][0], row["variance"][2]),
            "rotated_variance_ratio": row["variance"][0] / row["variance"][2],
        }
        for name, row in cases.items()
    }
    return values == {
        "carrier": {"isotropic_linear_frame": True,
                    "diagonal_and_mixed_variances": (Fraction(5, 48),) * 2,
                    "rotated_variance_ratio": 1},
        "fcc": {"isotropic_linear_frame": True,
                "diagonal_and_mixed_variances": (Fraction(1, 6), Fraction(1, 8)),
                "rotated_variance_ratio": Fraction(4, 3)},
    }, values


def _probe_central_plane_screening() -> tuple[bool, dict]:
    values = geometric_screening.central_plane_blanket_census()
    carrier, fcc = values["icosahedron"], values["cuboctahedron"]
    return (
        values["oph_isomorphism"] and carrier["inverse_matches_seams"] and fcc["inverse_matches_seams"]
        and carrier["plane_classes"] == ((4, 4, 4, 2, "10/9", 15),)
        and fcc["plane_classes"] == ((4, 4, 4, 0, "1", 3), (6, 3, 3, 0, "1", 4))
        and (carrier["screening_planes"], fcc["screening_planes"]) == (0, 7)
        and (carrier["minimum_separator_size"], fcc["minimum_separator_size"]) == (5, 4)
    ), values


def _probe_repair_centroid_memory() -> tuple[bool, dict]:
    values = centroid_dynamics.centroid_gaussian_census()
    expected = {
        "carrier": ("11/12 + 1/60 sqrt 5", "50/33 + 10/33 sqrt 5", 6, "1/180", "1/60"),
        "cuboctahedron": ("23/24", "40/33", 6, "1/72", "1/48"),
    }
    fields = ("retention", "one_axis_variance", "zero_fibre_size", "next_target_given_zero",
              "next_target_given_previous_target_and_zero")
    checks = ("graph_matches_reference", "uniform_stationary", "coordinate_mode_verified",
              "conditional_mean_verified", "stationary_variance_verified",
              "witness_matches_edge_moves", "witness_matches_reversible_formula",
              "stationary_centroid_non_markov")
    return all(
        tuple(row[field] for field in fields) == expected[name] and all(row[key] for key in checks)
        for name, row in values.items()
    ), values


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
        claim="The carrier Laplacian L = 5I - A defines the seam-disagreement energy and the mean port-load operator I - L / 60. In the Gaussian reading, L is an intrinsic precision and is positive definite on the zero-sum subspace. Its first positive eigenvalue is 5 - sqrt(5). The mean port-load operator has gap (5 - sqrt(5)) / 60 under uniform seam sampling and a fair coin.",
        oph="The Laplacian gives the rate of relaxation of mean port loads under the uniform attempt protocol.",
        fep="The scalar sheaf Laplacian gives the Gaussian quadratic form. The centred intrinsic model restricts that form to the zero-sum subspace.",
        synergetics="The icosahedron's graph: a twelve-vertex kissing arrangement, but not the A3 coordination shell, which is the cuboctahedron.",
        status="exact",
        source=_CARRIER_SOURCE,
        probe=_probe_spectral_gap,
    ),
    Correspondence(
        key="shell-counts",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="The total-two load shell contains 78 configurations, of which 66 are binary. Enumeration agrees with C(13, 2) = 78 by stars and bars and C(12, 2) = 66 by choice of occupied ports. These states define the native-rule transition census and its binary conditioning subset.",
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
        claim="On the carrier with a proper Gaussian record prior, harmonic extension, a separate KKT constrained solve and predictive-coding flow agree with the posterior interior mean. The variational free energy computed from the normalised joint density equals record surprisal at the posterior. With a nonquadratic asymmetric seam the energy minimiser and posterior mean separate.",
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
        claim="On the carrier with a Gaussian record prior, the update residual bounds Bayesian surprise through the extreme interior eigenvalues. The energy drop at fixed new records and the divergence between successive Gaussian posteriors both equal the residual quadratic. The normalised variation of record surprisal uses the Schur complement plus the record prior precision.",
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
        source=_SHEAF_SOURCE,
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
        claim="The grounded carrier Laplacian has determinant 5184000, the number of spanning trees by Kirchhoff's matrix-tree theorem. The product of the nonzero Laplacian eigenvalues divided by 12 gives the same count. The grounded Gaussian integral is proportional to the inverse square root of this determinant.",
        oph="The seam structure of the twelve-port carrier, whose Laplacian generates mean repair.",
        fep="As a grounded Gaussian precision, the matrix contributes the inverse square root of its determinant to the interface integral.",
        synergetics="The icosahedron's own graph, counted by its seam structure.",
        status="exact",
        source=_CARRIER_SOURCE,
        probe=_probe_tree_number,
    ),
    Correspondence(
        key="gmrf-independence",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="The proper Gaussian with precision L + I/n has the carrier's seam graph as its conditional-dependence graph. Each port's blanket has five neighbours and each seam's precision entry is minus one. Separately, the variance of a port's deviation from the mean in the centered intrinsic field is 7/36, the diagonal of L^+.",
        oph="Each seam touched by a local repair is one conditional-dependence edge in the proper Gaussian.",
        fep="For this proper Gaussian, a port's Markov blanket is its graph neighbourhood. Its precision comes from the scalar interface energy and the diagonal regularisation.",
        synergetics="Every vertex of the icosahedron is five-coordinated.",
        status="exact",
        source=f"{_CARRIER_SOURCE}; {_PAPER_SOURCE}, Section 7.2",
        probe=_probe_gmrf_independence,
    ),
    Correspondence(
        key="sheaf-cohomology",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="For the carrier's scalar sheaf, dim H^0 = 1 counts the constant global sections and dim H^1 = 19 gives the dimension of the space of seam assignments modulo coboundaries. The alternating sum is 12 - 30. The twenty face boundaries span the cycle space, so the triangulated surface has Betti numbers (1, 0, 1) and Euler characteristic 2. OPH's oriented faces sum to a fundamental cycle.",
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
        claim="The carrier's automorphism group has order 120 and acts transitively on the twelve ports, thirty seams and twenty faces. Its action on the oriented source faces distinguishes sixty rotations from sixty orientation-reversing symmetries.",
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
        claim="The carrier Laplacian has four eigenspaces of dimensions 1, 3, 5 and 3, corresponding to the eigenvalues 0, 5 - sqrt(5), 6 and 5 + sqrt(5). These spaces are the irreducible components of the carrier's multiplicity-free permutation representation. Any symmetric operator invariant under the group acts as a scalar on each component. Those scalars can coincide, in which case the corresponding eigenspaces combine.",
        oph="The dimensions of the irreducible components determine the multiplicities of the distinct Laplacian eigenvalues used in the mean repair certificate.",
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
        claim="Equal-covariance Gaussian KL is half the squared Mahalanobis distance. At unit covariance in the integer quadray embedding, the IVM's first shell is at squared distance eight, giving four nats of divergence. Displacements k m along a neighbour direction have squared norm 8 k^2 and divergence 4 k^2. The calculation compares endpoints under the chosen covariance.",
        oph="The native law evaluates a different quantity, the change in quadratic load energy.",
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
        claim="On the total-load-two shell the ungated native rule is confluent but does not terminate: mismatch-one moves preserve load energy. Gating at broken seams terminates but does not give unique normal forms.",
        oph="The native law has no normal form on this shell. The terminating, nonconfluent Lean demonstration is a separate rule. OPH's consensus theorems apply to the integer and transactional contracts.",
        fep="A load-energy plateau along a path is compatible with exact conditioning on the binary shell and with decreasing ensemble relative entropy. Under the uniform attempt schedule, the binary chain has a unique uniform stationary law and the port relaxation gap.",
        synergetics="No counterpart.",
        status="obstructed",
        source="OPH model/PREMISES.md:95 (PR-01) and Lean/ObserverPatchHolography/Primitives.lean:1198 (demoCarrier_not_confluent)",
        probe=_probe_pr01,
    ),
    Correspondence(
        key="free-energy-name-collision",
        kind=DISAGREEMENT,
        pair=OPH_FEP,
        claim="Minimum interface energy, the negative log of the unnormalised interior integral, and normalised record surprisal are distinct. The second includes the interior determinant; in the third it cancels against the joint normaliser. Two chains with equal Schur complements and the same record prior have the same normalised record density, yet their interior posterior divergences differ under the same contrast change.",
        oph="Minimum of the quadratic seam disagreement at fixed records.",
        fep="Variational free energy bounds normalised record surprisal, with equality at the posterior.",
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
        oph="Eliminating a mediator can induce a boundary term that cannot be represented by pairwise seam energies.",
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
        fep="The zero precision entry gives conditional independence at stationarity. The tested drift can still couple exterior and interior directly.",
        synergetics="No counterpart.",
        status="obstructed",
        source=f"{_PAPER_SOURCE}, Proposition 11; Da Costa et al., Bayesian Mechanics for Stationary Processes",
        probe=_probe_stationary_blanket,
    ),
    Correspondence(
        key="solenoidal-freedom",
        kind=DISAGREEMENT,
        pair=OPH_FEP,
        claim="One stationary density is shared by a family of laws: the diffusion dx = -(Gamma - Q) H x dt + sqrt(2 Gamma) dW has stationary law N(0, H^{-1}) for every antisymmetric Q, and its entropy production rate is tr(Q^T Gamma^{-1} Q H); in two dimensions with Gamma = H = I and Q = q [[0, 1], [-1, 0]] the spectrum is -1 +/- i q and the entropy production is 2 q^2.",
        oph="The Gaussian example motivates checking whether the same reference law can support different native transition kernels.",
        fep="A stationary law does not fix the law that preserves it; the solenoidal part is a free parameter.",
        synergetics="No counterpart.",
        status="obstructed",
        source=f"{_PAPER_SOURCE}, Proposition 9; Hwang, Hwang-Ma and Sheu (1993); Lelièvre, Nier and Pavliotis (2013)",
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
        oph="In the two-event native observer, the first signed seam current selects the second seam.",
        fep="An observation-only state can fail the response criterion while every solve converges.",
        synergetics="No counterpart.",
        status="obstructed",
        source=f"{_PAPER_SOURCE}, Proposition 4 and Section 4.2; OPH model/observer/native.py (record-controlled selector)",
        probe=_probe_reduction,
    ),
    Correspondence(
        key="existence-without-reachability",
        kind=DISAGREEMENT,
        pair=OPH_FEP,
        claim="Exact case analysis gives the origin as the only fixed point of the clipped two-node map. It is a source for the associated continuous relaxation, for which the probe estimates a radius range and period from a finite integration window. The synchronous map has an exact cycle on the clipping boundary and can therefore fail to reach the fixed point.",
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
        claim="For J = [[-1, 2], [0, -1]] and a state loss derivative c = (1, 0), the direct nudging response is (1, 0) while the adjoint is (1, 2). The adjoint gives the loss gradient with respect to an additive drive.",
        oph="This test concerns sensitivity of a directed equilibrium to a perturbation.",
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
        claim="Quadray coordinates give exact IVM geometry: twelve neighbour moves at squared distance eight, shell counts 10k^2 + 2, and the tetravolume table 1, 3, 4, 6, 20 and 96 from explicit vertex sets. The coordinate sum modulo four distinguishes the FCC sites from the two tetrahedral hole classes and the octahedral hole class.",
        oph="The icosahedral carrier realises a twelve-sphere kissing arrangement with five seams per port. The FCC neighbour shell supplies a separate geometric comparison.",
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
        claim="For the scalar energy in the probe, the contrast between parameter derivatives at the free and nudged equilibria equals (theta - target) / (1 + beta). It converges to the loss gradient as the nudge vanishes. The contrast illustrates how an energy-based learning method can obtain a parameter gradient from equilibria.",
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
        claim="For a symmetric payoff A, the replicator field is the Shahshahani gradient of F(x) = x^T A x / 2, with mobility M = diag(x) - x x^T. The probe checks the identity and the nonnegative potential rate. The tested asymmetric payoff breaks this quadratic-potential identity at every sampled point.",
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
        claim="Twelve is the kissing number in three dimensions and both arrangements attain it. With the twelve sphere centres at distance 2r from a central sphere of radius r, the IVM's cuboctahedral shell has edge equal to circumradius, so each outer sphere touches the centre and four other outer spheres. The icosahedral arrangement has squared edge over squared circumradius exactly 2 - 2 sqrt(5) / 5, so its outer spheres touch only the centre. The geometric icosahedral realisation is isomorphic to OPH's labelled carrier.",
        oph="The twelve-port carrier is the icosahedral kissing arrangement, the one with slack between neighbours.",
        fep="No fixed carrier.",
        synergetics="The vector equilibrium has contacts between outer spheres. The icosahedral arrangement has gaps between them.",
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
        claim="The IVM's Voronoi cell is a rhombic dodecahedron with fourteen vertices, twelve rhombic faces, twenty-four edges and Euler characteristic two. Its exact tetravolume is six. The inscribed ball has tetravolume pi sqrt 2, so the face-centred cubic packing density is pi / (3 sqrt 2).",
        oph="The IVM's cuboctahedral arrangement and the carrier's icosahedral arrangement share the kissing number. Their outer contacts differ.",
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
        claim="Under their stated conventions, the three calculations give distinct kinds of exact quantity. The unit-covariance divergence of an IVM neighbour displacement is four nats. The mean repair operator has gap (5 - sqrt 5) / 60. The synergetic tetravolume table is integral (1, 3, 4, 6, 20, 96), while the icosahedron's tetravolume is sqrt 2 (15/2 + (5/2) sqrt 5), with square 175 + 75 sqrt 5. The values describe a displacement, relaxation and volume respectively.",
        oph="Integer loads; spectral data in the quadratic field Q(sqrt 5).",
        fep="Real-valued information costs in nats.",
        synergetics="Tetravolumes of the prime vector: integral on the IVM, algebraic irrational on the icosahedron.",
        status="exact",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}; {_PAPER_SOURCE}, Section 6.2 (mean gap only); arithmetic comparison derived here",
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
        source=f"{_CARRIER_SOURCE}; resistance ladder derived here from the carrier Laplacian",
        probe=_probe_effective_resistance,
    ),
    Correspondence(
        key="two-record-resistance",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="On every pair of carrier ports the two-record Schur complement is the difference form divided by the effective resistance. The minimum interface energy therefore equals the squared record contrast divided by twice the resistance. A proper record prior adds its own quadratic and the resulting record log normaliser.",
        oph="The carrier seam network seen from a fixed pair of records.",
        fep="The Schur term that sets the contrast dependence of the record marginal.",
        synergetics="The graph distance classes select the three resistance values.",
        status="exact",
        source=f"{_CARRIER_SOURCE}; Schur complement of the interface precision and the grounded Laplacian definition of resistance",
        probe=_probe_two_record_resistance,
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
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}; conjugation and norm derived here",
        probe=_probe_conjugation_norm,
    ),
    Correspondence(
        key="repair-is-conditioning",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="On the binary shell each native repair, under its fair coin, is exact Bayesian conditioning: its law is the uniform reference conditioned on the datum the repair leaves visible (the seam's total load and every load off the seam), in all 3960 cases over the 66 configurations and both orientations of the thirty seams. That law is the unique minimiser of the variational free energy F[q] = E_q[ln q(s) - ln pi(s, o)] over the fibre of the datum, and the minimum is the surprisal -ln pi(o): log 33 on a seam that holds one unit, log 66 otherwise. A biased coin pays exactly its divergence from the fair one and a rule that always waits or always swaps pays log 2. Under a repair K the relative entropy to the reference falls by exactly D(p || pK). The identity belongs to the shell: with the uniform law on all 78 configurations as reference it fails in sixty state-and-seam cases, one for each doubled port and incident seam. The repair then moves one unit across with certainty.",
        oph="On the binary shell each native repair resamples the seam within the fibre of its repaired visible datum for the uniform reference, and each seam preserves that reference (PR-07).",
        fep="A variational step with the uniform reference as prior and the repaired datum as observation: the recognition density is the exact posterior, so the bound is tight and its value is the surprisal of the datum.",
        synergetics="No counterpart.",
        status="exact-on-domain",
        source=f"{_PAPER_SOURCE}, Sections 6.2 and 7.3; OPH model/PREMISES.md:100 (PR-07) and the synthesis paper's section on entropy and thermodynamic laws from repair",
        probe=_probe_repair_is_conditioning,
    ),
    Correspondence(
        key="binary-conditional-variance",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="For uniform seam choice and fair repair on the binary shell, the Dirichlet form of native repair is exactly the mean conditional variance over seam records. The two quadratic forms have identical coefficients on all 66 states.",
        oph="The random-seam native repair kernel on the shell of two occupied ports.",
        fep="The mean residual variance after conditioning the uniform reference on a seam record.",
        synergetics="The carrier's icosahedral seams index the resampling steps.",
        status="exact-on-domain",
        source=f"{_CARRIER_SOURCE}; Chlebicka, Łatuszyński and Miasojedow (2025), §2.2, for Gibbs projections; finite-shell comparison derived here",
        probe=_probe_binary_conditional_variance,
    ),
    Correspondence(
        key="carrier-spectral-chords",
        kind=ISOMORPHISM,
        pair=FEP_SYNERGETICS,
        claim="On the edge-length-two regular icosahedron, each squared chord equals 40 times the contribution of the first positive Laplacian eigenspace to the intrinsic unit-precision Gaussian's port contrast variance. The full variance is effective resistance and is not proportional to squared chord length. The squared-distance embedding ranks are 11 and 3, respectively.",
        oph="The labelled carrier graph supplies the seam Laplacian and its slowest nonconstant mean-load modes.",
        fep="Projection of the intrinsic Gaussian onto the first positive Laplacian eigenspace.",
        synergetics="Euclidean chords of the regular icosahedron at the stated edge scale.",
        status="exact",
        source=f"{_CARRIER_SOURCE}; Osting (2023) for the regular-polyhedron eigenspace; Klein and Randić (1993) for resistance; exact comparison derived here",
        probe=_probe_carrier_spectral_chords,
    ),
    Correspondence(
        key="repair-first-entry",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="Under uniform fair attempts on the load-two shell, the exit time from a doubled state is geometric with parameter 1/6 and independent of the first adjacent pair entered. A uniform doubled start gives a uniform law over seams at entry. Binary stationarity gives the same one-port marginals, with total variation distance 6/11 between the two joint laws.",
        oph="A doubled load splits deterministically when an incident seam is selected.",
        fep="First-entry and stationary ensembles give different joint distributions with identical single-port marginals.",
        synergetics="The icosahedral graph determines how many seams meet the initial port.",
        status="exact-on-domain",
        source=f"{_CARRIER_SOURCE}; Kemeny and Snell (1976), §3.2 for first passage; exact transient and entrance blocks derived here",
        probe=_probe_repair_first_entry,
    ),
    Correspondence(
        key="repair-distance-lumping",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="Occupied-port graph distance, with zero for a double, is a strong four-state lumping of the uniform fair repair chain on the load-two shell. The quotient also lifts class-uniform laws to class-uniform laws. On the binary block its nonconstant modes are (25 ± sqrt 3)/30.",
        oph="Distance between the two token positions is a record of the microstate.",
        fep="Exact Markov reduction for arbitrary initial laws. The class-uniform lift preserves relative entropy to binary stationarity.",
        synergetics="The carrier intersection array supplies the interclass attempt counts.",
        status="exact-on-domain",
        source=f"{_CARRIER_SOURCE}; Kemeny and Snell (1976), theorem 6.3.2 for strong lumpability; exact quotient derived here",
        probe=_probe_repair_distance_lumping,
    ),
    Correspondence(
        key="grounded-forest-precision",
        kind=ISOMORPHISM,
        pair=OPH_FEP,
        claim="With a nonempty set of clamped records in the unit Gaussian carrier, the determinant of the interior precision counts spanning forests with one record in each component. All pair and triple record layouts agree with an independent determinant of the anchored record-contrast covariance times the tree number.",
        oph="Record placement grounds a chosen set of carrier ports.",
        fep="The grounded determinant fixes conditional Gaussian entropy for a fixed number of records.",
        synergetics="Distances between the selected ports classify the tested layouts.",
        status="exact-on-domain",
        source=f"{_CARRIER_SOURCE}; Chaiken (1982) for the all-minors matrix-tree theorem; exact layout census derived here",
        probe=_probe_grounded_forest_precision,
    ),
    Correspondence(
        key="greedy-record-placement",
        kind=DISAGREEMENT,
        pair=OPH_FEP,
        claim="On the unit Gaussian carrier, extending any optimal pair of clamped records by a third gives determinant 743040. Every optimal triple has determinant 846720. The greedy conditional entropy exceeds the optimum by half the logarithm of 49/43.",
        oph="Antipodal pairs maximise the pair score. Optimal triples have pair distances (2, 2, 2).",
        fep="Greedy minimisation of conditional entropy by adding records can miss every optimum at the next budget.",
        synergetics="The optimal pair and triple layouts have distinct carrier geometries.",
        status="obstructed",
        source=f"{_CARRIER_SOURCE}; exhaustive grounded determinant census here",
        probe=_probe_greedy_record_placement,
    ),
    Correspondence(
        key="planar-dual-covariance",
        kind=ISOMORPHISM,
        pair=THREE_WAY,
        claim="With corresponding edge orientations and unit precision, the primal and dual intrinsic Gaussian edge covariances sum to the identity. Exact projectors are checked for the carrier and its dodecahedral dual, and for the FCC shell and its Voronoi dual. In each pair the spanning-tree counts agree.",
        oph="The oriented carrier faces induce a dodecahedral dual with one matched edge per seam.",
        fep="Gradient covariances are orthogonal projectors onto complementary subspaces in the common edge space.",
        synergetics="The same identity holds for the cuboctahedral shell and its rhombic dodecahedral dual.",
        status="exact-on-domain",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}; Furrer, Hungerbühler and Jantschgi (2020) for dual resistance; exact full matrices derived here",
        probe=_probe_planar_dual_covariance,
    ),
    Correspondence(
        key="cuboctahedral-spectral-chords",
        kind=ISOMORPHISM,
        pair=FEP_SYNERGETICS,
        claim="On the edge-length-two FCC first shell, squared chords equal 32 times the first positive Laplacian contribution to contrast variance in the centred unit-precision Gaussian. Polynomial projectors and grounded solves agree on the full variance of every pair.",
        oph="The carrier supplies a comparison at the same edge scale, with factor 40.",
        fep="The first nonzero covariance modes give an exact squared-distance embedding of the cuboctahedral chords.",
        synergetics="The first coordination shell has an explicit integer quadray realisation.",
        status="exact-on-domain",
        source=f"{_QUADRAY_SOURCE}; Klein and Randić (1993) for effective resistance; exact projector and chord comparison derived here",
        probe=_probe_cuboctahedral_spectral_chords,
    ),
    Correspondence(
        key="cuboctahedral-distance-split",
        kind=DISAGREEMENT,
        pair=FEP_SYNERGETICS,
        claim="On the FCC first shell, pairs two graph steps apart have squared chords 8 or 12 at edge length two. They have intrinsic Gaussian contrast variances 7/12 and 5/8, respectively. The full variance thus depends on more than graph distance on this shell.",
        oph="The icosahedral carrier gives one resistance per distance class.",
        fep="The Gaussian contrast variance separates two pair types at the same number of graph steps.",
        synergetics="The geometry of the cuboctahedron distinguishes the two types.",
        status="obstructed",
        source=f"{_QUADRAY_SOURCE}; all-pairs resistance comparison derived here",
        probe=_probe_cuboctahedral_distance_split,
    ),
    Correspondence(
        key="repair-contrast-quotient",
        kind=ISOMORPHISM,
        pair=THREE_WAY,
        claim="The intrinsic Gaussian contrast variance of the loaded pair is a strongly lumpable record for fair, uniform binary seam repair on both the OPH icosahedron and the cuboctahedral FCC first shell. Its classes are the squared-chord classes and have stationary masses proportional to their sizes.",
        oph="Fair binary repair is run with a uniform edge and a fair coin.",
        fep="Unit intrinsic contrast variance gives a Markov record with an exact stationary measure.",
        synergetics="The same local rule on the first IVM shell has four geometric pair classes.",
        status="exact-on-domain",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}; Kemeny and Snell (1976) for strong lumpability; native enumeration and hop comparison derived here",
        probe=_probe_repair_contrast_quotient,
    ),
    Correspondence(
        key="repair-distance-counterexample",
        kind=DISAGREEMENT,
        pair=THREE_WAY,
        claim="Graph distance is a Markov record for fair binary repair on the OPH carrier. On the FCC comparison graph the two distance-two classes have different fair repair probabilities of reaching adjacency in one attempt, 1/12 and 1/24, and different intrinsic contrast variances, 7/12 and 5/8.",
        oph="The native one-attempt moves follow the same fair edge protocol on both graphs.",
        fep="The two distinguishable contrast variances repair the failed distance-only reduction.",
        synergetics="Two kinds of cuboctahedral pairs are two steps apart and have different common-neighbour counts.",
        status="obstructed",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}; exact hop and variance enumeration derived here",
        probe=_probe_repair_distance_counterexample,
    ),
    Correspondence(
        key="reciprocal-dual-weights",
        kind=ISOMORPHISM,
        pair=THREE_WAY,
        claim="For both matched spherical graph pairs, the conductance-normalised intrinsic edge covariance projectors sum to identity exactly when the product of matched positive conductances is constant. Reciprocal and globally rescaled reciprocal weights pass a full exact covariance check.",
        oph="Opposite boundary incidences give signed edge coordinates for the carrier and its face dual.",
        fep="Weighting edge differences by the square root of conductance makes their Gaussian covariance an orthogonal projector.",
        synergetics="The first FCC shell and the rhombic dodecahedral Voronoi cell form the second matched graph pair.",
        status="exact-on-domain",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}; Furrer, Hungerbühler and Jantschgi (2020), weighted diagonal identity; full projector criterion and checks derived here",
        probe=_probe_reciprocal_dual_weights,
    ),
    Correspondence(
        key="copied-dual-weight-obstruction",
        kind=DISAGREEMENT,
        pair=THREE_WAY,
        claim="Copying a nonuniform positive conductance pattern to the dual destroys Gaussian edge covariance complementarity on both pairs. A single conductance changed to four suffices on each graph.",
        oph="A single carrier seam has increased conductance in the test.",
        fep="Copying that conductance onto its dual edge makes the sum of normalized diagonal variances exceed one.",
        synergetics="The same obstruction occurs for the cuboctahedral shell and its Voronoi dual.",
        status="obstructed",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}; nonuniform weight obstruction and exact witnesses derived here",
        probe=_probe_copied_dual_weight_obstruction,
    ),
    Correspondence(
        key="quadratic-even-field",
        kind=ISOMORPHISM,
        pair=THREE_WAY,
        claim="Evaluating trace-free quadratic forms at unit vertex directions identifies their five-dimensional coordinate space with the centred antipodally even fields on either shell. Pulling back the intrinsic Gaussian precision gives 48/5 on the carrier and separates diagonal and off-diagonal precisions of 6 and 8 on the FCC shell.",
        oph="The antipodal pairs of the native icosahedral carrier define the even port fields.",
        fep="The Laplacian quadratic pulls back to the Gaussian precision of the trace-free coefficients.",
        synergetics="Quadratic evaluations at the unit directions of the first IVM shell provide the comparison.",
        status="exact-on-domain",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}; Delsarte, Goethals and Seidel (1977) for spherical designs; exact frames and precisions derived here",
        probe=_probe_quadratic_even_field,
    ),
    Correspondence(
        key="quadratic-isotropy-obstruction",
        kind=DISAGREEMENT,
        pair=THREE_WAY,
        claim="The carrier and FCC shell both have isotropic linear evaluation frames. Quadratic Gaussian coefficient covariance remains isotropic on the carrier but is anisotropic on the FCC shell: rotating a unit diagonal form into an off-diagonal one changes its variance from 1/6 to 1/8.",
        oph="The carrier retains equal variance across all normalized trace-free forms.",
        fep="An isotropic linear frame does not fix rotational invariance for the covariance of quadratic coefficients.",
        synergetics="The fourth moments of the cuboctahedral directions distinguish Frobenius-unit forms related by a rotation.",
        status="obstructed",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}; exact quadratic covariance comparison derived here",
        probe=_probe_quadratic_isotropy_obstruction,
    ),
    Correspondence(
        key="central-plane-screening",
        kind=DISAGREEMENT,
        pair=THREE_WAY,
        claim="No plane through the centre of the OPH carrier screens its strict halves in the seam Gaussian. All seven vertex-spanned central planes do so on the FCC shell. The carrier's fifteen vertex-spanned planes leave conditional mutual information (1/2) log(10/9) in the clamped unit intrinsic field.",
        oph="The carrier's seam dependencies cross every vertex-spanned equator even when its ports are clamped.",
        fep="A plane screens exactly when the conditional precision has no cross-block entries. An exact determinant ratio measures residual dependence.",
        synergetics="Both coordinate and body-diagonal equators of the FCC first shell remove all cross-plane seams.",
        status="obstructed",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}; Dempster (1972) for Gaussian precision zeros; plane and determinant censuses derived here",
        probe=_probe_central_plane_screening,
    ),
    Correspondence(
        key="repair-centroid-memory",
        kind=DISAGREEMENT,
        pair=THREE_WAY,
        claim="On both twelve-vertex shells, the stationary binary repair centroid has the same full lag covariance sequence as an isotropic sampled Ornstein–Uhlenbeck process, but its next-step law depends on the previous observed centroid even after conditioning on the current one.",
        oph="The native transition counts give exact conditional mean retention and three-time event probabilities.",
        fep="The matching linear Gaussian process reproduces equilibrium second-order statistics while missing memory in the observed process.",
        synergetics="The geometric coordinate fields supply the slow modes on the carrier and the FCC shell. Opposite loaded pairs collapse to the same centroid.",
        status="obstructed",
        source=f"{_CARRIER_SOURCE}; {_QUADRAY_SOURCE}; Uhlenbeck and Ornstein (1930) for the Gaussian process; exact centroid comparison derived here",
        probe=_probe_repair_centroid_memory,
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

    An unknown filter value or an empty intersection raises ``ValueError`` so
    a verification run must exercise at least one probe.
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
    if not selected:
        raise ValueError("no correspondence probes match the requested filters")
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
