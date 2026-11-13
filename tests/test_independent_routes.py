"""Independent routes to claims that were previously checked against themselves.

Each test exercises a second computation that shares no code path with the
first: the surface cohomology beside the graph cohomology, the closed-form shell
gauge beside the breadth-first shells, an exact case analysis beside a sampled
orbit, a KKT solve beside an interior solve, energy minus entropy beside
surprisal plus divergence, and transition-level equivariance beside outcome
counts.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

from oph_fep_synergetics.core import controls, free_energy, invariants, sheaf, symmetry
from oph_fep_synergetics.core.carrier import OPH_ORIENTED_BASE_FACES, icosahedral_carrier
from oph_fep_synergetics.geometry import cuboctahedral_gauge, generate_shell, squared_norm


def test_the_carrier_surface_is_a_sphere():
    report = sheaf.surface_cohomology(icosahedral_carrier(), OPH_ORIENTED_BASE_FACES)
    assert report["cells"] == (12, 30, 20)
    assert report["delta0_rank"] == 11
    assert report["delta1_rank"] == 19
    assert report["delta1_after_delta0_vanishes"] is True
    assert report["betti"] == (1, 0, 1)
    assert report["euler_characteristic"] == report["cell_euler_characteristic"] == 2
    assert report["faces_span_graph_cycles"] is True
    assert report["oriented_faces_form_a_cycle"] is True


def test_unoriented_faces_are_not_a_cycle_but_have_the_same_cohomology():
    report = sheaf.surface_cohomology(icosahedral_carrier())
    assert report["betti"] == (1, 0, 1)
    assert report["oriented_faces_form_a_cycle"] is False


def test_sheaf_euler_characteristic_uses_stalk_dimensions():
    edges = tuple(
        sheaf.EdgeRestriction(u=i, v=i + 1, map_u=np.eye(2), map_v=np.eye(2)) for i in range(2)
    )
    chain = sheaf.CellularSheaf(stalk_dims=(2, 2, 2), edges=edges)
    report = sheaf.cohomology_dimensions(chain)
    assert report["h0"] == 2
    assert report["h1"] == 0
    assert report["euler_characteristic"] == report["h0_minus_h1"] == 2


def test_exact_rank_handles_an_empty_and_a_dependent_matrix():
    assert sheaf._exact_rank([]) == 0
    assert sheaf._exact_rank([[1, 2], [2, 4], [0, 0]]) == 1


@pytest.mark.parametrize("k", range(0, 9))
def test_the_cuboctahedral_gauge_is_the_shell_index(k):
    assert all(cuboctahedral_gauge(site) == k for site in generate_shell(k))


def test_the_free_energy_range_on_a_shell():
    for k in range(1, 9):
        energies = [squared_norm(site) // 2 for site in generate_shell(k)]
        assert max(energies) == 4 * k * k
        assert sum(1 for value in energies if value == 4 * k * k) == 12
        assert min(energies) == 2 * k * k + (2 if k % 2 else 0)


def test_the_gauge_rejects_a_non_lattice_point():
    with pytest.raises(ValueError):
        cuboctahedral_gauge((1, 0, 0, 0))


def test_the_clipped_map_has_the_origin_as_its_only_fixed_point():
    assert controls.two_node_fixed_points() == ((0, 0),)
    assert np.allclose(controls.two_node_map((0.0, 0.0)), (0.0, 0.0))


def test_the_kkt_minimiser_is_the_harmonic_extension():
    model = free_energy.demo_chain_model(4)
    for records in ((0.3, -0.2), (1.0, 2.0), (-0.7, 0.0)):
        assert np.allclose(
            model.constrained_minimiser(records), model.harmonic_extension(records), atol=1e-12
        )


def test_energy_minus_entropy_is_the_variational_free_energy():
    model = free_energy.demo_chain_model(4)
    records = (0.3, -0.2)
    mean, covariance = model.posterior(records)
    for shift, scale in ((0.0, 1.0), (0.2, 1.0), (0.0, 2.0), (-0.3, 0.7)):
        q_mean = mean + shift
        q_covariance = scale * covariance
        assert model.energy_minus_entropy(q_mean, q_covariance, records) == pytest.approx(
            model.variational_free_energy(q_mean, q_covariance, records), abs=1e-12
        )


def test_energy_minus_entropy_rejects_an_improper_covariance():
    model = free_energy.demo_chain_model(4)
    mean, covariance = model.posterior((0.3, -0.2))
    with pytest.raises(ValueError):
        model.energy_minus_entropy(mean, covariance @ np.diag([1.0, -1.0]), (0.3, -0.2))


def test_every_transition_commutes_with_every_automorphism():
    report = symmetry.transition_equivariance()
    assert report["equivariant"] is True
    assert report["checks"] == 120 * 78 * 60 * 2
    assert report["witness"] is None


def test_the_replicator_potential_needs_a_symmetric_payoff():
    report = controls.replicator_asymmetric_control()
    assert report["symmetric_points"] == 25
    assert report["symmetric_all_gradient"] is True
    assert report["asymmetric_not_gradient"] is True


def test_the_reciprocal_spectrum_sum_needs_no_spectrum():
    assert invariants.reciprocal_spectrum_sum() == Fraction(7, 3)
    assert invariants.marginal_variance() == Fraction(7, 36)
    pinv = np.linalg.pinv(icosahedral_carrier().laplacian())
    assert np.allclose(np.diag(pinv), 7 / 36)


def test_the_exact_inverse_inverts():
    matrix = [[Fraction(2), Fraction(1)], [Fraction(1), Fraction(1)]]
    assert invariants._exact_inverse(matrix) == [[1, -1], [-1, 2]]
