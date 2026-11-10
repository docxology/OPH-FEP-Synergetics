"""The free-energy dictionary on a linear-Gaussian cellular sheaf."""

from __future__ import annotations

import math

import numpy as np
import pytest

from oph_fep_synergetics.core import (
    GaussianSheafModel,
    SurpriseTriple,
    demo_chain_model,
    gaussian_kl,
    linear_chain_sheaf,
    three_surprises,
)


def test_gaussian_kl_is_zero_for_identical_gaussians():
    mean = np.array([1.0, -2.0])
    covariance = np.array([[2.0, 0.5], [0.5, 1.0]])
    assert abs(gaussian_kl(mean, covariance, mean, covariance)) < 1e-12


def test_gaussian_kl_matches_the_scalar_closed_form():
    for mean in (0.25, 1.0, -3.0):
        assert abs(gaussian_kl([0.0], [[1.0]], [mean], [[1.0]]) - 0.5 * mean * mean) < 1e-12
    ratio = 2.0
    expected = 0.5 * (ratio - 1.0 - math.log(ratio))
    assert abs(gaussian_kl([0.0], [[2.0]], [0.0], [[1.0]]) - expected) < 1e-12


def test_gaussian_kl_is_nonnegative_and_grows_with_distance():
    covariance = np.array([[1.0]])
    near = gaussian_kl([0.0], covariance, [0.5], covariance)
    far = gaussian_kl([0.0], covariance, [3.0], covariance)
    assert 0.0 < near < far


def test_gaussian_kl_rejects_bad_shapes_and_signs():
    with pytest.raises(ValueError):
        gaussian_kl([0.0], np.eye(1), [0.0, 1.0], np.eye(2))
    with pytest.raises(ValueError):
        gaussian_kl([0.0], np.eye(2), [0.0], np.eye(1))
    with pytest.raises(ValueError):
        gaussian_kl([0.0], np.array([[-1.0]]), [0.0], np.eye(1))


def test_demo_chain_needs_at_least_three_nodes():
    with pytest.raises(ValueError):
        demo_chain_model(2)


def test_demo_chain_shape(chain_model):
    assert chain_model.record_nodes == (0, 3)
    assert chain_model.interior_nodes == (1, 2)
    assert chain_model.record_dim == 2
    assert chain_model.interior_dim == 2


def test_blocks_are_copies_of_the_sheaf_laplacian_blocks(chain_model):
    blocks = chain_model.blocks
    assert set(blocks) == {"L_xx", "L_xo", "L_ox", "L_oo"}
    assert np.allclose(blocks["L_xx"], [[2.0, -1.0], [-1.0, 2.0]])
    assert np.allclose(blocks["L_xo"], [[-1.0, 0.0], [0.0, -1.0]])
    assert np.allclose(blocks["L_ox"], blocks["L_xo"].T)
    assert np.allclose(blocks["L_oo"], np.eye(2))
    blocks["L_xx"][0, 0] = 99.0
    assert chain_model.blocks["L_xx"][0, 0] == 2.0


def test_harmonic_extension_is_the_conditional_mean(chain_model):
    assert np.allclose(chain_model.harmonic_extension([3.0, 6.0]), [4.0, 5.0])
    assert np.allclose(chain_model.harmonic_extension([0.0, 0.0]), [0.0, 0.0])


def test_harmonic_extension_agrees_with_the_sheaf(chain_model):
    nodes = chain_model.sheaf.harmonic_extension({0: [0.5], 3: [-1.5]})
    assert np.allclose(
        chain_model.harmonic_extension([0.5, -1.5]),
        [float(nodes[1][0]), float(nodes[2][0])],
    )


def test_schur_complement_is_the_symmetric_marginal_precision(chain_model):
    marginal = chain_model.schur_complement()
    assert np.allclose(marginal, (1.0 / 3.0) * np.array([[1.0, -1.0], [-1.0, 1.0]]))
    assert np.allclose(marginal, marginal.T)


def test_posterior_matches_the_interior_block_inverse(chain_model):
    mean, covariance = chain_model.posterior([1.0, 2.0])
    assert np.allclose(mean, chain_model.harmonic_extension([1.0, 2.0]))
    assert np.allclose(covariance, (1.0 / 3.0) * np.array([[2.0, 1.0], [1.0, 2.0]]))


def test_surprisal_differences_are_exact(chain_model):
    marginal = chain_model.schur_complement()
    first, second = [1.0, 2.0], [-3.0, 0.5]
    quadratic = 0.5 * float(first @ marginal @ first) - 0.5 * float(second @ marginal @ second)
    difference = chain_model.surprisal(first) - chain_model.surprisal(second)
    assert abs(difference - quadratic) < 1e-12


def test_prior_quadratic_shifts_the_surprisal_exactly(chain_model):
    records = [0.25, 0.75]
    assert abs(
        chain_model.surprisal(records, prior_quadratic=2.5)
        - chain_model.surprisal(records)
        - 2.5
    ) < 1e-12


def test_variational_free_energy_is_minimised_at_the_posterior(chain_model):
    records = [2.0, -1.0]
    mean, covariance = chain_model.posterior(records)
    surprisal = chain_model.surprisal(records)
    assert abs(chain_model.variational_free_energy(mean, covariance, records) - surprisal) < 1e-12
    assert chain_model.excess_free_energy(mean, covariance, records) < 1e-12
    assert chain_model.excess_free_energy([5.0, 5.0], covariance, records) > 0.0
    assert chain_model.variational_free_energy([0.0, 0.0], covariance, records) > surprisal


def test_excess_free_energy_rejects_a_bad_covariance(chain_model):
    with pytest.raises(ValueError):
        chain_model.excess_free_energy([0.0, 0.0], np.eye(3), [1.0, 2.0])


def test_perceive_returns_the_harmonic_extension(chain_model):
    records = [1.0, 3.0]
    interior, residual, iterations = chain_model.perceive(records)
    assert np.allclose(interior, chain_model.harmonic_extension(records), atol=1e-8)
    assert residual <= 1e-10
    assert iterations > 0


def test_perceive_accepts_an_initial_state(chain_model):
    records = [1.0, 3.0]
    interior, _, iterations = chain_model.perceive(records, initial=[2.0, 2.0])
    assert np.allclose(interior, chain_model.harmonic_extension(records), atol=1e-8)
    assert iterations > 0


def test_perceive_rejects_a_bad_initial_state(chain_model):
    with pytest.raises(ValueError):
        chain_model.perceive([1.0, 3.0], initial=[1.0])


def test_surprises_vanish_for_unchanged_records(chain_model):
    records = [1.0, 2.0]
    assert chain_model.numerical_surprise(records, records) < 1e-12
    assert chain_model.bayesian_surprise(records, records) < 1e-12
    assert chain_model.numerical_surprise([0.0, 0.0], [1.0, 0.0]) > 0.0


def test_bayesian_surprise_is_the_residual_quadratic_form(chain_model):
    gradient = chain_model.blocks["L_xo"] @ np.array([1.0, 2.0])
    expected = 0.5 * float(
        gradient @ np.linalg.solve(chain_model.blocks["L_xx"], gradient)
    )
    assert abs(chain_model.bayesian_surprise([0.0, 0.0], [1.0, 2.0]) - expected) < 1e-12


def test_three_surprises_keeps_the_three_names_apart(chain_model):
    old, new = [0.0, 0.0], [1.0, 2.0]
    triple = three_surprises(chain_model, old, new)
    assert isinstance(triple, SurpriseTriple)
    assert triple.numerical == chain_model.numerical_surprise(old, new)
    assert triple.bayesian == chain_model.bayesian_surprise(old, new)
    assert triple.surprisal == chain_model.surprisal(new)
    assert triple.numerical > 0.0
    assert triple.bayesian > 0.0


def test_record_injectivity_delegates_to_the_sheaf(chain_model):
    assert chain_model.record_injectivity() is True


def test_model_rejects_bad_record_nodes():
    sheaf = linear_chain_sheaf(3)
    with pytest.raises(ValueError):
        GaussianSheafModel(sheaf, (0, 0))
    with pytest.raises(IndexError):
        GaussianSheafModel(sheaf, (0, 9))


def test_a_fully_recorded_model_has_no_interior():
    model = GaussianSheafModel(linear_chain_sheaf(2), (0, 1))
    assert model.interior_dim == 0
    assert model.harmonic_extension([1.0, 2.0]).shape == (0,)
    assert model.schur_complement().shape == (2, 2)
    mean, covariance = model.posterior([1.0, 2.0])
    assert mean.shape == (0,)
    assert covariance.shape == (0, 0)
    assert model.bayesian_surprise([1.0, 2.0], [0.0, 0.0]) == 0.0
    assert model.excess_free_energy([], np.zeros((0, 0)), [1.0, 2.0]) == 0.0
    interior, residual, iterations = model.perceive([1.0, 2.0])
    assert (interior.shape, residual, iterations) == ((0,), 0.0, 0)
