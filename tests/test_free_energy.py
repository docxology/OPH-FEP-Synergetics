"""The free-energy dictionary on a linear-Gaussian cellular sheaf."""

from __future__ import annotations

import math
from functools import cache
from itertools import product

import numpy as np
import pytest

from oph_fep_synergetics.core import (
    CellularSheaf,
    GaussianSheafModel,
    SurpriseTriple,
    carrier_sheaf,
    demo_chain_model,
    gaussian_kl,
    icosahedral_carrier,
    linear_chain_sheaf,
    three_surprises,
)
from oph_fep_synergetics.core.free_energy import carrier_grounded_forest_census


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
    with pytest.raises(ValueError, match="means must be vectors"):
        gaussian_kl(0.0, np.eye(1), 0.0, np.eye(1))
    with pytest.raises(ValueError, match="second mean must be finite"):
        gaussian_kl([0.0], np.eye(1), [math.nan], np.eye(1))
    with pytest.raises(ValueError, match="finite and symmetric"):
        gaussian_kl([0.0, 0.0], [[1.0, 0.5], [0.0, 1.0]], [0.0, 0.0], np.eye(2))
    with pytest.raises(ValueError, match="positive definite"):
        gaussian_kl([0.0, 0.0], np.eye(2), [0.0, 0.0], -np.eye(2))


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


def test_grounded_carrier_determinants_count_rooted_forests():
    census = carrier_grounded_forest_census()
    pair_layouts = census["pair_layouts"]
    triple_layouts = census["triple_layouts"]
    assert tuple(row[:3] for row in pair_layouts) == (
        ((1,), 30, 1900800),
        ((2,), 30, 2419200),
        ((3,), 6, 2592000),
    )
    assert tuple(row[:3] for row in triple_layouts) == (
        ((1, 1, 1), 20, 522720),
        ((1, 1, 2), 60, 604800),
        ((1, 2, 2), 60, 712800),
        ((1, 2, 3), 60, 743040),
        ((2, 2, 2), 20, 846720),
    )
    assert census["matches_contrast_covariance"] is True
    assert census["best_pair_example"] == (0, 3)
    assert census["best_triple_example"] == (0, 2, 8)
    assert census["worst_triple_example"] == (0, 1, 5)
    assert census["greedy_triple_example"] == (0, 1, 2)
    assert census["greedy_triple_determinants"] == ((743040, 60),)
    assert census["optimum_to_greedy_ratio"] == "49/43"

    carrier = icosahedral_carrier()
    sheaf = carrier_sheaf(carrier)
    neighbour_masks = tuple(
        sum(1 << other for other in carrier.neighbours(port)) for port in range(carrier.port_count)
    )

    @cache
    def tree_count(vertices: int) -> int:
        root = vertices & -vertices
        return branch_partitions(root.bit_length() - 1, vertices ^ root)

    @cache
    def branch_partitions(root: int, vertices: int) -> int:
        if not vertices:
            return 1
        # The first vertex determines which branch of the root to choose next.
        first = vertices & -vertices
        optional = vertices ^ first
        subset = optional
        count = 0
        while True:
            branch = first | subset
            attachments = (neighbour_masks[root] & branch).bit_count()
            if attachments:
                count += attachments * tree_count(branch) * branch_partitions(root, vertices ^ branch)
            if subset == 0:
                return count
            subset = (subset - 1) & optional

    # Independently count forests for one representative of each geometry, by
    # assigning every interior port to a root's component and adjoining trees.
    for _, _, determinant, records in pair_layouts + triple_layouts:
        model = GaussianSheafModel(sheaf, records, record_prior_precision=1.0)
        free_bits = tuple(1 << vertex for vertex in model.interior_nodes)
        count = 0
        for assignment in product(range(len(records)), repeat=len(free_bits)):
            components = [1 << root for root in records]
            for bit, slot in zip(free_bits, assignment, strict=True):
                components[slot] |= bit
            count += math.prod(tree_count(component) for component in components)
        assert count == determinant
        _, covariance = model.posterior([0.0] * len(records))
        sign, log_determinant = np.linalg.slogdet(covariance)
        assert sign == 1
        assert log_determinant == pytest.approx(-math.log(determinant))


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


def test_proper_record_prior_normalises_the_record_density():
    model = demo_chain_model(4, record_prior_precision=0.5)
    records = (0.3, -0.2)
    # Common and contrast coordinates have precision 1/2 and 7/6.
    normaliser = math.log(2.0 * math.pi) - 0.5 * math.log(7.0 / 12.0)
    common, contrast = [(records[0] + sign * records[1]) / math.sqrt(2) for sign in (1, -1)]
    expected = normaliser + 0.5 * (0.5 * common**2 + (7.0 / 6.0) * contrast**2)
    assert model.record_surprisal(records) == pytest.approx(expected, abs=1e-12)
    mean, covariance = model.posterior(records)
    assert model.variational_free_energy(mean, covariance, records) == pytest.approx(expected)
    assert model.energy_minus_entropy(mean, covariance, records) == pytest.approx(expected)
    with pytest.raises(ValueError, match="configured record prior"):
        model.surprisal(records, prior_quadratic=1.0)


def test_record_density_requires_a_proper_joint_precision(chain_model):
    with pytest.raises(ValueError, match="requires a record prior"):
        chain_model.record_surprisal((0.0, 0.0))
    for precision in (0.0, -0.5, math.inf, math.nan):
        with pytest.raises(ValueError, match="finite and positive"):
            demo_chain_model(4, record_prior_precision=precision)
    sheaf = CellularSheaf((1, 1), ())
    with pytest.raises(ValueError, match="does not ground"):
        GaussianSheafModel(sheaf, (0,), record_prior_precision=0.5)
    # A fully recorded model has zero interior dimension and a proper prior.
    model = GaussianSheafModel(sheaf, (0, 1), record_prior_precision=1.0)
    assert model.energy_minus_entropy([], np.zeros((0, 0)), [1.0, 2.0]) == pytest.approx(
        math.log(2.0 * math.pi) + 2.5
    )


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


def test_both_free_energy_routes_require_a_positive_definite_covariance(chain_model):
    for evaluate in (chain_model.variational_free_energy, chain_model.energy_minus_entropy):
        with pytest.raises(ValueError, match="positive definite"):
            evaluate([0.0, 0.0], -np.eye(2), [1.0, 2.0])


def test_perceive_returns_the_harmonic_extension(chain_model):
    records = [1.0, 3.0]
    interior, residual, iterations = chain_model.perceive(records)
    assert np.allclose(interior, chain_model.harmonic_extension(records), atol=1e-8)
    assert residual <= 1e-10
    assert iterations > 0


def test_perceive_does_not_require_a_direct_solve(chain_model, monkeypatch):
    def fail_solve(*args):
        pytest.fail("iterative perception called the direct solve")

    monkeypatch.setattr(np.linalg, "solve", fail_solve)
    interior, residual, _ = chain_model.perceive([1.0, 3.0])
    assert np.allclose(interior, [5 / 3, 7 / 3])
    assert residual <= 1e-10


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
