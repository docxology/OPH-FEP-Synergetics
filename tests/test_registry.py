"""SC5: one registry executes every correspondence and every disagreement."""

from __future__ import annotations

from fractions import Fraction

import pytest

from oph_fep_synergetics.analysis import (
    COMPLEMENTARITY,
    CORRESPONDENCES,
    DISAGREEMENT,
    FEP_SYNERGETICS,
    ISOMORPHISM,
    KINDS as REGISTRY_KINDS,
    OPH_FEP,
    PAIRS,
    SYNERGETICS_OPH,
    THREE_WAY,
    by_key,
    by_kind,
    by_pair,
    summary,
    verify,
)

STATUS_VOCABULARY = {"exact", "exact-on-domain", "obstructed", "not-reproduced"}
KINDS = (ISOMORPHISM, DISAGREEMENT, COMPLEMENTARITY)
assert KINDS == REGISTRY_KINDS, "the registry and its tests must agree on the kind vocabulary"


def test_registry_entries_are_well_formed():
    assert len(CORRESPONDENCES) == 42
    assert len({entry.key for entry in CORRESPONDENCES}) == 42
    for entry in CORRESPONDENCES:
        assert entry.kind in KINDS
        assert entry.pair in PAIRS
        assert entry.status in STATUS_VOCABULARY
        assert entry.claim and entry.source
        assert entry.oph and entry.fep and entry.synergetics
        assert callable(entry.probe)


def test_registry_covers_all_three_relations():
    counts = {kind: sum(1 for entry in CORRESPONDENCES if entry.kind == kind) for kind in KINDS}
    assert counts == {ISOMORPHISM: 17, DISAGREEMENT: 15, COMPLEMENTARITY: 10}


def test_registry_covers_every_pairwise_scope():
    assert PAIRS == (OPH_FEP, FEP_SYNERGETICS, SYNERGETICS_OPH, THREE_WAY)
    index = by_pair()
    assert set(index) == set(PAIRS)
    counts = {pair: len(entries) for pair, entries in index.items()}
    assert counts == {OPH_FEP: 24, FEP_SYNERGETICS: 2, SYNERGETICS_OPH: 14, THREE_WAY: 2}
    assert sum(counts.values()) == len(CORRESPONDENCES)
    for pair, entries in index.items():
        assert all(entry.pair == pair for entry in entries)


def test_by_kind_indexes_every_entry():
    index = by_kind()
    assert set(index) == set(KINDS)
    assert sum(len(entries) for entries in index.values()) == len(CORRESPONDENCES)
    for kind, entries in index.items():
        assert all(entry.kind == kind for entry in entries)


def test_by_key_indexes_every_entry():
    index = by_key()
    assert set(index) == {entry.key for entry in CORRESPONDENCES}
    assert index["carrier-port-chart"].kind == ISOMORPHISM
    assert index["pr01-obstruction"].kind == DISAGREEMENT
    assert index["quadray-ivm"].kind == COMPLEMENTARITY


def test_every_probe_passes_and_reports_values():
    results = verify()
    assert len(results) == 42
    for result in results:
        assert result.passed is True
        assert result.pair in PAIRS
        assert result.values
    report = summary(results)
    assert report["total"] == 42
    assert report["passed"] == 42
    assert report["failed"] == 0
    assert report["failures"] == ()
    assert report["by_kind"][ISOMORPHISM] == {"passed": 17, "failed": 0}
    assert report["by_kind"][DISAGREEMENT] == {"passed": 15, "failed": 0}
    assert report["by_kind"][COMPLEMENTARITY] == {"passed": 10, "failed": 0}
    assert report["pairs"] == PAIRS
    assert report["by_pair"][OPH_FEP] == {"passed": 24, "failed": 0}
    assert report["by_pair"][FEP_SYNERGETICS] == {"passed": 2, "failed": 0}
    assert report["by_pair"][SYNERGETICS_OPH] == {"passed": 14, "failed": 0}
    assert report["by_pair"][THREE_WAY] == {"passed": 2, "failed": 0}


def test_selection_limits_the_run():
    results = verify(keys=("carrier-port-chart", "pr01-obstruction"))
    assert [result.key for result in results] == ["carrier-port-chart", "pr01-obstruction"]


def test_pair_selection_limits_the_run_to_one_pairwise_scope():
    results = verify(pairs=(THREE_WAY,))
    assert [result.key for result in results] == ["tetravolume-arithmetic", "conjugation-norm"]
    for pair in PAIRS:
        selected = verify(pairs=(pair,))
        assert {result.key for result in selected} == {entry.key for entry in by_pair()[pair]}


def test_kind_selection_limits_the_run_to_one_kind():
    for kind in KINDS:
        selected = verify(kinds=(kind,))
        assert {result.key for result in selected} == {entry.key for entry in by_kind()[kind]}
        assert all(result.kind == kind for result in selected)
    both = verify(pairs=(SYNERGETICS_OPH,), kinds=(DISAGREEMENT,))
    assert {result.key for result in both} <= {entry.key for entry in by_pair()[SYNERGETICS_OPH]}
    assert all(result.kind == DISAGREEMENT for result in both)
    assert both, "the synergetics-oph scope has disagreements"


def test_unknown_kinds_are_rejected_instead_of_selecting_nothing():
    with pytest.raises(ValueError):
        verify(kinds=("no-such-kind",))


def test_unknown_pairs_are_rejected_instead_of_selecting_nothing():
    with pytest.raises(ValueError):
        verify(pairs=("no-such-pair",))


def test_unknown_keys_are_rejected_instead_of_selecting_nothing():
    with pytest.raises(ValueError):
        verify(keys=("no-such-correspondence",))


def test_probe_values_are_reportable_scalars():
    for entry in CORRESPONDENCES:
        for name, value in entry.run().values.items():
            assert isinstance(name, str)
            reportable = (
                isinstance(value, (str, int, float, bool, list, tuple, dict, Fraction))
                or value is None
                or hasattr(value, "tolist")
            )
            assert reportable, f"{entry.key}.{name} is not reportable"


def test_carrier_port_chart_values():
    values = by_key()["carrier-port-chart"].run().values
    assert (values["ports"], values["seams"], values["faces"]) == (12, 30, 20)
    assert set(values["degrees"]) == {5}
    assert values["derived_equals_declared"] is True
    assert values["faces_equal_graph_triangles"] is True


def test_spectral_gap_values():
    values = by_key()["spectral-gap"].run().values
    exact = 5.0 - 5.0**0.5
    assert values["expected_lambda_1"] == exact
    assert values["lambda_1"] == pytest.approx(exact, abs=1e-12)
    assert values["shell_gap"] == pytest.approx(exact / 60.0, abs=1e-15)
    assert values["laplacian_is_5I_minus_A"] is True


def test_shell_count_values():
    values = by_key()["shell-counts"].run().values
    assert values["total_two"] == 78
    assert values["binary"] == 66
    assert values["repair_module_total_two"] == 78
    assert values["repair_module_binary"] == 66


def test_native_repair_law_values():
    values = by_key()["native-repair-law"].run().values
    assert values["states"] == 78
    assert values["transitions"] == 4680
    assert "d^2" in values["closed_form"]


def test_variational_equivalence_values():
    values = by_key()["variational-equivalence"].run().values
    assert values["perceive_matches_harmonic_extension"] is True
    assert values["posterior_mean_is_harmonic_extension"] is True
    assert values["residual"] < 1e-9
    assert values["stationary_gradient_norm"] < 1e-9
    assert float(values["excess_free_energy_at_posterior"]) == 0.0
    assert values["min_perturbed_energy"] > values["energy_at_minimiser"]


def test_free_energy_dictionary_values():
    values = by_key()["free-energy-dictionary"].run().values
    assert values["record_injectivity"] is True
    assert values["bayesian_surprise"] == values["closed_form_half_g_linv_g"]
    assert values["sandwich_inside"] is True
    assert values["sandwich_lower"] <= values["bayesian_surprise"] <= values["sandwich_upper"]
    assert values["variational_free_energy"] == values["surprisal"]
    assert values["numerical_surprise"] > 0.0


def test_sheaf_agreement_values():
    values = by_key()["sheaf-agreement"].run().values
    assert (values["nodes"], values["edges"]) == (12, 30)
    assert values["laplacian_equals_carrier_laplacian"] is True
    assert values["global_section_dimension"] == 1
    assert values["chain_global_section_dimension"] == 1
    assert values["injectivity_all_ports"] is True
    assert values["chain_injectivity_one_record"] is True
    assert values["chain_injectivity_no_record"] is False


def test_certificate_species_values():
    values = by_key()["certificate-species"].run().values
    assert values["is_symmetric"] is True
    assert values["is_nonnegative"] is True
    assert list(values["row_sums"]) == [1.0]
    assert values["second_eigenvalue"] == pytest.approx(
        values["expected_second_eigenvalue"], abs=1e-12
    )
    assert values["induced_infinity_norm"] == 1.0


def test_pr01_obstruction_values():
    values = by_key()["pr01-obstruction"].run().values
    assert values["states"] == 78
    assert values["total_transitions"] == 4680
    assert (values["strict"], values["zero"], values["noop"]) == (120, 600, 3960)
    assert values["strict"] + values["zero"] + values["noop"] == 4680
    assert values["excursion_energy"] == 0.5
    assert values["excursion_states"] == values["excursion_distinct_states"] == 7
    assert values["coin_zero_current"] != values["coin_one_current"]
    assert values["coin_zero_decrease"] == values["coin_one_decrease"] == 0.0


def test_free_energy_name_collision_values():
    values = by_key()["free-energy-name-collision"].run().values
    assert values["is_positive_definite_interior"] is True
    assert values["interface_energy_minimised_over_interior"] > 0.0
    assert float(values["difference"]) == pytest.approx(
        values["expected_model_constant"], abs=1e-12
    )
    assert float(values["variational_free_energy_at_posterior"]) != float(values["difference"])


def test_mediator_obstruction_values():
    values = by_key()["mediator-obstruction"].run().values
    assert values["pairwise_representation_exists_iff_epsilon_ge_one_third"] is True
    assert values["margin_at_threshold"] == 0.0
    assert values["margin_below_threshold"] < 0.0
    assert values["decomposition_below_is_none"] is True
    assert values["decomposition_at_threshold_exists"] is True
    assert values["decomposition_is_positive_semidefinite"] is True
    assert values["decomposition_matches_target"] is True


def test_stationary_blanket_values():
    values = by_key()["stationary-blanket"].run().values
    assert values["conditional_independence"] is True
    assert tuple(values["expected_mu_drift"]) == (1.5, -0.25, -2.0)
    assert (
        values["mu_drift_external"],
        values["mu_drift_blanket"],
        values["mu_drift_internal"],
    ) == tuple(values["expected_mu_drift"])
    assert values["cross_drift_difference"] > 0.0


def test_solenoidal_freedom_values():
    values = by_key()["solenoidal-freedom"].run().values
    assert values["stationary_residuals"] == [0.0, 0.0, 0.0]
    assert set(values["spectral_real"]) == {-1.0}
    assert values["entropy_production"] == values["expected_entropy_production"]
    assert values["entropy_production"][0] == 0.0


def test_d_stability_values():
    values = by_key()["d-stability"].run().values
    assert values["routh_first_column_at_k4"] == (1.0, 6.0, -2.0, 24.0)
    assert values["charpoly_matches_closed_form"] is True
    assert values["P_matches_source"] is True
    assert values["Q_matches_source"] is True
    assert values["Q_is_antisymmetric"] is True
    assert values["scaled_is_unstable"] is True
    assert values["determinant"] == pytest.approx(115.0 / 27.0, abs=1e-12)
    assert values["lyapunov_residual"] < 1e-12
    assert values["logistic_embedding"] == {"f_zero_norm": 0.0, "jacobian_error": 0.0}


def test_schedule_path_values():
    values = by_key()["schedule-path"].run().values
    assert values["matches_source_rational"] is True
    assert Fraction(values["contrast_difference_exact"]) == Fraction(921633999, 40000000000)
    assert values["contrast_difference_float"] == pytest.approx(0.023040849975, abs=1e-18)
    assert values["full_iteration_steps"] == 12


def test_reduction_criterion_values():
    values = by_key()["reduction-criterion"].run().values
    assert values["retaining_the_record_is_exact"] is True
    assert values["same_observation"] == 1.0
    assert (values["first_after"], values["second_after"]) == (1.0, -1.0)
    assert values["response_distance"] == 2.0
    assert values["worst_case_error_lower_bound"] == 1.0


def test_existence_without_reachability_values():
    values = by_key()["existence-without-reachability"].run().values
    assert values["sync_map_min_radius"] >= 0.99
    assert values["relaxation_min_radius"] > values["source_relaxation_min_radius"] > 0.8
    assert abs(values["relaxation_period"] - values["source_relaxation_period"]) < 0.1
    assert values["relaxation_max_radius"] > values["relaxation_min_radius"]


def test_directed_gradient_values():
    values = by_key()["directed-gradient"].run().values
    assert tuple(values["direct"]) == (1.0, 0.0)
    assert tuple(values["adjoint"]) == (1.0, 2.0)
    assert values["direct"] == values["source_direct"]
    assert values["adjoint"] == values["source_adjoint"]


def test_four_cycle_values():
    values = by_key()["four-cycle-local-minimum"].run().values
    assert values["is_local_minimum"] is True
    assert values["trap_violations"] == 2
    assert set(values["flip_violations"]) == {2}
    assert set(values["compatible_violations"]) == {0}


def test_quadray_ivm_values():
    values = by_key()["quadray-ivm"].run().values
    assert values["neighbor_moves"] == values["distinct_moves"] == 12
    assert values["neighbor_squared_distance"] == 8
    assert values["min_pairwise_squared_distance"] == 8.0
    assert set(values["squared_norms"]) == {8.0}
    assert tuple(values["shell_counts"]) == (1, 12, 42, 92, 162)
    assert tuple(values["cumulative_counts"]) == (1, 13, 55, 147, 309)
    assert set(values["volume_ratios_present"]) == {1.0, 3.0, 4.0, 6.0, 20.0, 96.0}
    assert values["max_shell"] == 32


def test_equilibrium_propagation_values():
    values = by_key()["equilibrium-propagation"].run().values
    assert values["exact_gradient"] == -1.25
    assert values["closed_form_error"] < 1e-9
    assert values["distance_to_gradient"] == pytest.approx(
        values["expected_distance_to_gradient"], abs=1e-12
    )


def test_replicator_potential_values():
    values = by_key()["replicator-potential"].run().values
    assert values["field_equals_mobility_gradient"] is True
    assert values["rates_agree"] is True
    assert values["field_norm"] > 0.0


def test_kissing_twelve_values():
    values = by_key()["kissing-twelve"].run().values
    assert values["carrier_ports"] == 12
    assert set(values["carrier_degrees"]) == {5}
    assert values["ivm_moves"] == values["ivm_shell_one"] == 12
    assert values["move_squared_norm"] == 8.0
    assert values["geometric_adjacency_isomorphic"] is True
