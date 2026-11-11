"""SC2: the icosahedral symmetry group of the carrier and what it fixes."""

from __future__ import annotations

import pytest

from oph_fep_synergetics.core import (
    automorphism_group_is_closed,
    carrier_automorphisms,
    descent_outcome_multiset,
    eigenspace_irreducibility,
    icosahedral_carrier,
    law_equivariance,
    orbit_sizes,
    reflection_subgroup,
    rotation_subgroup,
)
from oph_fep_synergetics.core import carrier as carrier_module


def test_carrier_has_one_hundred_and_twenty_automorphisms():
    group = carrier_automorphisms()
    assert len(group) == 120
    assert len(set(group)) == 120


def test_every_automorphism_preserves_adjacency():
    neighbours = carrier_module.OPH_NEIGHBOURS
    for permutation in carrier_automorphisms():
        for port in range(12):
            assert sorted(permutation[j] for j in neighbours[port]) == list(neighbours[permutation[port]])


def test_the_identity_is_an_automorphism():
    assert tuple(range(12)) in carrier_automorphisms()


def test_automorphism_search_rejects_a_truncated_limit():
    with pytest.raises(AssertionError):
        carrier_automorphisms(limit=5)


def test_automorphism_group_is_closed():
    assert automorphism_group_is_closed() is True


def test_rotations_and_reflections_halve_the_group():
    rotations = rotation_subgroup()
    reflections = reflection_subgroup()
    group = set(carrier_automorphisms())
    assert len(rotations) == 60
    assert len(reflections) == 60
    assert set(rotations).isdisjoint(reflections)
    assert set(rotations) | set(reflections) == group


def test_composition_stays_inside_the_rotation_subgroup():
    rotations = set(rotation_subgroup())
    reflections = set(reflection_subgroup())
    size = 12
    for left in rotation_subgroup():
        for right in rotation_subgroup():
            assert tuple(left[right[i]] for i in range(size)) in rotations
    for left in reflection_subgroup():
        for right in reflection_subgroup():
            assert tuple(left[right[i]] for i in range(size)) in rotations
    assert len(rotations) == 60
    assert len(reflections) == 60


def test_orbit_sizes_are_one_on_ports_seams_and_faces():
    assert orbit_sizes() == {"ports": 1, "seams": 1, "faces": 1}


def test_eigenspaces_are_irreducible_under_the_group():
    report = eigenspace_irreducibility()
    assert len(report) == 4
    assert sorted(entry["dim"] for entry in report.values()) == [1, 3, 3, 5]
    for entry in report.values():
        assert entry["irreducible"] is True
        assert entry["orbit_rank"] == entry["dim"]
        if "min_orbit_rank" in entry:
            assert entry["min_orbit_rank"] <= entry["orbit_rank"]


def test_the_native_law_is_equivariant():
    assert law_equivariance() == {
        "invariant": True,
        "configurations": 78,
        "group_order": 120,
        "orbits": 4,
        "distinct_outcome_multisets": 3,
        "witness": None,
    }


def test_descent_outcome_multiset_covers_sixty_transitions():
    carrier = icosahedral_carrier()
    for configuration in carrier.shell_configurations(total=2):
        outcomes = descent_outcome_multiset(configuration)
        assert set(outcomes) == {"strict", "zero", "noop"}
        assert all(isinstance(value, int) and value >= 0 for value in outcomes.values())
        assert sum(outcomes.values()) == 60


def test_the_census_totals_seventy_eight_times_sixty():
    carrier = icosahedral_carrier()
    totals = {"strict": 0, "zero": 0, "noop": 0}
    for configuration in carrier.shell_configurations(total=2):
        for key, value in descent_outcome_multiset(configuration).items():
            totals[key] += value
    assert sum(totals.values()) == 78 * 60 == 4680
