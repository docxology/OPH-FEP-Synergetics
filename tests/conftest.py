"""Shared fixtures for the correspondence test suite.

The fixtures are session-scoped because the carrier, the sheaf models and the
exhaustive shell enumerations are deterministic: rebuilding them per test would
only spend time, never change a result.
"""

from __future__ import annotations

import pytest

from oph_fep_synergetics.core import GaussianSheafModel, demo_chain_model, icosahedral_carrier


@pytest.fixture(scope="session")
def carrier():
    """OPH's labelled twelve-port carrier."""
    return icosahedral_carrier()


@pytest.fixture(scope="session")
def chain_model() -> GaussianSheafModel:
    """A four-node scalar chain with both endpoints observed."""
    return demo_chain_model(4)
