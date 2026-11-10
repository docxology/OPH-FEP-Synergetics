"""The correspondence registry, the verification battery, and the provenance model.

``registry``
    The catalogue of isomorphisms, disagreements and complementarities, each
    with an executable probe and a source.
``provenance``
    Source pins, licence terms and reuse policy for OPH, FEP_Formal and QuadMath.
"""

from __future__ import annotations

from .provenance import (
    DEFAULT_SOURCES,
    SOURCES_ENV,
    CitedIdentifier,
    LicenceTerm,
    PinCheck,
    SourcePin,
    check_pins,
    head_of,
    load_sources,
)
from .registry import (
    COMPLEMENTARITY,
    CORRESPONDENCES,
    DISAGREEMENT,
    ISOMORPHISM,
    Correspondence,
    CorrespondenceResult,
    by_key,
    summary,
    verify,
)

__all__ = [
    "COMPLEMENTARITY",
    "CORRESPONDENCES",
    "DEFAULT_SOURCES",
    "DISAGREEMENT",
    "ISOMORPHISM",
    "SOURCES_ENV",
    "CitedIdentifier",
    "Correspondence",
    "CorrespondenceResult",
    "LicenceTerm",
    "PinCheck",
    "SourcePin",
    "by_key",
    "check_pins",
    "head_of",
    "load_sources",
    "summary",
    "verify",
]
