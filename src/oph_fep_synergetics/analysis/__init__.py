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
    FEP_SYNERGETICS,
    ISOMORPHISM,
    OPH_FEP,
    PAIRS,
    SYNERGETICS_OPH,
    THREE_WAY,
    Correspondence,
    CorrespondenceResult,
    by_key,
    by_kind,
    by_pair,
    summary,
    verify,
)

__all__ = [
    "COMPLEMENTARITY",
    "CORRESPONDENCES",
    "DEFAULT_SOURCES",
    "DISAGREEMENT",
    "FEP_SYNERGETICS",
    "ISOMORPHISM",
    "OPH_FEP",
    "PAIRS",
    "SOURCES_ENV",
    "SYNERGETICS_OPH",
    "THREE_WAY",
    "CitedIdentifier",
    "Correspondence",
    "CorrespondenceResult",
    "LicenceTerm",
    "PinCheck",
    "SourcePin",
    "by_key",
    "by_kind",
    "by_pair",
    "check_pins",
    "head_of",
    "load_sources",
    "summary",
    "verify",
]
