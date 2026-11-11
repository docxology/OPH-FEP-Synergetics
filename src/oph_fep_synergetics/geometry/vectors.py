"""Generic 3-vector arithmetic shared across the geometry modules.

Plain tuple-based helpers, deliberately free of ``numpy``: the geometry of
the quadrays, the IVM and the subdivision is small dense arithmetic over
floats and integers, and keeping the operations explicit keeps the
floating-point route -- and therefore every measured value -- reproducible
bit for bit.

One shared note on the tolerances defined here.  The geometry's natural
length scale is the IVM edge, ``sqrt(8) ~ 2.83``, so ``EDGE_TOLERANCE``
(1e-9) is an absolute slack of roughly 4e-10 of an edge and
``TIGHT_TOLERANCE`` (1e-12) roughly 4e-13; both absorb only accumulated
float noise, never real geometry.  ``PIVOT_TOLERANCE`` (1e-15) is smaller
by orders of magnitude because it guards the pivot of a Gaussian
elimination on matrices whose entries are exact small integers or simple
rational images of them, where a genuine nonzero pivot sits far above the
round-off floor and only a truly singular column can reach it.
"""

from __future__ import annotations

import math
from typing import Sequence

__all__ = [
    "EDGE_TOLERANCE",
    "PIVOT_TOLERANCE",
    "TIGHT_TOLERANCE",
    "Vector",
    "add",
    "cross",
    "dot",
    "norm",
    "scale",
    "subtract",
]

Vector = tuple[float, float, float]

#: Singular-pivot guard for the small dense integer-adjacent linear solves.
PIVOT_TOLERANCE: float = 1e-15
#: Tight absolute tolerance: duplicate-vertex floors and near-exactness checks.
TIGHT_TOLERANCE: float = 1e-12
#: Edge-scale absolute tolerance for shape and feasibility checks.
EDGE_TOLERANCE: float = 1e-9


def subtract(left: Sequence[float], right: Sequence[float]) -> Vector:
    return (left[0] - right[0], left[1] - right[1], left[2] - right[2])


def add(left: Sequence[float], right: Sequence[float]) -> Vector:
    return (left[0] + right[0], left[1] + right[1], left[2] + right[2])


def scale(vector: Sequence[float], factor: float) -> Vector:
    return (vector[0] * factor, vector[1] * factor, vector[2] * factor)


def dot(left: Sequence[float], right: Sequence[float]) -> float:
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2]


def cross(left: Sequence[float], right: Sequence[float]) -> Vector:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def norm(vector: Sequence[float]) -> float:
    return math.sqrt(dot(vector, vector))
