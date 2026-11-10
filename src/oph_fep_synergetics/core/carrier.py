"""The twelve-port icosahedral carrier.

The carrier is the graph of the regular icosahedron: twelve ports (vertices),
thirty seams (edges) and twenty faces.  It is the geometric object shared by
Observer Patch Holography's executable observer and the kissing-number-12
configuration of the face-centred-cubic / A3 lattice.  The module exposes the
graph, its Laplacian, and the finite shell of load configurations used by the
native repair law.

Everything here is exact geometry on the unit sphere up to floating-point
rounding; the combinatorial counts (12, 30, 20, 5-regular) are asserted by the
test suite.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from typing import Iterator

import numpy as np

PHI = (1.0 + math.sqrt(5.0)) / 2.0

#: OPH's twenty oriented faces of the twelve-port carrier, transcribed from
#: ``evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/
#: self_readback_repair_closure.py`` (``ORIENTED_BASE_FACES``) at commit
#: ``5517b0c77a6b11bc141ab663647bf7ec2b2d7661``.  OPH pins that archive by
#: SHA-256 in ``evidence/native_geometric_source_20260925/spec.json``.
#:
#: The adjacency below is *derived* from these faces rather than transcribed, so
#: this package uses OPH's own port chart and the derivation is checkable: the
#: triangles of the derived graph are exactly these twenty triples.  OPH's
#: antipodal involution is not ``i <-> 11 - i``; the graph-distance-three pairs
#: are ``0-3, 1-2, 4-7, 5-6, 8-11, 9-10``.
OPH_ORIENTED_BASE_FACES: tuple[tuple[int, int, int], ...] = (
    (0, 11, 5),
    (0, 5, 1),
    (0, 1, 7),
    (0, 7, 10),
    (0, 10, 11),
    (1, 5, 9),
    (5, 11, 4),
    (11, 10, 2),
    (10, 7, 6),
    (7, 1, 8),
    (3, 9, 4),
    (3, 4, 2),
    (3, 2, 6),
    (3, 6, 8),
    (3, 8, 9),
    (4, 9, 5),
    (2, 4, 11),
    (6, 2, 10),
    (8, 6, 7),
    (9, 8, 1),
)


def neighbours_from_faces(
    faces: tuple[tuple[int, int, int], ...], *, port_count: int = 12
) -> tuple[tuple[int, ...], ...]:
    """Derive sorted port adjacency from an oriented face list."""
    adjacency = [set() for _ in range(port_count)]
    for face in faces:
        if len(set(face)) != 3 or any(not 0 <= port < port_count for port in face):
            raise ValueError("a face must be three distinct in-range ports")
        for u, v in zip(face, face[1:] + face[:1]):
            adjacency[u].add(v)
            adjacency[v].add(u)
    return tuple(tuple(sorted(row)) for row in adjacency)


#: Port adjacency of OPH's executable carrier, derived from the oriented faces
#: above.  Every port has degree five.
OPH_NEIGHBOURS: tuple[tuple[int, ...], ...] = neighbours_from_faces(OPH_ORIENTED_BASE_FACES)


def _icosahedron_vertices() -> tuple[tuple[float, float, float], ...]:
    """Return the twelve vertices of a regular icosahedron on the unit sphere."""
    raw: list[tuple[float, float, float]] = []
    for s1 in (1.0, -1.0):
        for s2 in (1.0, -1.0):
            raw.extend(
                [
                    (0.0, s1 * 1.0, s2 * PHI),
                    (s1 * 1.0, s2 * PHI, 0.0),
                    (s2 * PHI, 0.0, s1 * 1.0),
                ]
            )
    norm = math.sqrt(1.0 + PHI * PHI)
    return tuple((x / norm, y / norm, z / norm) for x, y, z in raw)


def _squared_distance(
    a: tuple[float, float, float], b: tuple[float, float, float]
) -> float:
    return sum((x - y) ** 2 for x, y in zip(a, b))


@dataclass(frozen=True)
class IcosahedralCarrier:
    """A finite icosahedral carrier graph.

    Attributes
    ----------
    vertices:
        Twelve unit vectors.
    edges:
        Thirty unordered port pairs (seams), sorted.
    faces:
        Twenty unordered port triples, sorted.
    """

    vertices: tuple[tuple[float, float, float], ...]
    edges: tuple[tuple[int, int], ...]
    faces: tuple[tuple[int, int, int], ...]

    @property
    def port_count(self) -> int:
        return len(self.vertices)

    @property
    def seam_count(self) -> int:
        return len(self.edges)

    @property
    def face_count(self) -> int:
        return len(self.faces)

    def neighbours(self, port: int) -> tuple[int, ...]:
        """Return the ports joined to ``port`` by a seam."""
        return tuple(
            b if a == port else a
            for a, b in self.edges
            if a == port or b == port
        )

    def degree(self, port: int) -> int:
        return len(self.neighbours(port))

    def seam_index(self, a: int, b: int) -> int:
        """Return the index of the seam joining distinct ports ``a`` and ``b``."""
        if a == b:
            raise ValueError("a seam joins two distinct ports")
        pair = (min(a, b), max(a, b))
        try:
            return self.edges.index(pair)
        except ValueError as exc:  # pragma: no cover - defensive
            raise ValueError(f"ports {a} and {b} are not joined by a seam") from exc

    def adjacency_matrix(self) -> np.ndarray:
        matrix = np.zeros((self.port_count, self.port_count), dtype=float)
        for a, b in self.edges:
            matrix[a, b] = matrix[b, a] = 1.0
        return matrix

    def laplacian(self) -> np.ndarray:
        """Return the unweighted graph Laplacian ``L = D - A``."""
        adjacency = self.adjacency_matrix()
        return np.diag(adjacency.sum(axis=1)) - adjacency

    def laplacian_spectrum(self) -> np.ndarray:
        """Return the sorted Laplacian eigenvalues."""
        return np.sort(np.linalg.eigvalsh(self.laplacian()))

    def spectral_gap(self) -> float:
        """Return the first positive Laplacian eigenvalue."""
        spectrum = self.laplacian_spectrum()
        positive = spectrum[spectrum > 1e-9]
        if positive.size == 0:  # pragma: no cover - impossible for this graph
            raise ValueError("carrier graph is disconnected")
        return float(positive[0])

    def shell_configurations(self, total: int = 2) -> tuple[tuple[int, ...], ...]:
        """All nonnegative integer load vectors of the given total.

        For ``total=2`` there are ``C(13, 11) = 78`` configurations; the 66 with
        every load in ``{0, 1}`` are the binary shell used by OPH.
        """
        if total < 0:
            raise ValueError("total load must be nonnegative")
        return tuple(_compositions(self.port_count, total))

    def binary_shell_configurations(self) -> tuple[tuple[int, ...], ...]:
        """The 66 configurations with two unit loads in distinct ports."""
        states: list[tuple[int, ...]] = []
        for a, b in itertools.combinations(range(self.port_count), 2):
            state = [0] * self.port_count
            state[a] = state[b] = 1
            states.append(tuple(state))
        return tuple(states)


def _compositions(width: int, total: int) -> Iterator[tuple[int, ...]]:
    """Yield the weak compositions of ``total`` into ``width`` parts."""
    if width == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in _compositions(width - 1, total - first):
            yield (first, *tail)


def icosahedral_carrier() -> IcosahedralCarrier:
    """Construct OPH's labelled carrier from its canonical adjacency list.

    The faces are the triangles of the adjacency graph.  The result is the
    combinatorial carrier used by OPH; :func:`geometric_carrier` supplies the
    Euclidean realisation used for the kissing/geometry checks.
    """
    vertices = _icosahedron_vertices()
    edges = tuple(
        sorted(
            (a, b)
            for a, neighbours in enumerate(OPH_NEIGHBOURS)
            for b in neighbours
            if a < b
        )
    )
    faces = _triangles(OPH_NEIGHBOURS)
    return IcosahedralCarrier(vertices=vertices, edges=edges, faces=faces)


def _triangles(neighbours: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, int, int], ...]:
    adjacency = [set(row) for row in neighbours]
    return tuple(
        sorted(
            (a, b, c)
            for a, b, c in itertools.combinations(range(len(neighbours)), 3)
            if b in adjacency[a] and c in adjacency[a] and c in adjacency[b]
        )
    )


def geometric_carrier() -> IcosahedralCarrier:
    """Construct the carrier by measuring a regular icosahedron on the sphere."""
    vertices = _icosahedron_vertices()
    pairwise = [
        (a, b, math.sqrt(_squared_distance(vertices[a], vertices[b])))
        for a, b in itertools.combinations(range(len(vertices)), 2)
    ]
    edge_length = min(distance for _, _, distance in pairwise if distance > 1e-12)
    edges = tuple(
        sorted(
            (a, b)
            for a, b, distance in pairwise
            if abs(distance - edge_length) <= 1e-9
        )
    )
    adjacency = {vertex: set() for vertex in range(len(vertices))}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    faces = tuple(
        sorted(
            (a, b, c)
            for a, b, c in itertools.combinations(range(len(vertices)), 3)
            if b in adjacency[a] and c in adjacency[a] and c in adjacency[b]
        )
    )
    return IcosahedralCarrier(vertices=vertices, edges=edges, faces=faces)


def graph_signature(carrier: IcosahedralCarrier) -> tuple[tuple[int, ...], tuple[float, ...]]:
    """A degree/spectrum signature used to compare two constructions."""
    degrees = tuple(carrier.degree(port) for port in range(carrier.port_count))
    spectrum = tuple(round(float(x), 9) for x in carrier.laplacian_spectrum())
    return degrees, spectrum


def _ring_signature(adjacency: list[set[int]], vertex: int) -> tuple:
    """A two-ring degree signature used to order the isomorphism search."""
    first = sorted(len(adjacency[neighbour]) for neighbour in adjacency[vertex])
    second = sorted(
        len(adjacency[far])
        for neighbour in adjacency[vertex]
        for far in adjacency[neighbour]
    )
    return (len(adjacency[vertex]), tuple(first), tuple(second))


def carrier_isomorphism(
    left: IcosahedralCarrier,
    right: IcosahedralCarrier,
) -> tuple[int, ...] | None:
    """Return a port permutation carrying ``left``'s graph onto ``right``'s.

    The returned tuple ``p`` maps port ``a`` of ``left`` to port ``p[a]`` of
    ``right`` and preserves adjacency exactly.  ``None`` means the two graphs are
    not isomorphic.  The search orders ports by a two-ring degree signature and
    backtracks, which decides isomorphism at this size; it is a witness, so a
    returned permutation can be checked directly.
    """
    count = left.port_count
    if count != right.port_count:
        return None
    if left.seam_count != right.seam_count:
        return None
    left_adjacency = [set(left.neighbours(port)) for port in range(count)]
    right_adjacency = [set(right.neighbours(port)) for port in range(count)]
    left_signature = {port: _ring_signature(left_adjacency, port) for port in range(count)}
    right_signature = {port: _ring_signature(right_adjacency, port) for port in range(count)}
    if sorted(left_signature.values()) != sorted(right_signature.values()):
        return None
    order = sorted(range(count), key=lambda port: (-len(left_adjacency[port]), left_signature[port]))
    mapping: dict[int, int] = {}
    used: set[int] = set()

    def consistent(source: int, target: int) -> bool:
        if left_signature[source] != right_signature[target]:
            return False
        return all(
            (other in left_adjacency[source]) == (image in right_adjacency[target])
            for other, image in mapping.items()
        )

    def extend(index: int) -> bool:
        if index == count:
            return True
        source = order[index]
        for target in range(count):
            if target in used or not consistent(source, target):
                continue
            mapping[source] = target
            used.add(target)
            if extend(index + 1):
                return True
            del mapping[source]
            used.discard(target)
        return False

    if not extend(0):
        return None
    return tuple(mapping[port] for port in range(count))
