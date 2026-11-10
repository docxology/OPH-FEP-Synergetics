"""Cellular sheaves over a finite graph.

A cellular sheaf assigns a real vector space (a stalk) to every node and a
linear restriction map to every edge endpoint.  A global section is an
assignment on which every edge's two restrictions agree.  The sheaf Laplacian
``L = delta^T delta`` turns agreement into an energy, and harmonic extension
with clamped boundary values is exactly the linear-Gaussian perception flow of
the free-energy dictionary in :mod:`oph_fep_synergetics.core.free_energy`.

The implementation is deliberately small and explicit: dense matrices,
finite-dimensional stalks, and no topology beyond a graph.  That is enough to
state the three-way correspondence between OPH seams, FEP blankets, and
Synergetic cell complexes at the level of linear algebra.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .carrier import IcosahedralCarrier


@dataclass(frozen=True)
class EdgeRestriction:
    """The two linear restriction maps of one edge.

    An edge joins nodes ``u`` and ``v``.  ``map_u`` has shape
    ``(edge_dim, stalk_dim[u])`` and ``map_v`` has shape
    ``(edge_dim, stalk_dim[v])``.  Agreement on the edge means
    ``map_u z_u == map_v z_v``.
    """

    u: int
    v: int
    map_u: np.ndarray
    map_v: np.ndarray


class CellularSheaf:
    """A finite cellular sheaf on a graph."""

    def __init__(
        self,
        stalk_dims: Sequence[int],
        edges: Sequence[EdgeRestriction],
        *,
        node_names: Sequence[str] | None = None,
    ) -> None:
        self.stalk_dims = tuple(int(d) for d in stalk_dims)
        if not self.stalk_dims or any(d <= 0 for d in self.stalk_dims):
            raise ValueError("stalk dimensions must be positive")
        self.node_names = tuple(node_names) if node_names is not None else tuple(
            f"n{i}" for i in range(len(self.stalk_dims))
        )
        if len(self.node_names) != len(self.stalk_dims):
            raise ValueError("node_names must match the number of stalks")
        normalized: list[EdgeRestriction] = []
        for edge in edges:
            if edge.u == edge.v:
                raise ValueError("an edge joins two distinct nodes")
            if not (0 <= edge.u < len(self.stalk_dims) and 0 <= edge.v < len(self.stalk_dims)):
                raise IndexError("edge endpoint outside the node set")
            map_u = np.atleast_2d(np.asarray(edge.map_u, dtype=float))
            map_v = np.atleast_2d(np.asarray(edge.map_v, dtype=float))
            if map_u.shape[0] != map_v.shape[0]:
                raise ValueError("edge restriction maps must share an edge stalk")
            if map_u.shape[1] != self.stalk_dims[edge.u] or map_v.shape[1] != self.stalk_dims[edge.v]:
                raise ValueError("edge restriction map shape does not match its stalk")
            normalized.append(EdgeRestriction(edge.u, edge.v, map_u, map_v))
        self.edges = tuple(normalized)

    # -- identity and sizes -------------------------------------------------
    @property
    def node_count(self) -> int:
        return len(self.stalk_dims)

    @property
    def total_dim(self) -> int:
        return int(sum(self.stalk_dims))

    def offsets(self) -> tuple[int, ...]:
        offsets = [0]
        for dim in self.stalk_dims:
            offsets.append(offsets[-1] + dim)
        return tuple(offsets)

    def pack(self, vectors: Sequence[np.ndarray]) -> np.ndarray:
        if len(vectors) != self.node_count:
            raise ValueError("one vector per node is required")
        return np.concatenate([np.asarray(v, dtype=float) for v in vectors])

    def unpack(self, vector: np.ndarray) -> tuple[np.ndarray, ...]:
        vector = np.asarray(vector, dtype=float)
        if vector.shape != (self.total_dim,):
            raise ValueError("packed vector has the wrong total dimension")
        offsets = self.offsets()
        return tuple(vector[offsets[i] : offsets[i + 1]] for i in range(self.node_count))

    # -- sheaf operators ----------------------------------------------------
    def coboundary(self, vector: np.ndarray) -> np.ndarray:
        """Return the stacked edge disagreements ``(delta z)_a``."""
        nodes = self.unpack(vector)
        pieces = [
            edge.map_v @ nodes[edge.v] - edge.map_u @ nodes[edge.u]
            for edge in self.edges
        ]
        return np.concatenate(pieces) if pieces else np.zeros(0)

    def laplacian(self) -> np.ndarray:
        """Return the sheaf Laplacian ``delta^T delta`` as a dense matrix."""
        matrix = np.zeros((self.total_dim, self.total_dim), dtype=float)
        offsets = self.offsets()
        for edge in self.edges:
            u, v = edge.u, edge.v
            fu, fv = edge.map_u, edge.map_v
            su = slice(offsets[u], offsets[u + 1])
            sv = slice(offsets[v], offsets[v + 1])
            matrix[su, su] += fu.T @ fu
            matrix[sv, sv] += fv.T @ fv
            matrix[su, sv] -= fu.T @ fv
            matrix[sv, su] -= fv.T @ fu
        return matrix

    def energy(self, vector: np.ndarray) -> float:
        """Return ``0.5 * ||delta z||^2``."""
        disagreement = self.coboundary(vector)
        return 0.5 * float(disagreement @ disagreement)

    def global_section_basis(self, tolerance: float = 1e-9) -> np.ndarray:
        """Return an orthonormal basis of global sections as columns."""
        eigenvalues, eigenvectors = np.linalg.eigh(self.laplacian())
        return eigenvectors[:, eigenvalues <= tolerance]

    def is_global_section(self, vector: np.ndarray, tolerance: float = 1e-9) -> bool:
        return float(np.linalg.norm(self.coboundary(vector))) <= tolerance

    # -- clamped boundary problems -----------------------------------------
    def harmonic_extension(
        self,
        observations: Mapping[int, Sequence[float]],
        *,
        tolerance: float = 1e-12,
    ) -> tuple[np.ndarray, ...]:
        """Solve for free stalks with the observed stalks clamped.

        This is the unique minimiser of the sheaf energy over assignments that
        agree with ``observations`` when the free block is positive definite.
        """
        observed = {int(node): np.asarray(value, dtype=float) for node, value in observations.items()}
        for node, value in observed.items():
            if not 0 <= node < self.node_count:
                raise IndexError(f"observed node {node} is outside the sheaf")
            if value.shape != (self.stalk_dims[node],):
                raise ValueError(f"observation at node {node} has the wrong dimension")
        free = [node for node in range(self.node_count) if node not in observed]
        laplacian = self.laplacian()
        offsets = self.offsets()
        solution = np.zeros(self.total_dim, dtype=float)
        for node, value in observed.items():
            solution[offsets[node] : offsets[node + 1]] = value
        if not free:
            nodes = self.unpack(solution)
            if tolerance is not None and not self.is_global_section(solution, tolerance=tolerance):
                raise ValueError("the clamped assignment is not a global section")
            return nodes
        free_indices = [
            index
            for node in free
            for index in range(offsets[node], offsets[node + 1])
        ]
        observed_indices = [
            index
            for node in sorted(observed)
            for index in range(offsets[node], offsets[node + 1])
        ]
        l_ff = laplacian[np.ix_(free_indices, free_indices)]
        l_fo = laplacian[np.ix_(free_indices, observed_indices)]
        free_vector = np.linalg.solve(l_ff, -(l_fo @ solution[observed_indices]))
        solution[free_indices] = free_vector
        return self.unpack(solution)

    def record_injectivity(self, record_nodes: Sequence[int], *, tolerance: float = 1e-9) -> bool:
        """Whether records separate global sections.

        Equivalently, the interior block of the Laplacian is positive definite.
        """
        records = set(int(node) for node in record_nodes)
        if any(not 0 <= node < self.node_count for node in records):
            raise IndexError("record node outside the sheaf")
        interior = [node for node in range(self.node_count) if node not in records]
        if not interior:
            return True
        offsets = self.offsets()
        indices = [
            index
            for node in interior
            for index in range(offsets[node], offsets[node + 1])
        ]
        block = self.laplacian()[np.ix_(indices, indices)]
        return bool(np.linalg.eigvalsh(0.5 * (block + block.T))[0] > tolerance)

    def diffuse(
        self,
        initial: np.ndarray,
        observations: Mapping[int, Sequence[float]],
        *,
        steps: int = 200,
        step_size: float = 0.1,
    ) -> np.ndarray:
        """Clamped sheaf diffusion: ``z <- z - eta L z`` with boundary values held."""
        if steps < 0:
            raise ValueError("steps must be nonnegative")
        state = np.asarray(initial, dtype=float).copy()
        if state.shape != (self.total_dim,):
            raise ValueError("initial vector has the wrong total dimension")
        offsets = self.offsets()
        laplacian = self.laplacian()
        for _ in range(steps):
            state = state - step_size * (laplacian @ state)
            for node, value in observations.items():
                state[offsets[node] : offsets[node + 1]] = np.asarray(value, dtype=float)
        return state


def linear_chain_sheaf(node_count: int, *, edge_weight: float = 1.0) -> CellularSheaf:
    """A scalar sheaf on a path with identity restrictions and equal weights."""
    if node_count < 1:
        raise ValueError("a chain needs at least one node")
    edges = tuple(
        EdgeRestriction(
            u=i,
            v=i + 1,
            map_u=np.array([[np.sqrt(edge_weight)]]),
            map_v=np.array([[np.sqrt(edge_weight)]]),
        )
        for i in range(node_count - 1)
    )
    return CellularSheaf(stalk_dims=(1,) * node_count, edges=edges)


def carrier_sheaf(carrier: IcosahedralCarrier) -> CellularSheaf:
    """The scalar identity sheaf whose Laplacian is the carrier Laplacian."""
    edges = tuple(
        EdgeRestriction(u=a, v=b, map_u=np.array([[1.0]]), map_v=np.array([[1.0]]))
        for a, b in carrier.edges
    )
    return CellularSheaf(stalk_dims=(1,) * carrier.port_count, edges=edges)
