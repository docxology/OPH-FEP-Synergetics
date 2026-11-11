"""The icosahedral symmetry group of the carrier, and what it fixes.

The carrier's adjacency has exactly ``120`` automorphisms.  They form a closed
group, act transitively on the twelve ports, on the thirty seams and on the
twenty faces, and split into the sixty rotations and the sixty reflections.  The
rotation subgroup is detected without coordinates: an automorphism is a rotation
when it maps OPH's *oriented* generating faces onto the oriented face set, and a
reflection when it reverses every orientation.

Two further results belong to the same object:

* the Laplacian's eigenspaces of dimensions ``3``, ``5`` and ``3`` are
  irreducible under the group, so the spectrum is an irreducible-group
  invariant rather than a numerical accident;
* the native repair law is equivariant with respect to the whole group: the
  multiset of descent outcomes out of a configuration is unchanged by every
  relabelling, so the census is a group invariant of the carrier.
"""

from __future__ import annotations

from itertools import product

import numpy as np

from .carrier import IcosahedralCarrier, OPH_ORIENTED_BASE_FACES, icosahedral_carrier

__all__ = [
    "automorphism_group_is_closed",
    "carrier_automorphisms",
    "descent_outcome_multiset",
    "eigenspace_irreducibility",
    "law_equivariance",
    "orbit_sizes",
    "reflection_subgroup",
    "rotation_subgroup",
]


def _orientation(faces: tuple[tuple[int, int, int], ...]) -> frozenset[tuple[int, int, int]]:
    """Return the canonical form of an oriented face list.

    A face's orientation is its cyclic order, so ``(1, 7, 8)`` and ``(7, 8, 1)``
    are the same orientation while ``(1, 8, 7)`` is the opposite one.  Each face
    is therefore rotated so that its smallest port leads, which makes set
    comparison decide orientation.  (The directed-edge set cannot: on a closed
    surface it already contains both directions of every edge, so reversing every
    face leaves it unchanged.)
    """
    canonical: list[tuple[int, int, int]] = []
    for face in faces:
        start = face.index(min(face))
        canonical.append(face[start:] + face[:start])
    return frozenset(canonical)


def carrier_automorphisms(
    carrier: IcosahedralCarrier | None = None,
    *,
    limit: int = 400,
) -> tuple[tuple[int, ...], ...]:
    """Return every adjacency-preserving permutation of the ports.

    The search is a backtracking extension that only places a port where its
    neighbourhood relation to the already-placed ports still matches, so the
    traversal is bounded by the group itself rather than by ``12!``.
    """
    carrier = carrier or icosahedral_carrier()
    neighbours = [set(carrier.neighbours(port)) for port in range(carrier.port_count)]
    found: list[tuple[int, ...]] = []

    def extend(images: list[int | None]) -> None:
        if len(found) > limit:
            return
        index = next((i for i, value in enumerate(images) if value is None), None)
        if index is None:
            found.append(tuple(int(value) for value in images))
            return
        used = {value for value in images if value is not None}
        for candidate in range(carrier.port_count):
            if candidate in used:
                continue
            if any(
                value is not None
                and (candidate in neighbours[value]) != (index in neighbours[i])
                for i, value in enumerate(images)
            ):
                continue
            images[index] = candidate
            extend(images)
            images[index] = None

    extend([None] * carrier.port_count)
    if len(found) > limit:
        raise AssertionError("the automorphism search did not terminate in the group")
    return tuple(found)


def automorphism_group_is_closed(carrier: IcosahedralCarrier | None = None) -> bool:
    """Whether the automorphisms are closed under composition."""
    carrier = carrier or icosahedral_carrier()
    group = carrier_automorphisms(carrier)
    members = set(group)
    size = carrier.port_count
    return all(
        tuple(left[right[i]] for i in range(size)) in members for left in group for right in group
    )


def orbit_sizes(carrier: IcosahedralCarrier | None = None) -> dict[str, int]:
    """Return the number of group orbits on ports, seams and faces.

    An orbit is the set of images of one cell under the whole group, so the count
    is how many distinct such sets exist.  All three counts are one here: the
    group is transitive on the twelve ports, on the thirty seams and on the
    twenty faces.
    """
    carrier = carrier or icosahedral_carrier()
    group = carrier_automorphisms(carrier)
    ports = {
        frozenset(permutation[i] for permutation in group) for i in range(carrier.port_count)
    }
    seams = {
        frozenset(frozenset((permutation[a], permutation[b])) for permutation in group)
        for a, b in carrier.edges
    }
    faces = {
        frozenset(frozenset(permutation[index] for index in face) for permutation in group)
        for face in carrier.faces
    }
    return {"ports": len(ports), "seams": len(seams), "faces": len(faces)}


def rotation_subgroup(carrier: IcosahedralCarrier | None = None) -> tuple[tuple[int, ...], ...]:
    """Return the orientation-preserving automorphisms.

    OPH's generating faces are oriented, so an automorphism either carries the
    oriented face set to itself or mirrors every face.
    """
    carrier = carrier or icosahedral_carrier()
    oriented = _orientation(OPH_ORIENTED_BASE_FACES)
    rotations: list[tuple[int, ...]] = []
    for permutation in carrier_automorphisms(carrier):
        image = _orientation(
            tuple(
                (permutation[face[0]], permutation[face[1]], permutation[face[2]])
                for face in OPH_ORIENTED_BASE_FACES
            )
        )
        if image == oriented:
            rotations.append(permutation)
    return tuple(rotations)


def reflection_subgroup(carrier: IcosahedralCarrier | None = None) -> tuple[tuple[int, ...], ...]:
    """Return the orientation-reversing automorphisms, the mirror half of the group."""
    carrier = carrier or icosahedral_carrier()
    mirrored = _orientation(
        tuple((face[0], face[2], face[1]) for face in OPH_ORIENTED_BASE_FACES)
    )
    reflections: list[tuple[int, ...]] = []
    for permutation in carrier_automorphisms(carrier):
        image = _orientation(
            tuple(
                (permutation[face[0]], permutation[face[1]], permutation[face[2]])
                for face in OPH_ORIENTED_BASE_FACES
            )
        )
        if image == mirrored:
            reflections.append(permutation)
    return tuple(reflections)


def eigenspace_irreducibility(
    carrier: IcosahedralCarrier | None = None,
    *,
    probes: int = 4,
) -> dict[str, dict[str, int | bool]]:
    """Return, per eigenvalue, the dimension and whether the group acts irreducibly.

    Irreducibility is decided by Burnside's criterion in the constructive form:
    the orbit of a single generic vector must span the whole eigenspace.  Several
    deterministic probe vectors are used so that one unlucky choice cannot pass,
    and ``probes`` bounds how many of them are tried per eigenspace.
    """
    carrier = carrier or icosahedral_carrier()
    group = carrier_automorphisms(carrier)
    laplacian = carrier.laplacian()
    eigenvalues, eigenvectors = np.linalg.eigh(laplacian)
    groups: dict[float, list[int]] = {}
    for index, value in enumerate(eigenvalues):
        groups.setdefault(round(float(value), 8), []).append(index)
    report: dict[str, dict[str, int | bool]] = {}
    for value, indices in sorted(groups.items()):
        dimension = len(indices)
        if dimension <= 1:
            report[f"{value:.6f}"] = {"dim": dimension, "orbit_rank": dimension, "irreducible": True}
            continue
        basis = eigenvectors[:, indices]
        ranks: list[int] = []
        for weights in product((1.0, -0.5, 2.0), repeat=min(dimension, 3)):
            seed = np.array(list(weights) + [0.0] * (dimension - len(weights)))[:dimension]
            if not np.any(seed):
                seed[0] = 1.0
            vector = basis @ seed
            orbit = np.array([vector[list(permutation)] for permutation in group])
            ranks.append(int(np.linalg.matrix_rank(orbit, tol=1e-7)))
        report[f"{value:.6f}"] = {
            "dim": dimension,
            "orbit_rank": max(ranks),
            "min_orbit_rank": min(ranks),
            "irreducible": bool(max(ranks) == dimension),
        }
    return report


def descent_outcome_multiset(configuration: tuple[int, ...]) -> dict[str, int]:
    """Return how the sixty native transitions out of a configuration divide.

    Strict descents lower the seam energy, zero descents leave it unchanged with
    a changed configuration, and no-ops leave the configuration alone.  The
    classification is computed from the law itself, not from the census.
    """
    from .carrier import icosahedral_carrier as _carrier
    from .repair import energy, native_repair

    carrier = _carrier()
    outcomes = {"strict": 0, "zero": 0, "noop": 0}
    before_energy = energy(configuration)
    for a, b in carrier.edges:
        for coin in (0, 1):
            step = native_repair(configuration, (a, b), coin)
            if step.after == step.before:
                outcomes["noop"] += 1
            elif energy(step.after) < before_energy - 1e-12:
                outcomes["strict"] += 1
            else:
                outcomes["zero"] += 1
    return outcomes


def law_equivariance(carrier: IcosahedralCarrier | None = None) -> dict[str, object]:
    """Whether the native law commutes with the whole automorphism group.

    Every configuration's outcome multiset is compared against the multiset of
    every relabelled configuration, so a single symmetry-breaking witness would
    make the result false.
    """
    carrier = carrier or icosahedral_carrier()
    group = carrier_automorphisms(carrier)
    configurations = list(carrier.shell_configurations(total=2))
    outcomes = {configuration: descent_outcome_multiset(configuration) for configuration in configurations}
    invariant = True
    witness: tuple[tuple[int, ...], tuple[int, ...]] | None = None
    for permutation in group:
        for configuration in configurations:
            moved = tuple(configuration[index] for index in permutation)
            if outcomes[configuration] != outcomes[moved]:
                invariant = False
                witness = (configuration, moved)
                break
        if not invariant:
            break
    orbits = {frozenset(tuple(configuration[index] for index in permutation) for permutation in group) for configuration in configurations}
    return {
        "invariant": invariant,
        "configurations": len(configurations),
        "group_order": len(group),
        "orbits": len(orbits),
        "distinct_outcome_multisets": len({tuple(sorted(outcome.items())) for outcome in outcomes.values()}),
        "witness": witness,
    }
