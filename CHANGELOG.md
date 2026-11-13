# Changelog

## Unreleased

Corrections (a claim was wrong):

- The coordination shell `10k^2 + 2` is the set of lattice points at graph
  distance `k` (the surface of the frequency-`k` vector equilibrium), not the
  Euclidean sphere of squared norm `8k^2`; it now has a closed form, the
  cuboctahedral gauge, checked against every breadth-first shell through 12.
- `shell-vs-free-energy` no longer claims that a bounded free energy fails to
  bound the shell (it does bound it: the minimum on shell `k` is `2k^2`); the
  obstruction is that neither determines the other.
- The carrier Laplacian is the generator of the mean repair operator and the
  seam-disagreement energy, not the native repair energy `sum x^2 / 2`.
- Conditional *dependence* lives on the seams; independence is off them.
- 5184000 counts spanning trees, not acyclic seam sets.
- The proposition on the carrier graph says it is five-regular (it said twelve).
- The Lyapunov certificate solves `A P + P A^T = -2I` for the unscaled matrix;
  the antisymmetric `Q` is the Helmholtz solenoidal part, not the right-hand
  side. D-stability is now defined (Arrow and McManus).
- At the mediator threshold the pairwise blocks become rank one; the Schur
  complement is not singular there.
- Only the adjoint `-J^{-T} c` is the gradient; the direct response is what a
  forward nudge measures.
- The rhombic dodecahedron is not four-regular; the icosahedron's twelve
  induced five-cycles are the dual dodecahedron's faces, not its vertices.
- The carrier is the icosahedral kissing arrangement (with slack), not the
  densest packing's cuboctahedral one; the normalised gap `(5 - sqrt 5) / 60`
  is a quadratic surd, not rational; the port variance `7/36` is the diagonal
  of `L^+`, the variance of a port's deviation from the carrier mean.
- Two bibliography DOIs were wrong (Biehl et al. 2021; Da Costa et al., now
  dated 2020) and are corrected.

Additions:

- `geometry.exact`: exact arithmetic in `Q(sqrt 5)`, exact convex hulls, the
  six-entry tetravolume table derived from explicit integer vertex sets, the
  icosahedron's tetravolume as `sqrt 2 (15/2 + (5/2) sqrt 5)`, and an exact
  comparison of the two twelve-sphere kissing arrangements.
- Surface cohomology of the carrier's two-complex: Betti numbers `(1, 0, 1)`,
  Euler characteristic 2, every graph obstruction a face boundary, and OPH's
  oriented faces a fundamental cycle.
- Independent routes for previously self-referential probes: a KKT constrained
  minimiser, the free energy as energy minus entropy, closed-form shell counts,
  transition-level equivariance over all 1123200 transitions, an exact
  fixed-point enumeration of the clipped map, an exact rational inverse for the
  reciprocal spectrum sum, and an asymmetric-payoff negative control for the
  replicator potential.
- Twenty-one verified primary references (Kirchhoff 1847, Gauss 1840, Hales
  2005 and 2017, Conway and Sloane 1997, Caspar and Klug 1962, Arrow and McManus
  1958, Dempster 1972, Rue and Held 2005, Fields et al. 2022 and others), and
  DOIs for every entry that has one.
- Source pins re-read at current upstream for all three sources, with corrected
  line ranges and GitHub casing.

Maintenance:

- Bound values cast every probe-derived scalar to a Python number, so the
  manuscript no longer renders `np.float64(...)` wrappers.
- `scripts/build.py` prints the interpreter and numpy versions of each build.

## 0.1.0 — initial public edition

- 42 registered correspondences (17 isomorphisms, 15 disagreements, 10
  complementarities) with executable probes; `oph-fep-synergetics verify`
  reports pass, fail, and exit-code discipline for typos as well as failures.
- A 347-test suite with no mocks, no network, and a 90% coverage gate enforced
  by `pyproject.toml`.
- 22 figures with an accessibility registry, and a render contract that reads
  the delivered PDF back: front-matter page placement, measured margins,
  registry keys in the printed text, and a clean LaTeX log.
- Eleven body sections plus the graphical abstract, in which every number is a
  bound token resolved from the package that produces it.
- The record of changes between editions is the git history; the editorial log
  in [docs/PROVENANCE.md](docs/PROVENANCE.md) covers source-pin corrections
  only.
