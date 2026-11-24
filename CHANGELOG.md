# Changelog

## Unreleased

- Require Python 3.14 or newer and use Python 3.14 for local CI. The pre-push
  hook requires a Python 3.14 pass receipt for the exact commit sent to `main`.
- Follow Fuller's fixed-bar Jitterbug through octahedral joint coalescence.
  Five probes distinguish nearest-only disconnection, conductance-policy
  obstructions, a hard Gaussian constraint, finite-weight pair-load memory
  and the averaged repair generator on an accelerated clock.
- Provide the Jitterbug graphical abstract on page 6 and generate a separate
  vector poster with the figure suite.
- Add Bernhard Mueller as second author with affiliations to Pragma Research
  and the Active Inference Institute. Record his ORCID in the manuscript,
  citation and Zenodo metadata.
- Add a structured ontology concordance and render its term comparisons,
  domains and evidence references as a manuscript supplement.
- Expand the cover into a graphical abstract and revise the manuscript figures
  and README diagram. Compare squared chords with Gaussian variance on the
  cover using adjacent port pairs as the common reference.
- Track `uv.lock` for reproducible installation of the runtime, development and
  figure dependencies.
- Add a README overview with one executable example for each pair of programs.
  Cite *Agreement and Surprise* at the results it supplies and distinguish
  those results from derivations and examples introduced here.
- Verify finite linear-Gaussian inference on the carrier, with a chain control
  and a nonquadratic counterexample. A proper record prior separates interface
  energy, unnormalised record potential and normalised surprisal. The two-record
  Schur complement recovers the exact resistance ladder.
- Distinguish ungated repair confluence from termination on the total-load-two
  shell. The broken-seam gate terminates with multiple reachable normal forms.
  The binary-chain probe separates pathwise load energy from ensemble relative
  entropy and checks the port relaxation gap.
- Replace the project GitHub Actions workflow with the local CI gate in
  `scripts/check_local_ci.py`.
- Round measured values in the generated tables for portability. The verifier
  continues to return its raw values.
- Distinguish the proper Gaussian's seam Markov blankets from the intrinsic
  zero-sum field's covariance and its constrained conditional laws.
- Include the pinned source register in regular installations so the
  `sources` command works in a wheel-only environment.
- Check the template's project lane before rendering so a manuscript build
  targets the requested checkout.
- Validate finite positive-definite Gaussian covariances. Clamp sheaf records
  before diffusion and use the same rank tolerance for sections and cohomology.
- Keep native repair energies in exact rationals and the cuboctahedral gauge
  in integers, including at scales beyond floating-point precision.
- Reject empty correspondence selections, distinguish absent pins from
  unavailable HEADs, and check for unresolved references in the rendered PDF.
- Reuse exact tetravolume derivations within a run and report tests, probes,
  and generated documents as distinct CI steps.

Corrections (a claim was wrong):

- The kissing-arrangement probe measures outer-sphere contacts. The
  `cuboctahedron_neighbours_touch` result replaces
  `cuboctahedron_is_rigid`, which remains a deprecated alias.
- The `13/15` and `1/60` transition probabilities and the `35/264` port-zero
  occupation of *Agreement and Surprise* were recorded `not-reproduced`. All
  three are determined by protocols in OPH's public model
  (`model/common_source/NATIVE_FIELDS.md`, `model/observer/native.py`) and are
  reproduced exactly: `pr01-obstruction` carries the two probabilities and
  `reduction-criterion` the occupation with both selector controls.
- *Agreement and Surprise* is a public preprint (PhilPapers record MUEAAS-2);
  the citation, the registry source strings and the provenance notes said
  "private draft". OPH is cited by its synthesis paper (PhilPapers record
  MUEFOC) beside the repository.
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
- The carrier's intersection array is `{5, 2, 1; 1, 2, 5}`. The manuscript
  token printed `(2, 1, 0; 1, 2, 5)`, the sequence `b_1, b_2, b_3` with `b_0`
  left out, against the registry claim and the standard form.
- The docstring of `path_multiplicity_census` assigned the self-loop
  multiplicities to the wrong classes: 52 belongs to the 30 binary states whose
  units share a seam, 50 to the other 36 binary states and to the 12 states
  with both units on one port.
- `docs/SCOPE.md` called OPH's patches "Planck-scale screen units". OPH defines
  a patch as a bounded local system and assigns it no length scale; the caution
  on scale separation says so and agrees with the manuscript's version.
- `.zenodo.json` related this record to the OPH repository as `isCitedBy`; the
  relation is `cites`.

Additions:

- `binary-conditional-variance` checks the Dirichlet form of uniform-seam
  repair against the mean conditional variance on all 66 binary states.
- `carrier-spectral-chords` compares exact icosahedral chords with the first
  positive Laplacian eigenspace contribution to Gaussian contrast variance.
  The probe reports both embedding ranks and the comparison with full
  effective resistance.
- Extend the equal-Schur chain control to interior posterior divergences and
  verify quadray residue classes for FCC sites and tetrahedral and octahedral
  interstices.
- `repair-is-conditioning`, a registered OPH-FEP identity between the native
  law and a free-energy functional: on the binary shell each native repair is
  the uniform reference conditioned on the datum the repair leaves visible, and
  that law minimises `F[q] = E_q[ln q - ln pi(s, o)]` over the fibre with
  minimum the surprisal of the datum. The probe asserts where the identity
  fails off the shell.
- `certificate-species` averages the native rule itself over the sixty
  seam-and-coin draws on all 78 configurations and finds `x - L x / 60` in exact
  rationals, so the mean repair operator is derived from the law.
  `core.repair.mean_native_step`, `core.controls.native_transition_probabilities`
  and `core.controls.record_feedback_occupation` are the new functions, and five
  OPH identifiers are added to `data/sources.json`.
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
