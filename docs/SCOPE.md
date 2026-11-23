# Scope and assumptions

The registry tests finite correspondences among OPH, FEP, and Synergetics.
Each verdict is specific to the model, parameters, and inputs described by
its probe.

## Tested here

| subject | how |
| --- | --- |
| OPH's twelve-port carrier | adjacency is derived from OPH's `ORIENTED_BASE_FACES`. The 20 graph triangles recover the face vertex sets. |
| OPH's native seam law | reimplemented on exact Python integers; conservation and the energy identity are checked on every one of the 4680 shell transitions |
| binary repair as inference | fair repair resamples the uniform binary reference conditioned on the preserved seam record. Uniformly sampled seams give a repair Dirichlet form equal to the mean conditional variance. |
| PR-01 | the obstruction is reproduced exhaustively: 120 strict descents, 600 zero-descent mismatch-one moves, 3960 no-ops |
| patch agreement as a sheaf condition | a cellular sheaf over the carrier whose Laplacian is asserted equal to the carrier Laplacian |
| the free-energy dictionary | on the carrier with a Gaussian record prior, and on a chain control: harmonic extension, a separate KKT solve, Schur complement, three surprises, and two routes to variational free energy |
| obstruction and dynamics controls | the D-stability family, the mediator reduction, the reduction criterion, the two-node orbits, the solenoidal family, the four-cycle trap, the replicator potential |
| QuadMath's conventions | the quadray gauge and embedding, twelve neighbour moves, coordination shells compared with the closed-form cuboctahedral gauge, quadray residue classes for the FCC sites and two kinds of interstices, and the tetravolume table derived from explicit vertex sets |
| icosahedral subdivision | geodesic and Goldberg counts, duality, and the equality with the IVM shell sequence |
| the carrier's counting identities | the reduced Laplacian determinant in exact rationals (the spanning-tree number), its cross-check against the spectral product, the graph sheaf's two cohomology dimensions, and the surface's Betti numbers with its Euler characteristic |
| the carrier's symmetry | the automorphism group by exhaustive search, closure under composition, the orbit decomposition of ports, seams and faces, classification by preservation of face orientation, eigenspace irreducibility, and equivariance of each shell transition |
| two Gaussian readings of the carrier | the proper precision `L + I/n` gives the conditional-independence pattern, blankets and pairwise precisions. The diagonal of `L^+` gives the variance in the centered intrinsic field |
| the information cost of a lattice displacement | equal-covariance Gaussian divergence at unit covariance gives 4 nats for one IVM neighbour displacement and `4k^2` nats from the origin along a neighbour direction. The value uses the declared embedding and covariance scale. |
| the carrier's solid and the twelve-around-one cell | the icosahedron's tetravolume computed exactly in `sqrt 2 * Q(sqrt 5)` and by floating-point measurement, the two kissing arrangements compared exactly, and the Voronoi cell's counts, exact tetravolume and packing density |
| the carrier as a network | exact effective resistances as rationals from a grounded solve, their distance-class ladder, the Kirchhoff index, the distance distribution and intersection array around a port, and the subgraph census |
| Gaussian port contrasts and icosahedral chords | on the edge-length-two carrier, squared chord length is 40 times the contrast variance contributed by the first positive Laplacian eigenspace at unit intrinsic precision. Full contrast variance has a higher embedding dimension. |
| the carrier's curvature | the angle defect at every vertex of the derived chart, summed against Gauss-Bonnet, plus the defect of the dual solid's face set |
| the twelve-around-one families | compare the cuboctahedron and icosahedron as solids with distinct adjacency and volume. The cube's line graph independently realises the cuboctahedral graph. |
| the quadratic field of the tetravolume | the Galois conjugation of $\sqrt{5}$ (fixing $\sqrt{2}$) as a probe, in exact field arithmetic: it fixes the counting invariants and moves the icosahedron's tetravolume, whose product with its conjugate is exactly 50 |
| the transfer of binary repair to the FCC first shell | replay the same fair coin and uniform seam sampling on the cuboctahedral graph. Contrast variance gives an exact quotient on both graphs. The FCC graph-distance and both centroid processes retain observable memory under stationarity. |
| geometric Gaussian comparisons | test reciprocal conductance duality on two spherical graph pairs, compare the pulled-back seam precision for unit-direction quadratics, and use scalar Gaussian precision and vertex cuts to test central-plane screening. |
| transfer to larger graphs and cellulations | compare two-token repair on a four-cube, complete bipartite graphs and periodic FCC. Test covariance closure by Laplacian invariance, toroidal harmonic covariance, quadratic precision on nearest-neighbour root shells and screening on larger FCC boundaries. |
| Fuller’s Jitterbug | track the fixed bars through the regular icosahedron to octahedral joint coincidences. Compare hull, bar and nearest-pair graphs, then test conductance policies, the constrained Gaussian and the accelerated limit of fair binary repair. |

## Cited but not implemented

- **OPH's quantum event algebra and its continuum limit.** Nothing here
  evaluates them, and `data/sources.json` records no identifier for either.
- **General FEP methods.** The Gaussian inference model is linear with
  quadratic energies. An asymmetric seam gives a non-Gaussian control and the
  binary shell gives a finite Bayesian example. General nonlinear inference,
  path integrals, and stochastic-thermodynamic formulations remain outside
  this scope.
- **Fuller's full synergetic apparatus.** Only the pieces the correspondence
  needs are here: quadray coordinates, the IVM, the vector-equilibrium shell,
  volume ratios, subdivision counts, and the fixed-bar Jitterbug branch
  through octahedral coalescence.
- **OPH's transactional and integer consensus contracts.** OPH's consensus
  theorems hold for those contracts, not for the native law. This package tests
  the native law, which is the one the carrier actually runs.
- **Equilibrium propagation as a training method.** It appears here as a
  complementarity from energy-based learning, which the FEP reaches through
  predictive coding: it supplies the parameter gradient that the native repair
  law does not define. No learning is implemented, and the probe is a scalar
  energy model.
- **Other Jitterbug motion branches.** Fuller and Stachel describe a separate
  continuation from the octahedron toward a tetrahedral configuration. The
  tests follow the branch from the open vector equilibrium to the
  octahedral limit and treat the joint positions as an ideal linkage.

## Three source quantities and their protocols

*Agreement and Surprise* states three numbers for the native law. Each is
determined by a protocol in OPH's public model, and each is reproduced exactly
by a probe:

| quantity | protocol | OPH source | probe |
| --- | --- | --- | --- |
| `13/15` and `1/60` | one attempt, a uniform seam and a fair coin (sixty equally likely draws), from `(1, 1, 0, ..., 0)`; `13/15` is the self-transition and `1/60` is one particular hop, both at zero quadratic decrease | `model/common_source/NATIVE_FIELDS.md` | `pr01-obstruction` |
| `35/264` | uniform over the 66 binary configurations; event one repairs `(0, 1)` and records the signed current; event two repairs `(0, 5)` when that current is positive and `(0, 1)` otherwise; the readout is the mean final load at port zero over the 264 paths; `1/6` with the selector disconnected or independent | `model/observer/native.py` | `reduction-criterion` |

The cited lines are recorded in `data/sources.json`. The status `not-reproduced`
stays in the vocabulary for a quantity whose protocol no source states; no entry
carries it.

## Four cautions that shape the results

1. **Disambiguate "free energy".** Minimum interface energy, the unnormalised
   record potential and normalised record surprisal have different meanings.
   The interior Gaussian determinant contributes to the raw potential and
   cancels in the normalised probability. The `free-energy-name-collision`
   probe reports all three values. Two chains with equal Schur complements
   and the same record prior give equal normalised record densities and
   different interior posterior divergences under the same record change.
2. **Respect OPH's premise discipline.** On the total-load-two shell the ungated
   native rule is confluent without termination. Gating at broken seams gives
   termination with several normal forms reachable from each doubled state.
   The separate Lean copy-repair example does not characterize the native
   rule. The `pr01-obstruction` probe distinguishes constant pathwise load
   energy from the binary ensemble's nonincreasing relative entropy.
3. **Specify the objects being counted.** The IVM shell count, Goldberg face
   count, and geodesic vertex count each give `10k^2 + 2` on the stated
   subdivision sequence. The equality is combinatorial. The `quadray-ivm`
   entry uses it as a complementarity.
4. **Keep scale separations explicit.** In OPH a patch is a bounded local
   system with a specified state, access boundary, records and operations. Its
   documents assign a patch no length scale, and a graph vertex is a model
   representation of one (`docs/INTERPRETATION_AND_TERMS.md` in the OPH
   repository). The twelve ports of this package are the boundary ports of one
   local carrier. An FEP blanket is a conditional-independence structure in a
   chosen probability model. The probes compare a finite carrier graph and
   precisions constructed from it. A correspondence between OPH patch
   boundaries and FEP blankets or a continuum extension requires further
   assumptions.

## Domain limits on the isomorphism claims

The registry records the domain of each conditional identity. The following
examples spell out the assumptions that make the Gaussian and repair identities
hold.

- `variational-equivalence` holds on the carrier with a Gaussian record prior.
  A clamped, nonquadratic seam gives a posterior mean different from the
  energy minimiser.
- `free-energy-dictionary` holds in the same corner, where the surprisal is
  quadratic and the Bayesian surprise is the residual quadratic form. The
  grounded asymmetric seam gives different surprisals to opposite records of
  equal magnitude, outside that centred quadratic formula. The
  `two-record-resistance` probe verifies the Schur complement for every carrier
  record pair against a rational network solve.
- `repair-is-conditioning` holds on the binary shell with the uniform reference
  and a fair coin on the selected seam. The probe verifies 60 failures across
  doubled-state and incident-seam pairs when the reference is uniform on all
  78 configurations.
- `binary-conditional-variance` uses the same binary reference and fair coin
  and samples seams uniformly. On this domain the Dirichlet form of the
  random-seam repair kernel is the seam average of the conditional variances.
- `repair-contrast-quotient` applies to the two specified static graphs.
  Its readout is full intrinsic Gaussian contrast variance at unit seam
  precision. The FCC first shell excludes the central site.
- `reciprocal-dual-weights` uses matched signed edge coordinates and the
  conductance-normalised covariances of intrinsic Gaussian edge differences
  on spherical dual pairs. The cut-and-cycle dimension argument uses the
  spherical Euler characteristic.
- `quadratic-even-field` uses unit vertex directions and the centred
  Laplacian Gaussian. In the plane calculation, `central-plane-screening`
  grounds the plane ports for the intrinsic comparison so the free half-field
  coordinates carry a Dirichlet density.
- `surface-harmonic-covariance` uses closed orientable cellulations with
  simple face boundaries, distinct incident faces per edge and connected
  primal and dual graphs. Positive reciprocal weights leave a harmonic
  part of dimension twice the genus.
- `linear-lag-closure` concerns stationary second moments of a fixed binary
  load under positive seam-sampling weights and fair coins. A Gaussian match
  of all lag covariances can coexist with observed non-Markov repair.
- `jitterbug-coordinate-closure` concerns stationary second moments on distinct
  joints along the contraction branch, with unit bars and a common tuned
  matching conductance. The unit-hull and inverse-square examples in
  `jitterbug-weight-policy-obstruction` show where other rules lose lag closure.
- `jitterbug-octahedral-gaussian-limit` uses common pair values for a centred
  scalar Gaussian with unit bars and infinite matching conductance. This
  static reduction adds two parallel bars per octahedral edge.
- `jitterbug-averaged-repair-limit` uses fair repair of two binary loads and
  accelerates the clock by total attempt weight. Fast matching repair mixes
  within each pair-count fibre. The finite-weight stationary obstruction is
  recorded separately in `jitterbug-pair-load-memory`.
- `centered-gaussian-screen-obstruction` compares a grounded Gaussian with
  the same free field conditioned on an exact or noisy global sum. Screening
  under ordinary clamping alone does not carry to the constrained field.

The general variational identity separates divergence from log evidence when
the expectations are defined and the evidence is positive. The quadratic
record formula and the agreement of the posterior mean and energy minimiser
depend on the Gaussian assumptions above. The `exact` entries also retain
their stated domains. For example, a finite enumeration certifies the named
carrier or shell, and a lattice metric comparison fixes the covariance scale.
