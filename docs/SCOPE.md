# Scope: what is reimplemented, what is cited, what is out

The point of this repository is to make three-way identifications *decidable*.
That requires saying clearly where each identification is being tested, what
domain it is being tested on, and what is deliberately left outside. This
document is part of the deliverable: a claim is only as good as its stated
domain.

## Tested here

| subject | how |
| --- | --- |
| OPH's twelve-port carrier | the port chart is *derived* from OPH's `ORIENTED_BASE_FACES`; the 20 triangles of the derived graph are asserted to be exactly those faces |
| OPH's native seam law | reimplemented on exact Python integers; conservation and the energy identity are checked on every one of the 4680 shell transitions |
| PR-01 | the obstruction is reproduced exhaustively: 120 strict descents, 600 zero-descent mismatch-one moves, 3960 no-ops |
| patch agreement as a sheaf condition | a cellular sheaf over the carrier whose Laplacian is asserted equal to the carrier Laplacian |
| the free-energy dictionary | the linear-Gaussian corner, exactly: harmonic extension, a separate KKT constrained minimiser, Schur complement, the three surprises, and the variational free energy computed both as surprisal plus divergence and directly as energy minus entropy, equal at the posterior and separated by the divergence away from it |
| obstruction and dynamics controls | the D-stability family, the mediator reduction, the reduction criterion, the two-node orbits, the solenoidal family, the four-cycle trap, the replicator potential |
| QuadMath's conventions | quadray gauge and embedding, the twelve neighbour moves, the coordination shells against their closed-form cuboctahedral gauge, and the tetravolume table derived exactly from explicit vertex sets |
| icosahedral subdivision | geodesic and Goldberg counts, duality, and the equality with the IVM shell sequence |
| the carrier's counting identities | the reduced Laplacian determinant in exact rationals (the spanning-tree number), its cross-check against the spectral product, the graph sheaf's two cohomology dimensions, and the surface's Betti numbers with its Euler characteristic |
| the carrier's symmetry | the automorphism group by exhaustive search, closure under composition, the orbit decomposition of ports, seams and faces, the rotation/reflection split from OPH's oriented faces, eigenspace irreducibility, and the equivariance of every single transition over the full group |
| the Gaussian reading of the same matrix | the precision's zeros as the conditional-independence pattern, the Markov blanket of every port, the pairwise precisions, and the exact variance of a port's deviation from the carrier mean (the diagonal of `L^+`) |
| the information cost of a lattice step | equal-covariance Gaussian divergence against the IVM's squared distances at unit covariance: one step in nats, the sequence `4k^2`, and the stated dependence on the unit of length |
| the carrier's solid and the twelve-around-one cell | the icosahedron's tetravolume computed exactly in `sqrt 2 * Q(sqrt 5)` and by floating-point measurement, the two kissing arrangements compared exactly, and the Voronoi cell's counts, exact tetravolume and packing density |
| the carrier as a network | exact effective resistances as rationals from a grounded solve, their distance-class ladder, the Kirchhoff index, the distance distribution and intersection array around a port, and the subgraph census |
| the carrier's curvature | the angle defect at every vertex of the derived chart, summed against Gauss-Bonnet, plus the defect of the dual solid's face set |
| the twelve-around-one families | the obstruction is exhibited rather than asserted: the cuboctahedron, the icosahedron and the cube's line graph share the twelve-around-one description and differ in counts and volumes |
| the quadratic field of the tetravolume | the Galois conjugation of $\sqrt{5}$ (fixing $\sqrt{2}$) as a probe, in exact field arithmetic: it fixes the counting invariants and moves the icosahedron's tetravolume, whose product with its conjugate is exactly 50 |

## Cited but not implemented

- **OPH's quantum event algebra and its continuum limit.** Nothing here
  evaluates them, and `data/sources.json` records no identifier for either.
- **FEP beyond the linear-Gaussian corner.** Nonlinear generative models,
  non-Gaussian densities, multi-modal posteriors, and the path-integral and
  stochastic-thermodynamic formulations are outside this package.
- **Fuller's full synergetic apparatus.** Only the pieces the correspondence
  needs are here: quadray coordinates, the IVM, the vector-equilibrium shell,
  the volume ratios, and the subdivision counts.
- **OPH's transactional and integer consensus contracts.** OPH's consensus
  theorems hold for those contracts, not for the native law. This package tests
  the native law, which is the one the carrier actually runs.
- **Equilibrium propagation as a training method.** It appears here as a
  complementarity from energy-based learning, which the FEP reaches through
  predictive coding: it supplies the parameter gradient that the native repair
  law does not define. No learning is implemented, and the probe is a scalar
  energy model.
- **Fuller's jitterbug.** The continuous motion that carries the vector
  equilibrium through the icosahedron is cited to explain why Synergetics
  distinguishes the two twelve-vertex solids; the kissing slack that permits it
  is measured, the motion itself is not simulated.

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

1. **Disambiguate "free energy".** The interface energy minimised over the
   interior and the variational free energy at the posterior are *different
   objects*; identifying them numerically is a gauge choice. This is registered
   as the `free-energy-name-collision` disagreement, and the probe reports both
   values so the constant difference is visible.
2. **Respect OPH's premise discipline.** The native law is not a descent method:
   it accepts repairs that lower nothing, and it is not confluent. This is
   registered as `pr01-obstruction`, and it is why the FEP-side reading of
   "repair" is stated as a *constrained relative-entropy minimisation on a
   domain* rather than as the native law's own guarantee.
3. **Do not let numerology stand in for a frozen prediction.** Agreement between
   `10k^2 + 2`, the IVM shell count, the Goldberg face count and the geodesic
   vertex count is a genuine identity of one quadratic form seen four ways, and
   it is tested as such. It is not evidence that any two of the three programs
   describe the same physics. The `quadray-ivm` entry is registered as a
   *complementarity*, not an isomorphism of the theories.
4. **Keep scale separations explicit.** In OPH a patch is a bounded local
   system with a specified state, access boundary, records and operations. Its
   documents assign a patch no length scale, and a graph vertex is a model
   representation of one (`docs/INTERPRETATION_AND_TERMS.md` in the OPH
   repository). The twelve ports of this package are the boundary ports of one
   local carrier. An FEP blanket is a statistical boundary of whatever system
   the model describes. Nothing in this repository maps one to the other, and
   no probe claims to. The linear-Gaussian results are exact at finite size and
   say nothing about a continuum limit.

## Domain limits on the isomorphism claims

The two `exact-on-domain` entries say so in their status:

- `variational-equivalence` holds for finite linear-Gaussian cellular sheaves
  over the carrier complex.
- `free-energy-dictionary` holds in the same corner, where the surprisal is
  quadratic and the Bayesian surprise is the residual quadratic form.

Outside that corner the correspondence is an open question, not a result. The
remaining `exact` entries are combinatorial or algebraic identities that hold
without a domain restriction.
