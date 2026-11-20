# Method: executable claims and their scope {#sec:method}

## The unit of work

The paper's unit is a **correspondence**: a claim about the same object made in
two or three vocabularies, together with executable evidence and an explicit
pairwise scope.

::: {#def:correspondence .definition}
A **correspondence** is a record $c = (\kappa, \rho, \sigma, \pi)$. The claim
$\kappa$ states the identification, obstruction or complementarity. The scope
$\rho$ names the programmes it relates, one of `oph-fep`, `fep-synergetics`,
`synergetics-oph` or `three-way`. The status $\sigma$ describes the claim's
domain or outcome. The executable probe $\pi$ computes named values and a
pass verdict from the checks specified for that entry.
:::

The same carrier supports comparisons in different scopes. Its icosahedral
symmetry connects OPH and Synergetics, while the conditional independence of
its Gaussian precision connects OPH and the FEP. Use
`oph-fep-synergetics verify --pair synergetics-oph` to select one scope, or
`oph-fep-synergetics verify --kind disagreement` to select the obstructions.

Named values let readers inspect the evidence behind a verdict and identify
which result moved when an implementation changes. The proof or derivation in
the text states the mathematical domain. The probe tests specific instances or
exhausts a declared finite domain. A passing obstruction probe reproduces the
counterexample to an identification.

| Status | Meaning |
|:-------|:--------|
| `exact` | the probe decides the claim on the declared finite object |
| `exact-on-domain` | the probe supports the identity on the explicitly named domain |
| `obstructed` | the probe exhibits the failure described by the claim |
| `not-reproduced` | the source states a result whose derivation this package has yet to reproduce |

: The status vocabulary. It is closed: a new status requires a new definition.
{#tbl:status}

## The registry, grouped by pairwise scope

[@tbl:registry] is the registry, reproduced from
`src/oph_fep_synergetics/analysis/registry.py`. The grouping is the paper's
argument in miniature: {{REGISTRY_ISOMORPHISMS}} claims survive as
identifications, {{REGISTRY_DISAGREEMENTS}} fail as identifications while
remaining informative, and {{REGISTRY_COMPLEMENTARITIES}} are complements in
which one program supplies what another lacks without either reducing to the
other. By scope, {{REGISTRY_OPH_FEP}} claims relate OPH to the FEP,
{{REGISTRY_SYNERGETICS_OPH}} relate Synergetics to OPH,
{{REGISTRY_FEP_SYNERGETICS}} relate the FEP to Synergetics, and
{{REGISTRY_THREE_WAY}} entries need all three vocabularies at once.

| Key | Kind | Status | Pair | What it tests |
|:-----------------------------------|:----------------|:-------------------|:-----------------|:-------------------------------|
| `carrier-port-chart` | isomorphism | `exact` | `synergetics-oph` | the triangles of the derived graph reproduce the vertex sets of the oriented faces |
| `spectral-gap` | isomorphism | `exact` | `oph-fep` | one matrix read two ways (mean-repair generator and precision): $\lambda_1 = 5-\sqrt{5}$, normalised gap $(5-\sqrt{5})/60$ |
| `shell-counts` | isomorphism | `exact` | `oph-fep` | the total-two load shell and its binary sub-shell, enumerated and checked against closed forms — the finite domain both oph-fep readings evaluate on |
| `native-repair-law` | isomorphism | `exact` | `oph-fep` | OPH's integer rule reimplemented, with its closed-form energy drop |
| `variational-equivalence` | isomorphism | `exact-on-domain` | `oph-fep` | harmonic extension, constrained (KKT) minimiser and perception flow coincide; free energy at the posterior equals the surprisal |
| `free-energy-dictionary` | isomorphism | `exact-on-domain` | `oph-fep` | Gaussian carrier prior gives normalised surprisal, and the residual quadratic equals the KL between posteriors |
| `sheaf-agreement` | isomorphism | `exact` | `oph-fep` | the cellular sheaf Laplacian is the carrier Laplacian |
| `certificate-species` | isomorphism | `exact` | `oph-fep` | the certificate species of the shell and its spectral gap |
| `tree-number` | isomorphism | `exact` | `oph-fep` | one determinant, two readings: {{CARRIER_TREE_NUMBER}} counts spanning trees and normalises the FEP's Gaussian interface integral |
| `gmrf-independence` | isomorphism | `exact` | `oph-fep` | off-diagonal entries of the proper precision give dependence on seams and blanket size {{GMRF_BLANKET_SIZE}}; centred intrinsic variance {{GMRF_MARGINAL_VARIANCE}} |
| `sheaf-cohomology` | isomorphism | `exact` | `oph-fep` | patch agreement is the global-section condition: $\dim H^0 = {{SHEAF_H0}}$, $\dim H^1 = {{SHEAF_H1}}$ on the graph; with the faces, Betti numbers {{SHEAF_SURFACE_BETTI}} |
| `icosahedral-symmetry` | isomorphism | `exact` | `synergetics-oph` | order {{SYMMETRY_GROUP_ORDER}}, transitive on ports, seams and faces, with {{SYMMETRY_ROTATIONS}} rotations and {{SYMMETRY_REFLECTIONS}} orientation-reversing symmetries |
| `eigenspace-irreducibility` | isomorphism | `exact` | `synergetics-oph` | the eigenspaces are the group's non-trivial real irreducibles: dimensions {{SYMMETRY_EIGENSPACE_DIMENSIONS}} forced by the group and established by the commutant test |
| `law-equivariance` | isomorphism | `exact` | `synergetics-oph` | every transition commutes with every automorphism; {{SYMMETRY_ORBITS}} orbits, {{SYMMETRY_OUTCOME_MULTISETS}} outcome multisets |
| `kl-lattice-metric` | isomorphism | `exact` | `fep-synergetics` | at unit covariance, the divergence of a neighbour displacement is four nats and a displacement of squared norm $8k^2$ gives $4k^2$ |
| `pr01-obstruction` | disagreement | `obstructed` | `oph-fep` | ungated repair does not terminate; the gate can give different normal forms |
| `free-energy-name-collision` | disagreement | `obstructed` | `oph-fep` | equal record marginals can carry different interior posterior responses |
| `mediator-obstruction` | disagreement | `obstructed` | `oph-fep` | pairwise boundary reduction needs a threshold |
| `stationary-blanket` | disagreement | `obstructed` | `oph-fep` | a stationary blanket can coexist with direct exterior-to-interior drift |
| `solenoidal-freedom` | disagreement | `obstructed` | `oph-fep` | stationary blankets leave entropy production free |
| `d-stability` | disagreement | `obstructed` | `oph-fep` | stability is not generic in the coupling |
| `schedule-path` | disagreement | `obstructed` | `oph-fep` | iteration schedules are not interchangeable |
| `reduction-criterion` | disagreement | `obstructed` | `oph-fep` | competing reduction criteria disagree |
| `existence-without-reachability` | disagreement | `obstructed` | `oph-fep` | existence of an equilibrium does not imply reachability |
| `directed-gradient` | disagreement | `obstructed` | `oph-fep` | direct and adjoint gradients differ |
| `four-cycle-local-minimum` | disagreement | `obstructed` | `oph-fep` | the four-cycle has spurious local minima |
| `cuboctahedron-vs-icosahedron` | disagreement | `obstructed` | `synergetics-oph` | twelve vertices is not one solid: tetravolume {{CUBOCTAHEDRON_TETRAVOLUME}} against ${{ICOSAHEDRON_TETRAVOLUME}}\ldots$ |
| `shell-vs-free-energy` | disagreement | `obstructed` | `fep-synergetics` | the coordination shell (a polyhedral norm) and the divergence (a Euclidean one) determine neither each other; levels {{SHELL_SHARED_FREE_ENERGY_LEVELS}} are shared |
| `quadray-ivm` | complementarity | `exact` | `synergetics-oph` | exact tetravolumes from vertex sets and tetrahedral and octahedral FCC hole classes from quadray residues |
| `equilibrium-propagation` | complementarity | `exact` | `oph-fep` | the energy-based learning identity holds exactly |
| `replicator-potential` | complementarity | `exact` | `oph-fep` | symmetric payoffs give a Shahshahani gradient of half the mean payoff; the tested asymmetric payoff breaks that identity |
| `kissing-twelve` | complementarity | `exact` | `synergetics-oph` | both twelve-around-one arrangements attain the kissing number, with outer contacts in the cuboctahedron and gaps in the icosahedron |
| `icosahedron-tetravolume` | complementarity | `exact` | `synergetics-oph` | the carrier's solid is ${{ICOSAHEDRON_TETRAVOLUME}}\ldots$, squared exactly $175 + 75\sqrt{5}$ |
| `packing-density` | complementarity | `exact` | `synergetics-oph` | Voronoi cell ${{VORONOI_VERTICES}}/{{VORONOI_FACES}}/{{VORONOI_EDGES}}$, tetravolume six, density $\pi/(3\sqrt{2})$ |
| `tetravolume-arithmetic` | complementarity | `exact` | `three-way` | integer nats, a surd gap, and an integral table beside a non-integral carrier |
| `effective-resistance` | complementarity | `exact` | `oph-fep` | effective resistance {{RESISTANCE_ADJACENT}}, {{RESISTANCE_DISTANCE_TWO}}, {{RESISTANCE_ANTIPODAL}} by distance class, Kirchhoff index {{RESISTANCE_KIRCHHOFF_INDEX}} twice over |
| `two-record-resistance` | isomorphism | `exact` | `oph-fep` | the inverse effective resistance sets the Schur energy of a pair's record contrast |
| `vertex-shells` | disagreement | `obstructed` | `synergetics-oph` | shells {{CARRIER_SHELL_SIZES}} against the lattice's {{LATTICE_SHELL_SIZES}}; intersection array {{CARRIER_INTERSECTION_ARRAY}} |
| `subgraph-census` | isomorphism | `exact` | `synergetics-oph` | {{CARRIER_TRIANGLES}} triangles, {{CARRIER_PENTAGONS}} induced five-cycles, {{CARRIER_PERFECT_MATCHINGS}} perfect matchings and {{CARRIER_HAMILTONIAN_CYCLES}} Hamiltonian cycles |
| `angle-defect` | isomorphism | `exact` | `synergetics-oph` | OPH's derived chart is the icosahedron's triangulation, so the angular-defect accounting applies unchanged: $\pi/3$ at each port, $4\pi$ in total, which is $2\pi\chi$ |
| `twelve-around-one-families` | disagreement | `obstructed` | `synergetics-oph` | the phrase covers solids with {{CARRIER_TREE_NUMBER}} and ${{CUBOCTAHEDRON_SPANNING_TREES}} = {{CUBOCTAHEDRON_SPANNING_TREES_ROOT}}^2$ spanning trees |
| `dual-cells` | complementarity | `exact` | `synergetics-oph` | the cell and the cuboctahedron, ${{VORONOI_VERTICES}}/{{VORONOI_FACES}}/{{VORONOI_EDGES}}$ against ${{CUBOCTAHEDRON_VERTICES}}/{{CUBOCTAHEDRON_FACES}}/{{CUBOCTAHEDRON_EDGES}}$, volumes {{PACKING_CELL_TETRAVOLUME}} and {{CUBOCTAHEDRON_TETRAVOLUME}} |
| `conjugation-norm` | complementarity | `exact` | `three-way` | conjugation fixes the integer counts and moves the volume, whose norm is {{TETRAVOLUME_GALOIS_NORM}} |
| `repair-is-conditioning` | isomorphism | `exact-on-domain` | `oph-fep` | binary repair conditions the uniform binary reference and minimises free energy over the fibre. Conditioning the uniform total-two reference fails in {{CONDITIONING_OFF_SHELL_FAILURES}} doubled-state/seam cases |
| `binary-conditional-variance` | isomorphism | `exact-on-domain` | `oph-fep` | the binary repair Dirichlet form equals mean conditional variance over seam records |
| `carrier-spectral-chords` | isomorphism | `exact` | `fep-synergetics` | squared chords equal {{CHORD_VARIANCE_FACTOR}} times the lowest-mode contrast variance, with different embedding ranks for the full variance and the chords |
| `repair-first-entry` | isomorphism | `exact-on-domain` | `oph-fep` | doubled starts have a geometric first binary entry, with hazard {{ENTRY_HAZARD}}. Entry and stationarity have the same one-port marginals |
| `repair-distance-lumping` | isomorphism | `exact-on-domain` | `oph-fep` | graph distance gives an exact Markov chain with class sizes {{DISTANCE_CLASS_SIZES}} |
| `grounded-forest-precision` | isomorphism | `exact-on-domain` | `oph-fep` | grounded precision determinants count rooted forests and agree with contrast covariance minors on all record pairs and triples |
| `greedy-record-placement` | disagreement | `obstructed` | `oph-fep` | the best record pair cannot extend to a best triple under the conditional entropy criterion |
| `planar-dual-covariance` | isomorphism | `exact-on-domain` | `three-way` | matched primal and dual Gaussian edge covariances sum to the identity on both dual pairs |
| `cuboctahedral-spectral-chords` | isomorphism | `exact-on-domain` | `fep-synergetics` | on the FCC shell, squared chords equal {{FCC_CHORD_FACTOR}} times the lowest-mode contribution |
| `cuboctahedral-distance-split` | disagreement | `obstructed` | `fep-synergetics` | two pair types at the same graph distance have different full Gaussian contrast variances |
| `repair-contrast-quotient` | isomorphism | `exact-on-domain` | `three-way` | Gaussian contrast classes give an exact binary repair quotient on both shells |
| `repair-distance-counterexample` | disagreement | `obstructed` | `three-way` | merging the FCC middle classes leaves different transition rates for the same distance record |
| `reciprocal-dual-weights` | isomorphism | `exact-on-domain` | `three-way` | weighted dual projectors sum to identity at constant matched weight product |
| `copied-dual-weight-obstruction` | disagreement | `obstructed` | `three-way` | copying a nonuniform conductance pattern onto a dual breaks complementarity |
| `quadratic-even-field` | isomorphism | `exact-on-domain` | `three-way` | trace-free quadratics parameterise centred even fields |
| `quadratic-isotropy-obstruction` | disagreement | `obstructed` | `three-way` | equal linear isotropy leaves different quadratic variances on the FCC shell |
| `central-plane-screening` | disagreement | `obstructed` | `three-way` | cuboctahedral equators can screen while icosahedral equators leave dependence |
| `repair-centroid-memory` | disagreement | `obstructed` | `three-way` | matching Gaussian lag covariances leaves non-Markov centroid probabilities |
| `distance-regular-repair` | isomorphism | `exact-on-domain` | `three-way` | distance-regular graphs give the exact binary distance quotient, with variance as an equivalent record |
| `bipartite-contrast-collision` | disagreement | `obstructed` | `three-way` | equal resistance on a complete bipartite family can hide different one-step repair rates |
| `periodic-fcc-contrast-collision` | disagreement | `obstructed` | `three-way` | resistance, distance and common-neighbour counts fail to determine a repair quotient on a regular FCC graph |
| `linear-lag-closure` | isomorphism | `exact-on-domain` | `three-way` | invariance of the linear observation space is necessary and sufficient for Gaussian Markov covariance closure |
| `single-seam-lag-obstruction` | disagreement | `obstructed` | `three-way` | changing a single seam weight on the twenty-four-cell yields a rank-one lag-two defect |
| `surface-harmonic-covariance` | isomorphism | `exact-on-domain` | `three-way` | reciprocal edge weights leave a harmonic projector of rank twice the genus |
| `copied-surface-weight-obstruction` | disagreement | `obstructed` | `three-way` | copied weights can preserve the torus trace and break positivity |
| `root-shell-quadratic-precision` | isomorphism | `exact-on-domain` | `three-way` | the quadratic covariance is isotropic only in the four-dimensional member of the root-shell family |
| `fcc-plane-screen-transfer` | isomorphism | `exact-on-domain` | `three-way` | parity gives separating planes on FCC boundaries at all frequencies |
| `centered-gaussian-screen-obstruction` | disagreement | `obstructed` | `three-way` | conditioning on a total-field sum couples the strict sides of a seam cut |

: The {{REGISTRY_ENTRIES}} correspondences, in registry order:
{{REGISTRY_ISOMORPHISMS}} isomorphisms, {{REGISTRY_DISAGREEMENTS}}
disagreements and {{REGISTRY_COMPLEMENTARITIES}} complementarities, scoped
`oph-fep` {{REGISTRY_OPH_FEP}}, `synergetics-oph`
{{REGISTRY_SYNERGETICS_OPH}}, `fep-synergetics`
{{REGISTRY_FEP_SYNERGETICS}}, `three-way` {{REGISTRY_THREE_WAY}}. Every row is
executed by `oph-fep-synergetics verify`.
{#tbl:registry}

The registry is the source of truth for the paper. The prose in
[@sec:correspondences], [@sec:disagreements] and [@sec:complementarities] says
what the probes found; it does not add claims the probes do not carry.
[@fig:scope-matrix] summarizes the same {{REGISTRY_ENTRIES}} entries by
scope, kind and evidence status. [@tbl:registry] lists individual entries.

![Correspondence counts by scope. Each row aligns counts by claim kind with
counts by evidence status and a shared scope total. The matrix labels zeros
explicitly and reads from the executable registry.](../output/figures/fig_scope_matrix.png){#fig:scope-matrix width=100%}

## Reading a source without copying it

The three source programs are read at pinned commits, and the pins are data, not
prose: `data/sources.json` records the commit, the branch, the licence terms and,
for each borrowed identifier, the file, the line range and the symbol. The
command `oph-fep-synergetics sources --vault-root <root>` re-reads the live
checkouts and reports `match`, `moved`, `absent` or `unavailable` for each pin.

OPH applies different licences to different material. Its software, including
the Lean development and model code, is Apache-2.0. Its papers, book,
documentation, figures, data and ledgers are CC BY-NC-SA 4.0 [@oph]. This
package reimplements the software behaviour with attribution. It derives the
carrier adjacency from the vertex triples in OPH's software and checks the
resulting triangles against those triples. The computations use the package's
own code and inputs.

FEP_Formal is CC BY 4.0 and is cited by DOI and by topic identifier
[@fepformal]. The catalogue is read through its machine-readable metadata
(`config/catalogue_metadata.yaml`), which is the schema-2 roster sealed at
`fep-001` through `fep-{{FEP_FORMAL_TOPICS}}`; the generated projection
`config/topics.yaml` is never authored or edited, because FEP_Formal declares it
a projection and regenerates it.

QuadMath is Apache-2.0 [@quadmath]. The constants this package needs are
reimplemented and the source is cited by symbol: the neighbour moves, the
embedding, the shell formula and the tetravolume convention. This paper quotes no
QuadMath prose.

The provenance record identifies the OPH software containing the oriented face
list and the QuadMath revision examined for the lattice derivations. The
recorded revision gives readers the version to consult if a later checkout
moves. See `docs/PROVENANCE.md`.

## What a matching value establishes

An equilibrium calculation and a trajectory can agree on the location of a
fixed point while differing on whether the trajectory reaches it
([@prop:existence]). The registry records the fixed point, the dynamical
evidence and the scope of the conclusion separately.
