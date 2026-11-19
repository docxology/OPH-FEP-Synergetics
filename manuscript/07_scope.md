# Scope and limits {#sec:scope}

## What is verified

`oph-fep-synergetics verify` runs {{REGISTRY_ENTRIES}} probes and exits non-zero
if any fails. The `--pair` option selects a scope. Each probe reports the
computed values behind its verdict.

**Domain.** The relation between variational free energy, posterior divergence,
and evidence holds under the conditions in [@eq:vfe]. The quadratic record
formula, the equality of posterior mean and energy minimiser, and the
perception calculation use the linear-Gaussian model in [@sec:dictionary]
and [@prop:ve]. The conditioning result uses the uniform binary-shell
reference [@prop:conditioning]. Counterexamples probe the boundaries of
those assumptions.

**Size.** The carrier has {{CARRIER_PORTS}} ports, and the repair census is
exhaustive over the {{SHELL_TOTAL_TWO}} configurations of total load two. Other
obstructions use explicit control families. Extension to larger carriers or a
continuum is a separate question.

**Provenance.** *Agreement and Surprise* supplies dynamical controls, and the
OPH synthesis paper supplies context for the native model
[@mueller2026agreement; @mueller2026observers]. The registry distinguishes
source statements from values computed here. Commit pins and source
identifiers record the versions used for software comparisons.

**Arithmetic.** The spanning-tree determinant and grounded resistance solves
use rational arithmetic. The surface cohomology check uses exact incidence
ranks, alongside a tolerance-based numerical rank check for the scalar
sheaf. The symmetry and repair censuses use finite enumeration. Volumes are
derived from explicit coordinates in the rationals or in the quadratic field
$\mathbb{Q}(\sqrt{5})$, with a separate floating-point measurement for the
icosahedron. The eigenvalue locus, orbit and Gaussian free-energy calculations
report floating-point values, tolerances or residuals. The conjugation probe
illustrates the difference between the integer counts and the algebraic
volume [@sec:conjugation].

## Three source quantities and their protocols

The source states three numbers for the native law [@mueller2026agreement]. Each
is determined by a protocol in OPH's public model [@oph], and the registry
reproduces all three exactly ([@tbl:reproduced]).

| Quantity | Protocol | Probe |
|:---------|:---------|:------|
| {{PR01_STAY_PROBABILITY}} and {{PR01_HOP_PROBABILITY}} | one attempt, a uniform seam and a fair coin, from the configuration with one unit on each of two adjacent ports: staying put, and one particular hop | `pr01-obstruction` |
| {{FEEDBACK_PORT_ZERO_OCCUPATION}} | uniform binary source, two repair events, the second seam chosen from the first recorded current: the mean final load at port zero | `reduction-criterion` |

: The source quantities and the protocols that determine them. [@sec:pr01] and
[@sec:mediator] report the values with their controls. {#tbl:reproduced}

The status `not-reproduced` stays in the vocabulary for a quantity whose protocol
no source states. No entry of the registry carries it.

## What the comparisons establish

The results concern relations between the specified mathematical models.
Empirical adequacy remains open. In particular, OPH's empirical status
document leaves the selection of a native record law by laboratory evidence
undetermined [@oph]. The contribution here is to bring established
methods into one executable comparison and report both agreements and
obstructions with their assumptions visible.

## Distinctions used in the comparisons

The source dossier states four cautions, and this paper adopts them
[@dossier2026].

**Free energy.** At fixed records, name the minimum interface
energy, the unnormalised record potential, or the normalised record surprisal
that bounds variational free energy [@sec:collision].

**Repair rules.** OPH's PR-01 distinguishes the native implementation from the
contracts with consensus theorems. The ungated and broken-seam relations
have different reachability properties [@prop:native-rewrite].

**Counting.** The cuboctahedron, rhombic dodecahedron, and icosahedron all have
twelve-fold features. Identifying a solid requires specifying the feature and
its incidence structure. Here, for example, the rhombic dodecahedron has
twelve faces while the carrier has twelve vertices [@sec:cuboctahedron].

**Scale.** The carrier's ports are the boundary coordinates of one local OPH
model. A Markov blanket specifies conditional independence in a distribution.
Extending either notion to a continuum or identifying the two requires
additional assumptions.

## How to challenge a result

A state-changing transition at $\lvert d\rvert = 1$ with nonzero energy drop
would contradict [@prop:census]. The termination and confluence claims have
their own test using the state-changing transition graph
[@prop:native-rewrite]. A feasible pairwise representation of the stated
mediator model below $\varepsilon = {{MEDIATOR_THRESHOLD}}$ would
contradict [@prop:mediator]. A drift in the stated solenoidal family whose
entropy production disagrees with the derived formula would contradict
[@prop:solenoid]. Stability of $D_k A$ at $k \geq \sqrt{10}$ would
contradict [@prop:dstability]. A lattice site whose breadth-first shell
differs from its cuboctahedral gauge would contradict [@def:coordination].
A repair transition that fails to commute with a carrier automorphism would
contradict [@prop:symmetry]. Each check refers to an explicit model, input and
comparison that a reader can recompute.
