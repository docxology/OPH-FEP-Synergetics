# Scope, limits and falsifiers of a correspondence battery {#sec:scope}

## What is verified

Everything in [@sec:correspondences], [@sec:disagreements] and
[@sec:complementarities] is executed by `oph-fep-synergetics verify`, which runs
{{REGISTRY_ENTRIES}} probes and exits non-zero if any of them fails. The command
can be narrowed to one pairwise scope with `--pair`, so a reader who only wants
the Synergetics–OPH ledger does not have to take the rest on trust. The values
quoted in the prose are the values the probes report; a change in any of them
changes the probe's output and, where the change is a failure, changes the exit
code.

Three limits of that verification should be stated, and a fourth applies to the
newer entries specifically.

The first is domain. The variational equivalence [@prop:ve], the dictionary
[@sec:dictionary] and the two FEP identities hold in the linear-Gaussian corner.
The results are exact there and are not claimed elsewhere. Where a probe can test
the boundary of the domain — the `exact-on-domain` entries — it does, and the
status vocabulary records the difference.

The second is size. The carrier has twelve ports, and the shell census is
exhaustive over a set of {{SHELL_TOTAL_TWO}} configurations. Both are small enough
to enumerate and too small to settle any question about scaling. Nothing in this
paper claims that the identities survive a continuum limit, and the obstruction
results in particular are statements about this carrier.

The third is provenance. The source for the dynamics is a private draft
[@mueller2026agreement], not a published paper. Its propositions are cited as the
source states them, the package's own measurements are reported alongside, and a
reader who disagrees with the draft can still use the measurements. The three
software programmes are pinned to commits, and the pins are point-in-time
readings rather than perpetual guarantees.

The fourth is arithmetic. The newer identifications fall into two groups, and the
probes treat them differently. The counting identities are exact: the
spanning-tree number is an integer determinant computed over the rationals, the
cohomology dimensions are ranks, the group order and its rotation/reflection
split are integers, and the equivariance claim is a finite search over
{{SYMMETRY_GROUP_ORDER}} group elements and {{SYMMETRY_CONFIGURATIONS}}
configurations. The geometric claims are not of that kind: tetravolumes and the
packing density are floating-point measurements compared against closed forms,
and the probes report the residual rather than only a verdict. The statuses are
`exact` in both groups because the residuals are at the level of double
precision, and a reader who wants a stronger statement should read the residual,
which is reported. The cleanest statement of the split is the conjugation of
$\sqrt{5}$, which [@sec:conjugation] turns into a probe: it fixes the integer
counts and moves the volume to {{ICOSAHEDRON_TETRAVOLUME_CONJUGATE}}, and the two
volumes multiply to the exact norm {{TETRAVOLUME_GALOIS_NORM}}. A counting claim
in this paper is therefore invariant under a symmetry that a geometric claim is
not, and comparing one against the other without saying so is the error the
section is written to prevent.

## What is not reproduced

Two quantities stated in the source literature are not reproduced here, and the
registry records them as `not-reproduced` rather than dropping them or inventing
a protocol that would produce them.

| Quantity | Source | Why it is not reproduced |
|:---------|:-------|:-------------------------|
| `35/264` port-zero occupation | [@mueller2026agreement] | the record and selector protocol that would determine the occupation is not specified in the paper's text |
| `13/15` and `1/60` move probabilities | [@mueller2026agreement] | the probability space and the conditioning event are not determined by the text |

: Quantities recorded as not reproduced, with the reason. Recording them keeps the
gap visible; inventing a protocol would have produced numbers that look like
verification and are not. {#tbl:notrepro}

Both quantities are plausibly derivable from a fuller specification of the
record process. Neither is asserted anywhere in this paper, and neither appears
in a probe. If a later version of the source specifies the protocol, the
correspondence can be upgraded to `exact` by adding a probe; until then the honest
status is `not-reproduced`.

## What is not claimed

This paper does not claim that Observer Patch Holography is a free energy
principle. It does not claim that the Free Energy Principle is a holographic
theory. It does not claim that Synergetics predicts either. It does not claim that
any of the three is empirically adequate; OPH's own empirical status document is
explicit that no retained laboratory comparison selects a native record law, and
this paper inherits that position rather than improving on it [@oph].

It also does not claim novelty for the mathematics. The spectral theory of
cellular sheaves, the Routh criterion, the Markov-blanket critique, and the
Goldberg construction are all prior work; what is new here is the *conjunction* of
them on one finite object, and the discipline of reporting which conjunctions
fail.

## Four cautions for anyone continuing this work

The source dossier states four cautions, and this paper adopts them
[@dossier2026].

**Disambiguate "free energy".** At least two inequivalent functionals carry the
name, and they differ by a model constant [@sec:collision]. A sentence that does
not say which one it means is not yet a claim.

**Respect premise discipline.** OPH records its obstructed premises, including
PR-01, and a synthesis that assumes the premise it needs is doing the opposite of
what the programme asks. The census in [@prop:census] exists to make that
temptation visible.

**Distinguish numerology from prediction.** The twelve-fold and sixty-fold
coincidences in this paper are derived, and the derivation is part of the
probe. A numerical coincidence that is not derived is not evidence, and the
Synergetics literature contains both kinds. The two newest obstructions are the
sharpest examples: twelve vertices are shared by a cuboctahedron, a rhombic
dodecahedron's neighbour graph and an icosahedron, and only one of the three is
the carrier [@sec:cuboctahedron].

**Respect scale separation.** The linear-Gaussian results are exact at finite
size and say nothing about a continuum limit. A claim that survives at twelve
ports may fail at a thousand, and the paper's obstruction results are the ones
most likely to change with scale.

## What would change the conclusions

The conclusions are falsifiable in a specific sense, and it is worth saying how.

A single strict descent at $\lvert d\rvert = 1$ would disprove [@prop:census] and remove the
PR-01 obstruction. A pairwise representation below $\varepsilon = 1/3$ would
disprove [@prop:mediator]. A member of the solenoidal family with a different
entropy production would disprove [@prop:solenoid]. A stable $D_kA$ at
$k > \sqrt{10}$ would disprove [@prop:dstability]. Each of these is a computation,
not an interpretation, and each is executable against the package as it stands.

That is the property this paper was built for: the disagreements are as checkable
as the agreements, and either can be overturned by running a function.
