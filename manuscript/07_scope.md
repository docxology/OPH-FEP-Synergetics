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

Four limits of that verification should be stated.

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

The third is provenance. The source for the dynamical controls is a preprint
[@mueller2026agreement], and OPH's own results are stated in its synthesis paper
[@mueller2026observers]. The preprint's propositions are cited as the source
states them, the package's own measurements are reported alongside, and a reader
who disagrees with the preprint can still use the measurements. The three
software programmes are pinned to commits, and the pins are point-in-time
readings rather than perpetual guarantees.

The fourth is arithmetic. The counting identities are exact: the spanning-tree
number is an integer determinant computed over the rationals, the cohomology
dimensions are exact ranks, the reciprocal spectrum sum comes from an exact
rational inverse, the group order and its rotation/reflection split are
integers, and the equivariance claim is a finite search over
{{SYMMETRY_GROUP_ORDER}} group elements and {{SYMMETRY_CONFIGURATIONS}}
configurations. The geometric claims are exact too wherever the geometry allows:
the tetravolume table, the cell's volume and the icosahedron's tetravolume are
computed in the rationals and in $\mathbb{Q}(\sqrt{5})$ from exact vertex sets,
and the floating-point measurements are reported beside them as a second route.
What stays floating-point is what has to: the eigenvalue locus, the relaxation
orbit and the free-energy values of the linear-Gaussian model, whose probes
report residuals and use tolerances stated in the code. The cleanest statement of
the two kinds of claim is the conjugation of $\sqrt{5}$, which
[@sec:conjugation] turns into a probe: it fixes the integer counts and moves the
volume to {{ICOSAHEDRON_TETRAVOLUME_CONJUGATE}}, and the two volumes multiply to
exactly {{TETRAVOLUME_GALOIS_NORM}}. A counting claim in this paper is therefore
invariant under a symmetry that a geometric claim is not, and comparing one
against the other without saying so is the error the section is written to
prevent.

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
Synergetics literature contains both kinds. The counting obstructions are the
sharpest examples: twelve-fold features are shared by a cuboctahedron, a
rhombic dodecahedron (twelve faces) and an icosahedron, and only one of the three
is the carrier [@sec:cuboctahedron].

**Respect scale separation.** The linear-Gaussian results are exact at finite
size and say nothing about a continuum limit. A claim that survives at twelve
ports may fail at a thousand, and the paper's obstruction results are the ones
most likely to change with scale.

## What would change the conclusions

The conclusions are falsifiable in a specific sense, and it is worth saying how.

A single strict descent at $\lvert d\rvert = 1$ would disprove [@prop:census] and remove the
PR-01 obstruction. A pairwise representation below $\varepsilon = {{MEDIATOR_THRESHOLD}}$ would
disprove [@prop:mediator]. A member of the solenoidal family with a different
entropy production would disprove [@prop:solenoid]. A stable $D_kA$ at
$k > \sqrt{10}$ would disprove [@prop:dstability]. A lattice site whose
breadth-first shell differs from its cuboctahedral gauge would disprove
[@def:coordination], and a transition that fails to commute with an automorphism
would disprove [@prop:symmetry]. Each of these is a computation, not an
interpretation, and each is executable against the package as it stands.

That is the property this paper was built for: the disagreements are as checkable
as the agreements, and either can be overturned by running a function.
