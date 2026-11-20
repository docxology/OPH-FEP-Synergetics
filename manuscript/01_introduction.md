# Introduction {#sec:introduction}

## Three programmes, one comparison

Observer Patch Holography (OPH) investigates discrete reconstruction from local
patches. We use its finite carrier and seam repair rule on integer loads
[@oph; @mueller2026observers; @dossier2026]. The Free Energy Principle (FEP)
motivates a variational description of inference and a statistical account of
boundaries. We use a recognition density to bound record surprisal in a specified
generative model [@friston2010; @friston2013life; @parr2022active]. From
Fuller's Synergetics we use the close-packed isotropic vector matrix and
its tetrahedral unit of volume [@fuller1975synergetics]. We work with the
quadray notation documented by Urner and implemented in QuadMath to
count lattice shells and compare solid volumes in this unit
[@urner1997quadrays; @quadmath].

Where two or three programmes refer to comparable finite objects, we state
their readings in mathematical terms and compute the quantities each predicts.
Agreement means equality under the stated correspondence. A counterexample
locates the failure of a proposed identification. A complementarity records
how a construction from one programme helps analyse another.

## What a correspondence is

A correspondence states a claim, identifies the programmes involved, and links
the claim to a probe and status [@def:correspondence]. The probe reports the
computed values. The status records an exact identity, an identity under
explicit assumptions, an obstruction, or a quantity whose source leaves its
reproduction underdetermined [@sec:method]. The results distinguish exact
arithmetic from numerical corroboration at a stated tolerance.

## Contributions

1. **Carrier and repair.** We derive the port chart from the oriented faces
   and check the closed-form spectrum [@prop:carrier; @prop:spectrum]. We
   compute the carrier's symmetry and cohomology [@prop:symmetry;
   @prop:cohomology] and derive network and surface invariants
   [@prop:tree; @prop:resistance; @prop:subgraph-census; @prop:defect]. The
   repair census covers every transition on the shell of total load two.
   Its reachability analysis distinguishes the ungated law from the
   broken-seam gate [@prop:census; @prop:native-rewrite].
2. **Inference on the carrier.** We reproduce the linear-Gaussian
   equivalence from *Agreement and Surprise* by comparing harmonic extension
   and the limit of predictive-coding perception with a constrained KKT solve
   [@mueller2026agreement, theorem 7, items 2–3]. The comparison is stated in
   [@prop:ve]. Two Gaussian chains share the record density and have
   different interior responses to the same record change [@sec:collision].
   Native repair is Bayesian conditioning on the uniform binary-shell
   reference. Its random-seam Dirichlet form is the mean conditional
   variance [@prop:conditioning; @prop:binary-conditional-variance].
3. **Dynamical obstructions.** We check the mediator threshold, the solenoidal
   entropy production, and the stability onset in the stated families
   [@mueller2026agreement, propositions 8–10 and 15]. We report those checks
   in [@prop:mediator; @prop:solenoid; @prop:dstability]. We also enumerate
   the clipping regimes and measure the orbit in the two-node control
   [@prop:existence].
4. **Lattice and volume.** We derive the cuboctahedral coordination-shell
   gauge [@prop:ivm], tetrahedral and octahedral hole classes, and
   tetravolume table.
   We compare equal-covariance Gaussian divergence with lattice distance
   and relate the carrier's Euclidean chords to its lowest Gaussian
   modes [@prop:quadray-ivm; @prop:kl; @prop:carrier-spectral-chords].
   We also compute exact volumes for the carrier and the close-packing
   Voronoi cell [@prop:ico-volume; @prop:packing].
5. **Executable provenance.** The registry links each comparison to its
   probe. The source register records software commits and the identifiers
   used in each reimplementation [@sec:reproducibility].

## Scope of the comparison

The scope is mathematical verification of finite models motivated by the
integration proposals and the controls developed in *Agreement and Surprise*
[@dossier2026; @proposal2026; @mueller2026agreement]. We reproduce the source's
native-law quantities from the protocols in OPH's public model [@sec:scope].
The empirical adequacy of the programmes and the selection of a native record
law by laboratory data remain open [@oph].

## Related work and positioning

The FEP literature distinguishes statistical and dynamical readings of a
Markov blanket. The statistical reading is a conditional-independence structure
[@pearl1988probabilistic]. Dynamical interpretations require further assumptions
[@bruineberg2022emperor; @biehl2021technical]. The conditions under which
variational free energy becomes tight are also central
[@friston2007synthese; @sajid2021demystified]. Bayesian mechanics gives the
stationary-process reading of the same mathematics
[@dacosta2021bayesian; @parr2020blankets], and the discrete-state synthesis
supplies the notation for the active-inference corner
[@dacosta2020discrete].

One use of *screen* in the FEP literature interprets a Markov blanket as an
information boundary [@fields2022generic]. In OPH, a patch has a discrete access
boundary with records and local operations. We compare the carrier graph with
a Gaussian precision constructed on that graph. In [@sec:stationary], an
explicit control tests how much a stationary blanket constrains the drift
[@aguilera2022particular].

Edmondson describes Synergetics' geometric conventions
[@edmondson1987fuller]. The class-I geodesic subdivision of the icosahedron
and its Goldberg dual exchange vertex and face counts
[@goldberg1937multisymmetric; @coxeter1973regular]. The subdivision is the
frequency-squared case of the triangulations used by Caspar and Klug
[@caspar1962physical].

The isotropic vector matrix used here is the face-centred cubic lattice,
which is isometric up to scale to the root lattice $A_3$
[@conway1999sphere]. Its lattice kissing arrangement attains the
three-dimensional kissing number established by Schütte and van der Waerden
[@schutte1953dreizehn]. Its coordination sequence is $10k^{2}+2$
[@conway1997lowdim].

The tools this paper uses have long pedigrees, and it claims none of them: the
algebraic connectivity of a graph [@fiedler1973connectivity], the Routh–Hurwitz
criterion for a cubic [@routh1877stability; @hurwitz1895bedingungen], Lyapunov
stability [@khalil2002nonlinear], contraction analysis [@lohmiller1998contraction],
energy landscapes in neural networks [@hopfield1982neural], the information
measures the surprises are built from [@shannon1948communication;
@kullback1951information; @amari2016information], topological signal processing
[@ghrist2014elementary], the Markov-blanket literature [@kirchhoff2018blankets],
and the predictive-processing reading of the free energy principle
[@clark2013whatever; @friston2017process; @friston2019particular].

The registry places these established tools in a common executable setting.
The results distinguish the identities that hold, the assumptions needed for
domain-specific results, and the mechanisms that obstruct a proposed match.
