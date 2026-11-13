# Graphical abstract: three pairwise ledgers in four panels {#sec:graphical-abstract}

\setcounter{page}{2}

Every claim in this paper is one **correspondence**: a statement about a single
object, put in two or three vocabularies at once, executed by a probe that
returns values rather than a verdict. The registry holds
{{REGISTRY_ENTRIES}} such claims: {{REGISTRY_ISOMORPHISMS}} identifications,
{{REGISTRY_DISAGREEMENTS}} obstructions and {{REGISTRY_COMPLEMENTARITIES}}
complementarities, scoped to the pairs `oph-fep` ({{REGISTRY_OPH_FEP}}),
`synergetics-oph` ({{REGISTRY_SYNERGETICS_OPH}}), `fep-synergetics`
({{REGISTRY_FEP_SYNERGETICS}}) and `three-way` ({{REGISTRY_THREE_WAY}}).

The four panels are the paper in pictures, one page each: the ledger of the
three pairwise scopes; the carrier every scope is about; the dynamics that
separates the readings; and the accounting, where the programmes disagree about
currencies rather than facts.

## Panel 1 — the pairwise ledger: what each pair of programmes agrees on, obstructs, and supplies {#sec:ga-triangle}

Each edge carries one pair's ledger in the form the command line executes it:
`oph-fep-synergetics verify --pair oph-fep`, `--pair synergetics-oph`, `--pair
fep-synergetics`. The centre carries the one claim that needs all three
vocabularies at once. Most OPH–FEP work is obstruction: {{REGISTRY_OPH_FEP}}
claims, most of them failures of an identification proposed in good faith. Most
Synergetics–OPH work is identification:
{{REGISTRY_SYNERGETICS_OPH}} claims about symmetry, volume and counting. The
FEP–Synergetics edge is short by {{REGISTRY_FEP_SYNERGETICS}} claims and holds
two exact results about information, in one of which the programmes disagree
about the meaning of a number.

![The pairwise ledger triangle. Vertices are the three programmes, Observer Patch Holography, the Free Energy Principle and Fuller's Synergetics. Each edge is annotated with the counts of identifications, obstructions and complementarities for that pair, and the three load-bearing identifications are named on the edge where they sit. The centre carries the three-way ledger: the two claims that need all three vocabularies at once, `tetravolume-arithmetic` -- three exact quantities in three number systems -- and `conjugation-norm`. Every count is read from the executed registry at draw time rather than transcribed.](../output/figures/fig_ga_triangle.png){#fig:ga-triangle width=100%}

The scope of a claim is part of the claim: `carrier-port-chart` is a statement
about OPH's chart and the icosahedron's vertex set, `spectral-gap` is a
statement about one matrix read two ways — the generator of OPH's mean repair
operator and the FEP's precision — with the gap inherited from the shared object,
and
`kl-lattice-metric` is a statement about the FEP's divergence and Synergetics'
lattice metric. [@sec:correspondences] reports the three ledgers separately, and
[@tbl:registry] carries the scope of every row.

\clearpage

## Panel 2 — the carrier: one graph, its symmetry, its counting identities, and its Gaussian reading {#sec:ga-carrier}

Every scope is read on the same object, so the paper can compare vocabularies
without first having to reconcile subjects: twelve ports, thirty seams, twenty
triangular faces, every port of degree {{CARRIER_DEGREE}}, spectral gap
{{CARRIER_LAMBDA1}}. OPH supplies the repair law on the seams of that graph. The
FEP reads the same graph as a Gaussian precision, whose conditional-independence
graph is exactly the seam graph: two ports are dependent given the rest exactly
when they share a seam. The panel also carries two counting identities that
neither programme states but both inherit from the graph: the number of spanning
trees, and the dimensions of the agreement sheaf ([@prop:sheaf];
[@prop:cohomology]).

![The carrier and its exact invariants. Left: the twelve ports and thirty seams of the icosahedral carrier drawn from its derived adjacency. Right: the invariants read from the same graph at draw time, the spanning-tree number of the seam graph, the dimension of the global sections and of the obstructions of the scalar sheaf with their Euler characteristic, the automorphism group order split into rotations and reflections, the number of shell configurations and the orbits the group has on them, the Markov blanket size and pairwise precision of the Gaussian reading, and the marginal variance of a port.](../output/figures/fig_ga_carrier.png){#fig:ga-carrier width=100%}

The reduced Laplacian of the seam graph has determinant
{{CARRIER_TREE_NUMBER}}, so the carrier has exactly that many spanning trees, and
the same number is recovered from the spectrum as a product of eigenvalues. The
scalar sheaf on the seam graph has a {{SHEAF_H0}}-dimensional space of global
sections and {{SHEAF_H1}} independent obstructions, the coboundary has rank
{{SHEAF_COBOUNDARY_RANK}}, and the graph's Euler characteristic is
{{SHEAF_EULER_CHARACTERISTIC}}, the twelve nodes minus the thirty edges; with the
twenty faces added every obstruction is filled, and the surface's Euler
characteristic is {{SHEAF_SURFACE_EULER_CHARACTERISTIC}}.
The carrier's automorphism group has order {{SYMMETRY_GROUP_ORDER}}, splits into
{{SYMMETRY_ROTATIONS}} rotations and {{SYMMETRY_REFLECTIONS}} reflections when
the orientation of OPH's own faces is used as the test, and acts on the
{{SYMMETRY_CONFIGURATIONS}} total-two configurations in {{SYMMETRY_ORBITS}}
orbits with {{SYMMETRY_OUTCOME_MULTISETS}} distinct descent-outcome multisets.
As a Gaussian precision the same matrix puts every seam's pairwise precision at
{{GMRF_PAIRWISE_PRECISION}}, every port's Markov blanket at
{{GMRF_BLANKET_SIZE}} neighbours, and the variance of a port's deviation from the
carrier mean at {{GMRF_MARGINAL_VARIANCE}}, exactly {{GMRF_MARGINAL_VARIANCE_EXACT}}.

\clearpage

## Panel 3 — the dynamics: one consensus rule, its census, and the free energy it is not a function of {#sec:ga-descent}

The repair law is a total function on shell configurations, and its census is
exact: of {{SHELL_TOTAL_TWO}} states, {{SHELL_BINARY}} are binary,
{{SHELL_TRANSITIONS}} transitions occur, {{SHELL_STRICT_DESCENT}} lower the
quadratic energy strictly, {{SHELL_ZERO_DESCENT}} lower it not at all, and
{{SHELL_NOOP}} change nothing. That census is a group invariant rather than a
numerical accident, which is what the equivariance result in the previous panel
buys. The FEP's currency is a different function of the same configuration: half
the squared norm, whose level sets cut across the shell index instead of
refining it.

![The descent census and the free-energy level sets. Left: strict descents, zero descents and no-ops over every transition of the shell, annotated with the multiplicities of the equilibrium classes. Right: every distinct value of half the squared norm over the first five shells, one column per shell, so that a level carried by more than one shell appears as a repeated height.](../output/figures/fig_ga_descent.png){#fig:ga-descent width=100%}

One lattice step costs {{KL_ONE_STEP}} nats at unit covariance, because the
Gaussian divergence between equal-covariance neighbours is half the squared norm
and the first shell sits at squared distance {{KL_STEP_SQUARED_DISTANCE}}, so the
sites at squared norm $8k^2$ cost $4k^2$. The first five shells carry
{{SHELL_FREE_ENERGY_LEVEL_COUNTS}} distinct free-energy levels, and the levels
{{SHELL_SHARED_FREE_ENERGY_LEVELS}} are each carried by two shells: the shell
index does not fix the free energy, and the free energy does not fix the shell.
The shell is a level set of a polyhedral norm and the free energy of a Euclidean
one, and the obstruction is that neither refines the other.

\clearpage

## Panel 4 — the accounting: tetravolumes, one packing, and three number systems {#sec:ga-synergetics}

Synergetics counts volume in tetrahedra of the prime vector and asks which
solids come out integral. The IVM's own solids do:
{{TETRAVOLUME_INTEGRAL_SOLIDS}}. The carrier's own solid does not: it enters at
{{ICOSAHEDRON_TETRAVOLUME}} tetravolumes, which is
$(5/2)(3\sqrt{2} + \sqrt{10})$ exactly. The same accounting settles the
twelve-around-one question completely, because the cell of the twelve neighbours
is a rhombic dodecahedron whose volume can be derived rather than looked up, and
the ball inscribed in it gives the densest packing density in three dimensions
[@hales2005kepler].

![Tetravolume accounting and the twelve-around-one cell. Left: the tetravolumes of the IVM's integral solids beside the regular icosahedron's measured tetravolume on a logarithmic axis, showing that the carrier's own solid is the one entry an integer table cannot express. Right: the Voronoi cell of the twelve neighbours with its vertex, face and edge counts, its derived tetravolume from half the squared norm at each of its vertices, and the packing density of the inscribed ball against the closed form.](../output/figures/fig_ga_synergetics.png){#fig:ga-synergetics width=100%}

The cell of the twelve neighbours has {{VORONOI_VERTICES}} vertices,
{{VORONOI_FACES}} faces and {{VORONOI_EDGES}} edges, tetravolume
{{PACKING_CELL_TETRAVOLUME}} derived from the halved squared norm at each of its
vertices, and an inscribed ball of {{PACKING_BALL_TETRAVOLUME}} tetravolumes
which fills density {{PACKING_DENSITY}}, exactly $\pi\sqrt{2}/6 = \pi/(3\sqrt{2})$
because the cell's volume is rational and the ball's tetravolume is $\pi\sqrt{2}$.
The carrier's icosahedron enters at {{ICOSAHEDRON_TETRAVOLUME}}, computed exactly
as $\sqrt{2}$ times an element of $\mathbb{Q}(\sqrt{5})$, with squared value
{{ICOSAHEDRON_TETRAVOLUME_SQUARED}}, exactly $175 + 75\sqrt{5}$, and is not
integral: the tetravolume table and the carrier's own solid live in different
arithmetics even though they are the same geometry. The three programmes' exact
quantities live in three number systems: an integer {{KL_ONE_STEP}} nats for the
FEP's step, the quadratic surd {{CARRIER_SHELL_GAP}} for OPH's relaxation rate,
and integral tetravolumes beside one biquadratic one in Synergetics'.

\floatplacement{figure}{tbp}

\clearpage

\ophContentsPage
