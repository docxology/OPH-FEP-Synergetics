# Graphical abstract: claim scope and model evidence {#sec:graphical-abstract}

\setcounter{page}{2}

The registry compares executable models of Observer Patch Holography (OPH),
the Free Energy Principle (FEP) and Fuller's Synergetics. It holds
{{REGISTRY_ENTRIES}} claims, each paired with a probe that reports values
and a verdict. The claims include {{REGISTRY_ISOMORPHISMS}} identifications,
{{REGISTRY_DISAGREEMENTS}} obstructions and {{REGISTRY_COMPLEMENTARITIES}}
complementarities across `oph-fep` ({{REGISTRY_OPH_FEP}}),
`synergetics-oph` ({{REGISTRY_SYNERGETICS_OPH}}), `fep-synergetics`
({{REGISTRY_FEP_SYNERGETICS}}) and `three-way` ({{REGISTRY_THREE_WAY}}).

The panels map the claims, show the carrier and its Gaussian reading,
compare native repair attempts with Gaussian costs on lattice shells, and
explain how tetravolumes are measured.

## Panel 1 — scope of the executable claims {#sec:ga-triangle}

Each edge shows the total for one pair. The ledger gives its claim kinds
and a representative result. Run
`oph-fep-synergetics verify --pair` with the scope name to check its entries.

![Claim totals on the framework triangle, with representative results and counts by kind. The centre counts claims that use all three frameworks.](../output/figures/fig_ga_triangle.png){#fig:ga-triangle width=100%}

`carrier-port-chart` links OPH's port chart to the icosahedron's vertex set.
`spectral-gap` links the mean repair operator to a Gaussian precision built
from the same Laplacian. `kl-lattice-metric` expresses equal-covariance
Gaussian divergence as a lattice cost. [@tbl:registry] gives the scope of
each claim. [@sec:transfer-criteria] tests which results transfer to
other finite graphs, lattice shells and cellulations.

\clearpage

## Panel 2 — the carrier graph and its Gaussian reading {#sec:ga-carrier}

The icosahedral carrier has degree {{CARRIER_DEGREE}} and Laplacian spectral
gap {{CARRIER_LAMBDA1}}. The OPH model repairs the graph's seams, and a diagonal
shift of the Laplacian gives a proper Gaussian precision. In that Gaussian,
two ports are conditionally dependent given the rest exactly when they share a
seam. The graph also determines the spanning-tree count and the agreement
sheaf's cohomology dimensions ([@prop:sheaf]; [@prop:cohomology]).

![The icosahedral carrier, with solid foreground seams and dashed rear seams. The ledger groups graph and sheaf invariants, symmetries, and the intrinsic and proper Gaussian quantities.](../output/figures/fig_ga_carrier.png){#fig:ga-carrier width=100%}

The reduced Laplacian has determinant {{CARRIER_TREE_NUMBER}}. By the
matrix-tree theorem, it equals the carrier's spanning-tree count. The count
also equals the product of the nonzero Laplacian eigenvalues divided by the
number of ports. The scalar agreement sheaf has a {{SHEAF_H0}}-dimensional
space of global sections and {{SHEAF_H1}} independent cycle classes.
The coboundary has rank {{SHEAF_COBOUNDARY_RANK}} and the graph's Euler
characteristic is {{SHEAF_EULER_CHARACTERISTIC}}. Adding the oriented faces
fills the cycle space and gives the surface Euler characteristic
{{SHEAF_SURFACE_EULER_CHARACTERISTIC}}.

The carrier's automorphism group has order {{SYMMETRY_GROUP_ORDER}}, with
{{SYMMETRY_ROTATIONS}} rotations and {{SYMMETRY_REFLECTIONS}} orientation-reversing
symmetries. It acts on the {{SYMMETRY_CONFIGURATIONS}} total-two configurations
in {{SYMMETRY_ORBITS}} orbits, with {{SYMMETRY_OUTCOME_MULTISETS}} distinct
descent-outcome multisets. The proper Gaussian has seam precision
{{GMRF_PAIRWISE_PRECISION}} and a blanket of {{GMRF_BLANKET_SIZE}} neighbours
per port. The centred intrinsic field gives port variance
{{GMRF_MARGINAL_VARIANCE}}, exactly {{GMRF_MARGINAL_VARIANCE_EXACT}}.

\clearpage

## Panel 3 — native repair and Gaussian costs on lattice shells {#sec:ga-descent}

Native repair acts on a configuration, a seam and a coin outcome. Of the
{{SHELL_TOTAL_TWO}} states in the census, {{SHELL_BINARY}} are binary. Over
{{SHELL_TRANSITIONS}} attempts, {{SHELL_STRICT_DESCENT}} lower the load energy,
{{SHELL_ZERO_DESCENT}} move with the energy unchanged, and {{SHELL_NOOP}}
leave the state unchanged. The repair rule commutes with the carrier
symmetries, so the census is invariant under their action.

![Left: counts and shares of native repair attempts, with the path multiplicities below. Right: distinct unit-covariance Gaussian costs, their ranges on each coordination shell, and the dashed curve of shell maxima. Open diamonds connected by dotted guides mark costs shared by shells.](../output/figures/fig_ga_descent.png){#fig:ga-descent width=100%}

Gaussian divergence at unit covariance is half the squared distance between
the means. One lattice step has squared length {{KL_STEP_SQUARED_DISTANCE}}
and costs {{KL_ONE_STEP}} nats. The shell maxima occur at squared norm $8k^2$
and cost $4k^2$. The plotted shells contain {{SHELL_FREE_ENERGY_LEVEL_COUNTS}}
distinct costs respectively, with shared levels {{SHELL_SHARED_FREE_ENERGY_LEVELS}}.
A coordination shell follows a polyhedral norm, while a level of Gaussian
divergence follows the Euclidean norm. A coordination shell can contain
several divergence values, and a single divergence value can occur on
several shells.

\clearpage

## Panel 4 — tetravolumes and packing {#sec:ga-synergetics}

The IVM table measures volume in regular tetrahedra whose edge is the
nearest-neighbour distance. Its selected solids have integer tetravolumes
{{TETRAVOLUME_INTEGRAL_SOLIDS}}. The carrier's regular icosahedron, measured in
tetrahedra of its own edge length, has tetravolume
$(5/2)(3\sqrt{2} + \sqrt{10}) \approx {{ICOSAHEDRON_TETRAVOLUME}}$.

![Left: exact IVM tetravolumes and the carrier icosahedron's approximate tetravolume on a logarithmic axis. Right: the IVM Voronoi cell with its incidence counts, exact tetravolume and packing density.](../output/figures/fig_ga_synergetics.png){#fig:ga-synergetics width=100%}

The IVM Voronoi cell is a rhombic dodecahedron with {{VORONOI_VERTICES}}
vertices, {{VORONOI_FACES}} faces and {{VORONOI_EDGES}} edges. Triangulating
its faces and joining them to the centre gives volume
{{PACKING_CELL_TETRAVOLUME}} in the IVM tetrahedral unit. Its inscribed ball
occupies {{PACKING_BALL_TETRAVOLUME}} tetravolumes, so the packing density is
{{PACKING_DENSITY}}, exactly $\pi/(3\sqrt{2})$. This packing attains the
three-dimensional maximum established by Hales [@hales2005kepler].
The carrier's tetravolume is $\sqrt{2}$ times an element of
$\mathbb{Q}(\sqrt{5})$, with square $175 + 75\sqrt{5} \approx
{{ICOSAHEDRON_TETRAVOLUME_SQUARED}}$. Together with the integer
{{KL_ONE_STEP}} nats for a lattice step and the quadratic surd
{{CARRIER_SHELL_GAP}} for the normalised repair gap, it illustrates distinct
units and number systems in the three comparisons.

\floatplacement{figure}{tbp}

\clearpage

\ophContentsPage
