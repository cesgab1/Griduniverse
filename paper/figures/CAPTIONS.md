# Figures for the paper: captions, analogies, and where each analogy stops working

Regenerate all figures with `python3 make_figures.py` (PNG at 200 dpi and vector PDF).
Each panel is tagged in small grey italics: **computed** means the numbers come from a calculation in this repository;
**picture** means a drawing that shows an idea and carries no data. Each analogy below is followed by the place where it
breaks, so a reader or auditor can see what the picture claims and what it does not.

---

## Figure 0 (overview / cover). Layers of the mosaic grid, the fluid between them, and a mass pulling on them
**Caption.** Seven layers of the mosaic grid (random cells joined by links that carry tension), with the fluid between
them (dark matter, violet), which fills every gap and has pooled around the mass. The mass pulls on the grid from every
side: the layer level with it is pulled inward sideways, the layers above are pulled down toward it, and the layers
below are pulled up toward it. The pull weakens with distance (softened inverse square, the shape computed in Figure 2).
The pattern is the same above and below because gravity pulls equally in all directions; the grid itself is random, not
symmetric. Gold: a circular orbit, the straightest path available in the pulled grid. Pale blue: light passing above and
below the mass bends toward it. Dotted lines show where it would go without the mass. Bending is computed by 3-D ray
tracing and exaggerated.

**Open point for the text.** In the relativistic version (Khronon) the layers are read as successive moments of time;
the figure draws them stacked in space to show how the pull acts above, beside and below a mass. Which reading the paper
uses is to be settled between the authors.

**Why not the rubber sheet.** The usual stretched-blanket picture shows space dipping *downward*. That needs an extra
direction for the dip and uses 'downhill' (gravity) to explain gravity. In 3-D there is no 'down': every part of the grid
is drawn toward the mass, from above, below and the sides alike. The draw-in is the drawing's way of showing where the grid
ticks slower (Figure 3). That slowing is what bends paths, not a slope.

**Analogy: a sponge squeezed toward a point inside it.** Push a pin into the middle of a block of sponge and pull it
inward: the sponge around it is drawn in from all sides. Layers above move down, layers below move up, the layer at the
pin's height moves sideways, and the effect fades with distance. Anything sliding through the sponge gets steered toward
the pin. That steering is gravity.

**Where it stops working.** A real sponge is pulled by the pin; in the model nothing pulls the grid with a rope. The
draw-in is how the grid's tension settles around the mass, and its observable effect is the slower ticking. Everything is
hugely exaggerated. Clean version: `fig0_grid3d.png`.

## Figure 1. Anatomy of the grid
**Caption.** Space is modelled as a stack of layers. Each layer is the whole of space at one moment: a random mosaic of
cells joined by links that carry tension. Moments follow one another layer by layer. A fluid (the 'ocean') fills the
gaps between layers; it has mass but does not shine or collide, and it plays the role of dark matter. The overall stretch
of the net is the dark energy.

**Analogy: a flip-book printed on fishing net.** Each page of a flip-book is one moment, and flipping the pages runs the
film. Here every page is a sheet of fishing net pulled taut on its frame. Pulling on a knot is felt by the neighbouring
knots, then by theirs, and so on outward: that sharing of pull is gravity (Figure 2). The net is knotted at random, not
woven like a window screen. A woven screen has 'grain' directions, and a ball rolling on it would notice them. A random
net has no grain, so every direction is the same. (We checked: a cube-shaped lattice leaves a 20% imprint of its axes;
a random one does not.) Between the pages, like ink soaked into the paper, sits the fluid.

**Where it stops working.** Real layers are not separated by any distance you could measure with a ruler. 'Between the
layers' means 'part of the grid that light does not interact with', not a physical gap. The spacing of the net is far
too small to see: gamma-ray-burst timing requires it to be below 6e-28 m.

## Figure 2. Gravity is the grid's tension, shared out link by link
**Caption.** (a) A single mass pulls on one knot of a random 2-D net that carries tension only. The equilibrium dip,
computed by solving the network's balance equations, spreads out smoothly even though the net is irregular. (b) On a
random 3-D net of 24,000 knots, the pull (the slope of the dip) falls with distance with a measured slope of -1.96;
Newton's inverse-square law gives -2. (c) The reason: the same pull is shared by all the links crossing any shell around
the mass, and a shell twice as far out has four times as many links.

**Analogy: a crowd passing along a rope pull.** Imagine one person in the middle of a field tugging on everyone around
them through ropes. The people next to them feel the tug strongly. Those people pass a share of it on to everyone in the
next ring out, and so on. Because each ring holds more people than the last (the area of a sphere grows as distance
squared), the tug each person feels thins out as 1/distance². Nobody has to 'know' the law: it comes from counting how
many hands the pull is split between.

**Where it stops working.** This explains the *shape* of Newton's law, not its strength (Newton's G), which is set by
the stiffness of the links and is put in by hand. In the network, gravity is the equilibrium of the net. In the full
theory the net is the relativistic Khronon field, which reproduces Einstein's gravity in the Solar System and
gravitational waves travelling at the speed of light.

## Figure 3. What gravity actually is: slower ticking near mass
**Caption.** (a) Each row of a marching band steps forward at a pace set by the local 'tick rate' of the grid. Where the
grid ticks slower (near a mass, shaded), the ranks there take shorter steps, so every row turns toward the slow side.
Nothing pulls on the marchers; the straightest available path itself bends. This toy is computed with step length equal
to the local tick rate. (b) Light crossing the same tick-rate field, computed by ray tracing: rays bend toward the mass,
and rays that pass closer bend more. The mass is exaggerated about 10^5 times; for the Sun the measured bending is 1.75
arcseconds.

**Analogy: a marching band crossing into mud (the refraction picture).** Every marcher walks straight ahead and keeps
time with the drum. When the left end of a row reaches mud, those marchers cover less ground per beat while the right end
keeps its stride, so the whole row swings left. Ask the marchers and each will say they walked straight. This is what
'falling' is in the grid picture: an apple does not get a rope-pull from the Earth. The grid ticks a little slower lower
down, and the apple's straightest path through space *and time* curves toward slower ticking. Clocks really do run
slower lower down: this has been measured in a tower (Pound and Rebka, 1959), in aircraft, and every day in GPS satellites,
whose clocks must be corrected for it.

**How Figures 2 and 3 fit together.** They are the same thing seen two ways. Figure 2 shows how much the net is pulled
(the 'shape' of gravity). Figure 3 shows what that pull does to motion: the dip in the net *is* a region of slower
ticking. The rubber-sheet picture in textbooks gets this wrong. It uses a ball rolling *downhill* on the sheet, which
assumes gravity to explain gravity. The marching band needs no 'downhill': only the tick-rate difference.

**Where it stops working.** In a real orbit the bending of time does almost all the work at low speeds, and the bending
of space adds an equal amount for light. That is why light bends twice as much as a naive 'falling light' estimate. The
toy shows the mechanism, not the factor of two; the ray tracing in (b) includes it.

## Figure 4. Dark energy is the grid's own tension, and it 'gives' when the stretching speeds up
**Caption.** (a) In the grid law the dark-energy density goes as one over the square root of the rate at which the
universe is stretching (ρ_DE ∝ ȧ^-1/2), so it is not constant (dashed line: Einstein's cosmological constant Λ). It
grew while the expansion was slowing down, peaked when the expansion began to speed up, and fades as it speeds up
further. (b) Its fingerprint in the equation of state w: w < -1 in the past, w > -1 today, crossing -1 exactly when
acceleration begins (z = 0.73 for Ω_m = 0.31). DESI's two-number form approximates this curve as (w0, wa) = (-0.80, -0.49).
(c) Where three combinations of real data put the best-fitting free (w0, wa) (our refits: DESI DR2 BAO + Planck + each
supernova sample). The grid law, which has no dark-energy parameter to adjust, sits among them. It beats Λ in all three
combinations (Δχ² -2.7, -5.1, -10.0).

**Analogy: a shear-thinning fluid, like ketchup or paint.** Ketchup resists when you tilt the bottle gently and gives
way when you shake it hard: its resistance drops as the rate of stirring goes up. The grid's tension behaves the same way
with respect to the *rate* of stretching. While gravity was braking the expansion, the slow stretching kept the tension
high and rising. Once the expansion began to speed up, the tension started to 'give'. The turnover sits exactly at the
moment the braking stops, because that is the moment the stretching rate stops falling. That is a sharp, testable
coincidence that ordinary dark-energy models have no reason to share.

**Where it stops working.** Ketchup thins because of its molecules; we do not have a microscopic mechanism for the
grid's -1/2 power. It is a law fitted in form and confirmed by data, not derived. The exponent fitted freely is
0.63 ± 0.22, consistent with 1/2.

## Figure 5. The fluid between the layers: pulled by gravity, invisible, and it passes through itself
**Caption.** (a) A galaxy (orange) sits in a pool of fluid (blue) about five to six times its own mass. (b) Why the pool
is needed: the measured rotation of NGC 3198 (SPARC data) stays flat far out, while stars and gas alone would make it
fall (orange). Adding a fluid pool (a cored halo fitted to these data; blue) matches it. (c) The Bullet Cluster: two
clusters after a collision. Their hot gas (red, seen in X-rays) collided and stopped in the middle, while most of the
mass, found by gravitational lensing, passed through. That is what the fluid does.

**Analogy: two swarms of bees and two water balloons.** Throw two water balloons at each other and they splash and stop
in the middle: that is the cluster gas. Fly two loose swarms of bees through each other and they come out the other side
almost untouched, because the bees are far apart and do not bump: that is the fluid. Light shows us where the balloons
are; only the bending of background light (lensing) shows where the swarms went.

**Where it stops working, and what is still open.** In the best-fit model the fluid behaves as standard cold dark matter.
That passes the CMB, galaxy clustering, clusters, the Bullet Cluster and lensing. It does not, by itself, explain why
galaxies follow tight rules such as the Tully-Fisher relation (slope 2.65 predicted vs 3.59 observed). This is the same
open problem standard cosmology has; see `no_mond/README.md`.

## Figure 6. The history of the universe in the grid picture
**Caption.** From a bounce (the grid has a maximum density, so a collapse rebounds instead of forming a singularity),
through the release of the CMB and the dark ages, to galaxies, the start of acceleration, today and the far future.
Each stage carries its status: green = tested against data and passing, orange = a prediction waiting for data, red =
tested and failed, grey = picture or consequence. The time axis is not to scale.

**Analogy: a spring-loaded ruler.** The universe is like a measuring tape pulled out of its case against a spring. Early
on the matter brakes it. The spring's pull (the grid tension) builds while the braking wins, then eases once the tape
starts running away on its own. The bounce at the start is the tape hitting its stop inside the case: it can be wound in
only so far.

**Where it stops working.** A bounce alone cannot replace inflation: it gives the wrong ripple spectrum (~30σ from
Planck), so inflation is still needed. The figure marks this in red.

## Figure 7. Matter as patterns in the grid: the electron's two halves on two layers
**Caption.** (a) In a lattice model of the electron as a grid pattern (the domain-wall construction), one electron splits
into a right-moving half on one layer and a left-moving half on another, separated by a slab. Computed: 98% of each half
sits on its own wall. (b) The electron's mass comes from the two halves 'leaking' into each other through the slab. It
halves with every extra layer between them, so a light electron needs only a modest gap. (c) Analogy.

**Analogy: a rope thrown across a river.** Two people stand on opposite banks holding a rope. The tug they share is the
electron's mass. Across a narrow river the tug is strong (a heavy particle); the wider the river, the weaker the tug, and
it weakens by the same factor for every extra metre (here: halving per layer). This is how the grid can produce
particles far lighter than its natural scale without fine-tuning. The two banks also explain why the weak force acts on
only one 'handedness': a force that touches only one bank feels only one half.

**Where it stops working.** This is a 2-D lattice model, borrowed from lattice field theory (Kaplan 1992), not a
derivation of the electron's measured mass. On a random grid the simplest version gives 8-13 fake copies of the electron.
The grid's own correction term (Wilson) or its exact version (overlap) leaves exactly one, as verified in
`grid_models/fermions_*.txt`.

## Figure 8. Black holes: the grid cannot be compressed without limit
**Caption.** (a) The density at the core of a collapsing star, computed for dust collapse with a density cap of the kind
loop quantum cosmology gives. In Einstein's theory (dashed) the density runs to infinity; with the grid's cap it peaks
at the maximum and rebounds. (b) Analogy.

**Analogy: a full lift.** People can be squeezed into a lift only until they are shoulder to shoulder. Push harder and
the crowd pushes back, so the squeeze turns into a rebound. A collapsing star's core does the same at the grid's maximum
density, forming a 'Planck star' hidden inside the black hole's horizon.

**Where it stops working.** From outside, nothing differs from an ordinary black hole for more than 10^26 years, so this
is consistent but not testable with current data.

## Figure 9. Was the electron heavier when the CMB was released?
**Caption.** Fit quality (Δχ², lower is better) against the electron mass at recombination, relative to today's value,
from Planck + ACT DR6 + SPT-3G D1 + DESI DR2 BAO + DES-Dovekie supernovae. All other cosmological parameters are
re-fitted at each point. Best fit: m_e = 1.0099 ± 0.0049, 2σ from today's value. The shaded band is the grid prediction
(0.4-1.1%), derived from the energy budget that links the electron to dark energy, not fitted to these data.

**Analogy: a guitar string tuned by the frame's tension.** The constants of nature are set by how tightly the grid is
strung, the way a string's pitch is set by its tension. Early on the grid was strung slightly differently, so the
electron 'note' was about 1% higher. It then relaxed to today's value, and the energy released went into the grid's
overall tension: the dark energy of Figure 4.

**Where it stops working.** A 2σ preference is a hint, not a detection. The decisive tests are the full SPT-3G survey
and the Simons Observatory, which should measure this at about the 0.2% level.

## Figure 10. The grid's building block: one cell, its links, and how big it may be
**Caption.** (a) A real cell from a random 3-D mosaic (computed). It is a many-sided polyhedron: this one has 17 faces, and
the average over all cells is 15.5 (a cube has 6). Each face is a wall shared with one neighbour. Through each face runs
one link (gold) from this cell's node (white) to the neighbour's node. Tension is carried along the links, and each link's
strength is set by the area of the face it crosses; this is the same rule used in our network calculations. (b) Cells are
not all the same size. Over 4,096 computed cells the volume spreads by ±43% around the average, and 98% lie between 0.26
and 2.2 times the average. Only the average size matters for physics, and this randomness is what keeps every direction
equal (a grid of identical cubes imprints its axes, about 20%). (c) The size window for a cell. MAX: below 5.7e-28 m.
This is measured: light from distant gamma-ray bursts (LHAASO, GRB 221009A) shows no sign of a grid, while a coarser grid
would make different colours arrive at different times. MIN: the Planck length, 1.6e-35 m. This is not measured; it is
where the model's density cap sits (the cap that turns a collapse into a bounce, ~0.4 of the Planck density), and below
it distances lose their meaning in any quantum theory of gravity. The window spans 7.5 powers of ten.

**Analogy: a sponge seen up close.** From across a room a sponge looks smooth; up close it is a tangle of cells of every
size and shape, joined wall to wall. No two cells are alike, yet a block of sponge squeezes the same way in every
direction because the cells are random. Space in the model is the same. Even the finest probe we have, the LHC, sees down
to only about 1e-19 m, and the cells are at least 100 million times smaller than that, so space looks perfectly smooth to
every experiment.

**Where it stops working.** A sponge's cells are made of material in empty space; the grid's cells *are* space, with
nothing around or between them. The minimum size is a property of the model, not a measurement. The maximum is a firm
observational limit.
