# What survives from iterations 1-69 that could solve dark energy's SIZE without the grid (Oct 4 2026)

## On the speed of light (Coalesce)
Light slowed in a lab travels through a medium (e.g. an ultracold atom cloud); the vacuum speed limit c is unchanged. The early
universe was such a medium (light trapped in the hot plasma until recombination). A change of the vacuum c itself was tested:
data prefer light 0.31% +/- 0.13% slower at recombination (2.4 sigma; partly the same signal as the heavier early electron);
quasars and atomic clocks forbid any change at late times (> ~1e-6 since z ~ 4). Kept as a live, testable piece.

## Pieces that do not need the grid
| piece | where from | status |
|---|---|---|
| uniform vacuum energy does not gravitate (unimodular structure) -> removes the 10^120 | it. 36, 68; GR track | standard physics; passes MICROSCOPE |
| uncertainty principle: dark energy x 4-volume ~ hbar/2 (needs only a smallest length, not a full grid) | it. 53, 68 | gives the right ballpark (0.04-35 x) |
| dark energy follows the cosmic 'now' (Claim 1), read as a tally across the horizon | it. 4E, 66; light-cone L5 | favoured by today's data (-5 to -7 in chi2) |
| open curvature of the 'now' surfaces; 'a' = physical curvature radius | light-cone model | testable (Omega_k > 0, now ~2 sigma) |
| slightly slower early light / heavier early electron (~alpha = 0.73%) | Claim 2, it. 57 | fits CMB + BBN + white dwarfs |
| neutrino link (lightest = 2.24 meV) | it. 67 | passes, weakly |
## Pieces that failed (not to reuse)
electron pays for dark energy (x10-25 overshoot, it. 64); excess into light (it. 65); 24% claims (it. 63); tracking tallies (it. 37, 66);
deriving C from curvature or age alone (light-cone derive_C_options).

## Most promising non-grid synthesis
Light-cone model + unimodular structure + uncertainty principle + one smallest length:
    rho_DE ~ hbar c / (2 l^2 sqrt(V4))
The light-cone model removes one convention by itself: everything began at one event, so the 4-volume we can know about is the
region between that event and us (our past light cone inside the first event's future cone -- the 'causal diamond'). From
iteration 68 (past light cone column): with the reduced Planck length (the one in Einstein's equations, 8 pi G) and Heisenberg's
1/2 -> 0.69 x measured; with factor 1 -> 1.39 x; with the plain Planck length -> 17 x / 35 x.
Remaining choices: the length unit (Planck vs reduced Planck, a factor 8 pi) and 1/2 vs 1; and the 'why now' problem (the formula
must be frozen, not re-set, to avoid tracking).

## Iteration 70: which smallest length? (PREREG_70.md; results iter70_smallest_length.txt)
- Black-hole entropy (Hawking, GR + hbar) singles out l = 2 l_P (1 unit of entropy per element): dark energy = 4.4 x measured.
- Einstein's 8 pi coupling gives 0.69 x but would need 2 pi units of entropy per element for black holes.
- The two natural requirements disagree; GR + hbar do not pin the smallest length. Size from this route: 0.35 - 4.4 x measured.
- Lesson: the route gets the size right to within a factor of a few (vs 10^120 naively) but is not a prediction yet. Next:
  'why now' brainstorm (Coalesce), then whether one principle can satisfy both black holes and dark energy.

## Iteration 71: dark-energy horizon treated exactly like a black-hole horizon (PREREG_71.md; iter71_*.txt)
- Rules (black-hole entropy length l = 2 l_P; the horizon's own 4-volume) turn the formula into holographic dark energy with a
  DERIVED coefficient c_h = 0.715.
- H1: with today's horizon (5.11 Gpc): 0.57 x measured. H2: but a self-consistent horizon-filling dark energy allows ANY present
  share -> size not pinned.
- H3 (today's data): its history is the reverse of Claim 1 (w -0.78 at z = 1, -1.10 today; dark energy grows like (1+z)^2 early,
  upsetting the CMB distance). Delta chi2 vs Lambda: +46 / +47 / +44 (Pantheon+ / DES / Union3) -> EXCLUDED. (Expectation H3
  'never above -1' was wrong: it starts above -1 and dives below.)
Lesson: dark energy shares horizon thermodynamics with black holes, but its size is NOT set like a black hole's (energy filling
its own horizon). Remaining non-grid size route: 0.35 - 4.4 x (iteration 70), still unpinned.
