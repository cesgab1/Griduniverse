# PREREG 123: "IT is space" -- testing the three ways out of the bit-count puzzle

Written and committed BEFORE collecting the measurements below and before any code is run.

## The puzzle (ledger, iteration 122 follow-up)
If space is made of bits filling volume, a 1 m ball holds >= 8.2e80 bits (cells <= 1.7e-27 m). The black-hole
area rule allows only ~1.7e70. Three ways out, each tested here.

## T1 -- LOCKED BITS (bits fill volume but most move together)
Model: distance measurements over a length L carry a jitter dL = l^(1-b) * L^b, with l the cell size (metres).
 b = 1/2 : bits independent (random walk)                        -> "free crowd"
 b = 1/3 : bits locked so independent count ~ area (holographic) -> "crowd doing the wave"  (Ng-van Dam form)
 b = 0   : no accumulation (jitter = one cell)
Measurements (published, to be collected): (a) energy-independent arrival-time spread of GRB 090510 spikes (ms
scale), our own calculation: dL < c * spike width; (b) energy-dependent random spread of photon speeds (Fermi,
Vasileiou et al. 2015) -> random-grid cell size; (c) Holometer correlated jitter (Chou et al. 2016/2017);
(d) image-phase blurring of distant sources (Perlman et al. 2015) -- CONTESTED in the literature; reported, not
counted in the headline.
Output: upper limit on l in metres for each b. Expectation (written now): locked bits (b = 1/3) are the hardest to
see; limits for b = 1/3 will be weak (>> 1e-27 m) from timing; b = 1/2 much stronger. Free choices: spike width
used (1); which distance (comoving) (1).

## T2 -- SURFACE BITS (bits live on surfaces; a black hole's surface grows in tiles)
Model (Bekenstein-Mukhanov): horizon area changes in steps of a tile a = alpha * hbar G / c^3 (measured constants
only; alpha is the unknown, a in m^2 is reported). A spinning black hole can then only ring at
   omega_n * M = (kappa M) * alpha * n / (8 pi) + 2 (Omega M)       (G = c = 1, l = m = 2 mode, n integer)
Test: take published remnant mass and spin (from the inspiral, independent of the ringdown) and the measured
ringdown frequency for the loudest events (GW150914, GW250114 if available); find which alpha put a line within the
measurement uncertainty. Named values checked: alpha = 4 ln 2 (one bit per tile, Coalesce's idea), 4 ln 3 (Hod),
8 pi (Bekenstein's original). POPULATION RULE: one alpha must fit all events.
Expectation (written now): small alpha (one bit) gives a comb of lines finer than current spin uncertainty ->
NOT testable yet (neither excluded nor supported); only alpha >~ 10-30 can be tested. Also collect published
results (Laghi et al. 2021; echo searches, LVK) and compare.
Free choices: line-match criterion (measurement uncertainty, headline; QNM half-width, alternative) (1); which
events (1).

## T3 -- THE AREA RULE BREAKS AT SMALL SCALES
Correction written before testing: the black-hole area rule S = A c^3/(4 G hbar) uses G at the HORIZON scale (km
for stellar black holes, measured there by gravitational waves), not at 1e-35 m. Its real small-scale extrapolation
is inside Hawking's derivation: outgoing light near the horizon is traced back to wavelengths far below any grid.
Tests (published): (a) how far down G is measured (short-range gravity, Lee et al. 2020); (b) whether real black-hole
area never decreases (Isi et al. 2021; LVK GW250114); (c) whether a medium with KNOWN graininess still produces
Hawking's thermal radiation (analog black hole, Steinhauer et al. 2019): compute the grain-to-horizon ratio there and
compare with the largest allowed ratio for a real black hole (cell 1.7e-27 m vs horizon).
Expectation (written now): no evidence the rule breaks; analog grains are relatively ~1e30 times coarser than any
allowed real grid and the thermal law still held -> avenue 3 loses support but is not excluded (no tiny-scale test).

## Rules
Results are limits, not detections; no search over many models for a signal -> no look-elsewhere correction needed.
If any T2 alpha is "supported", require it in BOTH events (independent confirmation) before saying so.
Planck length appears only as a comparison, never as an input.
