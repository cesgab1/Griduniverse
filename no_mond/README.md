# Re-assessment without MOND as a requirement (Coalesce's rule, Oct 2026)
Rule: keep only what the DATA demand ('absolute truths'). Anything that failed only because it disagreed with MOND is re-judged
without MOND. MOND-motivated requirements are then taken one at a time, and kept only if they are observational facts.

## 1. Absolute truths (observed, model-independent) vs MOND inferences
Observed facts any model must match:
- T1 Rotation curves stay flat far out; the visible matter predicts them closely (radial acceleration relation, RAR).
- T2 Baryonic Tully-Fisher: v_flat^4 proportional to baryonic mass, slope ~3.85 +/- 0.09, small scatter.
- T3 One acceleration scale (~1.2e-10 m/s^2) appears in all galaxies.
- T4 Galaxy weak lensing out to ~1 Mpc (KiDS); fits both MOND and LCDM.
- T5 Clusters weigh 6-10x their visible matter; the Bullet Cluster's lensing mass is offset from its gas.
- T6 CMB, BAO, Lyman-alpha forest, sigma8, supernovae.
- T7 Galaxy lensing does not drop between z ~ 0.9 and 0.3; high-z Tully-Fisher shows mild evolution.
- T8 Solar System and gravitational waves: GR to high precision.
Contested (not absolute): external-field-effect detections, wide binaries (Chae vs Banik), dwarf spheroidal details.
MOND inferences (NOT facts): 'galaxies have no dark halo', 'the extra pull is a modification of gravity', 'fluid must stay
out of galaxies', and any cluster 'shortfall' measured relative to MOND.

## 2. Failed list, re-judged
Failed on MOND-independent data, still failed: five late Hubble fixes, electron switch at recombination, everpresent-Lambda
random walk, Wang-Unruh, ended black holes for H0, bounce instead of inflation, fluid clumping only in clusters (Lyman-alpha).
Existed only to serve MOND, now moot (neither pass nor fail): slack / mosaic / gap / channel derivations of a0; superfluid DM;
AeST variants; Khronon mass term; rate-dependent stiffness; any-single-switch; density-pressure window; DBI 'galaxies
halo-free'; shaken ocean; isolated wells / dark fountain; Khronon one-fluid sloshing; light fermions / hot dark component /
missing baryons for the MOND cluster shortfall; frame-locked 'current'.
Failed only against MOND-based targets, re-tested without MOND:
| idea | without MOND | verdict |
|---|---|---|
| 'Ocean' hybrid (fluid halo at the cosmic ratio) | = cold-dark-matter-like halos; clusters 0.7-1.4 of measured (already computed), CMB/Lyman-alpha/lensing-vs-z fine | BEST FIT |
| DBI fluid (dust-like in wells) | with MOND off it is a CDM-like fluid in collapsed objects; same as above | merges with the best fit |
| 'Salt' (fluid = 5.4 x retained baryons, hot-gas profile) | SPARC Newtonian chi2/pt 22.4 (no free numbers), wrong Tully-Fisher slope (v^2 ~ M) | still FAILS |
| 'Current' (fluid streams past galaxies) | galaxies need dark mass in Newtonian gravity, so a fluid that avoids them fails T1 | still FAILS |
Best fit without MOND: the ocean is cold dark matter (clumps everywhere, collisionless at the Bullet Cluster), plus the grid
dark-energy law (which beat Lambda on its own), plus the grid picture (Planck stars, no singularity, fate).

## 3. MOND-inferred requirements, one at a time
R1 (T1) Rotation curves from the visible matter (sparc_without_mond.*; 171 galaxies, 3375 points, 5% error floor, Newtonian
halos vs MOND):
| model | free numbers per galaxy | chi2 / point |
|---|---|---|
| MOND, fixed M/L | 0 | 6.9 |
| CDM-like NFW halo from abundance matching (Moster+13, Dutton & Maccio c-M) | 0 | 17.6 |
| MOND, disk M/L free | 1 | 2.65 |
| NFW, halo mass free | 1 | 3.04 |
| cored halo, mass and core free | 2 | 1.15 (best BIC) |
Reading: rotation-curve SHAPES are fit fine by dark halos given freedom (cored halos are best). What MOND captures and a naive
halo doesn't is PREDICTIVITY: with no free numbers, MOND is 2.6x better. That is a real fact (T1/T3) that the 'ocean = CDM'
model must earn through galaxy-formation physics (feedback making halos track baryons). That is the standard LCDM position:
plausible, debated, not a refutation.
Next (one at a time): R2 Tully-Fisher slope and scatter (T2); R3 the single acceleration scale (T3).

## 4. Revisit: the mosaic with 'electron twins', and fluid/gas between the grid layers (layers_and_twins.*)
Question: did these fail only because of MOND?
| idea | earlier verdict and why | without MOND |
|---|---|---|
| Electron twins = the fake copies a grid gives every electron (doublers) | never judged on MOND. Re-run (`grid_models/fermions_*.txt`): on the mosaic the plain operator gives 8-13 copies (square grid 4); the grid's own Wilson term leaves exactly 1; the overlap operator gives exactly 1 with exact handedness symmetry (residual 1e-13) | PASSES (unchanged). The removed copies get grid-scale mass: no relic, and as copies of the electron they would be charged, so not dark matter |
| Two-layer electron (each wall carries one handed half; mass = leak between layers) | never judged on MOND | PASSES (unchanged) |
| Mosaic as the gravity grid | 'snaps all at once, no square-root law': a MOND criterion | moot; static links on the mosaic give plain Newton (`network/README.md`), which is all that is needed now |
| Light electron-like twin (2 eV) as the cluster's extra mass | failed: the amount needed for the MOND cluster shortfall is 20-30x what the CMB allows | as ALL the dark matter: Pauli limit from SPARC cores needs >= 200 eV (median galaxy 52 eV); Lyman-alpha needs > ~5.7 keV. A heavy twin is then cold-dark-matter-like and merges with the best fit; its galaxy/cluster split (the reason it was liked) was a MOND-era virtue and is gone. Mass not derived |
| Squeezed layers (a0 grows with depth) | failed against SPARC/KiDS: a MOND criterion | moot |
| Fluid in channels between the layers (exact 1-D gas law, P ~ rho^2 dense, rho^3 dilute) | failed as a MOND superfluid (its push bypasses lensing; a0 not derived): MOND criteria | re-tested with its own pressure law, no MOND: FAILS on two MOND-free facts. (i) A single pressure law fixes how core size goes with core density: slope 0 (rho^2) or +0.5 (rho^3); SPARC Newtonian cored-halo fits (110 galaxies) give -0.53 +/- 0.04 (12 and 24 sigma off). (ii) A pressure strong enough to make the ~3.7 kpc cores gives a Jeans length of 7 Mpc at recombination and 0.7 Mpc at z = 100 (must be < ~0.1 Mpc): CMB and Lyman-alpha excluded. With the pressure turned down it is plain cold dark matter living between the layers = the best fit again |
Reading: none of the two pictures was killed by MOND alone in a way that now revives them. The electron-twin work never depended on
MOND and still passes (it is particle physics, not the dark sector). The layer fluid fails on its own pressure law, so 'between the
layers' survives only as a picture of WHERE cold dark matter sits, not as new physics. The best fit stays: ocean = cold dark
matter + grid dark-energy law + grid picture.
