# Brainstorm tests: 'salt' (fluid shares the gas's fate) and 'current' (fluid pools where cosmic flows converge)

## Salt rule: fluid = 5.4 x retained baryons, spread like hot gas (gas_retention_test.*) -> FAILS, but informative
Inputs (no free numbers): Popesso et al. 2024 gas fractions (arXiv:2411.16555), approximate stellar fractions (Gonzalez+13),
Brouwer+21 KiDS hot-gas model (hot = M*, rho ~ r^-2 to 143 kpc) and cold-gas formula.
A. Groups/clusters at R500, predicted/observed mass:
   | M500 | MOND alone | MOND + 5.4x retained baryons |
   | 1e13 | 0.99 | 1.52 |  | 1e13.5 | 0.76 | 1.21 |  | 1e14 | 0.66 | 1.15 |  | 1e14.5 | 0.64 | 1.27 |  | 1e15 | 0.66 | 1.57 |
B. KiDS (chi2 full | 0.1 dex floor | reliable points only): MOND stars+cold 202 | 68 | 34;  MOND + hot gas 103 | 52 | 31;
   MOND + hot gas + 5.4x fluid 343 | 146 | 119.
C. SPARC: MOND alone 4568; + fluid spread like hot gas 13464.
What the data actually say: galaxies (with their hot gas included) and groups up to ~1e13.5 need NO extra dark mass under MOND;
clusters (>= 1e14) need extra ~0.35 x M500, i.e. ~2-4 x their baryons (not 5.4). Extra mass appears only in the most massive systems.

## Implication for the 'current' picture
The cosmological fluid (Omega ~ 0.26, needed for the CMB) must by today sit almost entirely OUTSIDE galaxies and groups:
pooled in massive nodes (clusters, ~2-4 x their baryons) and spread through the cosmic web (filaments, sheets) as moving streams.
Next test: can streams outrun capture by galaxies when they form (z ~ 1-3)? Then: filament lensing and the large-radius
galaxy-lensing signal should show fluid in the web, not in bound galaxy halos.

## Current timing check (current_timing.*)
(1) Gravity-driven currents: matter inside a halo's turnaround radius is bound by construction, and filament streams feed galaxies
    directly at z ~ 2 (cold streams). A purely gravitational current traps fluid in halos of every mass. FAILS.
(2) Fluid locked to the cosmic rest frame (Khronon's preferred slicing), halos moving through it at their peculiar speeds
    (~310-350 km/s at z ~ 0.4-1.1): qualitatively right (galaxies V_esc 120-240 pass, groups and clusters trap), dividing mass
    ~5e12 Msun with no free numbers. Quantitatively wrong with a Maxwellian spread of speeds:
    | log M500 | trapped | needed by data |
    | 13.0 | 0.96 | 0.02 |  | 13.5 | 1.00 | 0.54 |  | 14.0 | 1.00 | 0.71 |  | 14.5 | 1.00 | 0.57 |  | 15.0 | 1.00 | 0.38 |
    (a Milky-Way-mass galaxy would still trap ~35%). And a fully frame-locked fluid could not cluster in the early universe (CMB).
What the data need, robustly: no extra mass up to ~1e13 Msun, roughly half the cosmic share at 1e13.5-1e14.5, less at 1e15
(subject to hydrostatic-mass bias, cluster baryon censuses and the MOND function at cluster accelerations).

## Three avenues for the cluster shortfall (three_avenues.*)
Target: fraction of the cosmic fluid share each system must hold at R500 under MOND (from the salt test above).
1. Partially locked fluid (drag toward the grid's rest frame; fluid is trapped if its speed relative to the halo is below k x V_esc,
   Maxwellian speeds, 1 fitted number k):
   - as measured: k = 0.29, mismatch chi2 = 24.9; predicted 0.02/0.14/0.34/0.61/0.91/1.00 vs needed 0/0.02/0.54/0.71/0.57/0.38
     (log M 12, 13, 13.5, 14, 14.5, 15). Right rise from groups to clusters; wrong fall-off at 1e14.5-1e15.
   - with 20% hydrostatic-mass bias: needed 0/0.50/1.09/1.23/0.97/0.65; k = 0.56, chi2 = 6.8, predicted 0.12/0.55/0.87/0.99/1.00/1.00.
     Still can't give more than 1 x the cosmic share (needed 1.1-1.2 at 1e13.5-1e14) or the drop at 1e15.
   - cost: the implied drag time (0.4 Gyr, tau*H0 = 0.03) would damp the fluid's growth on all scales (sigma8). It survives only if
     the drag acts on fast relative motion inside collapsing halos and not on large-scale flows -- but bulk flows are also
     ~300 km/s, so a plain speed threshold doesn't separate them. Open.
2. Hot dark component (e.g. massive neutrinos): at most 0.06 of the cosmic share can be hot (Planck+BAO); needed 0.38-0.71. FAILS (6-12x short).
3. Missing ordinary matter: clusters would need baryon fractions 0.33-0.51, i.e. 2.1-3.2 x the cosmic ratio, gathered from ~3x the
   volume and hidden from X-ray/SZ/stellar censuses. Not strictly excluded, implausible.
Best natural fit: avenue 1 (shape right, one number), with an unsolved cost to cosmological growth.

## Tesla valve + grid vibrations (valve_vibration.*)
A. Tesla valve: rectifies only at high Reynolds number (fast flow), symmetric when slow. So it is a speed switch that ADDS
   trapping (fast infall can't bounce back out); it can't stop slow gravitational capture by young galaxies. Alone: no.
   With vibrations, valve + oscillation = valveless pump (a way vibrations could drive net flow). Kept for later.
B. 'Shaken ocean': grid vibrations heat the dark fluid LATE, galaxies/groups lose their halos, deep clusters keep theirs.
   - Timing derived from our dark-energy law (rho_DE ~ adot^-1/2): rho_DE stops rising and starts falling when acceleration
     begins, z = 0.65; by today 6.5% of rho_DE (0.045 rho_crit c^2) is released. Heating needs only 1e-5 of it.
     (The gravitational-wave background is 400-1000x too weak, so the vibration must be of the grid/dark-energy field itself.)
   - Escapes the 'galaxies in the middle' no-go by timing: the early universe/Lyman-alpha (z 2-5) see a cold fluid;
     galaxies formed early with halos and later lose them; clusters form late and deep.
   - Shape (one number, sigma_h = the fluid's speed spread): as measured sigma_h = 587 km/s, chi2 25.5 (same as the locked fluid);
     with 20% hydrostatic bias sigma_h = 301 km/s, chi2 6.7 (MW-size galaxies keep 10%, need ~0; 1e15 still overshoots).
   - Cost (two-fluid linear growth from z = 0.65): sigma8 / unheated = 0.96 (300 km/s), 0.88 (600), 0.80 (800), 0.69 (1100).
     At 300 km/s it lowers sigma8 by ~4% -- in the direction of the weak-lensing S8 tension, not against it.
   - Open: fluid near halo centres is bound ~2.5x deeper than at R500; if it must be cleared there too, sigma_h -> ~750 km/s
     and sigma8 drops ~18% (excluded). Needs a proper heated-halo equilibrium calculation (next).
   - Predictions if it survives: galaxy lensing at lens redshift z > ~0.7 shows full dark halos (KiDS z ~ 0.2-0.4 shows none);
     groups at z > 0.7 hold the full cosmic share; disks expand as their halos leave (size growth since z ~ 0.7).

## Heated-halo equilibrium (heated_halo.*) -- the 'centres bound 2.5x deeper' worry
Heated fluid settles isothermally, rho_f = rho_bg exp(DeltaPhi/sigma_h^2), in the MOND well of the observed baryons plus the fluid's
own (Newtonian) gravity; MOND's log well is cut off by the external field g_e = 0.01-0.05 a0. (Same profile for a collisionless
gas with an isotropic speed spread, so it does not require the fluid to collide -- the Bullet Cluster needs that.)
- Galaxies: a Milky-Way-like galaxy keeps essentially nothing. Fluid/baryons < 0.01 inside 30 kpc and 0.01-0.05 inside 300 kpc for
  sigma_h >= 400 km/s. The surroundings are so dilute that even the deeper centre holds a negligible amount. Worry resolved.
- Groups/clusters: once the well is deep compared with sigma_h, no equilibrium exists -- the fluid stays a bound, self-gravitating
  halo, i.e. the system keeps its full share (it can't gather more than its supply).
- Averaged over surroundings (g_e 0.01-0.05 a0, 1-3x mean density), fraction of cosmic share kept:
  | sigma_h | 1e12 | 1e13 | 1e13.5 | 1e14 | 1e14.5 | 1e15 | chi2 measured | chi2 20% bias |
  |   400   | 0.00 | 0.62 | 1.00 | 1.00 | 1.00 | 1.00 | 54 | 6.1 |
  |   450   | 0.00 | 0.48 | 0.94 | 1.00 | 1.00 | 1.00 | 46 | 5.7 |
  |   500   | 0.00 | 0.35 | 0.81 | 1.00 | 1.00 | 1.00 | 37 | 8.1 |
  Best: sigma_h ~ 450 km/s with hydrostatic bias; almost all the remaining misfit is the 1e15 point (needs 0.65, gets 1).
  Without the bias, clusters keep too much (data say 0.4-0.7 of the share).
- Cost at 450 km/s: sigma8 x 0.925 (-7.5%) -- about the size of the weak-lensing S8 deficit relative to Planck (not fitted to it).
- Still not derived: sigma_h itself (needs the coupling between grid vibrations and the fluid); the fluid's gravity was taken unboosted.

## Lensing across redshift (lensing_z_test.*) -- the shaken ocean FAILS
Prediction: at fixed stellar mass, galaxy lensing (DeltaSigma at 50-300 kpc) at z ~ 0.37 (after heating) should be
0.37-0.59 of its value at z ~ 0.88 (before heating) -- roughly a factor 2 drop (Milky-Way-mass: 0.45-0.59).
Data: Leauthaud+2012 (COSMOS, z 0.2-1; arXiv:1104.0928): pivot Mh/M* constant (~27), low-mass scaling 'does not evolve
significantly'. Hudson+2015 (CFHTLenS, z 0.2-0.8; arXiv:1310.6784): peak M*/Mh falls 4.5% -> 3.4% toward today, i.e. halos
get relatively heavier, not lighter (blue galaxies constant). No factor-2 drop: excluded.
Escape route: finish the heating before z ~ 1 (start z >~ 1.2-1.5). Then the derived timing (z = 0.65) is lost, and the cost
grows: sigma8 / unheated at sigma_h = 450 km/s = 0.92 (from z 0.65), 0.88 (1.0), 0.83 (1.5), 0.80 (2.0) -- excluded.
General lesson: whatever keeps dark fluid out of galaxies must already be in place by z ~ 1 (lensing looks the same at
z 0.3 and 0.9), while the fluid must still be cold and clumpy at z 2-5 (Lyman-alpha) and keep growing (sigma8).

## Isolated wells (user idea): heat gained inside a well stays there; outside fluid stays cold (isolated_well.*)
Translation: heat the fluid locally on infall (fixed heat per mass, e.g. latent heat of a state change, triggered by infall faster
than v_th) instead of globally. Good: nothing changes with redshift (passes the lensing test by construction), large scales stay
cold, and it is two conditions (fell into a well AND the well is shallow), so it is not a single switch.
Catch: fraction of the cosmic fluid that has passed through wells deeper than v_th (Press-Schechter):
  | v_th | z=5 | z=3 | z=2 | z=1 | z=0 |
  |  20  | 0.16 | 0.32 | 0.43 | 0.57 | 0.72 |
  |  40  | 0.09 | 0.23 | 0.34 | 0.49 | 0.65 |
  |  80  | 0.03 | 0.13 | 0.23 | 0.38 | 0.57 |
  | 150  | 0.01 | 0.06 | 0.13 | 0.26 | 0.46 |
If the ejected fluid stayed hot (~450 km/s), 6-32% hot fluid at z = 3 removes ~12-54% of small-scale power (Lyman-alpha allows ~2%).
So the user's second half is required, not optional: escaped fluid must give its heat back (cool) quickly once outside the well,
i.e. heat is only ever held inside wells. That makes galaxies a 'dark fountain': in, heated, out, cooled, back in.
Next test: steady-state fluid held by a galaxy in the fountain (inflow x residence time) vs KiDS/SPARC, and whether cooled
fluid can re-enter clusters/groups at the needed rate. Unknown numbers: heat per mass (sigma ~ 450 km/s) and cooling time.

## Dark fountain (fountain.*) -- FAILS
Heated fluid must (1) cool fast outside wells so the cosmic fluid stays cold (Lyman-alpha z = 3: hot fraction <~ 1%), but
(2) get beyond the galaxy's reach before it cools, else it falls back and the galaxy keeps it all as a cycling cloud.
  | v_th | t_c max (Lyman-alpha) | heat needed to escape (turnaround / MOND reach) | clusters allow |
  |  20  | 0.08 Gyr | 2170 / 3090 km/s | <~ 500-600 |
  |  80  | 0.10 Gyr | 1650 / 2350 km/s | <~ 500-600 |
  | 150  | 0.16 Gyr | 1040 / 1480 km/s | <~ 500-600 |
Even the most lenient corner is 2x off, and with v_th = 150 km/s most z = 3 galaxies (V < 150) would never heat and keep halos.
Heat that big would also empty groups and clusters. Recycling (ignored) only tightens (1).
Lesson: heating can't do the separation. Lyman-alpha forces any heat outside wells to vanish within ~0.1 Gyr, i.e. within
~50-200 kpc of where it was made, which is inside the galaxy's own reach. 'Shallow wells lose fluid' by heat and
'the ocean stays cold' cannot both hold. Same family: shaken ocean (global heat) failed on lensing vs redshift.
