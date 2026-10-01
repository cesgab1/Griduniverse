# Clusters and where the dark fluid settles (Khronon base)

1. Quasi-static mass term (khronon_mass_term.*): fails. Enough extra mass for cluster cores (1/mu ~ 1 Mpc) multiplies isolated-galaxy
   lensing at 0.3-1 Mpc by 5-500x; for 1/mu <~ 0.5 Mpc no static solution exists (the infrared instability). MOND alone gives ~1/3 of
   the measured cluster pull.

2. Fluid that clumps only in clusters? (condensate_jeans.*): fails against Lyman-alpha.
   Khronon's fluid stops self-clumping below a crossover length lambda_* = 2 pi c_ad / sqrt(4 pi G rho) (pressure cancels self-gravity).
   | c_ad^2 | crossover today | P(k)/P_LCDM at z=3, k = 2 / 5 / 10 per Mpc |
   |---|---|---|
   | 1e-11 | 0.14 Mpc | 1.000 / 0.998 / 0.992 |
   | 1e-9  | 1.4 Mpc  | 0.971 / 0.843 / 0.581 |
   | 1e-8 (AeST paper's Cosh value, used in all our fits) | 4.4 Mpc | 0.772 / 0.365 / 0.134 |
   | 1e-7  | 14 Mpc   | 0.283 / 0.063 / 0.018 |
   A cluster-only window (crossover 1-4 Mpc) removes 15-60% of the small-scale power the Lyman-alpha forest sees at z ~ 3; the
   warm-dark-matter bound (m > ~5 keV) allows only ~2% at k = 5/Mpc. So c_ad^2 <~ 1e-10: the fluid clumps like cold dark matter
   down to ~0.3 Mpc.

3. Consequence for our fits: the paper's Z0 = 1e-9 is Lyman-alpha-excluded. Re-evaluating at Z0 = 1e-12 (same parameters):
   all chi2 rise by +1.7 to +2.2, but the comparisons hold: beta = 1/2 vs Lambda -6.7 (no exchange), -3.8 (exchange).

4. Open problem this exposes: with Lyman-alpha-safe settings the fluid should gather around galaxies like cold dark matter, yet our
   galaxy fits (and MOND) need no halos. What an irrotational fluid (it cannot virialise like particles) does nonlinearly is not
   worked out for Khronon or AeST in the literature. This is now the key question for the model.

5. Can the fluid's own pressure keep it out of galaxy halos but let it into clusters? (jeans_window.*)  No.
   For any pressure law P = kappa rho^Gamma (Gamma = 2: dense channel gas = Khronon's quadratic K; Gamma = 3: dilute channel gas),
   pressure grows toward the past, so the early-universe bound (w <= 0.016 at a = 10^-4.5) caps kappa. The largest allowed Jeans mass
   at halo densities is 3e4-3e5 Msun (Gamma = 2) or ~1e-9 Msun (Gamma = 3), far below galaxy masses (1e12). Pressure that depends
   only on density cannot do it; the fluid falls into galaxy halos too.

6. Literature (summary of a dedicated search): no simulation or analytic model of halo formation exists for AeST, Khronon, mimetic DM
   or "dust of dark energy". A single-valued (irrotational) flow forms caustics in finite time (Babichev & Ramazanov 1704.03367);
   proposed completions (higher derivatives, complex/wave field, pressure) all change the theory. Khronon's MOND limit is derived for
   stationary systems only (2404.06584); how infalling fluid ends up in that stationary state is unknown.

Status: the deciding question is a dynamical one. Either the fluid builds halos around galaxies (then halos + MOND double-count and
the galaxy successes fail), or Khronon's nonlinear dynamics turn infalling fluid into the stationary MOND configuration (then the
fluid's energy must go somewhere). Needs a spherical-collapse solution of the full Khronon equations (weak-field: Newtonian gravity +
the khronon's Hamilton-Jacobi flow with K(Q) and J(Y)).

7. Khronon's DBI fluid (dbi_cosmology.*, dbi_largelam.*, dbi_window.*, dbi_phase_change.*, dbi_scan2d.*, table: dbi_scan_table.md)
   - The paper's stated setting (1/mu = 22.3 Mpc, lambda_D ~ 1) does NOT give LCDM in our solver (sigma8 0.22, CMB off by >100%);
     our implementation reproduces the paper's own c_ad^2 formula and its quadratic-K bound (1/mu <~ 0.2 kpc). LCDM-like behaviour
     needs lambda_D ~ 1e7-1e12 (depending on mu), which keeps the fluid on its pressureless branch until late.
   - Identity for any K(Q) (w << 1): c_ad^2 mu_eff^2 = K'/2Q ~ 4 pi G rho, so the cosmological crossover length today = 2 pi / mu_eff,
     the same number that sets the extra pull around galaxies. Cosmology (sigma8, Lyman-alpha) and galaxy lensing pull in opposite
     directions; a narrow window remains (1/mu_eff ~ 15-25 Mpc).
   - Phase change: in a static well the khronon is pushed toward the DBI limit, x = x0 + |phi|/c^2 -> x_lim. Past the critical depth
     v_crit = c sqrt(x_lim - x0) the fluid's pressure saturates and it becomes collapsing dust (dark-matter-like); shallower wells keep
     it smooth (no halo, MOND only). Inside the cosmology window v_crit ~ 500-700 km/s: Milky-Way-like galaxies (central depth ~500
     km/s with our MOND potentials) stay fluid; groups (~1200) and clusters (~2000-2800) collapse. Massive spirals (~670) are borderline.
   - Best point found: 1/mu = 300 Mpc, lambda_D = 9e6: v_crit 694 km/s, sigma8 0.790, Lyman-alpha 0.979 (at the ~2% limit), galaxy
     lensing +5% at 1 Mpc.
   Status: PLAUSIBLE, not established. Caveats: sits on the Lyman-alpha limit; well depths depend on where the potential is zeroed
   (3 Mpc here; MOND potentials grow logarithmically); the runaway past x_lim is inferred, not simulated; two free parameters (mu,
   lambda_D) are chosen, not derived. It is the first mechanism in this project that gives "smooth in galaxies, dark matter in
   clusters" from one fluid, and it refines lesson 5: density-only power laws fail, but a saturating (DBI) pressure can work.

8. Time-dependent check (collapse_sim.*): dark fluid as Lagrangian shells with the DBI pressure, in static MOND wells of baryons,
   starting at the cosmic mean density at rest, 10 Gyr.
   | well | 1/mu = 300 Mpc, lambda_D 9.2e6 (v_crit 694) | 1/mu = 22.3 Mpc, lambda_D 1e10 (v_crit 514) |
   |---|---|---|
   | Milky-Way-like galaxy (~500 km/s deep) | settles smooth, centre 1.3x mean (no halo) | runaway after 1.9 Gyr |
   | massive spiral (~670) | runaway after 1.3 Gyr | runaway after 0.1 Gyr |
   | group, clusters (1200-2800) | runaway within 0.1-0.15 Gyr | runaway within 0.1 Gyr |
   The dynamics follow the static critical depth: the phase change is real in time evolution, not just a static argument.
   Limits of this test: baryons are a fixed, already-formed well; the fluid reservoir is today's mean density within 3 Mpc
   (~4e12 Msun), far less than a cluster gathers from its ~10 Mpc formation region, so cluster masses can't be read off yet; the
   run stops at the runaway. Needed next: cosmological spherical collapse (expanding background, baryons and fluid together).

9. Epoch problem (vcrit_vs_z.*, cosmo_collapse.*)
   The critical depth is not fixed: the background fluid was denser in the past, i.e. closer to the DBI limit, so
   v_crit(z) (km/s), 1/mu = 300 Mpc / lambda_D 9.2e6:  z=0: 694 | 0.5: 208 | 1: 88 | 2: 26 | 3: 11 | 6: 2.
   Galaxy halos assemble at z ~ 1-3, when v_crit was 10-90 km/s: the fluid then falls into every galaxy as pressureless dust, like
   cold dark matter. The cosmological shell run (exploratory; purely radial orbits, so baryon radii are unphysical) agrees: galaxy-mass
   objects end up holding dark mass at or above the cosmic share.
   The mechanism can only give MOND-only galaxies today if halos captured early "re-melt": as v_crit rises after z ~ 0.5, a halo
   whose well is shallower than v_crit loses pressure balance and should expand back out. Whether that happens is a field-theory
   question (the fluid has formed caustics by then) that a shell code with permanent conversion cannot answer.
   If re-melting happens, it is a sharp prediction: galaxies carried dark-matter halos until z ~ 0.5 and lost them since, while
   groups and clusters kept theirs. Galaxy-galaxy lensing vs redshift and high-z rotation curves can test it.
   Status: downgraded from 'plausible' to 'plausible only with re-melting (unproven)'.

10. Re-melting test with KiDS-1000 lensing (remelt_kids_test.*) -> FAILS.
   At the KiDS lens epoch (z ~ 0.25) v_crit = 359 km/s. A captured halo deepens its own well: abundance-matched halos of the two
   massive bins have depths 461 and 649 km/s, so they are predicted to still be there (and to stay until v_crit exceeds their own
   depth, i.e. essentially until today). With those halos plus MOND (what Khronon gives), the predicted pull is ~5x the measured one:
   chi2 10,715 and 32,553 (15 points each) vs 46 and 61 for MOND alone. Over all four bins: MOND-only 202, MOND+halo 50,514
   (crude LCDM reference: 720; not a fair LCDM test: point-mass baryons, no two-halo term).
   Lesson: halos are self-sustaining, since a captured halo keeps the well deep enough to hold itself. With cosmology-allowed parameters
   (v_crit today <~ 700 km/s) massive galaxies would keep halos to the present, contradicting both KiDS and SPARC.
   Consequence: the DBI phase change cannot make galaxies halo-free. More generally, any relativistic-MOND theory whose dark-matter-like
   cosmological fluid is cold when galaxies form (AeST, Khronon) must explain why that fluid does not build galaxy halos, which
   (with MOND acting on top) KiDS excludes by a wide margin. This is now the central open problem for this class of theories.
