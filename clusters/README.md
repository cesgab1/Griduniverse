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
