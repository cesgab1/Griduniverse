# The 'ocean' hybrid: MOND from baryons (through the geometry) + an unboosted fluid halo of the cosmic ratio

Rule (one rule for all systems): g_obs = nu(g_bar/a0) g_bar  +  Newtonian pull of a fluid halo of mass R x M_baryons,
R = Omega_fluid/Omega_b = 5.4 (cosmic ratio, NOT fitted), halo with a superfluid-like core. Lensing sees both terms
(the MOND push acts through the geometry light follows); the halo is NOT boosted by MOND (only baryons source the push).

Results (hybrid_test.*, hybrid_cosmic_ratio.*):
| data | MOND only | + cosmic-ratio halo, core 100 kpc | + core 30 kpc | + NFW |
|---|---|---|---|---|
| KiDS-1000 lensing, 4 mass bins, chi2 / 60 pts | 202 | 129 | 153 | 175 |
| KiDS best-fit R (2 sigma) | - | 11.7 (8.9-14.6) | 4.9 (3.3-6.6) | 3.5 (2.2-5.0) |
| Cluster pull / measured at 0.1, 0.3, 0.6, 1 Mpc | 0.35 0.35 0.34 0.33 | 0.92 1.13 0.94 0.79 | 1.47 1.15 0.94 0.80 | 1.50 1.15 0.94 0.80 |
| SPARC rotation curves, chi2 / 2788 pts | 4568 | 4594 (+26) | 4994 (+426) | - |

With a ~100 kpc core, one cosmic-ratio rule keeps SPARC galaxies essentially unchanged, improves KiDS lensing, and supplies the missing
cluster mass. Contrast: the same halo with MOND acting on it too (Khronon/AeST as they stand) overshoots KiDS 5x (clusters/README §10);
a full CDM-size galaxy halo (abundance matching) is excluded (chi2 755-1795).

Caveats
- Core size: 3 values tried (30, 100 kpc, none); 100 kpc preferred. Superfluid cores (polytrope) give core radii ~ (K M)^(1/5); to derive.
- 'Fluid = 5.4 x present baryons' in galaxies means galaxy halos are ~5x lighter than CDM abundance matching: fluid and baryons are
  retained together. Mechanism not derived.
- Inputs: crude KiDS baryons (median stellar mass per bin, rough cold-gas fractions), one model cluster, SPARC M/L 0.5.
- Theory: needs a relativistic theory where visible matter and light share a metric that carries the MOND push sourced by baryons
  only, while the dark fluid couples to the plain metric (so it pulls Newtonianly and does not source MOND). TeVeS-like
  two-metric structure inside a Khronon/AeST-safe framework (GW speed, Solar System). Not built yet.
