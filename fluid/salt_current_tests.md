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
