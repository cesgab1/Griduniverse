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
