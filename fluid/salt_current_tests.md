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
