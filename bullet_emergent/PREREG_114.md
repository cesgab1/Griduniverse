# PREREG 114 -- emergent gravity faces the Bullet Cluster (non-spherical version)
Committed BEFORE running iter114_bullet_eg.py.
Why: Verlinde's formula is spherical-only; published tests call the Bullet 'not a valid test case' (Tamosiunas+2019,
arXiv:1901.05505). Coalesce: build the non-spherical version and test it.
Version A (published root, no new choice): Hossenfelder 2017 covariant emergent gravity (arXiv:1703.01415). Its static
weak-field limit has an 'imposter' field with kinetic term chi^(3/2) -> extra potential phi_D obeying
   div( |grad phi_D| grad phi_D ) = 4 pi G a_V rho_visible,   a_V = c H0 / 6 = 1.13e-10 m/s^2,
total pull g = g_Newton + g_D (additive). A1: curl-free solution g_D = sqrt(a_V g_N) along g_N (fine grid, 20 kpc).
A2: full non-linear solution (convex energy minimised on a 40 kpc grid) vs A1 on the same grid -> size of the curl term.
Version B (our construction, 1 extra free choice): Verlinde's spread term made local,
   g_D^2 = a_V ( g_N + 4 pi G rho_visible r ),  r = distance from the system's centre of visible mass, along g_N.
Data and visible-matter model: exactly iteration 102 (Clowe+2006 Table 2; Plummer blobs fitted to aperture masses;
Sigma_crit for sources at z = 1). Same scoring: kappa(galaxy) - kappa(gas) at main and sub, and the location of the tip.
Measured differences: +0.31 +/- 0.085 (main), +0.18 +/- 0.078 (sub).
Expectations:
- A behaves like the galaxy rule (iteration 102: tip on the gas, combined miss 5.1 sigma): tip within 100 kpc of a gas peak,
  combined miss 4-6 sigma. Curl term changes kappa by < 20%.
- B: worse, because its extra term follows visible DENSITY (gas) times distance -> tip on the gas, miss > 5 sigma.
- If either version puts the tip on the galaxies with a miss < 2 sigma, emergent gravity passes its hardest test.
Free choices: version B's centre (1); A has none beyond iteration 102's. Look-elsewhere: 2 versions.
