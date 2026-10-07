# PREREG 132: does the extra-pull threshold grow with the universe's expansion rate?
Committed before data and code.
Coalesce's compression picture: gravity = slope of space's compression; the baseline is set by the whole universe,
so departures from Newton appear where the local slope ~ the cosmic background. If the background tracks the
expansion rate H, the threshold acceleration should scale as a_t(z) = a_t0 * H(z)/H0.
Models compared:
  C (compression): a_t(z) = a_t0 * H(z)/H0   (x2.97 at z = 2 for Om = 0.3)
  K (constant):    a_t(z) = a_t0              (no evolution)
  a_t0 = 1.2e-10 m/s^2 (local galaxies, SPARC; same curve family as our boost map).
Data (to collect): per-galaxy rotation measurements at z ~ 0.6-2.6 with baryonic mass (stars + gas) and circular
velocity at a stated radius (e.g. RC100, Nestor Shachar et al. 2023; KMOS3D / Genzel et al.). Also any published
measurement of a_t or of the baryonic Tully-Fisher zero point vs redshift.
Method: g_bar = G M_b(<R)/R^2 (enclosed baryons as published; if only totals, half inside R_e -- stated),
g_obs = V_c^2/R. Curve g_obs = g_bar / (1 - exp(-sqrt(g_bar/a))). Per galaxy fit a; median log a in redshift bins
with bootstrap errors; compare to C and K by chi^2. Prefer a model only if delta chi^2 > 4; if both miss by > 2 sigma,
report 'neither'. Published BTFR zero-point evolution compared as a second, independent check
(C predicts at fixed V: baryonic mass lower by H(z)/H0, i.e. -0.47 dex at z = 2).
Conventions/free choices (3): stellar-mass IMF (Chabrier vs Salpeter, x1.7), gas-mass recipe, radius used.
Expectation (written now): the two models differ by ~0.1-0.2 dex in g_obs at the relevant g_bar; per-galaxy scatter
~0.2-0.3 dex and mass conventions ~0.2 dex -> likely INCONCLUSIVE or weak (< 3 sigma) preference; high-z galaxies
are compact (g_bar >~ a_t0), limiting sensitivity.
