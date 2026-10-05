# PREREG 92 -- The universe without dark energy: where are its gaps, and does our model close them?
Committed BEFORE running iter92_gaps.py.

## Question (user's framing)
Take the dark-energy-free universe. Wherever it disagrees with the real universe, list the gap
and size it. Then check which gaps our model (Claim 1 law, beta = 1/2) closes.

## Dark-energy-free universes tested (2 free choices, fixed in advance)
- EdS : flat, matter only (Einstein-de Sitter). Free: h, omega_b.
- OPEN: matter + negative curvature, no dark energy. Free: Om, h, omega_b.
References: LCDM and LAW (flat), same data, same code (fit_law.py machinery).

## Observables (fixed in advance; 5)
1. Age of the universe vs oldest stars: 13.3 +/- 0.5 Gyr (globular clusters, Valcin et al. 2020, stat+sys).
2. Supernova distance SHAPE (Pantheon+, magnitude marginalised), fitted alone.
3. BAO (DESI DR2).
4. CMB distance priors (Planck 2018 compressed).
5. Growth of structure: sigma8 today for the same early-universe amplitude, linear growth only.
   Scaled from Planck LCDM sigma8 = 0.811 by the ratio of growth factors (transfer-function shape
   change ignored -- approximation disclosed). Compared to measured S8 ~ 0.78-0.83.

## Expectations
- E1 Age: EdS gives ~9-10 Gyr, a ~3-4 Gyr gap (>5 sigma). OPEN somewhat older but still short when other data included.
- E2 SN shape alone: EdS worse than LCDM by delta chi2 > 100. OPEN closes part of it but stays worse by > 25.
- E3/E4 BAO+CMB combined: no dark-energy-free model fits; delta chi2 > 500 vs LCDM.
- E5 Growth: EdS overshoots sigma8 by 20-30% (no late slowdown of clustering).
- E6 LAW closes every one of these gaps to the same level as LCDM (within delta chi2 ~ 10 on combined data),
     with its known small leftovers (H0, S8 slightly worse) unchanged.
Look-elsewhere: 5 observables x 2 variants; nothing here can PROVE our law, only check it closes the gaps.
