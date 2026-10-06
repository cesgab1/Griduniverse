# PREREG 111 -- gamma-ray-burst distances: our law vs constant dark energy (and no-dark-energy universes)
Committed BEFORE running iter111_grb.py.
Data: A118 long-GRB sample (Khadka et al. 2021, arXiv:2105.12692, Table 7 = ar5iv Table 3), z = 0.34-8.2. Read from the paper's
HTML in a browser; all five column sums of grb/A118.csv match the page exactly (118 rows, z 285.3842, Ep 121079.538,
Ep_err 16353.12, S 787.947, S_err 47.585).
Method (same as the paper, no low-z calibration, so no circularity): Amati relation log Ep = a + b log(Eiso/1e52 erg),
Eiso = 4 pi D_L^2 S_bolo / (1+z); scatter^2 = s_int^2 + (sig log Ep)^2 + b^2 (sig log Eiso)^2. Fit a, b, s_int jointly with the
cosmology; H0 is absorbed by a. Report -2 ln L for each model at its best fit, and Om.
Models: LCDM (Om free), LAW (beta = 1/2, Om free) -- main comparison. Secondary: EdS (Om = 1, no dark energy) and
BARE-OPEN (atoms only, Om = 0.05, curvature, no dark energy). 4 comparisons -> look-elsewhere x4.
Free choices: Amati form (1), A118 sample (1), no calibration (1).
Expectations (from grb/power.py and memory of the paper):
- transcription check: b ~ 1.1-1.3, s_int ~ 0.35-0.45 dex. Outside -> suspect data, stop.
- LAW vs LCDM: |delta(-2 ln L)| < 1 -> GRBs cannot tell them apart (the power check said ~10-25x too weak).
- Om from GRBs alone: loose (tens of percent), consistent with 0.3 within ~1-2 sigma.
- EdS and BARE-OPEN: guess, disfavoured by less than ~3 sigma -- GRBs are weak; NOT a rescue of no-dark-energy if so, since
  supernovae + BAO + CMB already exclude them (iterations 92-93).
