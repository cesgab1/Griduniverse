# What the grid predicts that standard cosmology does not (checked Oct 2026)
1. Dark-energy law rho_DE ~ adot^-1/2 (law_signature.*, fit_*.txt)
- Shape, no free number: w < -1 while the expansion slows, w > -1 once it speeds up; the crossing is EXACTLY at the start of
  acceleration (z = 0.70-0.76 for Om 0.32-0.30). Early on w -> -1 - 1/6 = -1.167 (matter era); today w0 = -0.82.
  w(z), Om 0.31: z 0: -0.82 | 0.3: -0.90 | 0.5: -0.95 | 1: -1.05 | 2: -1.13 | 3: -1.15.
- Mimicked by DESI's two-number form as (w0, wa) = (-0.80, -0.49).
- Refit to DESI DR2 BAO + Planck distance priors + each supernova set (our reproductions of the free (w0, wa) fit):
| supernovae | free (w0, wa) best fit | law vs LCDM dchi2 | law vs free (w0, wa): dchi2 / dAIC / dBIC |
|---|---|---|---|
| Pantheon+ | -0.86, -0.50 | -2.7 | +4.4 / +0.4 / -10.3 |
| DES-Dovekie | -0.82, -0.62 | -5.1 | +5.2 / +1.2 / -9.8 |
| Union3 | -0.69, -0.96 | -10.0 | +3.2 / -0.8 / -4.1 |
  The law sits where the data put dynamical dark energy, with no dark-energy parameters; by information criteria it ties
  (AIC) or beats (BIC) the two-parameter fit. Growth: sigma8 1.4% below LCDM with the same early universe (Om 0.31).
- New since our fits: DES Y6 (Jan 2026, arXiv:2601.14559): 3x2pt S8 = 0.789 +/- 0.012, 2.6 sigma below the CMB in LCDM;
  joint wCDM w = -0.981 +/- 0.022 (constant w only, not a test of the crossing). The law's 1.4% lower sigma8 goes the right
  way but covers about a third of that gap. No public DES Y6 likelihood used here.
- Decisive future test: the crossing redshift equals the acceleration onset. DESI DR3 / Euclid can measure w(z) in bins at
  z 0.3-1.2; a crossing well away from z ~ 0.7 (e.g. at 0.4) would rule the law out.
2. 21-cm step from the electron switch (step21.*)
- Switch at z 100-200 with a 0.4-1.1% heavier early electron: the line shifts by 2 delta, so the dark-ages spectrum jumps at
  7-14 MHz. Size: |jump| <= 0.7 mK, about 2% of the dark-ages signal (trough ~ -41 mK in this simple model); near z ~ 100
  it nearly vanishes. Feature shift 60-310 kHz.
- Reach: LuSEE-Night (lunar far side, launch planned 2026) is a pathfinder aiming at the dark-ages signal itself; a sub-mK
  step needs a future far-side array. Distinct prediction, not near-term.
3. Early electron 0.4-1.1% heavier (at recombination), tested with SPT-3G D1 (spa/, cosmology_fits/run_spa.py)
- Data: 'SPA' = Planck (low-l + cut high-l) + ACT DR6 + SPT-3G D1 T&E (SPTlite, candl 2.2, official spt_candl_data)
  + DESI DR2 BAO + DES-Dovekie SN. Electron mass at recombination fixed at 1.000/1.004/1.008/1.012/1.016, everything else
  re-minimised (2-3 starts per point; minimiser noise ~1-2 in chi2).
| m_e / today | best chi2 - LCDM | H0 |
|---|---|---|
| 1.000 | 0 | 68.2 |
| 1.004 | -2.0 | 68.6 |
| 1.008 | -3.8 | 69.3 |
| 1.012 | -3.8 | 69.8 |
| 1.016 | -2.2 | 70.4 |
  Parabola: m_e = 1.0099 +/- 0.0049 (2.0 sigma from today's value), Delta chi2 -3.8. Before SPT (P-ACT): 1.0078 +/- 0.0044.
- Reading: the prediction (0.4-1.1%) SURVIVES and the centre sits inside it; adding SPT-3G D1 (4% of the sky) moves the
  centre up slightly but does not tighten it. Not decisive: LCDM is 2 sigma away. Decisive: SPT-3G full depth (~25% of sky)
  and Simons Observatory. H0 rises to ~69.5: eases but does not solve the Hubble tension (SH0ES 73).
