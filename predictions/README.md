# What the grid predicts that standard cosmology does not (checked Oct 2026)
1. Dark-energy law rho_DE ~ adot^-1/2 (law_signature.*, fit_*_beta0.5.txt)
CORRECTION (2 Oct 2026): the first version of this section used d ln rho_DE/d ln a = q, which is rho_DE ~ adot^-1 (beta = 1),
not the law (beta = 1/2). The CAMB/CLASS pipelines behind the main README always used beta = 1/2; only fit_law.py,
law_signature.py and Figure 4 had the error. Fixed; the beta = 1 numbers are kept in fit_*_beta1.txt for the record.
- Shape, no free number: w < -1 while the expansion slows, w > -1 once it speeds up; the crossing is EXACTLY at the start of
  acceleration (z = 0.66-0.72 for Om 0.32-0.30). Early on w -> -1 - 1/12 = -1.083 (matter era); today w0 = -0.91.
  w(z), Om 0.31: z 0: -0.91 | 0.3: -0.95 | 0.5: -0.98 | 1: -1.03 | 2: -1.06 | 3: -1.08.
- Mimicked by DESI's two-number form as (w0, wa) = (-0.90, -0.25): between Lambda and the free fits below.
- Refit to DESI DR2 BAO + Planck distance priors + each supernova set (beta = 1/2):
| supernovae | free (w0, wa) best fit | law vs LCDM dchi2 | law vs free (w0, wa): dchi2 / dAIC / dBIC |
|---|---|---|---|
| Pantheon+ | -0.86, -0.50 | -5.4 | +1.7 / -2.3 / -13.1 |
| DES-Dovekie | -0.82, -0.62 | -7.2 | +3.1 / -0.9 / -12.0 |
| Union3 | -0.69, -0.96 | -6.8 | +6.4 / +2.4 / -0.9 |
  With no dark-energy parameter, the law beats LCDM by 5.4-7.2 in chi2 and beats the two-parameter fit by information
  criteria in 2 of 3 sets (AIC) and all 3 (BIC). (For comparison beta = 1 gave -2.7 / -5.1 / -10.0.) Growth: sigma8 0.6%
  below LCDM with the same early universe.
  Simplified vs exact deceleration term (EQUATIONS.md E4; fit_*_beta0.5_exactq.txt): exact form gives -5.4 / -6.9 / -6.3,
  i.e. the simplification shifts results by 0.1-0.55 in chi2; the crossing redshift is unchanged.
  Oct 2026 fix (independent check): the sound horizon was integrated from the first grid point above z* instead of z*;
  fixed and l_A recalibrated to CAMB. Claim-1 numbers move by <= 0.1 (fit_*_rsfix.txt: -5.4 / -7.2 / -6.8; exact -5.4 / -6.9 / -6.3).
- New since our fits: DES Y6 (Jan 2026, arXiv:2601.14559): 3x2pt S8 = 0.789 +/- 0.012, 2.6 sigma below the CMB in LCDM;
  joint wCDM w = -0.981 +/- 0.022 (constant w only, not a test of the crossing). The law's 0.6% lower sigma8 goes the right
  way but covers only a small part of that gap. No public DES Y6 likelihood used here.
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

4. Slower light in the less-stretched early grid (Coalesce, 3 Oct 2026) -> the fine-structure constant alpha = e^2/(4 pi eps0 hbar c)
   would be LARGER early (alpha ~ 1/c if only c changes). Same data as item 3 (Planck + ACT DR6 + SPT-3G D1 + DESI DR2 +
   DES-Dovekie), alpha at recombination fixed at 0.992-1.008, everything else re-minimised (spa/spa_al*.json):
| alpha_rec / alpha_0 | 0.992 | 0.996 | 1.000 | 1.004 | 1.008 |
|---|---|---|---|---|---|
| dchi2 | +65.4 | +25.3 | 0 | -5.4 | +8.7 |
   Parabola: alpha_rec/alpha_0 = 1.0031 +/- 0.0013 (2.4 sigma), dchi2 -5.2  ->  light 0.31% +/- 0.13% slower at recombination
   (if c is what changed). The data prefer the direction the idea predicts. alpha and the electron mass both set atomic
   energies (~alpha^2 m_e), so this is partly the same signal as the early-electron shift (item 3), not independent evidence.
   Late-time limits: quasar absorption lines |d alpha/alpha| < ~1e-6 at z ~ 1-4; atomic clocks < ~1e-17 per year today
   -> any change must have finished long before z ~ 4 (as for the electron switch at z ~ 100-200), so cells cannot simply
   stretch with the expansion (same conclusion as Fig. 12d).
