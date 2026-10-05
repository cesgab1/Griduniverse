# PREREG 102 -- where is the tip of the Bullet Cluster's curvature, from visible matter only? (Coalesce)
Committed BEFORE running iter102_tip.py. No dark matter, no dark energy, nothing invented.
Data (Clowe et al. 2006, ApJ 648 L109, Table 2; 100 kpc apertures; z = 0.296, 4.413 kpc/arcsec):
| aperture | RA | Dec | gas M_X (1e12) | stars M* (1e12) | lensing kappa |
| main BCG | 06:58:35.3 | -55:56:56.3 | 5.5+/-0.6 | 0.54+/-0.08 | 0.36+/-0.06 |
| main plasma | 06:58:30.2 | -55:56:35.9 | 6.6+/-0.7 | 0.23+/-0.02 | 0.05+/-0.06 |
| sub BCG | 06:58:16.0 | -55:56:35.1 | 2.7+/-0.3 | 0.58+/-0.09 | 0.20+/-0.05 |
| sub plasma | 06:58:21.2 | -55:56:30.0 | 5.8+/-0.6 | 0.12+/-0.01 | 0.02+/-0.06 |
Transcription check passed: positions give BCG-BCG 722 kpc (published ~720) and main BCG-plasma 47" (published ~47").
Model: gas = 2 Plummer blobs at the plasma peaks, stars = 2 Plummer blobs at the BCGs; masses and sizes fitted to the 4 aperture
values of each (gas 4 numbers -> 4 unknowns; stars same). Collision assumed in the sky plane.
Gravity A (Newton/Einstein): kappa = Sigma / Sigma_crit (sources z ~ 1; Sigma_crit computed, disclosed).
Gravity B (galaxy rule, QUMOND on a 3-D grid, a0 = 1.2e-10): extra 'phantom' curvature computed from the visible matter, projected.
Note on lensing kappa: maps carry an unknown constant offset (mass-sheet); compare DIFFERENCES between apertures (offset cancels).
Expectations:
- Visible matter puts the tip at the GAS (main plasma highest), not at the galaxies. Measured tip is at the main BCG.
- Key difference kappa(main BCG) - kappa(main plasma): measured +0.31 +/- 0.085; Newton-visible ~0 (|..| < 0.03) -> > 3.5 sigma off.
- Galaxy rule: tip stays near the gas/baryon centre (within ~100 kpc of the main plasma), difference < 0.10 -> still > 2.5 sigma off.
