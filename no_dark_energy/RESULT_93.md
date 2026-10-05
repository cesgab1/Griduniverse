# RESULT 93 -- the bare universe (atoms + light + neutrinos; no dark energy, no dark matter)
Scripts: iter93_bare.py, iter93_at_h.py. Outputs: iter93_bare_{PANTHEON,UNION3,DESDOVEKIE}.txt, iter93_at_h.txt.

## Non-dark lessons plugged in
Cells widen (iter 84), no light fading (iter 91), constant G (iter 87). None of them changes the expansion
history at all -- every effect our grid has on cosmology so far runs through dark energy. Flatness (iter 81) used as a check.

## Results (BARE-OPEN = real atom count, curvature fills the rest; BARE-FLAT = atoms forced to 100%)
| Check | BARE-OPEN | BARE-FLAT | Data |
|---|---|---|---|
| Atoms (omega_b) | 0.0224 (matches BBN) | 0.45 (20x BBN, 863 sigma) | 0.0222 +/- 0.0005 |
| Flatness Omega_k | 0.95 (729 sigma) | 0 | 0.002 +/- 0.0013 |
| Age at H0 67.4 / 73 | 13.5 / 12.6 Gyr -- FITS | 9.7 / 8.9 Gyr (7-9 sigma) | stars 13.3 +/- 0.5 |
| Supernova shape, delta chi2 vs LCDM (Union3/Pantheon+/DES) | +10 / +50 / +36 | +269 / +697 / +926 | |
| BAO (13 pts) chi2 at H0 67.4 | ~11,700 | ~162,000 | LCDM ~12-31 |
| CMB acoustic angle l_A at H0 67.4 | 1881 (spots ~6x too small) | 488 | 301.46 +/- 0.09 |
| Lump growth since recombination | x78 -> lumps today ~1e-3 | x1000 -> ~0.01-0.1 | need ~1 (galaxies exist) |
Not computed, known to need dark matter: flat rotation curves, cluster lensing, CMB third peak.
Caveats: compressed CMB/BAO formulas were built for near-LCDM universes; sizes of the failures are indicative,
signs are not in doubt. BARE-OPEN's combined fit runs to the h = 0.30 edge; SN-only fits run to h = 1.0 (emptier preferred).

## Reading
- The bare universe gets the atom count and the age right, and nothing else.
- It has no galaxies: without dark matter, lumps grow ~100x too little. Without dark energy, the supernova
  shape and BAO distances are wrong. Space would be strongly curved, not the near-flat we measure.

## Scorecard vs PREREG_93: all six expectations MET (SN 10-60: +10..+50; BAO > 10 sigma; l_A > 100 sigma;
age closes; flatness > 500 sigma and BARE-FLAT 20x atoms; growth x78 in 30-100).
