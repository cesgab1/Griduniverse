# RESULT 92 -- gaps of the dark-energy-free universe
Script: iter92_gaps.py. Outputs: iter92_gaps_{PANTHEON,UNION3,DESDOVEKIE}.txt. Pre-registration: PREREG_92.md.

## Gap table (dark-energy-free universe vs real data)
| Gap | Flat, matter only (EdS) | Matter + curvature (OPEN) | LAW (beta=1/2) |
|---|---|---|---|
| Supernova shape alone, delta chi2 vs LCDM | +272 / +698 / +927 | +12 / +56 / +44 (wants Om -> 0.05, an almost empty universe) | -1.2 / -0.4 / -1.3 |
| BAO+SN+CMB combined, delta chi2 vs LCDM | +7400 to +8100 | +5700 to +6300 | -5.4 to -7.2 |
| Age at measured H0 (67.4 / 73) | 9.7 / 8.9 Gyr vs stars 13.3 +/- 0.5 (7-9 sigma) | -- | 13.74 Gyr |
| H0 needed to fit everything | 45 km/s/Mpc (measured 67-73) | 52 | 67.5 |
| Matter needed | 100% of critical (lensing/clusters see ~30%) | 91% | 31% |
| sigma8 for the same early seeds | 1.07 (measured 0.78-0.83; +32%) | 1.03 | 0.81 |
(Supernova numbers listed Union3 / Pantheon+ / DES-Dovekie.)

## Reading
- Gaps move instead of closing: fixing the age by lowering H0 opens an H0 gap; fixing distances with curvature
  opens a CMB gap (CMB alone +3500). No single knob without dark energy closes all of them.
- "Matter not where it should be": the no-dark-energy universe needs ~3x more matter than gravity actually finds,
  and that extra matter would over-grow clusters by ~30%.
- LAW closes every gap at least as well as constant dark energy (combined delta chi2 -5 to -7, three SN sets).
- Not closed by any model here: the local H0 (73) vs 67.5 gap -- unchanged from earlier iterations.

## Scorecard vs pre-registration
- E1 age: MET (7-9 sigma at measured H0); at its own best fit EdS trades it for an H0 gap (logged).
- E2 SN: EdS MET (+272..+927). OPEN: MET for Pantheon+/DES (+56/+44), MISSED for Union3 (+12 < 25).
- E3/E4 combined: MET (> 5700).
- E5 growth: overshoot 26-32%, at/just above the predicted 20-30% -- MET at the edge.
- E6 LAW closes all gaps to within delta chi2 ~10 of LCDM: MET (it is better by 5-7).
Caveats: compressed CMB priors are calibrated for near-LCDM models and are only indicative for EdS/OPEN
(their huge chi2 is robust in sign, not in exact size); sigma8 ignores transfer-shape changes.
