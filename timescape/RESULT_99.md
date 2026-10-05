# RESULT 99 -- uneven expansion (timescape) vs Lambda vs our law, late-universe data
Outputs: iter99_{PANTHEON,UNION3,DESDOVEKIE}.txt. Self-check passed (low-z distance = cz/H0 to 0.4%, as expected).

| | Pantheon+ | Union3 | DES-Dovekie |
|---|---|---|---|
| Timescape - LCDM, supernovae only | +5.1 | -1.4 | +2.0 |
| Timescape - LCDM, BAO shape only | +4.6 | +4.6 | +4.6 |
| Timescape - LCDM, SN + BAO | +28.9 | +2.5 | +31.6 |
| Our law - LCDM, SN + BAO | -4.9 | -6.1 | -6.4 |
Best void fraction: supernovae 0.86 / 0.77 / 0.85, BAO 0.69. The two data sets want different amounts of void, which is what
costs timescape in the combination.

Scorecard vs PREREG_99:
- E1 SN within |5| of LCDM: MET for Union3 and DES; Pantheon+ +5.1, at the line (marginal MISS).
- E2 f_v0 0.70-0.85: Union3 MET; Pantheon+/DES 0.856/0.853 just above (marginal MISS).
- E3 BAO worse by 5-30: +4.6, slightly better than predicted (MISS, small).
- E4 combined worse than LCDM and LAW by > 5: MET for Pantheon+ and DES (+29, +32); MISSED for Union3 (+2.5 vs LCDM; +8.6 vs LAW).
Reading: uneven expansion fits supernovae about as well as Lambda (as its authors report), but supernovae and BAO disagree on how
much void there must be; combined, it is disfavoured by ~29-32 in chi2 on two of three sets and roughly tied on Union3.
Our law is the best of the three on every combination.
Caveats: BAO D_H = dD_M/dz is an approximation for timescape; no CMB used; d_L formula from arXiv 1107.5596 eq 7, 1+z(t) from
memory (verified at t0 and by the low-z check); published timescape analyses (e.g. Lane et al. 2024) use different cuts.
