# PREREG 105 -- hemispheres toward / away from the candidate centre direction (l,b) = (260, +12) (Coalesce)
Committed BEFORE running iter105_hemi.py.
Data: Pantheon+ (zHD, STAT+SYS covariance). Main cut z > 0.05 (beyond the local flow found in iteration 103); check with z > 0.01.
Each hemisphere gets its own magnitude offset (so only the SHAPE of distance vs redshift is compared).
Significance: the same statistic for 300 random axes (handles uneven sky coverage) -> p-value for the candidate axis.
PART A (first, per Coalesce): NO dark energy, NO dark matter -- ordinary matter only (Omega_b = 0.049, curvature fills the rest;
no free shape parameter). Per hemisphere: penalty vs the best LCDM fit of that hemisphere, per supernova.
Statistic S_A = penalty/SN toward - penalty/SN away.
PART B: dark energy -- LCDM matter share Om fitted per hemisphere (dark-energy share = 1 - Om), and our law likewise.
Statistic S_B = Om(toward) - Om(away), with its own error, and vs random axes.
Expectations:
- A: the bare universe fails in BOTH hemispheres (needs a speed-up everywhere); S_A not unusual vs random axes (p > 0.05).
- B: Om toward and away agree within 2 sigma; p > 0.05 vs random axes; same for the law.
A centre in that direction would show as a hemisphere difference that is rare among random axes (p < 0.05) in A and/or B.
