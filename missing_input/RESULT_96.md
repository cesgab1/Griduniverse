# RESULT 96 -- missing measured input linking galaxies and clusters?
Outputs: iter96_missing_input.txt (pre-registered), iter96_robust.txt (POST-HOC, added after seeing the result).

## Pre-registered analysis (least squares over 165 SPARC galaxies)
| input | SPARC slope | needed x5.3 / x17 | tension x5.3 / x17 |
|---|---|---|---|
| baryonic mass | +0.13 +/- 0.04 | 0.16 / 0.27 | 0.8 / 3.3 sigma |
| well depth v^2 | +0.40 +/- 0.06 | 0.33 / 0.56 | -1.1 / 2.5 sigma |
| size | +0.36 +/- 0.10 | 0.32 / 0.55 | -0.4 / 2.0 sigma |
| gas fraction | -0.01 +/- 0.12 | 3.7 / 6.3 | 30 / 51 sigma |
| surface density | -0.10 +/- 0.11 | 3.9 / 6.6 | 38 / 64 sigma |
EXPECTATION MISSED for mass/depth/size: I predicted slopes ~0 and > 4 sigma exclusion; instead the pre-registered fit is
COMPATIBLE with the x5.3 cluster target (and 2-3.3 sigma short of x17). Gas fraction and surface density: excluded as expected
(clusters are not unusual galaxies in these inputs).

## Post-hoc robustness (not pre-registered; weigh accordingly)
- Quartile medians of a0 vs well depth: 0.57 / 1.17 / 1.08 / 1.21 e-10. The whole trend is the SMALLEST quarter of galaxies;
  the upper three quarters are flat.
- Theil-Sen (outlier-resistant): depth slope +0.24 +/- 0.06 -> clusters x3.4 (1.5 sigma short of x5.3, 5.5 sigma of x17).
- Upper 75% only: slopes +0.01 to +0.03 +/- 0.05-0.10 -> clusters x1.1 for every input.
Reading: the signal is a low-a0 tail in the smallest, gas-dominated galaxies (where fixed M/L, pressure support and distance
errors matter most), not a smooth power law that grows toward clusters. Among normal-to-large galaxies a0 does not rise with
mass, depth or size, so extrapolating to clusters gives x1.1, not x5-17.

## Verdict
No measured input tested here (5) tunes galaxies and clusters with one equation, on the robust reading. The pre-registered
reading leaves "a0 grows with well depth" (EMOND, Zhao & Famaey 2012) formally open at the x5.3 level. Clean way to settle it:
galaxy GROUPS (depth between galaxies and clusters) -- a depth law predicts a0 x2-3 there; the flat reading predicts x1.
