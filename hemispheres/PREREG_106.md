# PREREG 106 -- is the iteration-105 hemisphere split real or a survey effect? (Coalesce: 'we hope to learn, not to be right')
Committed BEFORE running iter106_checks.py.
Check 1 (Pantheon+, same method as 105, z > 0.05): survey make-up of each hemisphere; then drop each survey holding > 5% of the
sample, one at a time, and recompute the dark-energy-share difference and its random-axis p (300 axes).
  Expectation: the split is driven by one or two surveys -- dropping one of them takes p above 0.2. (If the split survives every
  drop with p < 0.1, survey calibration is an unlikely cause.)
Check 2 (Union3 per-supernova inputs, independent light-curve calibration, RA/Dec available): simplified standardisation
  mu = mB + alpha x1 - beta c - M (alpha, beta, M fitted globally; diagonal errors from the mB/x1/c covariance + 0.10 mag scatter;
  one offset per sample), z_CMB > 0.05. Hemisphere matter shares and random-axis p (300 axes).
  Expectation: difference in the same sense as Pantheon+ but weaker, p > 0.1 (independent data should not reproduce a 2-sigma
  calibration fluke at full strength). Simplified method -> numbers indicative only.
DES-Dovekie: no sky positions in the public distance file -> cannot be tested (disclosed).
