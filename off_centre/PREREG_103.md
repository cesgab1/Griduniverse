# PREREG 103 -- are we off-centre in a dispersal? Direction test with supernovae (Coalesce)
Committed BEFORE running iter103_dipole.py.
Note first (thinking, not data): in a pure 'debris sorted by speed' dispersal (speed = distance/time from the centre), every piece
sees the same even Hubble law, so SPEED alone cannot reveal position. Position shows up through what differs with direction:
density (the hollow centre vs the dense shell) and therefore how fast the expansion slows or speeds up in each direction.
So the test is: does the supernova distance-redshift relation differ by direction (a dipole), and does it point the same way
at all distances and toward the quasar-count anomaly?
Data: Pantheon+ (1590 SNe, z > 0.01, STAT+SYS covariance), redshift zCMB (our motion removed using the CMB dipole);
repeat with zHD (also corrected for modelled local flows) as a check. Baseline shape: flat LCDM fitted Om (shape only).
Model per distance shell: residual = M + dipole; dipole changes distance by fraction A cos(angle to direction d).
Shells: z 0.01-0.05, 0.05-0.15, > 0.15; and all together. Generalised least squares with the full covariance.
Reference directions: CMB dipole (l 264.0, b 48.3); quasar-count dipole (Secrest+2021, l 238.2, b 28.8; from memory).
Expectations:
- z > 0.05 shells: amplitude A consistent with zero within 2 sigma; 95% limit < ~1.5%.
- z < 0.05: possible small dipole (local bulk flows), 1-3%, ~2 sigma, with zCMB; reduced with zHD.
- No consistent direction across shells. Off-centre picture predicts a significant, same-direction dipole in every shell.
Caveat: different surveys cover different parts of the sky, so calibration differences between surveys can fake a dipole.
