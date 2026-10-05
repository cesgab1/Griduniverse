# PREREG 99 -- uneven expansion (Wiltshire 'timescape') as dark energy, vs Lambda and our law
Committed BEFORE running iter99_timescape.py.
Model: timescape tracker solution, one parameter f_v0 (today's void volume fraction). Formulas:
  dressed d_L (Wiltshire, via Smale & Wiltshire / arXiv 1107.5596 eq 7): Hbar0 d_L = (1+z)^2 y^2 [2y + (b/6) ln((y+b)^2/(y^2-by+b^2))
     + (b/sqrt3) atan((2y-b)/(sqrt3 b))]_y^{y0},  y^3 = Hbar0 t,  b^3 = 2(1-f_v0)(2+f_v0)/(9 f_v0);
  1+z = 2^(4/3) y (y^3+b^3) / [f_v0^(1/3) y^3 (2y^3+3b^3)^(4/3)]  (from memory; verified = 1 at t0 analytically);
  y0^3 = (2+f_v0)/3;  dressed H0 = (4f^2+f+4)/(2(2+f)) Hbar0.
Self-check required before use: low-z d_L -> c z / H0(dressed) within 1%.
Data (late universe only, so no early-universe calibration is needed by any model):
  (1) supernova shape, 3 sets (absolute magnitude marginalised). Free: f_v0 (timescape), Om (LCDM, LAW).
  (2) DESI DR2 BAO with the ruler r_d*H0 FREE for every model (shape only). D_H = dD_M/dz (approximation for timescape, disclosed).
  CMB not used: timescape needs its own early-universe treatment (disclosed limitation).
Free-parameter count equal for all three models (1 shape parameter; +1 ruler in BAO).
Expectations:
  E1 SN: timescape within |delta chi2| < 5 of LCDM on all three sets (it was built to fit supernovae).
  E2 best f_v0 0.70-0.85.
  E3 BAO shape: timescape WORSE than LCDM by 5-30.
  E4 SN+BAO combined: timescape worse than LCDM and LAW by > 5.
