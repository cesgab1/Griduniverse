# PREREG 108 -- a hot quasi-static phase held by the grid: seeds, predictions, and whether the grid can hold it (Coalesce)
Committed BEFORE running iter108_grid_static.py. Mechanism borrowed from string-gas cosmology (Nayeri-Brandenberger-Vafa 2006;
Brandenberger et al. 2007/2014) -- formulas from memory, flagged: P_Phi ~ (l_P/l_s)^4 / (1 - T/T_s); tensor tilt n_t = 1 - n_s.
Grid mapping: l_s -> grain size 1.67 l_P; saturation temperature T_s -> grain temperature.
Q1 (no free number): seed amplitude with grain-set constants. Expect >= 0.1 vs measured 2.1e-9 -> ~1e8 too large: FAILS.
Q2 (one free number, the phase's saturation temperature): amplitude fitted by construction. Predictions: n_t = +0.035 (blue), vs
   inflation n_t = -r/8 < 0. With r < 0.036 (BICEP/Keck 2021), compute the gravitational-wave background at LIGO frequencies
   (blue tilt helps): expect still ~1e7 below LIGO's limit -> not testable soon except via CMB B-modes (LiteBIRD / CMB-S4).
Q3 Can the grid hold it? Pressure of that hot phase vs the grid's measured tension today (5.3e-10 Pa, elastic_grid/) and the
   grid's natural Planck stiffness (~Planck pressure). Expect: needs ~1e95-1e96 Pa -> requires Planck-level stiffness then and
   dark-energy-level now: a ~1e105+ relaxation = the unsolved dark-energy SIZE problem in another form.
