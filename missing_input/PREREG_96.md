# PREREG 96 -- is there a missing MEASURED input that tunes galaxies and clusters with one equation?
Committed BEFORE running iter96_missing_input.py.
Rule (Coalesce): no new ingredients (no dark matter/energy). Allowed: let the galaxy law's scale a0 depend on a quantity every
system already has measured: a0_eff = a0 * (X / X_ref)^n, one universal exponent n.
Targets: clusters behave as if a0_eff = 5.3x the galaxy value (iteration 94, Coma vs Fornax dwarf; Virgo 2.7-16x) or ~17x
(cluster radial-acceleration relation, Tian+2024 CLASH, g-dagger 2.0e-9; from galaxies_lensing/cluster_review).
Data: 165 SPARC galaxies, per-galaxy best-fit a0 (fixed M/L, method of no_mond/a0_without_mond.py).
Candidate inputs X (5): baryonic mass; well depth v^2; size (last measured radius); gas fraction; baryon surface density.
Cluster values: M_b 3e14 Msun, v 1400 km/s (isothermal sqrt(2) x 1000), R 2000 kpc, gas fraction 0.83, Sigma = M_b/(pi R^2).
For each X: slope n measured inside SPARC (with error) vs slope n_req needed to reach the cluster target.
Expectations:
- Mass, well depth, size: n_req ~ 0.1-0.6; SPARC slope ~0 +/- 0.03 -> each EXCLUDED at > 4 sigma (galaxies span 4 decades in
  mass and keep one a0 to ~15%).
- Gas fraction and surface density: clusters are NOT outside the galaxy range (gas-rich dwarfs ~0.9; Sigma ~ 20 Msun/pc^2 is
  galaxy-typical) -> cannot separate clusters at all (n_req undefined/huge). EXCLUDED.
- Overall expectation: no single measured input works; what remains is a threshold put in by hand (e.g. hot X-ray gas) or
  extra mass. Known-matter options already excluded earlier (neutrinos by DESI; hidden ordinary gas by the cosmic baryon share).
Look-elsewhere: 5 inputs x 2 targets; a pass at ~2 sigma on one would be weak.
