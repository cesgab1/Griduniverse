# Does Khronon turn infalling fluid into the MOND halo? (infall.*)

## 1. The theory already says "the fluid IS the MOND halo" (verified independently)
Blanchet & Skordis (arXiv:2404.06584) eqs 3.13 and 3.16: one density
  rho_tau = -(1/4 pi G)[div(J_Y grad Xi) + mu^2 Xi]  = MOND phantom part + dust (mu^2) part,
obeys a continuity equation, d rho_tau/dt + div(rho_tau v) = 0, and gravity is Newtonian in (baryons + rho_tau).
So the MOND halo is made of the conserved fluid, not added on top of it. This is the paper's result, not ours.
(Independent reduction: charge along the khronon flow q = -S.n = -8 pi G (rho_K + rho_ph); same sign, same normalisation.)

## 2. The fluid's equation of motion, Dv/Dt = -grad phi + grad Xi, gives a saturation mechanism
- Static equilibrium (v = 0) means Xi = phi, which is exactly MOND, with the fluid excess = MOND phantom mass (Xi = phi + const).
- It is stable. Too much fluid gives a net outward push, too little gives a net inward pull. The repulsion factor f/(1-f) reaches ~4.
- In the MOND regime the fluid's own gravity is cancelled (over-cancelled). It is pulled in only by the baryons' Newtonian pull.
- Consequence for earlier tests: the hybrid model and the DBI collapse runs treated the fluid as extra Newtonian mass on top of the
  baryons' MOND pull. In Khronon, gravity is Newtonian in (baryons + fluid), and MOND comes from how the fluid arranges itself.
  Those runs double-counted, so their verdicts on galaxy halos need redoing with the correct bookkeeping.
  (The static mass-term test, clusters/khronon_mass_term.*, solved the paper's eq 3.21 directly and is not affected.)

## 3. But with no self-gravity, the fluid can't be gathered (mu -> 0 on galaxy scales)
Analytic bound (reviewer): with the baryons' pull only, from z = 30, the enclosed fluid excess is <~ (Omega_c/Omega_m) F M_b,
F = 4.8 -> <~ 4 M_b at z = 0 (3.2 M_b at z = 0.9), independent of radius. The simulation gives 1.5-2 M_b, consistent with this.
MOND needs phantom ~ M_b r/r_M (r_M = sqrt(G M_b/a0) = 10.8 kpc for 1e11, 3.4 kpc for 1e10).
-> MOND can be assembled only inside ~4 r_M (~40 kpc for 1e11, ~14 kpc for 1e10). At 300 kpc it is 7x (1e11) to 20x (1e10)
short, so outer rotation curves and KiDS lensing at 0.1-1 Mpc would fail.
Simulation caveats (reviewer): the inner <~10 kpc is invalid. The weak-field equation has no solution there, shells cross,
radial orbits are used, and only ~10 shells cover 30-300 kpc. The 'old model' numbers in infall.txt are radial-collapse
artifacts. A minor a(t) bug at z ~ 0 does not affect r < 300 kpc. The conclusion rests on the analytic bound, not the shell numbers.

## 4. What this leaves
The fluid must self-gravitate (CDM-like) on scales larger than the ~3 Mpc comoving region a galaxy has to drain, which
cosmology and Lyman-alpha already demand. It must sit in its MOND equilibrium on galaxy scales. In Khronon both are set
by one number, the crossover 1/mu (density-dependent for the DBI form). Open question: is there a crossover that does both, with
the correct one-fluid bookkeeping and the constant in Xi = phi + C fixed by how much fluid actually arrived (not set to 0)?
