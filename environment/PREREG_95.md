# PREREG 95 -- does living near/inside a cluster change how a galaxy behaves? (Milky Way + Virgo, and SPARC)
Committed BEFORE running iter95_environment.py.

Question (Coalesce): compare a lone galaxy with its nearest cluster. Three parts, all with the galaxy law v^4 = G M_b a0
(a0 = 1.2e-10) and the cluster constant kappa from iteration 94.

A. Virgo, the nearest cluster: kappa_Virgo = sigma^4/(G M_b a0), sigma = 638 km/s (Kashibadze+2020).
   Baryons: low 1.65e13 Msun (gas ~0.138 + stars inside r200 = 974 kpc, M200 = 1.05e14, Simionescu+2017);
   high 1.0e14 Msun (cosmic 0.157 x virial mass 6.3e14 inside 1.7 Mpc). Aperture mismatch disclosed: the big uncertainty.
   Expect: kappa_Virgo above the Fornax-dwarf value (0.03-0.06) by x2.5-16 -> Virgo shows the same cluster gap as Coma.
B. Galaxies inside a cluster environment vs field: SPARC galaxies with flat outer curves; residual
   r = log10( v_flat^4 / (G M_b a0) ). Ursa Major members (the 28 SPARC galaxies at the cluster distance 18.0 Mpc;
   a loose, spiral-rich cluster) vs all others. Expect: mean residuals agree within 2 sigma (no environment effect detectable).
   Caveat: UMa members share one distance, so their residuals are cleaner than field galaxies' (smaller scatter expected).
C. Milky Way (field galaxy, 16.5 Mpc from Virgo): M_b 5-7e10 Msun -> predicted v_flat 176-183 km/s vs outer rotation
   ~200 km/s at ~25 kpc (Eilers+2019, still slowly declining). Expect: within ~15%.
Virgo's pull on the Milky Way (G M/d^2) also computed, as a fraction of a0 -- expected < 0.01 a0 (too weak to matter).
