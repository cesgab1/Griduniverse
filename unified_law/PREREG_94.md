# PREREG 94 -- one law for galaxies and clusters (Coalesce's proposal, voice session Oct 5 2026)
Committed BEFORE running iter94_unified.py.

Proposal: v_tot^2 = v_ord^2 + v_rand^2;  v_ord^4 = G M a0;  v_rand^4 = k G M a0;  k = random/ordered kinetic energy = (v_rand/v_ord)^2.
Inputs: NGC 3198 M = 3e10 Msun, v_ord 150, random 10 km/s. Coma: visible mass, random 1000 km/s, ordered ~0. a0 = 1.2e-10.

Checks and expectations:
C1 Self-consistency: if k = (v_rand/v_ord)^2 and v_rand^4 = k v_ord^4 then k^2 = k, so k = 0 or 1 only. Expect: the
   formula as written is not self-consistent, and for Coma (v_ord = 0) k is infinite -> undefined. (Algebra; certain.)
C2 Well-posed version: v_rand^4 = kappa G M a0 with ONE universal kappa (this is MOND's mass-velocity relation for
   random-motion systems; isothermal theory value 4/81 = 0.049). Fit kappa separately from
   (a) Fornax dwarf spheroidal (random-motion galaxy): sigma 11.7 km/s, M_b 2-4e7 Msun;
   (b) Coma: sigma 1000 km/s, M_b = stars 0.5e14 + gas 2.5e14 = 3e14 Msun (galaxiesbook.org; Hughes 1989); also 1.5e14 and the proposal's 2e13.
   Expect kappa_Fornax ~ 0.03-0.06; kappa_Coma ~ 0.2-0.4 -> Coma needs ~4-8x the dwarf's kappa (= the known MOND
   cluster gap, ~x2.5-3.4 in galaxies_lensing/cluster_review.txt by a different method). PREDICTED: FAIL.
C3 NGC 3198 alone: v_ord = (G M a0)^(1/4) ~ 148 km/s. Expected PASS (the rotation leg is MOND, already verified).
Decision rule (Coalesce's): if one constant fits both random-motion systems -> unified; if Coma needs a different constant -> fails honestly.
