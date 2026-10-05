# RESULT 94 -- unified galaxy-cluster law: FAILS on clusters (as pre-registered)
Output: iter94_unified.txt
- C3 PASS: NGC 3198 rotation from (G M a0)^1/4 = 147.9 km/s (observed ~150).
- C1 FAIL (algebra): with k = (v_rand/v_ord)^2, the law v_rand^4 = k v_ord^4 forces k^2 = k. For NGC 3198 it returns 38 km/s of
  random motion from an input of 10 km/s; for Coma (no rotation) k is infinite. The formula as written is not self-consistent.
- C2 FAIL: the well-posed version (one universal kappa: v_rand^4 = kappa G M a0) gives kappa = 0.03-0.06 for the Fornax dwarf
  (matches MOND theory 4/81 = 0.049) but 0.21-0.42 for Coma (3.1 with the proposal's 2e13). Coma needs ~5x the dwarf's constant.
  With the dwarf's constant, Coma's speeds imply 1.6e15 Msun vs 3e14 Msun visible -> the known cluster gap (x5 here; x2.5-3.4 in
  galaxies_lensing/cluster_review.txt with radial profiles). Decision rule: Coma needs a different constant -> fails honestly.
Sources: Coma sigma ~1000 km/s, stars 0.5e14, gas 2.5e14 Msun (galaxiesbook.org ch. 5, Hughes 1989); Fornax sigma 11.7 km/s
(Walker et al. 2009), L_V 2e7 Lsun (McConnachie 2012).
