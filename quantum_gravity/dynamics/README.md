# Running the equations on the grid; leakage; electron sector (iterations 47-49, Oct 2026)
Expectations: PREREG_47_49.md (committed before running).

## 47a  Whole history from the grid's own action  -> PASS (D1, D2)
The computer was given only the discrete grid action (K read from link-length changes, Claim 1 entering only via K) and found
each time step by making the action stationary. Result: H(z) equals the Friedmann + Claim 1 solution to 2e-9; crossing
z = 0.665; w0 = -0.926; age since z = 1100 = 0.9456/H0 (13.7 Gyr at H0 = 67.5). The acceleration equation, never imposed,
holds automatically: residual 1.1e-4 -> 2.8e-5 -> 6.9e-6 as steps halve (~ step^2). Caveat: for a uniform universe the grid
reader is exact, so this tests action -> equations and Claim 1's consistency, not the reader.

## 47b  Gravitational ripples on 64,000 random points  -> G1, G2, G4 PASS; G3 not met at 1%
Speed of plane waves = smoothed-continuum prediction to 0.03%, same in all five directions (spread 0.03%): no preferred
direction. Travelling packet: first two speed measurements were biased by scattered energy (recorded in the script; the
measurement was revised twice); the coherent wave moves at 0.976 x the continuum value, but the follow-up (47c) shows this
number scatters by ~2% between wavelengths (+2.3% at wavelength 7), so it is consistent with no slowing at the ~2% level only.
NEW FINDING: a random grid scatters waves -- 29% of a wavelength-10 wave's coherent energy in 16 time units.

## 47c  Scattering vs wavelength (follow-up, expectations committed first)
Coherent loss per unit distance 0.063 / 0.025 / 0.0094 at wavelengths 7 / 10 / 14: slope 2.7 (expected 4; between the 2 and
4 cases). Real universe (cell <= 5.7e-28 m, LIGO waves ~3000 km, 40 Mpc): loss ~ 1e-40 with slope 2.7; even slope 2 gives
~2e-16. Negligible -- but an extrapolation over ~33 powers of ten from a factor-2 range; slope measured roughly.

## 48  Leakage factor c (quantum_gravity/lambda/iter48_leakage.*)  -> L1, L2, L3 as expected
Tree level c = 0 (VCDM class, BORROWED: GW speed exactly 1, G_eff = G below the horizon, stars as in GR). The frame shifts
the K^2 weight only for horizon-sized ripples: Delta lambda = Omega_DE/6 (0.11 today, 1e-35 at nucleosynthesis) -- needs a
perturbation calculation (open). Claim 1 must read K AVERAGED over > ~0.03 mm (K^(-1/2) is undefined where cell-scale ripples
make K cross zero); with that, loop leakage is suppressed to ~0. COST: the graviton-speed bound no longer pins the cell size
to the gamma-ray-burst limit; window back to 4e-32 .. 5.7e-28 m. A near-testable number is lost.

## 49  Electron sector (quantum_gravity/electron/)  -> E1 pass; E2 FAIL of the fitted switch shape, fixable
No Lorentz violation if the electron mass depends on the tension field as a scalar (requirement: no coupling to the frame's
direction). Time variation: the tanh switch used in every fit (width dz = 30) leaves a tail excluded by methanol at z = 0.89
(Delta mu/mu = (-1.0 +/- 1.3)e-7) for z_s = 100-150. Fix: sharper switch, z_s/dz >~ 5 (dz < 19 for z_s = 100). The CMB fit
is insensitive to the width. Width is a free choice (accommodation).
