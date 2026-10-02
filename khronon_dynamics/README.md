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

## 5. Crossover window with one-fluid bookkeeping (crossover.py, crossover_results.txt)
Full weak-field equation with the mass term: div((1-f) grad Xi) - mu^2 Xi = 4 pi G drho. The fluid self-gravitates above the
crossover lambda = 2 pi/mu (CDM-like) and is in the MOND regime below it. Schedule: comoving lambda = 100-200 kpc at z = 3
(constant sound speed before), growing as a power of a to lambda0 today.
Solver check: an early version differentiated Xi numerically and blew up where shells pile up. It now takes Xi' from the
integrated form, and in the small-lambda limit it reproduces pure Newtonian dust (control: 5.8-7.3e11 vs 6.6e11 within 300 kpc).
A. Linear: sigma8 unchanged (>= 0.999) for every schedule (lambda0 up to 10 Mpc); Lyman-alpha power 0.995-0.996 for
   100 kpc at z = 3 (0.982 for 200 kpc). The window exists in linear theory.
B. One galaxy (Mb = 1e11 in a 1.5e12 Msun region, radial shells), fluid excess vs MOND halo:
   - ~30 kpc (rotation-curve scale): CDM-like control 5.8e11 (2.5x MOND); with a crossover growing to >= 1 Mpc today
     1.2-1.9e11 (MOND 2.3e11). The excess IS expelled dynamically: the saturation works where SPARC measures.
   - 100-300 kpc (lensing scale): 2.7-6.3e11 vs MOND 8.8e11-2.7e12, i.e. 2-8x short at z = 0. With correlated surroundings
     (amplitude too high, so read z = 0.9): 4.6-6.1e11 at 100 kpc and 1.2-1.3e12 at 300 kpc, ~2x short of MOND.
   - Seed only (no surrounding overdensity): 1.1-1.9e11 at 30 kpc, 2.3-2.7e11 at 300 kpc.
   Note: MOND's halo grows without limit (~ r) out to the external-field radius. A CDM-like region saturates, so supply at
   >= 100 kpc depends on the surroundings. LCDM matches KiDS there with the halo plus correlated neighbours; that is the
   supply this picture would also need.
Limits: 1-D, radial orbits, baryons fixed from z = 30, crude environment profile. These are trends, not precision numbers.
Status: PLAUSIBLE/OPEN. Inner galaxy: works. Outer halo: short by ~2x with surroundings, more without. Next: realistic
environment (3-D or a peak-profile initial condition) and the external-field radius, which limits how far MOND's halo must extend.

## 6. Realistic surroundings + external field (peak_env.py, peak_env_results.txt)
- Surroundings: conditional mean density around a 1.27-sigma peak (1.5e12 Msun region, collapsing at z = 1) in the LCDM
  linear field (BBKS spectrum, sigma8 = 0.81). This replaces the crude R^-1.2 profile.
- MOND target with the external field (1-D estimate), Mb = 1e11, phantom mass within 30/100/300/1000 kpc:
  isolated 2.3e11 / 8.8e11 / 2.7e12 / 9.2e12; g_e = 0.025 a0: 2.1e11 / 6.9e11 / 1.4e12 / 1.9e12;
  g_e = 0.05 a0: 2.0e11 / 5.5e11 / 8.7e11 / 9.6e11. The external field caps the halo near ~1-2e12.
- Fluid excess at z = 0, crossover growing from 100 kpc (z = 3) to 1-10 Mpc today: 1.4-2.2e11 / 2.9-3.4e11 / 6.0-7.1e11 /
  1.2-1.3e12. CDM-like control: 4.9e11 / 5.6e11 / 7.5e11 / 1.2e12.
- Ratio to MOND for g_e = 0.025-0.05 a0: 30 kpc 0.7-1.1; 100 kpc 0.4-0.6; 300 kpc 0.4-0.8; 1 Mpc 0.6-1.3.
Reading: with realistic surroundings and the external field, the gap shrinks from 2-8x to ~1.2-2.5x and is worst at ~100 kpc.
At 30 kpc the CDM-like excess (2.3x MOND) is expelled to about the MOND amount. The 1-D radial-orbit model can't settle a gap
of this size (radial shells oscillate; no angular momentum or relaxation). Deciding it needs a 3-D run. Status: PLAUSIBLE, not shown.

## 7. 3-D particle-mesh run with the DBI fluid law (pm3d.py, static_check.*, remelt_check.*, run_*.log)
No hand-picked crossover: the fluid law is the DBI K(Q) with the parameters that passed our CLASS fits (1/mu = 300 Mpc,
lambda_D = 9.2e6). 6 Mpc periodic box, 128^3 particles and grid (47 kpc cells), conditional-mean peak (1.5e12 Msun, z_c = 1) plus a
Gaussian random field, fixed 1e11 Msun baryons at the centre, CDM control with the same initial conditions.
Solver: Newton/Picard with CG for the dust-density variable (cancellation-free DBI algebra). Xi is switched on at z = 7; an
independent check found the Xi force is 0.2-0.8% of gravity at z = 7-9, so this is harmless.
Results (independently reviewed; corrections included):
1. HOLDS. Until z ~ 0.5-0.7 the fluid inside halos is dust-like, so halos form exactly as in CDM. The Khronon run equals the CDM
   control at z = 2, and at z = 1.8 the Xi push is 0-10% of the fluid's gravity. The cause is DBI saturation: the Xi well depth a
   halo can hold is capped at ~c^2 (x_max - x_bg), growing as a^6. That is ~(30 km/s)^2 at z = 1.8, (90)^2 at z = 1, (210)^2 at
   z = 0.5 and (700)^2 at z = 0. So galaxy-size wells (150-200 km/s) can't be in the MOND regime before z ~ 0.5-0.7.
   Consequence: galaxies at z >~ 1 have CDM-like halos, and MOND phenomenology can only develop afterwards.
2. NOT SHOWN: "the CDM-built halo is already in MOND equilibrium". My first reading (net force 1-7% of gravity) was wrong.
   In the MOND regime the Xi push automatically cancels the fluid's own gravity, so the informative number is the residual
   relative to the baryons' pull. For the z = 0 CDM halo that is 55-85% (fluid still pulled in). The code's own equilibrium (no
   external field) needs ~1.5x the CDM excess at 300 kpc and ~2.5x at 800 kpc. That matches the 1-D finding.
   What does hold: the CDM excess (6.3e11 / 1.3e12 / 1.8e12 / 3.1e12 within 100 / 200 / 300 / 1000 kpc) lies between the 1-D MOND
   targets with external field 0.025 a0 and the isolated targets. The periodic box has no external field, so it can't decide
   between the two.
3. NOT SIMULATED: the approach to equilibrium after z ~ 0.7. The solver stalls below z ~ 2 (CG at its limit, line-search steps
   2^-9 to 2^-12), so the run was stopped. The a = 0.5 and 0.77 rows of remelt_check rescale the z = 0 halo and are toys.
Further caveats: the initial conditions add an unconstrained random field, so the halo is ~2x the designed 1.5e12. Results
inside ~200 kpc depend on the smoothing (unresolved). Baryons are fixed from z = 30 (3% of the halo).
Status: the DBI epoch behaviour is confirmed in 3-D. Whether halos relax to MOND after z ~ 0.7, and whether the external
field closes the 1.5-2.5x gap, needs a better solver (full Newton including the MOND-weight derivative, or multigrid), a larger
box or explicit external field, and a zoom for < 100 kpc.
Testable consequence if this picture holds: galaxy rotation curves and lensing at z >~ 1 should look CDM-like (cuspy inner
halos), with MOND-like behaviour appearing only at z <~ 0.5. Massive z ~ 2 disks with falling rotation curves (Genzel+2017)
are a first check.

## 8. The late universe: does the CDM-built halo turn into a MOND halo? (late_runs_results.txt, endstate_check.py)
Solver rewritten as a convex energy minimisation with a barrier at the DBI limit and a per-cell fraction-to-boundary step.
It converges in 1-3 Newton steps (residual ~1e-4) from z ~ 0.3 on. It still stalls around z ~ 0.7-2, where halo cells sit at
the DBI limit. Since the fluid is dust-like there anyway (section 7), the CDM run is used up to z = 0.3 and Khronon from there.
Results at z = 0 (fluid excess around the 1e11 Msun galaxy):
| | 100 kpc | 300 kpc | 1 Mpc |
|---|---|---|---|
| CDM at z = 0.3 (start) | 6.8e11 | 1.80e12 | 2.93e12 |
| CDM at z = 0 | 6.3e11 | 1.82e12 | 3.12e12 |
| Khronon z = 0.3 -> 0, no external field | 1.5e11 | 5.4e11 | 2.92e12 |
| Khronon, external field 0.025 a0 | 1.5e11 | 5.4e11 | 2.92e12 |
| MOND needs (external field 0.025 a0 / isolated) | 6.8e11 / 8.7e11 | 1.4e12 / 2.7e12 | 1.9e12 / 9.2e12 |
Once the fluid is in the MOND regime its self-gravity is cancelled, so the halo, which has the random (virial) motions it
gathered as dust, is held only by the baryons' Newtonian pull. It unbinds: within ~3 Gyr the fluid inside 300 kpc drops to
~30% of its CDM value and ~40% of what MOND needs. The material moves to 300 kpc - 1 Mpc. The force balance at z = 0 confirms
it is far from the MOND equilibrium: net inward force / baryon pull = 0.9-0.98 (0 at equilibrium). The external field makes
almost no difference. Galaxies would end up with neither CDM halos nor MOND halos (rotation curves near baryons-only at
100-300 kpc), so this FAILS.
Caveat that matters: Khronon's fluid is a single-stream, irrotational flow. Random (multi-stream) motions are outside the
theory, which forms caustics there; the particles here are the standard CDM-style stand-in. The real problem is therefore
deeper: no mechanism in Khronon turns a dust-assembled, moving halo into the static MOND configuration (the kinetic energy
has nowhere to go). Any fix needs either dissipation in the fluid or MOND active before halos assemble (which the
cosmology-passing DBI setting forbids, section 7).

## 9. Can dissipation turn the moving halo into the static MOND halo? (dissipation_results.txt; TAU option in pm3d.py)
Test: from z = 0.3, in the MOND regime, the fluid's random motions relax toward the local mean flow on a time TAU, as if the
energy were radiated away (a 'cooling' fluid). External field 0.025 a0. Fluid excess at z = 0 / 1-D MOND target:
| | 30 kpc* | 100 kpc | 200 kpc | 300 kpc | 500 kpc | 1 Mpc |
|---|---|---|---|---|---|---|
| no dissipation | 0.2 | 0.22 | 0.26 | 0.38 | 0.72 | 1.5 |
| TAU = 1 Gyr | 1.1 | 0.61 | 0.59 | 0.72 | 1.18 | 1.6 |
| TAU = 0.3 Gyr | 2.5 | 0.98 | 0.95 | 1.06 | 1.3 | 1.6 |
(*30 kpc is below the grid resolution.) Force balance (net inward / baryon pull; 0 = equilibrium) at 150-800 kpc:
0.87-0.95 without dissipation, 0.31-0.84 (TAU = 1 Gyr), 0.27-0.66 (TAU = 0.3 Gyr). Moving toward the MOND state, not there yet.
Reading:
- With fast dissipation (0.3 Gyr) the fluid at 100-300 kpc settles to the MOND amount (within ~5% of the 1-D external-field target).
  The kinetic-energy problem of section 8 is the right diagnosis, and removing that energy is a cure in principle.
- Costs: (1) the cooling time is a picked number; (2) the inner region over-collects (2.5x at 30 kpc, unresolved), as cooling
  flows do; (3) the dissipation must act only in the MOND regime. Acting during the dust-like era (z >~ 0.7) it would make
  dissipative dark matter (dark disks, cored or collapsed halos), which is constrained. (4) The energy must go somewhere dark
  (e.g. into grid vibrations), about M sigma^2 ~ 1e52-1e53 J per galaxy.
Status: PLAUSIBLE, with a specific requirement. A UV completion of the khronon fluid must dissipate random motions on
<~ 0.3-1 Gyr in the MOND regime but not in the dust regime. Nothing derives that yet.

## 10. Looking for the dissipation channel, and the redshift test (dissipation_channel.*)
Gate: free. The fluid couples to the Xi field with strength ~[k^2/(k^2 + mu_eff^2)]^2, which vanishes in the DBI dust era and is ~1
in the MOND era. Any channel that works through Xi switches on exactly when needed.
Channel: radiation of Xi waves by the halo's moving mass would remove the random motions in 0.3 Gyr only if those waves
travel at ~190 km/s. Gravitationally coupled aether/khronon modes slower than cosmic rays are excluded by the absence of
vacuum Cherenkov losses (Elliott, Moore & Stoica 2005), so the modes must travel at ~c. At c the dissipation time is ~1e21 Gyr.
FAILS. No grid-level dissipation channel found. The requirement (lose random motions in <~1 Gyr, MOND era only) stays open,
and it is now the sharpest single thing a completion of the theory must supply.
Redshift test (data, no new fit): if MOND switches on only at z <~ 0.5-0.7, galaxies at z ~ 1-2 are CDM-like and the baryonic
Tully-Fisher relation should differ from today's. MOND-at-all-times predicts no change (v^4 = G M a0). Observed:
KMOS3D (Ubler+2017) finds a negative BTFR zero-point change from z = 0 to 0.9 and a positive one from 0.9 to 2.3. Sharma+2024
(0.6 < z < 2.5) find slope 3.2 +/- 0.3 and a 'subtle deviation' from local (local slope ~3.85 +/- 0.09). Mild support for
'not MOND at z ~ 1-2', which our epoch result predicts. These measurements use velocities at ~2-3 disk scale lengths with
pressure corrections, so systematics are large; it is a pointer, not a confirmation.
