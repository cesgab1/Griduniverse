# Ideas ledger: every idea, grouped by family, so nothing is re-proposed under a new name
Check here before suggesting anything. A new idea must say which family it belongs to and what is different from the closed tests.

## Family 1. MOND from link mechanics (links stretch, slacken, squish, change length/width, take up slack)
Closed. Uniform loose links (Saturn drift 30-40x too big), connected mosaic (snaps all at once, no square-root law), gap strengths
from geometry (slope 0.34-0.38 vs 0.56), channel geometry for a0 (needs coupling 40-900). `grid_models/`, `fluid/a0_from_channels.py`.
Includes "deeper wells stretch links more, so links respond more strongly" (EMOND-like). Do not re-propose without a new
mechanism for WHY the response changes. The math of EMOND (a0 depends on well depth) is untested, but its grid story is this family.

Spider-web geometry (Coalesce, Oct 2026): Family 1, but new in one respect. A sheet-like (2-D) web spreads the pull over a
circle, not a sphere, so it falls as 1/r -- exactly the deep-MOND law, which no link model so far produced. Requirements:
(a) the switch from 3-D to web-like must happen where the pull drops to a0 (the MOND radius sqrt(GM/a0)); a web of fixed size
gives v^2 ~ M, but galaxies show v^4 ~ M (slope 3.85 +/- 0.09); (b) it must look the same in every direction (round lensing,
polar-ring galaxies), so 'thread density falling as 1/r' rather than one flat sheet. Bonus: a strong outside pull would tighten
the web back to 3-D, which is MOND's external-field effect. Still put in by hand: why the switch is at a0 (same as slack).

Fluid-tube network that thickens with flow (fungus/slime mould; Coalesce, Oct 2026): Family 1, NEW mechanism (adaptive growth,
not static slack). Tubes thickening as sqrt(flux) and saturating give the MOND field equation (AQUAL), the a0 switch radius
and v^4 ~ M automatically; 3-D lattice test passes (slope 1.01, isotropic to 5-10%). Put in: exponent 1/2 and a0. `network/`
Same-fashion geometry comparison (network/README.md): static links give Newton in every geometry. With the flow rule:
mosaic and fungal pass (tie), cubic imprints its axes (~20%, excluded), spider web never Newtonian and centred (excluded).
Earlier mosaic failures were the slack rule's, not the geometry's.
Fluid = contents of the thickened tubes (as the MOND halo): NO (general no-go: tube contents are local in g, the MOND halo
density ~ g/r is not; mosaic test mismatch >= 0.33 in ln). The network models the force law only.
Planks made of fluid with dead space between: the same tubes. (Earlier '10-um channels as light's grid' was a different claim.)

Revisit without MOND (no_mond/README.md section 4): electron twins (doublers) on the mosaic pass, never MOND-dependent (Wilson /
overlap leave exactly one electron); a light twin fermion as ALL dark matter must be >= 200 eV (Pauli, SPARC) and > ~5.7 keV
(Lyman-alpha) -> cold-DM-like, merges with the best fit. Fluid in channels between the layers with its own pressure law fails
without MOND (core-size slope -0.53 vs 0 or +0.5; Jeans length 7 Mpc at recombination); with pressure off it is just CDM.

## Family 2. Relativistic base
AeST: closed (fails Solar System by 1e3-1e4; early instability). Khronon: passes (current base). `lorentz/`, `aest_upgrade/`.

## Family 3. Dark-energy law rho_DE ~ adot^-1/2
Passes (beats Lambda, Delta chi2 -7). Consequences: fate (`future/`), z = 0.65 release. Hubble tension: not fixed by late physics.

## Family 4. The CMB needs a cold dark component; galaxies must not end up with it; clusters need some
The core open problem. Root cause found every time: the fluid is cold and clumpy early (CMB, Lyman-alpha), galaxies form early
and capture it, and lensing shows no halos at z 0.3-0.9.
- Keep it out by pressure or stiffness: mass term, Jeans window, density pressure, DBI phase change, viscoelastic: closed.
- Keep it out by a single switch (rate, speed, stress, mass, depth): closed (galaxies in the middle).
- Add it next to MOND (hybrid, salt rule): closed.
- Move it with currents or drag (gravity currents, frame-locked, partially locked): closed or costly (growth damping).
- Other matter (hot DM, missing baryons): closed.
- Heat it out (shaken ocean, isolated wells, dark fountain): closed (lensing vs redshift; Lyman-alpha vs escape).
- The fluid in a galaxy IS its MOND halo (one thing, not added), the original ocean picture: budget checks PASS
  (supply ample; total MOND halo mass of the universe = 0.3-1.5 x the CMB fluid for g_e = 0.01-0.05 a0). Khronon already has it
  (paper eqs 3.13/3.16: one conserved density = dust + phantom) and its dynamics give a stable saturation at exactly MOND.
  Problem: in the MOND regime the fluid has no self-gravity, so it gathers <~4 Mb, too little beyond ~4 r_M. Needs CDM-like
  self-gravity on >~ Mpc scales (crossover 1/mu). Earlier hybrid/DBI collapse runs double-counted (fluid + MOND). Crossover window (one-fluid): linear OK (Lya 0.995, sigma8
  unchanged); a crossover growing from 0.1 Mpc (z=3) to >= 1 Mpc today expels the excess at ~30 kpc (galaxy ~ MOND amount), but
  at 100-300 kpc the fluid is short of MOND; with realistic peak surroundings + external field (0.025-0.05 a0) the gap is
  ~1.2-2.5x (worst at ~100 kpc). 3-D run with the DBI fluid law (no hand-picked schedule): halos are dust-like (CDM) until
  z ~ 0.5-0.7 (DBI well-depth cap ~ a^6), confirmed; at z = 0 the CDM-built halo is 1.5-2.5x short of the isolated MOND
  equilibrium, between the external-field and isolated targets; late run (z 0.3 -> 0): once MOND switches on, the halo's random motions unbind it (fluid inside 300 kpc -> ~0.4 of
  the MOND need; external field no help). FAILS without dissipation. With dissipation of random motions (cooling time
  0.3 Gyr, picked) the fluid at 100-300 kpc settles to the MOND amount (+-6%); needs dissipation only in the MOND regime
  (not in the dust era), inner region over-collects. Grid channel searched: the gate is free (Xi coupling vanishes
  in the dust era), but khronon-wave radiation needs ~190 km/s waves, excluded by cosmic-ray Cherenkov limits. Open.
  High-z BTFR (KMOS3D, Sharma+24) mildly favours 'no MOND at z~1-2'. `khronon_dynamics/` sections 8-10

## Family 5. Black holes
Planck stars: consistent, no observable effect for >1e26 yr. Ended holes can't make voids or fix H0: closed.

## Family 6. Origin: Big Bang singularity vs grid bounce (bounce/)
Grid density cap -> bounce (LQC-like) instead of a singularity: consistent, not testable alone. Our DE law forbids turnaround,
so at most one bounce from an always-contracting branch.
Bounce INSTEAD of inflation (matter/LambdaCDM bounce, with Lambda or our DE law in the contracting branch): n_s = 0.965 forces a
running of +0.21 to +0.23, ~30 sigma from Planck. FAILS. Bounce + inflation: survives, but needs an inflaton (new ingredient).
Primordial black holes: not predicted either way (needs large early fluctuations; our early universe is the CMB-fitted one).
If Planck stars are real, primordial holes ending today weigh ~1e23 kg (asteroid mass), not ~1e11 kg (black_holes/).

## Chladni picture (Coalesce, Oct 2026): sand on a vibrating plate collects at the still lines (nodes)
Physics of the real effect: grains are kicked where the plate moves and come to rest where it doesn't, because every bounce
loses energy (inelastic). Without that energy loss there is no pattern. So it needs the same missing ingredient as the
dissipation requirement (khronon_dynamics section 9); it doesn't supply it. Fine powder goes the other way (to the moving
regions, carried by air currents).
Cosmic version: a standing vibration of the grid would gather matter on its nodes, i.e. a preferred spacing/lattice in
matter. Data: the CMB and galaxy surveys look like a Gaussian random field with one known preferred scale (BAO, already
explained by early-universe sound waves); old claims of 128 Mpc periodicity were not confirmed. So any such pattern must be
weak, or random-phase across many frequencies (which looks like ordinary initial fluctuations).
Checked (grid_patterns/): oscillation searches in BOSS + Planck limit any regular pattern to < 2-3% of the matter power on
2-22 Mpc/h scales (A_lin < 0.022-0.031, no detection).
Possible use: inside a galaxy, a vibrating grid plus a fluid that loses energy on each 'bounce' would settle fluid where the
grid is quiet. That is a candidate form for the dissipation channel, only if the energy loss itself is supplied.

## Where does the sloshing energy go? / never fall in fast (Coalesce, Oct 2026)
Budget: a Milky-Way halo's fluid sloshing energy ~ M sigma^2 ~ 1e12 Msun x (200 km/s)^2 = 8e52 J; released in 1 Gyr = 2.5e43 erg/s,
in 10 Gyr = 2.5e42 erg/s (comparable to or above the galaxy's whole starlight).
- As NEW free electrons: forbidden alone (charge conservation). As e+e- pairs: 1.5e49 pairs/s vs the Milky Way's 511 keV
  annihilation line ~5e43 /s, so 1e5-1e6 too many. As electron+proton (new hydrogen): 4.5e5 Msun per Gyr per galaxy (needs
  baryon creation; not excluded by this estimate alone, but a new law).
- As HEAT in existing free electrons (hot halo gas): the halo would shine in X-rays at 2.5e42-2.5e43 erg/s vs the observed
  Milky-Way hot-halo ~1e39-1e40 erg/s, i.e. 300-25000x too bright. Excluded.
-> the energy must stay dark (grid / dark sector), or there must be no sloshing to begin with.
Never fall in fast: the cleanest route (nothing to dissipate). Requires the MOND regime (fluid self-gravity off, gentle seepage)
to be on while galaxies assemble (z ~ 1-3). Two conflicts: Lyman-alpha (z 2-5) needs the fluid clumpy at ~1 Mpc, and with
the DBI law shallow structures switch first (wrong order). And gentle seepage alone gathers <~4 x baryons (1-D result), so
the outer halo needs supply from structure that clumped earlier. Checked (khronon_dynamics section 11): IMPOSSIBLE in Khronon-type theories. Both
handles (fluid law K(Q): deep wells look like an earlier universe; MOND function J(Y): low acceleration = calm) make the
diffuse forest calm before galaxies. A third handle (e.g. shell crossing) would still give an accretion shock at the calm
boundary. Only dissipation into a dark sink remains.

## Re-assessment without MOND (no_mond/README.md)
Rule: keep only data-required tests. Best fit without MOND: ocean = cold dark matter + grid DE law + grid picture. Salt and
'current' still fail without MOND; most fluid-steering ideas were moot (existed only to serve MOND). Requirement R1 (rotation
curves from visible matter): halos fit shapes (cored best, BIC), but zero-parameter MOND is 2.6x better than zero-parameter
abundance-matched halos. Must be earned by galaxy-formation physics. Next: R2 Tully-Fisher, R3 one acceleration scale.
