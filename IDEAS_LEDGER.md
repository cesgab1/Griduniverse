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
R2 Tully-Fisher without MOND: abundance-matched CDM halos give slope 2.65 vs observed 3.59 and scatter 0.087 > observed total
0.064 dex; passes only if LCDM galaxy-formation physics fixes it (conditional, not derived). `no_mond/btfr_without_mond.*`
R3 one a0 without MOND: CDM halos give an effective a0 of 1.6e-10 (observed 0.98e-10), drifting 2x with galaxy size; partial pass.
Open issue for the best fit: galaxies are tighter/more regular than halos placed by abundance matching (same as LCDM).

## Family 2. Relativistic base
AeST: closed (fails Solar System by 1e3-1e4; early instability). Khronon: passes (current base). `lorentz/`, `aest_upgrade/`.

## Family 3. Dark-energy law rho_DE ~ adot^-1/2
Passes (beats Lambda, Delta chi2 -7). Consequences: fate (`future/`), z = 0.65 release. Signature (predictions/): phantom
crossing exactly at acceleration onset (z ~ 0.7), mimics (w0, wa) = (-0.90, -0.25) [corrected: first pass used beta=1], beats LCDM by dchi2 5.5-7.3; ties w0wa by
AIC, beats it by BIC. 21-cm step <= 0.7 mK (not near-term). Next decisive: binned w(z) crossing (DESI DR3/Euclid). Early electron with SPT-3G D1 added (SPA+DESI+DES):
m_e = 1.0099 +/- 0.0049 (prediction 0.4-1.1% survives, 2 sigma, not decisive). Hubble tension: not fixed by late physics.

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

## Family 7. Fluid in contact with visible matter (energy exchange; contact/)
Purpose: make halos track the baryons (the open R2/R3 issue without MOND). Literature: Famaey-Khoury-Penco 2018 (needs
sigma = C a0/(n_DM v^2)). Needs ~0.3-18 cm^2/g in galaxies; same law at recombination is 6-110x over Planck's DM-baryon bound
-> must switch on after recombination; strong coupling in the IGM at z = 3 and near-strong at z = 100 -> Lyman-alpha
temperature and 21-cm are the tests. NOT the same as Family 4's drag/currents (fluid moved by the grid) or the hot-halo
X-ray budget (that was heat dumped into halo gas).
CLOSED (early_tests.*): 21-cm fine for heavy DM, but Lyman-alpha gas locked to the fluid 25-43x too fast (no pressure smoothing),
CMB 6-110x, and in galaxies the contact drags disc gas to a stop in 0.3-0.55 Gyr. A momentum sink in the grid is worse.

## Gravity as backdraft / inflow (Coalesce, Oct 2026), brainstorm only, not tested
Three readings: (1) as a PICTURE of gravity: exact. GR's 'river model' (Painleve-Gullstrand; Hamilton & Lisle 2008) has space
flowing into a mass at v = sqrt(2GM/r); a horizon is where the inflow reaches c. Same content as GR, no new prediction.
(2) as a literal sink (a medium consumed by mass, like air into a vacuum cleaner): a conserved inflow gives v ~ 1/r^2 and the
wrong force law; Newton needs the flow to be non-conserved, so the medium must be destroyed along the way. Old aether-sink and
Le Sage 'pushing gravity' theories died on drag/heating, aberration (Laplace) and eclipse shielding limits.
(3) as a real flow with its own rest frame: a 'headwind' for moving bodies = preferred-frame effects, bounded by
alpha1 < ~1e-5, alpha2 < ~1e-9 (pulsars, solar spin). Khronon passes because its slicing does not flow at 1PN.
As a way to move the dark fluid: Family 4 (currents/infall), closed.

## Time as a stamp on every layer and cell (Coalesce, Oct 2026), brainstorm, consistent with the base
This is Khronon theory's own variable: a field tau(x) (the 'khronon') whose value is the time stamp; layers are the
surfaces of equal stamp. Two quantities must be kept apart: the STAMP (a global label, same step between layers
everywhere) and the TICK RATE (proper time per stamp step, the lapse N, slower near mass). Gravity = the tick rate varying
across a layer. Consequences: stamps only increase along any path a body can take (no time loops); a global stamp means a
hidden 'true simultaneity' that experiments must not detect: preferred-frame bounds alpha1 < ~1e-5, alpha2 < ~1e-9 and
gravitational-wave speed = c to 1e-15; Khronon passes these. Inside black holes khronon theories have a 'universal horizon'.
No new test beyond those; it would settle the OPEN 'time layer' meaning in favour of moments in time.

## Crunching the grid ahead of a ship (warp; Coalesce, Oct 2026), brainstorm only
= Alcubierre (1994) warp: contract space ahead, expand it behind. Needs negative energy (no known matter supplies it);
cost set by the grid's stiffness c^4/G ~ 1e44 N (rough: shortening a 1 km region by 10% ~ 1e45 J, several Jupiter masses
of energy, and of the wrong sign). Grid-specific twist: with global time stamps (Khronon foliation) faster-than-light
travel would not create time-travel paradoxes, since stamps still only increase. Density cap limits how far cells can be
crunched. Not testable; no claim.

## Light slower when the grid was less stretched (Coalesce, Oct 2026): TESTED, mild support
Grid speed limit c = cell/tick; early grid -> slower c -> larger alpha. SPA+DESI+DES: alpha_rec/alpha_0 = 1.0031 +/- 0.0013
(dchi2 -5.2). Overlaps the early-electron signal (both set atomic energies). Must have frozen by z ~ 4 (quasars, clocks).
`predictions/README.md` item 4.

## Hierarchy of forces; universe as one wave function (Page-Wootters) (Coalesce, Oct 2026), brainstorm
Hierarchy: gravity is weak because the grid is stiff (c^4/G); the real puzzle is why particles are so light vs the Planck
mass. Two-layer electron gives m ~ (1-M)^W: m_e/M_Planck = 4.2e-23 = 2^-74.3 -> ~74 layers at halving per layer
(Randall-Sundrum/domain-wall logic). Relocates the puzzle to a modest integer; not a derivation; untestable at the LHC.
Page-Wootters: the khronon stamp is a ready-made global clock, which removes the 'problem of time' gap; remaining gaps:
measurement/Born rule for one universe, initial state (past hypothesis), no complete quantum gravity, clock ambiguity.

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

## Quantum gravity: where the model's gaps are (Coalesce, Oct 2026), brainstorm + one test
The grid's gravity is CLASSICAL (Khronon metric theory; the tension-link field equation is its Newtonian limit). Gap map:
- HAVE (borrowed): a time variable (the stamp = Khronon foliation, fixes the problem of time); a UV cutoff (cell size
  l_P to 5.7e-28 m, GRB-tested); no singularity (density cap / LQC bounce); a fermion construction (two-layer electron).
- MISSING: (1) quantum rules for the grid itself (what is in superposition: cell shapes? link tensions?); (2) a spin-2
  graviton: scalar tension links give only spin-0; the metric is borrowed, not built from links; (3) Lorentz invariance
  beyond 1PN from a random lattice + stamp (bounds alpha1 < 1e-5, alpha2 < 1e-9, GW speed = c to 1e-15); (4) the
  cell-creation rule (space must add cells); (5) black-hole information; (6) the Bekenstein-Hawking 1/4 coefficient.
- TEST (quantum_gravity/area_law.py): quantised tension links on a random 3-D Mosaic (2744 cells). Entanglement entropy
  of a ball goes as R^2.00 and tracks cut WALL AREA best (R^2 0.991 vs volume 0.976): the area law holds on the grid.
  Coefficient 0.060 per cell-wall area vs the 0.25 per Planck area black holes need; cutoff- and species-dependent
  (Srednicki/Jacobson 'species problem'), so not a prediction. BORROWED result reproduced; gap (6) stays open.
- TEST gap (2), graviton (quantum_gravity/graviton_test.py, regge_gauge.py). Scalar links: helicity 0 only; node
  springs: 0, ±1; the graviton (±2) needs a tensor per cell = cell SHAPE = its link lengths. Derivative-based tensor on
  the random grid: FAILS (gauge residual 0.9, method artefact; cubic lattice exact). Link-length (Regge) action on the
  random grid: gauge symmetry EXACT (375/375 zero modes, residual 2.5e-13). So gravity on our grid must live in link
  LENGTHS (deficit angles), not tensions; E1 is its Newtonian shadow. 4-D graviton: Rocek-Williams 1981 (BORROWED).
  Quantum version = Causal Dynamical Triangulations (random simplices + time slicing = our stamp) / Horava gravity
  (IR limit = khronometric). BORROWED. Possible new question: does CDT/Horava give the Tension-Rate Law?

## Does CDT/Horava quantum gravity give the Tension-Rate Law? (Coalesce, Oct 2026): CALCULATED (quantum_gravity/law_from_grid.md)
- Action: a preferred-slicing minisuperspace term P = -(2/3) C (H a^s)^(-1/2) reproduces the law exactly (sympy, Noether
  identity holds); no extra degree of freedom (VCDM / extended-cuscuton class) so crossing w = -1 is ghost-free. But an
  analytic low-energy expansion in H can never give H^(-1/2): standard CDT/Horava EFT does NOT produce it. Needs a
  collective sqrt(N) effect + a physical link length.
- Counting audit: law needs INTENSIVE (vacuum-like) tension, and links that stretch (s = 1) rather than cells being added
  (s = 0 gives rho ∝ H^(-1/2)). Data fit of s (DESI DR2 + Planck priors + SN): s = 0.91-1.03 ± 0.15-0.20; s = 0 excluded
  at 5.0-6.7 sigma. [WITHDRAWN, see follow-up (7): the stretching length is comoving, not the cell.]
- Size of dark energy NOT explained (kick per crossing 6e-154 Planck densities). Law must be cut off near the Bounce.
- Independent check found 2 bugs (exact-q only valid at s=1; r_s integrated from grid point above z*): fixed, Claim 1
  moves <= 0.1 (-5.4 / -7.2 / -6.8). Open: VCDM perturbations instead of PPF in the CAMB fits (expected small).
- FOLLOW-UP (same day), quantum_gravity/law_from_grid.md sections 5-7:
  (5) joint beta-s fit: s robust at 0.72-1.06 for any beta; beta weakly fixed (Pantheon+/DES allow 1/2 at 68%/95%,
  Union3 prefers 1-3, Claim 1 point +6.6 there). 2-parameter family = w0wa fit quality.
  (6) PPF vs no dark-energy fluctuations (VCDM proxy): CMB cosmic-variance chi2 0.008 (TT), 0.03 (EE), sigma8 0.09%,
  lensing <= 0.28%: treatment irrelevant.
  (7) CORRECTION: Planck cells cannot stretch (G/c constancy: s_cell < 1e-3), so cells ARE added today. The counting length
  in Claim 1 must be a COMOVING length distinct from the cell (candidates: matter/Ocean separation, parent cells of the
  mosaic, comoving correlation length). Earlier line "cells are not being added today" is withdrawn.

## Grid toy, iterations 1-2 (Coalesce, Oct 2026): quantum_gravity/toy_grid/LESSONS.md
- It. 1 (memory): the grid tension may lag the stretch rate by its memory time 1/(kappa H). Data (DESI DR2 + Planck +
  3 SN sets) need memory <~ 1/2 Hubble time; 1 Hubble time barely beats Lambda, 2+ is worse. Signature: w = -1 crossing
  comes AFTER acceleration onset by an amount set by the memory (memory 1/4 -> z 0.49 vs q=0 at 0.68). Monte Carlo 5000
  domains vs equation agree 1-2%.
- It. 2 (quantum kicks): quantum VACUUM fluctuations give no random walk (log only -> Lambda-like). REAL thermal quanta of
  a CONFORMAL field (CMB photons) give Var = Theta_c/(4 pi kappa adot) ∝ 1/adot: Claim 1's counting derived, no comoving
  length needed (resolves law_from_grid section 7). Not valid for gravitons. Independent check PASS with caveats.
  Next: coupling of grid tension to photons vs existing limits; what sets the memory.
- It. 3 (who jostles): CMB photons CANNOT (gauge invariance: ∫E dt is a boundary term, variance saturates; energy-density
  coupling gives s = 5, excluded). Needed: a massless conformally coupled SCALAR relic ('grid radiation', the links' own
  quanta). Prediction: Delta N_eff ≈ 0.027-0.057 (decoupled above QCD; T today ~0.9 K); below-QCD decoupling (0.30) is
  excluded by N_eff = 2.990 ± 0.070 (dN < 0.107). Testable by SO / CMB-S4-class. New requirement: negligible coupling
  to ordinary matter (fifth force).
- It. 4: memory DERIVED: tension = velocity of a minimally coupled grid field -> Hubble friction 3H (kappa = 3; conformal
  would give ~0.9, disfavoured). Zero-parameter fits (kappa = 3): quadratic energy -7.1 / -9.2 / -8.5 vs Lambda
  (Claim 1 -5.4/-7.2/-6.8; w0wa -7.1/-10.3/-13.2); crossing z 0.46, w0 -0.90. Look-elsewhere: 4 discrete variants tried.
  Relic must be SELF-thermalised (decay-made relic gives no walk). Iteration-3 Delta N_eff 0.027-0.057 WEAKENED to
  0 < dN < 0.107 (inflation + reheating set its temperature; Higgs-portal contact impossible with m_s < 1e-32 eV).
  No tree-level fifth force with s -> -s symmetry. Naturalness of m_s open (as for quintessence).

## Where next after toy iterations 1-4 (Coalesce, Oct 2026), brainstorm
Open issues in plain terms:
(1) HEIGHT: we explain the shape of the dark-energy curve, not its height. The coupling would need to be ~1e-150.
(2) TWO INGREDIENTS: the grid field must be 'ordinary' (minimal) and the jostler 'light-like' (conformal), chosen to work.
(3) FEATHERWEIGHT: the jostler must be lighter than 1e-32 eV, with nothing protecting it.
(4) DATA: Claim 1 ~2.5 sigma, Claim 2 ~2 sigma; galaxy regularity (R2/R3) still open.
Ideas:
 a. ONE FIELD (Occam): can the grid field's own thermal quanta jostle its own long-wavelength motion via self-interaction?
    If yes, (2) disappears. Testable in the toy (next iteration).
 b. Height from counting cells (an amplitude ∝ 1/sqrt(number of cells in the horizon), Sorkin-like) - CAUTION numerology
    risk; only pursue if it comes out of the toy's own equations, not by matching numbers.
 c. Featherweight from a shift symmetry of the grid field (a field that only appears through differences is massless).
    But conformal coupling breaks shift symmetry: tension with (a). Check.
 d. Decisive data: DESI DR3 binned w(z): the toy's crossing at z ~ 0.46 vs Claim 1's at 0.68 vs Lambda (none).
    Simons Observatory: N_eff > 0 and the early-electron shift (Claim 2).
- It. 5 (Coalesce: 'do the vibrations tell the height?'): height NO (coupling and relic temperature free). But the
  vibrations make dark energy LUMPY with a coupling-free ratio: correlation length ~ Hubble/3, delta rho/rho ~ 1.4 ->
  large-angle CMB imprint ~0.16 vs 1e-5 observed (~2e4 too strong, order-of-magnitude estimate). The patch-by-patch toy
  fails as built; Claim 1 itself (global law) unaffected. Escapes: (a) ~1e9 independent components, (b) global
  response only (VCDM-like), (c) noise not dominated by long waves.
- It. 6 (Coalesce: 'each layer has its own jostle; stacked they give a harmonious wave'; (a) segues to (b)): CONFIRMED as
  the only surviving microscopic picture. Lumpiness of N stacked layers sqrt(2(1/N + c^2)); CMB needs N >~ 3e10
  independent layers with layer-to-layer jostle correlation c <~ 6e-6; 74 layers far too few. 'Harmonious' = averaging,
  not synchronising (locked layers bring lumps back). Each layer needs its own relic jostlers (allowed if T_s <~ 4 mK
  at N = 1e10). Escape (c) short-wave jostling CLOSED: the 1/adot law needs horizon-scale noise. Open: what are the
  ~3e10 layers?

## What are the ~3e10 independent layers? (Coalesce, Oct 2026), brainstorm
Requirements (iteration 6): N >~ 3e10 nearly independent copies of the tension field (jostle correlation < 6e-6), each
with its own cold relic jostlers; gravity and matter see their SUM.
Candidates:
 1. SPECIES (favoured): N internal copies of the tension field per cell ('layers' = species, not places). Dvali's species
    bound says gravity's true cutoff is M_P/sqrt(N), so the effective cell size is sqrt(N) l_P. With the GRB cell bound
    (5.7e-28 m) this gives N <~ 1.2e15. Window: 3e10 <~ N <~ 1e15, i.e. cell 2.6e-30 to 5.7e-28 m, gravity cutoff
    7.6e13 to 3.5e11 GeV. Bonus (BORROWED, Dvali): many species is one explanation of why gravity is weak. Prior art:
    Dvali species bound; N-naturalness (Arkani-Hamed et al. 2016: N copies, each colder, Delta N_eff signature).
 2. HIDDEN-DIRECTION STACK: N nearly decoupled sheets along a 4th spatial direction (deconstruction / many-brane). Needs
    inter-sheet coupling tiny (fits c < 6e-6) so their modes stay nearly massless; KK modes of a smooth extra
    dimension do NOT work (massive). Equivalent to 1 seen from inside.
 3. NOT: time layers (already the memory), spatial patches (correlated over horizon/3), the 74 hierarchy layers (too few
    unless each holds ~4e8 species).
Tests: lower edge from large-angle CMB (residual dark-energy lumps ∝ 1/sqrt(N): a predicted ISW-like contribution just
below current level if N is near 3e10); upper edge from Lorentz-violation timing (LHAASO/CTA: the cell bound tightening
toward 2.6e-30 m would squeeze the window shut); Delta N_eff > 0 (small).
- It. 7 (Coalesce: gaps between layers limit interference to neighbours): RESCUES the stack. Correlation that falls off
  across the stack only costs a factor (1 + 2 sum c_k^2): leakage 0.9 per gap -> 2.5e11 layers needed (still inside the
  species window <~1e15). The real requirement is NO jostle common to all layers (1% common part -> lumps 1.4e-2,
  excluded). Picture: cold Ocean in the gaps does not carry fast jostles.
- It. 8: proper large-angle CMB (late ISW, l = 2-30) bound on the layer stack: N_eff >= 1.7e9 - 4.8e9 (2 sigma; two
  super-horizon treatments). Window 2e9-5e9 <~ N <~ 1.2e15 OPEN. Lump signature concentrated at l = 2-3 (steep fall);
  the real quadrupole is low, so no hint. Not yet included: matter response, lensing / galaxy clustering (may tighten).
- It. 9 + CORRECTION: iteration 8 had a time-weighting bug (sqrt(dt) for dt) overstating the bound ~100x; found by
  independent code, machinery validated on LCDM late ISW (D_2 = 74 muK^2). Corrected CMB bound N_eff >= 1.1e7 - 3.2e7
  (incl. matter response, +10%). Window ~1e7 <~ N <~ 1.2e15. Galaxy clustering: induced power <= 7e-3 of LCDM at the
  largest scales: no tightening.

## Coalesce's answers (Oct 2026): (1) each layer vibrates on its own, vibration set at the Big Bang; (2) a layer is pinned
## (pulled from many directions) but free in some
(1) Self-vibrating layers, started hot: matches iteration 4's requirement (made at the hot start after inflation, then
    self-thermalised). One ingredient instead of two: the tension (stretch) of a layer is jostled by the same layer's own
    vibrations. Scaling check (Langevin + fluctuation-dissipation): if the layer's own vibrations damp the tension at a
    rate Gamma ∝ T^n, the tension variance ∝ Gamma T/H ∝ T^(1+n)/H -> link-stretch exponent s = 1 + n.
    Data (stretch_scan): s = 1 fits; s = 2 (n = 1: a scale-free self-interaction, Gamma ∝ T) costs +57/+53/+27 vs Lambda,
    EXCLUDED; s = 0 (n = -1) +37/+30/+23, EXCLUDED. So the coupling between a layer's stretch and its own vibrations must
    be 'Ohmic' (damping independent of temperature), like the linear coupling of iteration 2, and weak (Gamma << H;
    otherwise the tension simply reaches equipartition and loses the 1/H). A Kirchhoff-type coupling (tension ∝ wobble
    gradient squared, as in a drum skin) gives s = 5: EXCLUDED.
(2) Pinned but free in some directions = massless wobbles along the free directions (protected by that freedom), massive
    along pinned ones. Insight: the stretch must couple LINEARLY to a free-direction wobble, which slightly breaks the
    freedom and gives the wobble a small mass. That breaking is controlled by the same tiny coupling that sets the
    (small) height of dark energy: as the coupling -> 0 the freedom is restored ('technically natural'). Next calculation:
    is the coupling required for the observed height small enough to keep the wobble lighter than ~1e-32 eV?
- It. 10 (naturalness of Coalesce's pinned-but-free wobble): FAILS in the minimal model. Linear link L = g Phi h (the only
  direct push) + canonical tension energy: height needs g >= 3.5e-25 eV^2 (best case N = 1.2e15) but stability/lightness
  need g <= (3H0)^2 = 1.9e-65 eV^2: 1e40x conflict. Pinned wobble, derivative coupling, potential energy, stiffness
  rescaling all fail too. The jostling picture explains the SHAPE (1/adot, memory, smoothing) not the SIZE. Claim 1
  unaffected. Open: nonlinear non-mixing link; noise via gravity / time stamp.

## Coalesce (Oct 2026): drop the 'non-simple link' (agreed: it would need its own explanation). New idea: the jostling
## comes from FREE ELECTRONS repelling each other near the layers, generating waves/vibrations. Checked (no new code):
- Scaling: free electrons thin out with expansion, n_e ∝ a^-3. Shot/collision kick rates ∝ n_e v (IGM temperature roughly
  constant after reionization; Coulomb rate ∝ n_e T^-3/2; Thomson ∝ n_e): all give link-stretch exponent s ≈ 3. Data:
  s = 0.72-1.06; s = 2.5 already costs +77 to +139 vs Lambda (stretch_scan). EXCLUDED.
- Electrons repel through electric fields = photons; iteration 3: EM fields cannot random-walk the tension (gauge: ∫E dt
  is a boundary term). Only density-type couplings remain, which give the a^-3 scaling above.
- Ionisation history: free-electron fraction drops ~5000x at recombination (z 1100) and returns at reionisation (z ~7):
  the jostling would switch off and on, with a jump in dark energy at z ~7.
- A light field coupled to electrons -> fifth force on electrons; equivalence-principle tests are extremely tight.
- Positive link kept: Claim 2 already ties electrons to the grid tension, as a ONE-TIME energy transfer (z ~100-200),
  not ongoing jostling.

## Coalesce (Oct 2026, 'way out there'): an atom type that constantly collects and releases free electrons. Checked:
- GOOD part: catch-and-release events are discrete and independent = shot noise, exactly the sqrt(N) counting Claim 1
  needs.
- Problem: the rate. Events per volume = (number of atoms per volume) x (cycles per atom). Atoms thin out as a^-3
  (s ≈ 3, excluded). Real case: intergalactic hydrogen is constantly photo-ionised and recombines; cycling per volume
  = recombination rate ∝ n_e n_p ∝ a^-6 (s ≈ 6, excluded). Dark atoms (atomic dark matter, Kaplan et al. 2010;
  Cyr-Racine & Sigurdson 2013) thin out the same way. Plus: electrons push through photons (iteration 3) and matter
  coupling -> fifth force.
- Insight: the data need events per volume ∝ 1/a. Only two kinds of source thin out that slowly: (i) the temperature of
  a relic glow (iteration 2), (ii) something living on SHEETS that are not added as space grows (area per volume ∝ 1/a;
  a frozen domain-wall network behaves this way). If the catch-and-release happens ON comoving sheets, the rate is
  right. Seed: are those sheets the grid's own walls? (Constraint to check: frozen walls must be light enough not to
  dominate, the Zel'dovich bound.)

## Coalesce (Oct 2026): stretching cools, like a forged bar drawn out from end to end: more surface, so the void cools it
## faster. Applied to the relic glow (#1), combined with the sheets (#2):
- Unification: if the glow is each layer's OWN vibration (#1), and the layer is a sheet that stretches without
  multiplying (#2), stretching cools the vibrations as T ∝ 1/a (wavelengths stretched) -> exactly the rate the data need
  (s = 1). The 'cooling glow' and 'non-multiplying sheets' are the same object: vibrations of stretching sheets.
- Forged-bar test (EXTRA cooling by leaking heat into the void/gaps): T ∝ a^-(1+eps) gives stretch exponent s = 1 + eps.
  Data (joint beta-s scan): s = 0.72-1.06 (68%), best 0.82-0.86 -> eps <~ +0.06 (68%), <~ +0.2 (95%); leak rate
  <~ 0.2 H. If anything the best fit has the glow cooling slightly SLOWER than 1/a (heat flowing in, not out), but
  s = 1 (no leak) is fully consistent.
- Where leaked heat would go: the gaps (the Ocean). The glow's energy is tiny (Delta N_eff < 0.1, ~1e-5 of matter), so
  warming of dark matter is negligible: no conflict with cold dark matter.
- It. 11 (Coalesce: jostling via gravity, rubber sheet pulled in all directions): no free coupling (G known), so
  predictive: Phi per layer 1e-53 to 1e-59 -> predicted dark energy 1e-99 to 1e-102 of observed. Gravity route closed
  for the SIZE. Conclusion after it. 10-11: jostling explains shape/memory/smoothness, not size. Structure: large
  baseline tension (size unexplained, like Lambda) MODULATED by jostling.

## Coalesce (Oct 2026): ponder the bounce ('our universe the result of a bounce back') in light of toy iterations 1-11
Known from before: the density cap gives a bounce (LQC, borrowed); a bounce ALONE fails (~30 sigma, needs inflation);
bounce + inflation is allowed.
1. The bounce as the 'tuning' of the baseline tension (the SIZE problem). Vacuum-like tension is NOT diluted by inflation
   (relics are), so a baseline set at or before the bounce survives to today. Prior art addressing size this way:
   Steinhardt & Turok 2006 (Science 312, 1180): in a cyclic universe the cosmological constant relaxes a little each
   cycle, so most cycles have a tiny positive value; also 'relaxing the cosmological constant' (JHEP12(2016)022).
   Tension with our law: with Claim 1 the far future is a ∝ t^5 (no recollapse), so OUR expansion never turns around;
   the bounce would be a one-time past event (or our cycle the last), not an endless cycle. A ratchet over cycles would
   need the law to change in late times. OPEN.
2. The law at the bounce: rho ∝ |adot|^-1/2 spikes where adot = 0. The bounce is a natural RESET moment for the tension
   (the counting breaks down there anyway) - a place where the baseline could be set.
3. The layers' vibrations (Coalesce #1, 'started at the Big Bang'): in a bounce picture, the hot start IS the bounce
   (compression heats); the glow is then rebuilt at reheating after inflation (iteration 4).
4. A testable overlap: LQC bounce + inflation predicts SUPPRESSED power at the largest angles (Ashtekar, Agullo et al.;
   linked to the low CMB quadrupole and power asymmetry). Our layer lumps ADD power at l = 2-3. Both act on the same
   multipoles: a combined fit could constrain bounce and layer count together. Caution: two effects tuned to cancel
   would be a fudge; only a joint, parameter-counted fit is meaningful.
- Its. 12-15 (bounce/baseline/future/self-tuning, cross-applied): real Planck low-l likelihood: bounce suppression
  only dchi2 -1.4 (not evidence); layer bound N_eff >= 2.9e7 (matches it. 9). Jostled fraction f: data want f >= 1,
  mildly f > 1 (negative baseline) -> 'large positive baseline' (it. 11) disfavoured. Negative baseline never
  recollapses: jostled energy ∝ 1/H self-tunes it away -> coasting future, bounce must be one-time (no cycle). Self-
  tuning + stability bound: H would be >= 1e26x too large -> size problem moved, not solved.

## DROPPED for now (Coalesce, Oct 2026): the bounce. Reason: no improvement (real Planck low-l: dchi2 -1.4, AIC +0.6;
## one-time only in the toy). Can be revived if a later result needs it.

## Penrose: gravitational entropy and the Weyl curvature hypothesis (WCH) (Coalesce, Oct 2026), discussion seeds
WCH: the universe began with zero Weyl (tidal) curvature = perfectly smooth = minimal GRAVITATIONAL entropy; clumping
under gravity raises gravitational entropy (black holes maximal) and gives the arrow of time. Seeds for the grid:
 1. Grid reading of Weyl: Ricci = how much cells are squeezed overall (matter); Weyl = how cells are DISTORTED in shape
    (stretched one way, squeezed another) = deficit-angle pattern beyond uniform squeezing (Regge picture,
    quantum_gravity/regge_gauge.py). WCH = 'the mosaic started with no shape distortion'. Gravitational entropy ~
    variety of cell-shape distortions.
 2. Arrow of time from cell creation: cells are ADDED as space grows (constants test, law_from_grid sect. 7); the number of
    grid configurations grows -> entropy capacity grows in the direction cells are added. 'Time runs forward because the
    grid grows.' The time stamp (Khronon) gives a preferred slicing but no direction; WCH + cell creation supply it.
 3. Our jostling is DISSIPATIVE only in expansion: Hubble friction (-3H) damps the tension's motion -> irreversible ->
    arrow aligned with expansion. In contraction the same term AMPLIFIES (anti-friction): the dark-energy mechanism is
    intrinsically time-asymmetric. (Consistent with dropping the bounce.)
 4. Testable coincidence to check (not assume): gravitational clumping (Weyl growth) slows when acceleration starts.
    Claim 1 crosses w = -1 at q = 0 (z 0.68), the toy at z 0.46. Compute the history of gravitational-entropy production
    (linear growth + halo formation, e.g. Clifton-Ellis-Tavakol 2013 gravitational-entropy measure) and see where it
    peaks/turns. If dark energy's turnover tracks gravitational-entropy production, that is a physical link (caution:
    must be computed, not matched).
 5. Penrose's conformal cyclic cosmology (CCC) is cyclic -> set aside with the bounce.
- Its. 16-19 (Penrose/WCH, quantum_gravity/arrow_of_time/README.md): (16) gravitational-entropy production peaks scatter
  z ~ 0-1.4; horizon-entropy match is circular -> NO link to dark energy's turnover (idea dropped). (17) random mosaic
  irregularity averages away (not Weyl) -> compatible with smooth start; rebuilt cells erase shape distortion -> Weyl and
  gravitational entropy live in LINK LENGTHS (consistent with Regge gauge result). (18) generalised second law holds for
  Claim 1 and the toy; unlike LCDM (finite de Sitter ceiling) our horizon entropy grows forever. (19) in contraction the
  tension grows as a^-6 (stiff, same as shear): mechanism intrinsically time-asymmetric; Penrose's Bang/Crunch asymmetry
  built in; supports dropping the bounce.

## Coalesce (Oct 2026): a smooth start implies the grid EXISTED (smooth) before the Big Bang; the 'Big Bang' is when
## distortion of the grid started.
- Agrees with a real open problem: the smooth (low-entropy) start must be explained or assumed (the 'past hypothesis';
  Penrose's objection that inflation itself needs special initial conditions).
- Combined with iteration 19: the smooth state CANNOT have come from a collapse (contraction makes the tension and the
  tidal chaos grow as a^-6). So if the grid existed before, it was QUIET/STATIC, not contracting. Consistent with
  dropping the bounce.
- Prior art matching this picture: (1) 'Emergent universe' (Ellis & Maartens 2004): the universe sits in a static state
  in the infinite past, then starts expanding (and inflates); no singularity, no bounce. (2) 'Quantum graphity' /
  geometrogenesis (Konopka, Markopoulou & Smolin 2006): space is a network that is disordered when hot and condenses
  into a smooth, regular geometric phase as it cools; the 'Big Bang' is that phase transition. Both are BORROWED
  frameworks close to 'a smooth grid existed, then began to stretch and distort'.
- In our model: the time stamp can exist in the quiet phase (labels with nothing changing); the arrow of time starts
  when stretching and distortion start (iterations 17-19).
- Wording: 'Big Bang' (Hoyle's joke name) is not an explosion from a point but expansion everywhere; 'the moment the
  grid started stretching' is a fair restatement.
- Tests to consider: emergent/static-phase models leave a largest-scale imprint (power suppression, like a bounce;
  iteration 12 found no preference, dchi2 -1.4); the static phase must be STABLE (the classical Einstein static universe
  is not; some quantum-gravity versions are). Candidate calculation: is a static stretched grid stable in our toy (does
  the jostled tension destabilise it, given it grows when H -> 0, iteration 14)?
- It. 20 (Coalesce's quiet pre-existing grid): (A) CORRECTION: the toy's tension motion energy is a stiff fluid fed by the
  glow (decelerating, energy-limited to < 1.3e-6 of critical) -> the dark-energy fits are valid only in the VCDM reading
  (jostling sets the SHAPE of a vacuum-like tension, energy via the gravity sector, SIZE free) - explains its. 10/11/13 in
  one line. (B) a quiet grid cannot stay quiet: no expansion = no friction = unbounded tension build-up -> self-starting
  Big Bang; direction not chosen by the toy; expanding branch smooth, contracting branch chaotic (it. 19). (C) large-angle
  imprint allowed, not favoured (it. 12).

## Self-audit (Coalesce, Oct 2026: 'if all iterations answer correctly, isn't that an issue?')
Honest tally of the toy/quantum-grid work (iterations 1-20 + side checks):
FAILED (idea dropped or picture broke): photons as jostlers (it. 3); single-layer lumps 2e4x too strong (it. 5);
  short-wave jostling (it. 6); height/featherweight naturalness 1e40x (it. 10); gravity-mediated jostling 1e100x
  (it. 11); free electrons (s ~ 3); catch-and-release atoms (s ~ 3-6); self-tuning vs stability (it. 15); entropy-timing
  link (it. 16); bounce not favoured (it. 12) and no cycle (it. 14); literal (stiff) reading of the toy (it. 20).
RESCUED by adding or re-reading something (ACCOMMODATIONS, not successes): lumps -> many layers (N >= 3e7) + gaps;
  'jostlers' -> a new dark scalar; size -> declared a free normalisation; toy energy -> VCDM reading (it. 20);
  conformal vs minimal coupling chosen to work.
GENUINE RISKY PASSES (could have failed, did not): Claim 1 beats Lambda with zero parameters (3 SN sets); data measure
  s ≈ 1 (could have been 0 or 3); memory kappa = 3 DERIVED from Hubble friction and inside the data window; layer bound
  reproduced by two independent methods (3.2e7 vs 2.9e7); second law (weak test: nearly guaranteed).
Lesson: the microscopic toy has become flexible (ingredients added after failures = epicycle risk). Its explanatory
  value now rests on a few risky predictions only:
  (P1) w = -1 crossing redshift: 0.68 (Claim 1) or 0.46 (memory kappa = 3), vs never (Lambda), DESI DR3 / Euclid;
  (P2) Claim 2 early-electron shift (Simons Observatory, SPT-3G full);
  (P3) Delta N_eff > 0 (weak); (P4) coasting future (mild; f > 1).
Rules adopted from now: (a) every new ingredient must come with a NEW prediction or be marked 'accommodation';
  (b) predictions written down BEFORE new data (pre-registration file); (c) keep counting free choices.

## PRE-REGISTRATION frozen 2026-10-03 18:19 PDT, commit cb5154f: predictions/preregistration/ (PREREGISTRATION.md,
## predictions.json, prediction_wz.png). P1 crossing z 0.68 (Claim 1) / 0.46 (memory kappa = 3) with binned w and kill
## criteria; P2 m_e(rec)/m_e0 = 1.004-1.011; P3 0 < dN_eff < 0.107 (weak); P4 consequences. No re-tuning after new data.

## Back to the stuck point: the SIZE of dark energy - intuitive candidate (Oct 2026), brainstorm + magnitude check
'Votes that nearly cancel': every grid cell in our causal past contributes a random +/- to the grid tension; the net is
the leftover imbalance, ~ 1/sqrt(number of cells). Planck 4-cells in one Hubble 4-volume: 5.2e243 -> 1/sqrt = 1.4e-122
(Planck units); observed rho_DE = 1.1e-123: within a factor ~12 (votes over one Hubble time) or ~110 (votes within the
toy's memory, 1/3 Hubble time) - versus the usual 120 orders of magnitude. This is Sorkin's 'everpresent Lambda'
(causal sets), which PREDICTED Lambda ~ 1e-120 before its 1998 discovery: BORROWED, but the grid gives it a home
(cells exist and are added; Coalesce's 'averaging, not marching in step', it. 6-7).
Also explains 'why now' (dark energy ~ critical density at every epoch) and allows a NEGATIVE value (random sign;
cf. it. 13's mild preference for f > 1).
Problems to test before believing: (1) the simplest everpresent-Lambda history fluctuates and tracks H^2 and was in
tension with data (earlier test; literature); (2) its time dependence (variance ∝ H^2) differs from our law's shape
(∝ (aH)^-1/2 or 1/(aH)); a combined model must fit DESI+SN+CMB; (3) the factor 12-110 must come out, not be tuned.
- Terminology: 'votes' -> RANDOM NUDGES (fluctuations); no choice or intelligence implied (Coalesce, Oct 2026).
- Coalesce: the binary nudge could be SPIN (two-valued, known physics, no new bias). Check: electron spins are too few:
  ~1e80 electrons -> leftover 1/sqrt = 1e-40 (photons ~1e89 -> 1e-45), not the ~1e-122 needed; electrons also dilute as
  a^-3 (s ~ 3, excluded) and spin-gravity couplings are tightly tested. Refined version that keeps the idea: give every
  grid LINK/CELL its own two-valued spin-like state. That is exactly loop quantum gravity's spin networks (links carry
  spin labels; areas come in units set by them): BORROWED. Count = Planck cells (~1e244) -> the right size. Bonus to
  explore: the electron is a pattern on the grid (two-layer electron), so its spin may be the grid's own spin
  showing through. OPEN.
- It. 21 (spin-flip size + memory, quantum_gravity/size/README.md): SIZE right within ~100x (N = 6e241 remembered flips ->
  75x critical vs 0.69x). HISTORY FAILS: a size recomputed from the remembered horizon tracks the critical density
  (Omega_DE = const x random); CMB (<~0.02 at recombination) + BBN (<~0.1) vs 0.69 today: 0 of 20000 histories pass.
  Lesson: size and shape are linked; the size must be accumulated/frozen, not recomputed each moment. Candidates: frozen
  at launch of the quiet grid (needs ~1e244 cells then; loses 'why now'), or a ratchet of flips that are never undone.

## Coalesce (Oct 2026): stretchable things already have built-in tension (rubber wants to go back when pulled)
- Matches established physics: dark energy IS a tension. In gravity, tension = NEGATIVE pressure, and negative pressure
  makes expansion SPEED UP (pressure gravitates; for w = -1 the tension equals the energy density). The rubber
  'wanting to go back' locally and the universe accelerating are the same property seen two ways.
- Built-in = a MATERIAL property, set when the material formed (rubber's chemistry fixes its elasticity). This is the
  'frozen' option that iteration 21 said is needed: the size should not be recomputed from the horizon each moment.
- Grid version of rubber: as space stretches, the grid adds cells (new material) to relieve the strain. Dark energy =
  the RESIDUAL strain left because relief lags behind stretching (stick-slip: a cell is added when local strain reaches a
  threshold eps_c). Two cases to compute:
   (a) relief tied to the Planck tick: residual strain ~ H t_P -> rho ~ rho_P (H t_P)^2 ~ H^2 M_P^2: right size but
       TRACKS the critical density -> same failure as it. 21 (generic for any size tied to H).
   (b) relief at a fixed material yield threshold eps_c: residual strain cycles 0..eps_c -> rho ~ stiffness eps_c^2/3,
       CONSTANT in time (frozen, like rubber's yield point), shape from the jostling. Size needs eps_c ~ 3e-62: the open
       question becomes 'why is the grid's yield threshold so small' - must be derived, not set.
- It. 22 (Coalesce's built-in tension, stick-slip): prompt relief -> leftover strain CONSTANT (frozen size, behaves as
  Lambda); short relief delay -> w > -1 always (no crossing); long delay -> runaway. Size needs eps_c = 8e-62 = one cell of
  stretch per row of ~1.2e61 cells = 1.4x today's Hubble length: either tracks the horizon (it. 21 failure) or matches today
  by coincidence. Data: = Lambda (loses Claim 1's gain); as positive baseline disfavoured (it. 13). Pattern across 21-22:
  natural size mechanisms either track the horizon or need a 'why now' coincidence.

## Coalesce (Oct 2026): it is common sense that the size is tied to our horizon - we cannot see past it in any direction.
- Agreed in substance: everything we can measure lies inside our horizon, so a dark-energy size set by the horizon is
  natural. Correction of picture: the Big Bang has no location (it happened everywhere); our horizon is a sphere around
  US, set by how far light has travelled.
- Prior art taking exactly this view: holographic dark energy (Cohen, Kaplan & Nelson 1999: vacuum energy limited by the
  horizon size, rho <~ M_P^2/L^2 -> right size). Using the PAST/Hubble horizon fails (it tracks matter, iteration 21's
  failure); Li 2004 uses the FUTURE event horizon (how far we will EVER see): accelerates, w ≈ -1, and can cross -1 -
  but for c < 1 it crosses from w > -1 (past) to w < -1 (future): the OPPOSITE order to Claim 1 and to DESI's
  preference. Test to run: fit Li's model (1 parameter c) to the same data and compare with Claim 1.
- Coalesce's correction (kept in mind, parked): 'no location for the Big Bang' holds INSIDE our horizon (CMB looks the
  same in every direction to 1e-5: no centre visible). Beyond the horizon we cannot know: a larger structure with a
  'birthplace' is not excluded (e.g. bubble universes in eternal inflation, where each bubble has a centre in the larger
  space). Possible test if revived: a large-scale dipole or gradient at the edge of what we see (CMB large-angle
  anomalies) could hint at it.
- It. 23 (horizon-tied size: Li 2004 holographic DE, future event horizon, 1 parameter c): Delta chi2 vs Lambda +44.5 /
  +45.5 / +41.5 (Claim 1: -5.4 / -7.2 / -6.8). Best c ~0.7 -> w0 ≈ -1.10 (phantom NOW, crossing the wrong way); BAO, SN
  and CMB all object. Horizon-based size: right magnitude, holographic SHAPE strongly excluded; supports Claim 1's
  crossing order. (First run had a z <= 30 truncation bug; fixed, no material change.)

## Trampoline: Pools (halos) sag the sheets and add tension (Synthesis + Coalesce, Oct 2026) -- CLOSED for the size/history
Motivated by the sky checklist (`quantum_gravity/size/ocean_sky_checklist.md`): DM is constant to a few % since the CMB, so DE
cannot be drained from it; only its arrangement can matter. Pre-registered (commit 2abb7e2), then fitted: tension ~ pooled
fraction gives Delta chi^2 +32 to +50; tension ~ sag energy +282 to +596. Bare size 4.6e-7 of the need. The local (galaxy)
version is Family 1 (static links -> Newton). Lesson: the data read the expansion clock (peak at q = 0), not the structure
clock. Do not re-propose clumping-driven dark energy without a reason it would peak at the onset of acceleration.

## Sheet <-> gap crossing sets the dark/ordinary ratio (Coalesce question, Oct 2026) -- OPEN, borrowed (asymmetric DM)
Shared leftover charge until crossing shut at T_f > ~7 GeV -> Ocean particle ~8.5-13 GeV (simplest content); no annihilation
signal; nucleon cross-section <~1e-49 cm^2 (below the neutrino fog). Pre-registered at bda169c; consistent with LZ/XENONnT nulls
(weak test). Trades 5.36 for a mass; would become an explanation only if the grid fixes m ~ 10 m_p. `ocean_origin/`

## Anglerfish: the Ocean grows heavy by eating free electrons (Coalesce, Oct 2026) -- CLOSED (budget)
All electrons weigh 4.8e-4 of ordinary matter; the Ocean is 5.36x ordinary matter -> eating every electron supplies 1/11,000
of its mass. Also the sky: electrons are not missing (CMB recombination, charge neutrality), and Bullet-Cluster gas (full of
free electrons) is not touched by dark matter. The KEPT kernel: 'heavier because built from the same kind of pieces, bundled'.
Literature version: twin (mirror) sector -- a copy of the strong force in the gaps makes twin baryons ~5 Lambda'_QCD, mass
2.5-100 GeV, ratio ~5 without tuning the mass to the proton (Garcia Garcia, Lasenby & March-Russell, PRL 115, 121801, 2015;
Farina 1506.03520). It predicts dark radiation Delta N_eff ~0.075-0.16 (our P3 window: < 0.107). NEXT: test 'gaps hold a twin
copy' against the checklist and Delta N_eff before adopting (it is a new ingredient; must bring predictions).
  Iteration 26 result (ocean_origin/README.md): twin neutrino-only passes Delta N_eff (0.06-0.07 < 0.107); twin photon fails
  (0.13-0.16); self-collisions pass; weight consistency FAILS by 1.3-2.4x in the crude sharing (needs m_b' 8.6-9.5 vs 12.5-21 GeV).
  Parked, not adopted. Re-test only with the full sharing (charge neutrality + sphalerons) and the same pre-set window.
  Iteration 27 (full sharing: neutrality, lepton number, masses): needs m_b' 7.7-9.8 GeV, never in the pre-set 12.5-21 window.
  Fraternal twin Ocean CLOSED in this form. Do not re-propose without a new reason for lighter twin quarks plus a new prediction.

## The book: time has a fixed total of pages x; 'page n of x' (Coalesce, Oct 2026) -- OPEN (borrowed: vacuum-energy sequestering)
Size from the TOTAL length (fixed fee, not recomputed): rho_DE ~ rho_P/x^2 needs x ~ 3e61 Planck ticks = 52 Gyr total, 3.7x
today's age (factor ~1-12x for an O(1) constant). 'Why now' becomes 'a random page of a finite book'. Kaloper & Padilla 2014:
finite total lifetime required (recollapse). Conflicts with pre-registered P4 (coasting future): adopting it = new model.
Sheets -> book length (iteration 28, book/): sheets = pages fails by >40 orders; sheets setting the tick (species bound) makes
the book end 800x-1e7x too soon. Sheets are side by side, not pages; x is a separate number. NEXT: test the recollapse
version's w(z) today (pre-register before fitting).
  Iteration 29 (book/iter29_can_the_book_close.py): Claim 1 + any negative constant has exactly one positive H at every a
  (analytic + numeric): no turnaround, coasting forever (same as the toy, iteration 14). Our mechanism FORBIDS a crunch, so a
  finite book could only end by the grid stopping, with no precursor. Book PARKED: untestable within this model; P4 unchanged.
  Reopen only if data show dark energy heading below zero (which would also contradict Claim 1).

## Count every cell of the visible universe, whole past (Coalesce: 'what we can see is all there is', Oct 2026)
Proper 4-volume of today's visible region (comoving radius 46.5 Gly) over the whole past = 1.9e245 Planck cells; the counting
rule rho_DE ~ 1/sqrt(N) gives 2.3e-123 vs measured 1.1e-123 (factor 2.1). This is Sorkin's everpresent-Lambda estimate (predicted
before 1998). Counting each of N layers separately ruins it (2e3-2e7 too small), so if counting is right the layers must be
species inside one cell (ledger candidate 1), not separate places. SAME CATCH as iteration 21: the visible region grows, so the
count grows and the size tracks the horizon (fails the early-universe history) unless the count is fixed once ('why today'
returns). Not new; recorded as the best size estimate so far, not a solution.
  Coalesce (Oct 2026): 'time need not last as long as the universe; nothing needs to be counting' (tally is physical, no
  observer; say 'tally/record', not 'count'). Combined with the cell tally: if the size is set by the FINAL tally of the
  visible region (fixed fee, recorded once), 1/sqrt(N) matches the measured size when the tally reaches 4.4x today's, i.e. time
  stops ~6.5 Gyr from now (range ~0.3-14 Gyr for an O(1) rule constant 0.5-2). No crunch (consistent with iteration 29).
  CORRECTION (Coalesce clarified): time lasts exactly as long as the universe and does not need observers. So 'time stops,
  space remains' is NOT the idea: if the book ends, the WHOLE grid stops (space and time together). If instead the universe
  never ends (our coasting future), the final tally is infinite and this size route fails. In the grid, cells tick by
  themselves, so there is always a clock (Penrose's 'no clocks, no time' does not apply here). Untestable in advance in our model; 'why now' becomes typicality. Related literature: Bousso, Freivogel,
  Leichenauer & Rosenhaus 2010 ('time will end', arXiv:1009.4698); Penrose (no clocks without mass in the far future).
  Coalesce's fork (Oct 2026): Branch A, the Big Bang is the true beginning -> the universe (and time) also ends -> finite book,
  size from the final tally, end ~0.3-14 Gyr from now (untestable). Branch B, something existed before the Big Bang -> other
  big bangs could appear elsewhere -> no end -> size needs a memory mechanism. Caveats noted: a beginning does not imply an end
  in physics (it is an assumption, labelled as such); the size match was USED to set the end date, so it is a fit, not support;
  Branch A needs a second, independent number to count as supported. Coalesce's earlier 'smooth grid existed before the Big
  Bang' (arrow-of-time discussion) belongs to Branch B -- the two views must be reconciled.

## Watermark (primordial ripples) from the grid (iteration 30, arrow_of_time/) -- OPEN as an option, borrowed (Mukohyama 2009)
Cell jitter: n_s = 4, excluded. Horava z = 3: n_s = 1, excluded at 8.7 sigma. z = 3 - eta with eta = 0.013-0.017 fits (one number
for one number); predicts zero running (data consistent). Amplitude not predicted. Needs a derivation of eta to count.

## Time from position: 'horizontal and vertical cells give position; the combined position is time' (Coalesce, Oct 2026)
Two readings. (1) Literal: time = a combination of one point's coordinates (e.g. a diagonal) -> FAILS: one diagonal direction
would be special (space is the same in all directions to high precision) and walking backwards along it would be travel back in
time. (2) Global: time = a property of the WHOLE arrangement of cells (total cell count / overall size / how the arrangement has
changed). This is 'relational time' (Barbour, 'The End of Time'; York time = mean expansion of the slice, which is exactly the
K = 3H our preferred slicing uses; cosmologists already use the size of the universe as a clock). In our grid, cells are being
added, so the total cell count is a built-in clock that nobody has to read. Possible link to the book: the page number = total
cell count. Untested; no new prediction yet (relational-time models usually reproduce the same physics).

## Random cell addition driven by electrons popping in and out (Coalesce, Oct 2026) -- CLOSED (both readings)
Real free electrons: 99.98% were captured into atoms at recombination (z ~ 1100) and freed again at reionization (z ~ 7); the
expansion shows no stutter at either, and their density dilutes as a^-3 while the needed addition per volume tracks 3H -> no.
Virtual electron pairs (vacuum fluctuations, real: Lamb shift, Casimir): their total energy (m_e/M_P)^4 = 3e-90 is 3e33x too big
(the electron version of the cosmological-constant problem); the random leftover of pair events in our horizon history is 9e43x
too small; one pair would have to add ~1e29 cells; and their rate is the same at all times, so they act like a constant, not
the past history. The same leftover counted over PLANCK cells gives 6.9e-123 (6x measured; 2x with the whole visible past).
Lesson: if growth is random, the randomness lives at the cell level, not at the electron level.

## 'Time adds it': the random unevenness accumulates instead of being recomputed (Coalesce, Oct 2026) -- does not work alone
(1) Accumulating the random +/- of every cell ever added: leftover ~ sqrt(N_added) spread over a volume ~ N -> density ~ a^-1.5,
i.e. w = -1/2 (dilutes like thin matter): excluded by every w(z) measurement. (2) Accumulating the tension with long memory: the
memory-length fits (toy LESSONS iteration 1) give Delta chi2 worse than Lambda for >= 2 Hubble times; the shape needs ~1/3 Hubble
time (derived: Hubble friction). (3) A separate accumulated constant plus the short-memory jostle: iteration 13 says a POSITIVE
constant part must be small; only a slightly NEGATIVE one is liked (an accumulated random sum can have either sign -- noted).
Lesson: the shape needs forgetting, the size needs remembering; one process cannot do both. Two-part idea still needs a size
for the jostled part.
  Iteration 31 (quantum_gravity/size/): bucket fed by the same kicks, no new parameter: 18x the sponge, phantom at all z,
  Delta chi2 +28 / +35 / +20 vs Lambda -> EXCLUDED (pre-registered). A remembered size needs a separate channel, not the same
  kicks kept longer.

## Shaking makes it forget (Coalesce, Oct 2026) -- agreed for the SHAPE; points to a CONSERVED stamp for the size
Shaking scrambles arrangements and motions (non-conserved things) but cannot change conserved things (a shaken jar of marbles
forgets the arrangement, never the NUMBER of marbles). So: shape = non-conserved jostled tension (forgets, kappa = 3); size =
a conserved quantity stamped once. Borrowed homes: unimodular gravity (Lambda is an integration constant, conjugate to the
4-volume, i.e. to 'cosmic time'; Henneaux-Teitelboim 1989; Sorkin's 1/sqrt(V) is its uncertainty) and projectable Horava gravity
('dark matter as integration constant', Mukohyama 2009). Constraint from iteration 13: a separate POSITIVE constant must be
small, so the stamp should set the jostling's coupling strength A (how hard kicks count), not add a second term.
Untested; value of the stamp still unexplained.

## CORRECTION: the 'within 2x' of the visible-universe cell tally is AUTOMATIC, not a clue (checked Oct 2026)
1/sqrt(N4) of the region visible at ANY epoch equals the total density at that epoch within a factor 1.4-5.6 (a = 0.001 to 1;
analytic reason: N4 ~ (t/t_P)^4 and rho_total ~ 1/(G t^2), so 1/sqrt(N4) ~ rho_total in Planck units always). Today dark energy
is 69% of the total, so the match is guaranteed. It restates the 'why now' coincidence; it does not derive the size.
Using the edge of what we can EVER see (event horizon) = holographic dark energy (Li 2004): already excluded (iteration 23).
Using the region seen at one special moment (e.g. onset of acceleration) is again automatic (dark energy ~ total density then).
The size remains unexplained; the cell tally only says the counting rule's coefficient is O(1).

## Whole-workbook inversion (Coalesce, Oct 2026; iteration 32, book/README.md) -- OPEN, one testable half
Measured dark energy -> total tally of all space x all time. Closed sphere excluded (needs Omega_k <= -0.04..-0.10; data +0.002).
Flat wrap-around space allowed only in a narrow window: cubic torus side 27.5-38 Gpc (rule constant 1; up to ~60 for 2), just above
the CMB circle-search limit; time then ends within ~0-11 Gyr. Testable half: deeper CMB topology searches (COMPACT).

## Project the expansion forward, read dark energy at the end, validate with an independent equation (Coalesce, Oct 2026)
Claim 1 forward (Pantheon+ fit): no end. a ~ t^5 at late times, acceleration forever, dark energy fades slowly (0.80 of today's at
25 Gyr, 0.43 at 62 Gyr, 0.17 at 140 Gyr, 0.03 at 460 Gyr; ~ t^-2) but never runs out. 'How long' = forever in our model; an end
needs Branch A with an external rule (circular for the size). Validation chain that DOES exist: the same equation run backward
-> pre-registered binned w(z) (DESI DR3); and redshift drift (Sandage-Loeb, an equation taking H(z) as input): Claim 1 minus
LCDM = -0.5 to -0.7 cm/s over 20 years at z = 0.5-4 (LCDM drift +5.3 to -10.1 cm/s) -- below expected ELT precision (~cm/s).

## Published checks that take the dark-energy amount as INPUT (Coalesce asked, Oct 2026)
(1) Local dark energy from Local Group dynamics (Chernin, Karachentsev et al., A&A 507, 1271, 2009): local density 0.8-3.7x the
global value (20-30% accuracy). Bearing on our model: if the tension responded to LOCAL expansion (rho ~ adot^-1/2), it would
blow up inside bound regions (adot = 0) -> far above 3.7x -> excluded. Our theoretical home (VCDM / constant-mean-curvature
slicing: K = 3H uniform on each slice) makes dark energy uniform -> PASS. This fixes a reading: the law uses the slice's
global expansion, never local motions. (2) Maximum turnaround radius R = (3GM/(Lambda c^2))^(1/3) (Pavlidou & Tomaras 2014):
Local Group 1.4 Mpc, Virgo ~11.5 Mpc; Claim 1 vs LCDM differ by 1.4% -- far below measurement errors; consistent, not
discriminating. (3) ISW, growth, lensing, redshift drift: already in our fits/predictions or too small to separate.

## Voids: compare void sizes to read dark energy (Coalesce, Oct 2026; iteration 33, voids/) -- real method, weak for us
Void cosmology is established (BOSS DR12 void size function: w = -1.1 +/- 0.2; Euclid voids alone: w ~10%, FoM 17). Claim 1 vs
Lambda (same matter and early ripples): growth differs by <= 0.5% (lower after z ~ 1, a hair higher before); realistic void counts
amplify this 4-14x -> -1% to -3% fewer large voids at z = 0.5, +0.1-0.3% at z = 1. Distinctive sign change, but below void-only
precision. Useful only in combination. (A textbook-threshold first try gave spurious 5-40% ratios from an extreme tail; rejected.)

## Direct leftover check (Coalesce: speed + mass at each moment -> is the push constant?) -- direct_check/leftover.py
Friedmann leftover from real DESI DR2 expansion rates + Planck r_d and omega_m: omega_DE(z) = 0.378 +/- 0.034 (z 0.51), 0.388 +/-
0.038 (0.71), 0.300 +/- 0.032 (0.93), 0.279 +/- 0.067 (1.32), then large errors. One constant fits (chi2 5.4/5 dof; best 0.345 +/-
0.020 vs Planck-Lambda 0.311, ~1.7 sigma). Claim 1 (no fitting) fits equally (chi2 5.4/6). Its predicted peak is only ~1% high,
far below the ~10% errors. Method is right; precision must improve ~10x to see our shape this way.
  How much is left (direct_check/how_much_left.py): today 5.5-6.2e-27 kg/m^3 (~3.3-3.7 hydrogen-atom masses per m^3) across
  readings. Left at 25 / 40 / 60 Gyr: Lambda 1 / 1 / 1; Claim 1 0.81 / 0.61 / 0.44; w0wa formula (if extrapolated) 0.33-0.53 /
  0.07-0.15 / 0.01-0.04. Future part is projection only; w0wa is not a law.
  Wobble check (direct_check/tracer_split.py): leftover from LRG bins (z 0.51, 0.71) = 0.382 +/- 0.026 vs other tracers 0.297 +/-
  0.029: a 2.3 sigma step exactly where DESI switches galaxy type (z ~ 0.9). Claim 1 predicts only a ~1% difference, so even if
  real, the step is ~25x bigger than our model's peak; most likely a fluctuation or tracer effect. Not claimed as support.

## Size as a BALANCE between the grid's tension/strength and slack (Coalesce, Oct 2026) -- discussed, related to iteration 22
The grid's stiffness is Planck-scale (c^4/G; Planck density 5e96 kg/m^3); dark energy is 1e-123 of it, so a balance needs a tiny
ratio from somewhere. Natural ways to write it: (stiffness) x (l_P/L)^2 with L ~ horizon (tracks the horizon: the percentage-tip
failure; stick-slip threshold of iteration 22 came out as 'one cell per horizon-long row', same issue), or (stiffness) x (l_P/L)^4
with L = 88 micrometres (energy 2.2 meV). The second predicts gravity changes below ~0.1 mm if the balance length is physical;
Eot-Wash torsion balances see Newton's law hold to ~40-50 micrometres -> a simple 'slack length ~0.1 mm' is already squeezed.
Coincidence noted (not used): 2.2 meV is near the neutrino mass scale.

## Link tolerance: stretch to a limit -> elasticity gives (slack); break -> leak (Coalesce, Oct 2026) -- Family: iteration 22
Readings: (1) BREAK = a new cell is inserted: a natural answer to 'why cells are added' (expansion stretches links to tolerance,
they break, a cell fills in); cell size stays fixed on average (consistent with constant G, c). (2) Stored strain just below
tolerance acts like a constant energy (fixed fee): stiffness x strain^2 = dark energy needs a tolerance strain of 3.4e-62, about
0.3x cell size / horizon size today (1.2e-61): the 'why now' coincidence again (same as iteration 22's threshold). (3) A
tolerance that depends on how fast links are pulled (strain rate ~ H) gives rho ~ H^-1/2, the s = 0 law, excluded at 5-6 sigma
(stretch_scan). (4) LEAK: energy leaking from the tension into dark matter or radiation is limited by DM constancy (few %) and
Delta N_eff < 0.107. Reading (1) kept as a picture; (2)-(3) do not set the size.

## Thresholds: how much dark energy it takes to affect matter vs the grid (Coalesce, Oct 2026)
To make expansion speed up today: rho_DE > rho_m/2 = 1.3e-27 kg/m^3 (measured 5.9e-27 = 4.3x this). Effect on the grid tension:
1.1e-123 of the stiffness (negligible; the grid does not notice). Upper window (Weinberg 1987): dark energy must not take over
before galaxies assemble: rho_DE < ~rho_m(z_gal) = 7e-26 to 9e-25 (z 2-6); measured is 0.6-8% of that maximum. So the value sits
in a 'Goldilocks' window: big enough to matter now, small enough for galaxies. This EXPLAINS the size only with selection
(many regions with different random stamps; observers arise only where galaxies form) = Coalesce's Branch B (other big bangs
elsewhere) + a random per-region stamp. Borrowed: Weinberg 1987 (predicted a small non-zero value before 1998), Martel, Shapiro &
Weinberg 1998. Not adopted; it is selection, not a mechanism.

## Quantum-gravity steps run proactively (Oct 2026)
- Iteration 34 (growth/): link breaking: explains cell addition and sqrt(N); cannot give the 1/3-Hubble memory (needs links to
  stretch ~165%, storing ~Planck energy) or the size. Memory stays with Hubble friction.
- Iteration 35 (arrow_of_time/): watermark eta from layers saturating the species bound: natural size c/(16 pi^2) = 0.003-0.019
  brackets the needed 0.013-0.017, but coefficient and sign are uncomputed -> consistent in size, not a prediction.
- Iteration 36 (stamp/): conserved stamp built into the action (unimodular-type, ghost-free): C is exactly conserved (shaking
  cannot erase it), but the classical action allows any value; hbar uncertainty gives 1e-246 (far too small); cell discreteness
  gives 1/sqrt(N4) = the cell tally. All three QG routes converge on the same fork: past tally (fails history) or whole finite
  history (Branch A book; testable half = small flat wrap-around universe).
  Iteration 37 (direct_check/iter37_both_options.*): both size options vs the direct leftover. Tracking (option 1) EXCLUDED by
  low-z data alone (chi2 154/5; 37/3 without red-galaxy bins). Fixed stamp (option 2) fits (5.4/5; 0.14/3). The 'rise toward
  today' (p = -0.72, ~1.9 sigma) disappears without the red-galaxy bins (p = -0.10 +/- 0.8) -> not robust; no revision.
  Iteration 38 (book/): wrap-around cube with side 27.5-38 Gpc vs real Planck low-l TT: Delta chi2 -0.6 to -0.1 (quadrupole 3-14%
  lower) -> allowed, not decisive. Testable half survives; decisive test = repeated-pattern / anisotropy searches (COMPACT).

## QG: Einstein gravity from the grid at large scales (iteration 39/39b, quantum_gravity/emergence/)
Well-shaped cells: grid action -> Einstein action (ratio 1.12 at 384k cells, error ~1/m, converging to 1). Fully random Delaunay
grid: no convergence at 4k-64k cells (scatter ~1, sliver cells). Constraint: cells must be well-shaped (as in CDT) or gravity read
off a coarse-grained grid. First run failed partly by design (bump too narrow); one revision bug (corner points) caught by the
flat-space check and discarded.

## QG: black-hole entropy 1/4 (iteration 40, quantum_gravity/entropy/) -- consistent, borrowed theorem
With gravity induced by the N layers (Sakharov), entanglement entropy = A/(4G) automatically (Jacobson 1994; Susskind-Uglum
1994) for scalar fields like our tension field. The grid's measured entanglement (s0 = 0.0603 per cell area) then fixes
a = 0.49 sqrt(N) l_P (species-bound factor, previously unknown). Window: N 2.9e7-5.2e15, cell 4.3e-32 to 5.7e-28 m.
  Coalesce (Oct 2026): maybe tidy cells are needed only because of HOW curvature is defined. Agreed in principle: the
  tidy-cell requirement belongs to Regge's definition (deficit angles on single links). Definitions built for randomness are
  known to converge on random grids: Ollivier-Ricci curvature (optimal transport between neighbourhoods; van der Hoorn et al.,
  PRR 3, 013211, 2021) and the causal-set Benincasa-Dowker action (averages over many points; converges in the mean with a
  smearing scale). Revised constraint: either tidy cells + Regge, or random cells + a smeared/nonlocal curvature definition.
  NEXT: test Ollivier-Ricci on our random grid against the known curvature.
  Iteration 41 (emergence/): Ollivier-Ricci on the random grid: ratio to the predicted curvature 18 -> 13 -> 3.7 -> 3.0 +/- 0.9 for
  neighbourhoods of 300 -> 3700 points; right sign, clear trend, but the pre-set criterion (within 25%) is not met and the kill
  line (>50% off) is crossed at reachable sizes. Inconclusive; far better than Regge on random grids (no trend).
