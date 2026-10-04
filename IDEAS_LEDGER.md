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
