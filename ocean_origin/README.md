# Where the Ocean came from: crossing between sheets and gaps (iteration 25)
Question (Coalesce): can dark matter become ordinary matter or vice versa?
Sky: today, almost never (no annihilation/decay signal; both amounts constant to a few % since the CMB). Early on, possibly.
Test: if crossing shared one leftover charge until it shut at T_f, does the ratio 5.36 come out? (BORROWED: asymmetric DM.)
Predictions are in `iter25_predictions.txt` and were committed BEFORE checking detector limits.

## Results (computed and committed at bda169c, then compared with the literature)
- **Number ratio:** with the simplest gap content (one particle type carrying a whole unit of charge), the shared charge
  splits about evenly (number ratio ~0.6-1). So the ratio 5.36 becomes a mass: **the Ocean particle weighs ~8.5-13 GeV**
  (9-14 proton masses).
- **When crossing had to shut:** above T_f ~ 7 GeV (about a nanosecond after the start). If it shut later, the Ocean particle
  is already too slow to hold its share, and 5.36 cannot be reached for any mass. Content carrying 1/3 unit needs ~100 GeV.
- **Detectors (P25c):** contact coupling fixed by "crossing shuts at T_f" gives sigma <~ 1e-49 cm^2. LZ (Dec 2025) and
  XENONnT see no dark matter at 3-9+ GeV and have reached the solar-neutrino "fog" (boron-8 CEvNS seen at 4.5 sigma).
  Present limits near 10 GeV are several orders of magnitude above 1e-49. **PASS, but weak**: the prediction is "never
  seen in these detectors".
- **Annihilation (P25b):** no anti-Ocean is left, so no signal is expected. Consistent with null searches. Watch item: the
  Fermi Galactic-centre GeV excess (often fit by ~40-50 GeV annihilating DM; pulsars are the favoured explanation). If it were
  ever shown to be dark-matter annihilation, P25b would fail.
- **Honest verdict:** this does NOT explain 5.36. It trades it for one number (the Ocean particle mass ~9-13 GeV) and a
  window (crossing shut above ~7 GeV). It becomes an explanation only if the grid fixes the Ocean particle's mass at
  ~10 proton masses. It is consistent with every sky row in `../quantum_gravity/size/ocean_sky_checklist.md`
  (heavy enough for Lyman-alpha; collisionless allowed).
Sources: LBNL news 8 Dec 2025 (LZ results); Kaplan, Luty & Zurek 2009 (asymmetric DM).

## Iteration 26: twin copy of matter in the gaps (checks committed at 7043860, then compared)
Data: N_eff = 2.990 +/- 0.070, Delta N_eff < 0.107 (95%) from CMB + BAO + primordial abundances (arXiv:2603.13226);
ACT DR6: N_eff = 2.86 +/- 0.13 (arXiv:2503.14454).
- **A. Dark radiation:** twin neutrino only: Delta N_eff = 0.062-0.074 -> PASS (about 1.8 sigma above the measured central value).
  With a twin photon (own light force): 0.133-0.159 -> FAIL. So: no twin light force, which also removes the dark-disk risk (D).
  Shared budget: the jostling radiation (P3) now has at most ~0.03-0.045 left.
- **B. Weight consistency:** the sharing of iteration 25 needs twin b' quarks of 8.6-9.5 GeV (dark matter particle ~26-29 GeV);
  twin Higgs expects 12.5-21 GeV (f/v = 3-5). **FAIL by a factor 1.3-2.4** (it would need f/v ~ 2.1-2.3, below the usually
  quoted f/v >~ 3). Caveat stated now: the sharing ignored charge neutrality and sphaleron effects, which are of this same size,
  so a full calculation could move it either way. A re-test is allowed only with this SAME window (12.5-21 GeV) fixed in advance.
- **C. Self-collisions:** sigma/m ~ 2e-6 to 1e-4 cm^2/g -> PASS (limit ~1 cm^2/g).
- **Verdict:** a twin world with no light force passes the sky checklist and the radiation limit, but the weight does not come out
  naturally in our crude sharing. Not adopted; parked as the best-motivated Ocean candidate, pending the full sharing calculation.

## Iteration 27: full sharing re-test of the twin weight (window fixed beforehand at 7043860)
`iter27_full_sharing.py`: chemical equilibrium with electric neutrality, separate lepton number (L = 0, or the standard sphaleron
leftover L = -(51/28) B), all masses (top, W, b, c, tau) and W bosons. Needed twin b' mass: 7.65-9.75 GeV (f/v = 1.8-2.3)
across T_f = 7-100 GeV and both lepton choices; the crude estimate gave 8.4-9.5. The extra effects moved it AWAY from the
window. **FAIL: never inside 12.5-21 GeV.**
Verdict: the fraternal twin Ocean (dark matter = three-b' twin baryon, sharing through baryon-number crossing) is CLOSED in
this form. What would still be open (each a new free choice, so an accommodation unless it brings a prediction): a twin with
lighter twin quarks (non-mirror Yukawas), a different crossing charge, or binding effects. Not pursued.

## Iteration 59: two-level rule for the Ocean particle mass (pre-registered: PREREG_59.md)
- GR level: gravity sees only density -> its twin is 'cold and collisionless' (passed). No number to compare.
- QG level, 11 pre-declared grid quantities vs 8.5-13 GeV: two hits -- m_e/alpha^2 = 9.60 GeV (blind: I had misjudged it as MeV
  before computing) and m_p/sqrt(alpha) = 10.98 GeV (not blind). Layer doubling 2^k misses (7.5 / 15 GeV).
- Look-elsewhere: ~6% chance of >= 2 hits by luck. NOT significant. Recorded as two leads; m_e/alpha^2 sits at the top of the
  simplest-content range (8.6-9.5 GeV) and echoes the electron's pop ladder (alpha m_e = what each electron pays dark energy).
- Learned: the dark-matter mass is invisible to GR; only the quantum side can fix it, and right now nothing derives it.
