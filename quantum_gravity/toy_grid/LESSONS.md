# Grid toy: what each iteration taught us

## Iteration 1: does the grid tension respond instantly? (iter1_memory.py, iter1_wz.txt, iter1_fits_*.json)

**Rules.**
- Fixed Planck cells, added by division; daughters inherit their parent's tension.
- Comoving domains carry an intensive tension, kicked ± at their walls at a rate ∝ 1/(physical domain size).
- The tension relaxes over a memory time of 1/(κH).
- Energy is either LINEAR in |tension| (string-like; the instant limit is Claim 1) or QUADRATIC (spring-like; the instant
  limit is β = 1).

**Check.** Because division copies the tension, the whole grid reduces exactly to one random variable per domain. A Monte
Carlo of 5000 domains agrees with the one-variable equation to 1–2% (iter1_partA.txt).

**Fits (DESI DR2 + Planck priors + supernovae), Δχ² vs Λ:**

| energy | memory (Hubble times) | Pantheon+ | DES-Dovekie | Union3 |
|---|---|---|---|---|
| linear | 4 | +4.1 | +5.6 | +3.1 |
| linear | 2 | +1.8 | +2.7 | +1.3 |
| linear | 1 | −1.6 | −1.7 | −1.6 |
| linear | 0.5 | −4.1 | −5.0 | −4.1 |
| linear | 0.25 | −5.2 | −6.5 | −5.4 |
| linear | 0.06 | −5.4 | −6.9 | −6.1 |
| quadratic | 1 | −1.9 | −2.0 | −1.9 |
| quadratic | 0.5 | −6.5 | −8.2 | −7.0 |
| quadratic | 0.25 | −7.0 | −9.2 | −9.0 |

**Lessons.**
1. **Memory must be short:** ≲ ½ a Hubble time for the linear case, ≲ ¼–½ for the quadratic. With a memory of one Hubble
   time (the "natural" value in the old counting story), the law barely beats Λ. With 2 or more Hubble times it does
   worse than Λ.
2. **Why:** a lagging tension crosses w = −1 later than the acceleration onset. Crossing redshift at q = 0 (z ≈ 0.68):
   memory 1 → z = 0.16; ½ → 0.35; ¼ → 0.49; 1/16 → 0.63. The instant law crosses exactly at q = 0.
3. **New signature:** any real memory makes the crossing come AFTER acceleration starts, never before. The gap between the
   two redshifts measures the grid's memory time. A binned w(z) from DESI DR3 or Euclid can bound it.
4. Linear (√N) and quadratic (N) energies both work if the memory is short. The quadratic one fits slightly better on
   DES and Union3, in line with the joint β–s scan.
5. The old counting story ("N crossings per stretch time", memory = one Hubble time, instant response) was internally
   inconsistent. A memory of a Hubble time implies a Hubble-time lag, and the data reject that lag.

## Iteration 2: where do the kicks come from, quantum mechanically? (iter2_quantum_kicks.py/.txt)

**Setup.** Treat the link field as a free massless quantum field: the long-wavelength limit of the quantised tension
links. A domain's tension is the field averaged over the domain and accumulated with memory rate γ. Compute its
quantum variance in two states: the pure vacuum, and a relic thermal background at temperature Θ.

**Results** (numerical, matching the closed forms):
1. **The vacuum alone gives NO random walk.** Its variance grows only as ln(memory)/(4π²); the growth exponent falls
   0.33 → 0.08 instead of 1. The fluctuations of a massless field are anti-correlated in time. A vacuum-only grid gives
   a near-constant dark energy (Λ-like), not Claim 1.
2. **Real thermal quanta give an honest random walk:** Var = Θ/(4πγ), independent of domain size (checked for R = 0.5–4,
   to better than 0.1% once memory ≫ R and ≫ 1/Θ).
3. **Expanding universe, conformal field:** Var = Θ_c/(4π κ ȧ) ∝ 1/ȧ. This is the counting of Claim 1 with NO comoving
   length put in. The physical kick rate Θ(t)/2π ∝ 1/a is exactly iteration 1's rate, now derived. The puzzle of
   law_from_grid.md section 7 ("which comoving length?") goes away: the 1/a is the redshifting temperature of relic
   quanta.
4. **For CMB photons today** (κ = 1, Planck window), the thermal variance exceeds the vacuum variance by about 4×10²⁷.

**Limits** (independent check, all items PASS):
- Valid only for a CONFORMAL field (photons, or a scalar with ξ = 1/6). NOT for relic gravitons or minimally coupled
  fields, which feel the expansion at horizon scale, exactly where the relevant modes sit. Photons are a vector, so the
  scalar treatment is an analogy.
- The formula "∝ 1/ȧ exactly" needs κ ≫ |q|. Otherwise the full memory equation of iteration 1 applies, and iteration 1
  shows the data want κ ≳ 2–4 anyway.
- Still put in by hand: the memory rate κH, linear vs quadratic energy, and the size of the coupling (so the amount of
  dark energy is still unexplained).

**What this changes.** Claim 1's √N counting now has a concrete physical reading. The grid tension random-walks because
it is jostled by REAL relic quanta of a conformal field, presumably the CMB photons. The amount of jostling over a
memory time ∝ 1/H scales as Θ/H ∝ 1/ȧ.

**New questions for iteration 3:**
- (i) If the jostlers are CMB photons, the grid tension couples to photons. Which existing limits (CMB spectral
  distortions, photon propagation) bound that coupling?
- (ii) What sets a memory of ≲ ½ a Hubble time?
- (iii) Before e⁺e⁻ annihilation and neutrino decoupling, the photon temperature did not fall exactly as 1/a. That is
  irrelevant for late dark energy, but it is a check.

Prior-art relatives, none the same: "Thermal dark energy" (Hardy & Parameswaran, PRD 101, 023503: a hidden-sector
finite-temperature potential); spacetime-diffusion models (PRD 7whh-9j22); stochastic dark energy from inflationary
fluctuations (EPJC 78, 5862).

## Iteration 3: can CMB photons be the jostlers? (iter3_who_jostles.py/.txt)

**Results.**
1. **Photons cannot do it through a linear coupling.** Gauge invariance lets the tension feel only E or B, and ∫E dt = −ΔA
   is a boundary term. The filtered variance SATURATES (2.27×10⁻² for memories 10³–10⁵) while the scalar's grows in
   proportion to the memory. There is no random walk; B is worse (an extra factor k).
2. **Photons through their energy density also fail.** By scaling (not computed in detail), the noise ∝ Θ⁵/R² gives a
   stretch exponent s = 5. The data allow s = 0.72–1.06, so this is excluded.
3. **What works is a massless, conformally coupled SCALAR relic** ("grid radiation"), presumably the quanta of the
   tension links themselves. It must have a thermal population, and that population counts as extra radiation:
   - Shared a temperature with ordinary matter above the electroweak scale (or born with the grid at the SM temperature):
     ΔN_eff = 0.027, temperature today 0.90 K.
   - Last in contact between the electroweak and QCD transitions: ΔN_eff ≈ 0.056.
   - Last in contact below QCD: ΔN_eff = 0.30. EXCLUDED: measured ΔN_eff < 0.107 at 95% (N_eff = 2.990 ± 0.070,
     arXiv:2603.13226).

**Lessons.**
- The "jostled by the CMB" idea from iteration 2 is dead. The jostlers must be a NEW relic: a scalar, massless and
  conformally coupled.
- **[Weakened in iteration 4: holds only if s shared a temperature with ordinary matter; in general 0 < ΔN_eff < 0.107.]**
  **First testable prediction of the quantum-grid toy:** extra radiation ΔN_eff ≈ 0.027–0.057. That is allowed today; the
  upper half is within reach of Simons Observatory (σ ≈ 0.05) and CMB-S4-class surveys (σ ≈ 0.03). A firm ΔN_eff = 0
  (σ ~ 0.01, beyond planned surveys) would kill the thermal-jostling mechanism, unless the scalar was born colder than
  ordinary matter.
- **New requirements:**
  - The scalar must barely couple to ordinary matter today: a massless scalar coupled to matter would give a fifth force
    (equivalence-principle and Solar-System bounds).
  - Conformal coupling (ξ = 1/6) is needed for its thermal spectrum to stay exact. Whether the grid produces that
    coupling is an open question.

**Next (iteration 4):**
- What sets the memory (≲ ½ Hubble time)?
- Whether the energy is linear or quadratic in the tension (√N vs N).
- How to estimate fifth-force limits on the scalar.

## Iteration 4: what sets the memory, and what must the jostling relic be? (iter4_memory_and_source.py/.txt)

**Results.**
1. **The memory is derived, not chosen.** Treat the jostled tension as the VELOCITY of a grid field Φ, so that it is
   jostled and also slowed by the expansion (Hubble friction).
   - For a minimally coupled Φ, the friction rate is exactly 3H (checked numerically): κ = 3, memory ⅓ Hubble time.
     That is inside the range iteration 1 found the data need (≲ ½).
   - A conformally coupled Φ would give about 0.9H, which the data disfavour. (I had guessed 2H; the numerical check
     corrected that.)
2. **Fits with κ = 3** (no free dark-energy parameter), Δχ² vs Λ:

   | energy | Pantheon+ | DES-Dovekie | Union3 |
   |---|---|---|---|
   | linear in tension (√N) | −4.9 | −6.0 | −5.0 |
   | quadratic (N) | −7.1 | −9.2 | −8.5 |
   | for scale: instant Claim 1 | −5.4 | −7.2 | −6.8 |
   | for scale: w0wa (2 free parameters) | −7.1 | −10.3 | −13.2 |

   Quadratic with κ = 3 crosses w = −1 at z = 0.46 (acceleration starts at 0.71); w0 = −0.90.
   **Look-elsewhere caveat:** four discrete variants were tried (linear/quadratic × κ = 2 or 3), so the best one gets a
   small selection bonus.
3. **The relic needs a true thermal tail.** A relic made by decays without self-interaction gives no random walk
   (its variance grows only as a log: 0.09 → 0.27 against 7.8 → 7957 for thermal). The relic must be self-thermalised:
   a conformal s⁴ self-coupling will do.
4. **Correction to iteration 3.** Inflation (which the model needs) dilutes any earlier relic, and a conformal field is
   not produced by inflation. So s must come from reheating and thermalise itself.
   - Its temperature is set by the inflaton's branching ratio, so the iteration-3 value ΔN_eff = 0.027–0.057 holds only
     if s also shared a temperature with ordinary matter.
   - That contact cannot have gone through the Higgs: a Higgs portal strong enough to thermalise s would give it a mass
     ~10⁸⁷ times too large.
   - **Weakened prediction:** 0 < ΔN_eff < 0.107. The amount is not fixed; it is absorbed into the (already free) size
     of dark energy.
5. **Fifth force.** A symmetry s → −s leaves no linear coupling to matter, so there is no long-range fifth force at tree
   level.
6. **Naturalness.** s must be lighter than about 10⁻³² eV, which is as unprotected as any quintessence field (an open
   problem, not ours alone).

**Picture after four iterations:**
- A minimally coupled grid field Φ, whose velocity is the tension.
- It is jostled by a self-thermalised dark radiation s (conformal, massless, with an s → −s symmetry, made at reheating).
- It is damped by Hubble friction (κ = 3).
- The dark energy is the tension energy (quadratic fits best).
- This gives a zero-free-parameter dark energy with Δχ² −7.1 / −9.2 / −8.5 vs Λ, crossing w = −1 at z ≈ 0.46.
- Open: the size of dark energy, the s coupling, and why Φ is minimally coupled and s conformally.

## Iteration 5 (Coalesce's question): do the jostling's vibrations tell us the height? (iter5_vibrations.py/.txt)

**Height: no.** The height is (coupling)² × (energy scale) × Θ/(κH). Both the coupling and Θ are free (0 < ΔN_eff < 0.107).
The vibrations fix how the height changes over time, not its value.

**Lumpiness: yes, and it is a problem.** The vibrations also make dark energy lumpy, and the lump-to-mean ratio does NOT
depend on the coupling.
- The tension's spatial correlation is C(d)/C(0) = (1 − e^(−d))/d, with d in units of c/(κH). That is a correlation
  length of about a third of the Hubble length, roughly 1500 Mpc.
- With quadratic energy, each such patch has δρ/ρ ~ √2.
- Rough estimate: potential Φ ~ 0.16, so the large-angle CMB imprint would be ~0.16, against the ~10⁻⁵ observed: about
  2×10⁴ times too strong.
- This is an order-of-magnitude estimate, not a Boltzmann calculation. The margin is large.

**Lesson.** The patch-by-patch stochastic picture (iterations 1–4) is in serious trouble as it stands. Escape routes:
- (a) ~10⁸⁻⁹ independent jostled components, which shrink the lumps by 1/√N;
- (b) jostling that acts on the universe as a whole: a global response, as in the VCDM embedding of Claim 1, which has
  no dark-energy lumps;
- (c) noise that is not dominated by long waves.

Claim 1 itself (a law depending on the global expansion) is NOT affected; this hits only the microscopic toy.
**Next:** decide between (a), (b) and (c), each with its own test.

## Iteration 6 (Coalesce's feedback: each layer has its own jostle; stacked they give a smoother, harmonious wave) (iter6_layers.py/.txt)

**A. Stacked layers turn escape (a) into escape (b).**
- Lumpiness of N stacked layers: δρ/ρ = √(2(1/N + c²(N−1)/N)), where c is the layer-to-layer correlation of the jostles.
  Monte Carlo agrees to 4 digits. (My first formula had c instead of c²; the Monte Carlo caught it.)
- As N grows, the stack's total becomes smooth and deterministic (law of large numbers), so the stack responds as one
  whole. That is the global response of escape (b).
- **Requirement:** CMB imprint < 10⁻⁶ needs N ≳ 3×10¹⁰ independent layers, AND correlation c ≲ 6×10⁻⁶.
- 74 layers (the count hinted at by the electron-mass hierarchy) give an imprint of 0.02, far too lumpy, if those are the
  same layers.
- **Subtlety:** "harmonious" must mean averaging, not synchronising. Layers that lock in step (strongly coupled jostles)
  bring the lumps straight back.

**B. Each layer needs its own jostlers.** That means N relic species, with ΔN_eff = N (4/7)(T_s/T_ν)⁴ < 0.107. For N = 10¹⁰
the relic must be colder than ~4 mK. That is allowed; the amplitude goes into the free height.

**C. Escape (c), short-wave jostling, is CLOSED.** For noise weight ∝ kᵖ:

| p | Var ∝ γ^−x | correlation length (c/γ = 1000) |
|---|---|---|
| 0 | x = 1.00 | 1597 |
| 0.5 | x = 0.51 | 529 |
| 1 | x = 0.13 | 56 |
| 2 | x = 0.00 | 3.5 |

The 1/ȧ law (x = 1) comes only from waves of about horizon size. A short correlation length destroys the law and leaves a
Λ-like constant.

**Lessons.**
- The ONLY surviving microscopic picture is Coalesce's stack: a very large number (≳ 3×10¹⁰) of nearly independent
  layers, each with its own relic jostlers, whose sum behaves as one smooth global tension.
- The layers must be stacked but not jostled together.
- **New question:** what are these layers? Not the 74 of the hierarchy, and not time layers (the time accumulation is
  already the memory).

## Iteration 7 (Coalesce: the gaps between layers are big enough that layers barely interfere with their neighbours) (iter7_gaps.py/.txt)

**Result.** Iteration 6's requirement (c ≲ 6×10⁻⁶) assumed a jostle SHARED by all layers. With gaps, any correlation falls
off with distance in the stack. Then lumps = √(2[1 + 2Σc_k²]/N): only the total leakage to neighbours counts, and it just
reduces the effective number of independent layers. Monte Carlo agrees to ~1% (N = 1000–2000; leakage per gap r = 0–0.99).

| leakage per gap r | cost factor on lumps | layers needed (imprint < 10⁻⁶) |
|---|---|---|
| 0 | 1 | 2.6×10¹⁰ |
| 0.5 | 1.3 | 4.3×10¹⁰ |
| 0.9 | 3.1 | 2.5×10¹¹ |
| 0.99 | 10 | 2.6×10¹² |

All of these sit inside the species window (≤ ~10¹⁵), even with heavy leakage.

**Lessons.**
- Coalesce's gaps rescue the picture. Neighbour interference is harmless; it only costs a modest factor.
- **The real requirement is NO COMMON JOSTLE:** nothing that shakes all layers together (a single shared relic radiation,
  or one long-range mode spanning the whole stack). Even a 1% common component leaves lumps of 1.4×10⁻², which is excluded.
- **Picture (no separate test):** the gaps hold the Ocean, which is cold and non-relativistic, so it would not carry fast
  jostles between layers. That is consistent, but only a picture.
- **Gravity is common to all layers**, but it is a smooth drive (the expansion), not noise. That is allowed.

## Iteration 8: leftover lumps vs the real large-angle CMB (iter8_cmb_lumps.py/.txt, _variant.txt)

**[CORRECTED in iteration 9:** the time integral was weighted with √dt instead of dt, which overstated the bound about
100×. Corrected numbers: N_eff ≥ 1.1×10⁷ – 2.9×10⁷. The numbers below are the buggy first version (iter8_*_BUGGY.txt).**]**

**Method.**
- Lumps are Gaussian (many layers), with spatial correlation ((1 − e^(−x))/x)² on the scale c/(3aH) and decorrelation
  rate 6H.
- Their gravitational dip enters through comoving Poisson. Their CMB imprint is the late ISW effect, computed as a full
  double time integral for ℓ = 2–30.
- The added power is compared with a Planck-like ΛCDM sky at cosmic variance.

**Result.**
- One layer adds D₂ ~ 2–6×10¹² μK², against 1022 observed. The lump power falls steeply with ℓ (D₉ is about 2% of D₂), so
  the quadrupole and octopole dominate.
- 2σ bound: **N_eff ≥ 4.8×10⁹** with plain Poisson, or **≥ 1.7×10⁹** with the super-horizon suppression k²/(k² + 3a²H²).
  This is consistent with iteration 6's rough 2.6×10¹⁰ (that was stricter by 5–15×).
- With neighbour leakage r per gap, multiply by (1 + 2r²/(1 − r²)).

**Window still OPEN:** about 2×10⁹ – 5×10⁹ ≲ N ≲ 1.2×10¹⁵ (species bound + gamma-ray-burst timing).

**Signature if N is near the lower edge:** extra power concentrated in the lowest multipoles, ℓ = 2–3, falling steeply with
ℓ. The real sky's quadrupole is LOW (the known anomaly), so this signature is disfavoured rather than hinted at.

**Not included (caveats):**
- matter falling into the lumps;
- CMB lensing and large-scale galaxy clustering from the lumps (could tighten the bound);
- a proper relativistic treatment beyond the horizon (bracketed by the two variants).

## Iteration 9: matter falls into the lumps; a bug found and fixed (iter9_matter_response.py, iter9_isw_check.py)

**Bug.** Iteration 8's time integral used √dt weights instead of dt, overstating the bound about 100×.
- Found because iteration 9's independent linear-algebra code disagreed with it.
- The method is validated on ordinary ΛCDM: it gives a late-ISW D₂ = 74 μK², the right order for the known
  ~10% contribution.
- The validation also caught a dropped early boundary term, which matters for matter (whose potential is nonzero early)
  but is small for the lumps. It is now included everywhere.

**Corrected results.**
- One layer adds D₂ ≈ 3.5×10¹⁰ μK², against 1022 observed. The power falls steeply with ℓ.
- CMB bound (2σ, ℓ = 2–30): **N_eff ≥ 2.9×10⁷** (plain Poisson), **1.1×10⁷** (super-horizon suppressed), and
  **3.2×10⁷** including matter falling into the lumps (+10% in χ²).
- Iteration 6's rough 2.6×10¹⁰ was far too strict. The lumps decorrelate in about 1/(6H), so their ISW imprint largely
  cancels along each line of sight.
- **Window now about 10⁷ ≲ N ≲ 1.2×10¹⁵: open and wider.**

**Galaxy clustering does not tighten it.** At N_eff = 3×10⁷ the induced matter power today is 7×10⁻³ of ΛCDM at
k = 0.001 h/Mpc, 10⁻³ at 0.002, and < 10⁻⁵ above 0.006. That is far below survey errors (10–50% at those scales).
The large-angle CMB stays the strongest test.

**Lessons.**
- Cross-checking with independent code caught a 100× error. Keep doing it.
- The stack needs at least ~10 million layers (not billions). With gaps that leak, multiply by (1 + 2r²/(1 − r²)).

## Iteration 10: can one small coupling explain both the height and the featherweight wobble? (iter10_naturalness.py/.txt)

**Model.**
- A stretch field Φ per layer, with tension τ = Φ̇ and energy ½τ².
- A wobble h along a free direction, linked by L = gΦh. This is the only way to get the direct push τ̇ ∝ h that
  iterations 2–4 need.
- Height: ρ = N g² T_s/(24πH), with T_s capped by ΔN_eff.
- Stability and featherweight: g ≤ m_Φ m_h ≤ (3H₀)² = 1.9×10⁻⁶⁵ eV².

**Result.** The height needs g ≥ 3.5×10⁻²⁵ eV² even at the most favourable N = 1.2×10¹⁵: about **10⁴⁰ times too strong**.
With that g the system is unstable within milliseconds.

**Every route checked fails:**
- coupling to a pinned (heavy) wobble: no zero-frequency noise;
- derivative coupling: the push becomes ḣ, which has no zero-frequency noise;
- potential instead of motion energy: smaller, not larger;
- a large "stiffness": equivalent to rescaling, so no change.

**Lessons.**
- The hope that height and featherweight merge into one "technically natural" small number FAILS in the minimal model.
  The stochastic-jostling mechanism, made fully quantitative, cannot give the observed height without wrecking the
  wobble's lightness.
- This is the hardest hit the microscopic toy has taken. The jostling picture explains the SHAPE of dark energy (the
  1/ȧ law, the memory, the smoothing by layers), but not its SIZE.
- Claim 1 itself is untouched: it is a law fitted to data, independent of this microscopic story.
- **Open escapes (each a new assumption):**
  - a nonlinear link between the wobbles and the tension energy that is not a mass mixing;
  - the noise reaching the tension through gravity or the time stamp rather than through a direct coupling.
