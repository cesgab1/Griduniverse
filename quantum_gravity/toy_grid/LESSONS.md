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
- **First testable prediction of the quantum-grid toy:** extra radiation ΔN_eff ≈ 0.027–0.057. That is allowed today; the
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
