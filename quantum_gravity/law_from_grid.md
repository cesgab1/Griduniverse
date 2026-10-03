# Does the grid (CDT/Hořava-type quantum gravity) give the Tension-Rate Law? Due-diligence calculations, Oct 2026

Question: can Claim 1 (ρ_DE ∝ ȧ^(−1/2)) be derived from a quantum-gravity grid with a preferred time slicing
(Causal Dynamical Triangulations / Hořava gravity, whose low-energy limit is our Khronon base)? Seven calculations.

## 1. Can it come from an action with a preferred time slicing? Yes, but not as an ordinary low-energy term (theory/minisuperspace.py)
Minisuperspace, flat FRW, lapse N, scale factor a. The expansion of the time slices is K = 3H. The grid ingredient is the
physical link length ℓ ∝ a^s.

    S = ∫ dt N a³ [ −3M²H² + P(H, a) ] − matter

- Varying N gives the Friedmann constraint, with ρ_DE = H ∂P/∂H − P. Choosing P = −(2/3) C (H a^s)^(−1/2) gives
  ρ_DE = C (H a^s)^(−1/2). This is unique up to a total derivative. For s = 1 it is Claim 1.
- The Noether (Bianchi) identity holds: the constraint is preserved by the a-equation, so the theory is consistent.
  w = −1 − (q + 1 − s)/6.
- An ordinary low-energy expansion P = Σ cₙ Hⁿ gives only ρ ∝ H⁰, H², H³, … It can never give H^(−1/2).
  **Hořava/CDT as a standard effective theory does not produce the law.** It needs:
  - (i) a non-analytic dependence on the expansion rate, which grows as expansion slows (an infrared, collective
    effect, as in √N counting statistics, not a local term);
  - (ii) a physical length that tracks stretching: the link length, or equivalently a conserved cell-number density.
- No extra degree of freedom is added: P depends only on the slicing's expansion and on ℓ. Theories of this type have only
  the 2 graviton polarisations, and they cross w = −1 without ghosts. These are VCDM / type-II minimally modified gravity
  (De Felice, Mukohyama et al.; arXiv:2004.12549, 2011.04188, 2508.03784) and "extended cuscuton" (Iyonaga, Takahashi & Kobayashi, JCAP 12 (2018) 002). In VCDM, perturbations
  equal ΛCDM's apart from the momentum constraint. This gives Claim 1 a ghost-free home, which a single scalar
  (quintessence) cannot provide. BORROWED framework.
- Perturbation treatment: checked in section 6 below. It does not matter.

## 2. Audit of the microscopic counting: it decides whether cells are being added today (stretch_scan/)
Each link is crossed about N = c/(2 H ℓ_phys) times per stretch time, and the net tension is ∝ √N.
- Links stretch with space (ℓ_phys ∝ a, s = 1) → N ∝ 1/ȧ → ρ ∝ ȧ^(−1/2), which is Claim 1.
- Cells added, link length fixed (s = 0) → N ∝ 1/H → ρ ∝ H^(−1/2), a different law that never crosses w = −1 at q = 0.

Requirement found: the tension must be INTENSIVE (vacuum-like: each crossing changes the energy density, not the energy
per link). If tension were energy per link, it would dilute (d ln ρ/d ln a = q/2 − 3), and that is excluded outright.

General law: d ln ρ_DE/d ln a = (q + 1 − s)/2. Fit to DESI DR2 BAO + Planck distance priors + each supernova set,
with no dark-energy parameter at fixed s (stretch_scan/stretch_scan.txt):

| supernovae | best s (simple q) | best s (exact q) | s = 0 (cells added now) excluded at |
|---|---|---|---|
| Pantheon+ | 0.91 ± 0.15 | 0.93 ± 0.16 | 6.0σ / 5.4σ |
| DES-Dovekie | 0.94 ± 0.14 | 0.96 ± 0.15 | 6.7σ / 6.1σ |
| Union3 | 1.01 ± 0.18 | 1.03 ± 0.20 | 5.4σ / 5.0σ |

**The data measure the stretch exponent at s ≈ 1 (0.9–1.0 ± 0.15):** the length the counting runs over stretches with
space. CORRECTION (section 7): this length cannot be the Planck cell itself. An earlier version of this file said "cells
are not being added today"; that was wrong.

## 3. Size: the counting fixes the shape of dark energy, not its amount (theory/magnitude.py)
With Planck-length links, N ≈ 4×10⁶⁰ and √N ≈ 2×10³⁰. The observed ρ_DE = 1.1×10⁻¹²³ Planck densities then needs a kick
of 6×10⁻¹⁵⁴ Planck densities per crossing (3×10⁻¹⁵⁰ for the largest allowed cells). Nothing in the model sets that
number. The cosmological-constant size problem remains as unsolved as in every other model.

## 4. Where the law must stop
ρ ∝ ȧ^(−1/2) diverges where ȧ = 0, at the Bounce. There the stretch time is infinite and the counting is meaningless, so
the law is a late-universe law that must be cut off near the Bounce. In the future ȧ grows (a ∝ t⁵), so ρ_DE falls and
there is no problem.

## 5. Is the √N exponent (β = 1/2) itself picked by the data? Joint fit of β and s (joint_scan/joint_scan.txt)
Both β and s left free in d ln ρ_DE/d ln a = β (q + 1 − s). Grid β 0.25–3, s 0.25–1.75 (simplified q):

| supernovae | best β, s | Δχ² vs Λ (2 params) | 68% range of s | 68% range of β | Claim 1 (β ½, s 1) above best |
|---|---|---|---|---|---|
| Pantheon+ | 1.08, 0.83 | −7.5 (w0wa −7.1) | 0.72–1.02 | 0.42–1.79 | +2.1 (inside 68%) |
| DES-Dovekie | 1.29, 0.82 | −10.3 (w0wa −10.3) | 0.74–0.97 | 0.59–2.02 | +3.1 (inside 95%) |
| Union3 | 2.19, 0.86 | −13.4 (w0wa −13.2) | 0.80–0.97 | 1.09–3.0 (edge) | +6.6 (just outside 95%) |

Lessons:
- **s is robust.** Whatever β is, the stretch exponent comes out between 0.72 and 1.06. The "stretching length" statement
  does not depend on the √N assumption.
- **β is weakly fixed.** The √N value ½ is allowed by two of three supernova sets. Union3 prefers a stronger dependence
  (β ≈ 1–3), and the 2-parameter family then fits exactly as well as DESI's w0wa.
- So the zero-parameter Claim 1 survives, but the data do not single out β = ½. If future data move β toward 1, that
  points to coherent rather than independent crossings (β = 1 is the "fully coherent" case in beta_derivation.py).

## 6. Does the treatment of dark-energy fluctuations matter? No (perturbations/ppf_vs_smooth.txt)
Same Claim-1 background in CAMB, computed two ways: PPF (used in all our CMB fits) and dark-energy fluctuations switched
off entirely. The second is a proxy for the ghost-free VCDM home, where dark energy has no fluctuations of its own.
- CMB TT: cosmic-variance χ² = 0.008. EE: 0.03. Even a perfect full-sky experiment could not tell them apart.
- σ8 differs by 0.09%.
- The lensing power differs by at most 0.28% over the measured range, against a measurement precision of about 2%.

The CMB results (Claim 2, slower early light) do not depend on this choice. Caveat: switching fluctuations off is a
proxy, not VCDM's exact equations.

## 7. Correction: the stretching length is NOT the grid cell (theory/cells_vs_constants.txt)
In the model, the cell size is the Planck length (c = ℓ_P/t_P; stiffness c⁴/G). If the cells stretched as a^s_cell,
then ℓ_P² = ħG/c³ would grow. Measured limits on that:

| test | limit on s_cell |
|---|---|
| lunar laser ranging (G today) | < 1.6×10⁻³ |
| nucleosynthesis (G in the first minutes) | < 2.4×10⁻³ |
| our CMB constant test (c via α) | < 8×10⁻⁴ |

The dark-energy law needs s ≈ 1. Therefore:
- **The Planck cells do NOT stretch: the grid IS adding cells today.** This agrees with the Stretch limit (E6), and it
  reverses what I said in the previous summary.
- **The counting in Claim 1 must run over a COMOVING length:** something carried along by the expansion, with fixed
  Planck cells forming and dividing underneath it. The memory of a comoving length is what the data detect (s ≈ 1).
- **Candidates (OPEN, untested):**
  - the separation of matter or Ocean particles (comoving by definition);
  - a coarse comoving structure of the mosaic that survives cell division ("parent cells");
  - the comoving Hubble-scale correlation length of the exchanges.

  Each must keep the tension intensive.
- **Consequence for the size problem (section 3):** with a comoving length L, N = c/(2ȧL). L is unknown, so the per-crossing
  kick needed is unknown too. The size problem is unchanged, only relabelled.

## Verification
An independent agent re-derived items 1 and 2 from scratch.
- Its own Union3 fit agrees with ours to about 0.2 in Δχ².
- It found two errors, both fixed and rerun:
  - (a) the exact-q option was self-consistent only at s = 1;
  - (b) fit_law.py integrated the sound horizon from the first grid point above z* rather than from z*. This made the
    CMB χ² surface jagged; l_A was then recalibrated to CAMB.
- Effect on Claim 1: ≤ 0.1 in Δχ² (now −5.4 / −7.2 / −6.8; exact −5.4 / −6.9 / −6.3).
- Effect on the s = 0 exclusion: it drops from 5.7–7.2σ (before the fixes) to 5.0–6.7σ.

## Bottom line
- The existing quantum-gravity frameworks with a time slicing (CDT, Hořava) do NOT produce the law as a standard
  low-energy term.
- The law fits into their ghost-free extensions (VCDM / extended cuscuton) only if two ingredients are added: a
  collective √N effect and a physical link length.
- Through the law, the data say that the length behind the counting stretches with space (s ≈ 0.7–1.05, robust to β).
  A fixed length (s = 0) is excluded at ≥ 5σ.
- Measured constancy of G and c says the Planck cells cannot be that length (s_cell < 10⁻³). The cells are being added,
  and the counting runs over a comoving length that is not yet identified. Identifying it is now the central open problem
  for deriving Claim 1.
- β = ½ (√N) is allowed but not singled out. The 2-parameter family fits as well as w0wa.
- The perturbation treatment does not matter.
- The amount of dark energy is not explained.
