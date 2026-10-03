# Does the grid (CDT/Hořava-type quantum gravity) give the Tension-Rate Law? Due-diligence calculations, Oct 2026

Question: can Claim 1 (ρ_DE ∝ ȧ^(−1/2)) be derived from a quantum-gravity grid with a preferred time slicing
(Causal Dynamical Triangulations / Hořava gravity, whose low-energy limit is our Khronon base)? Four calculations.

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
- Not done yet: our CAMB fits use PPF dark-energy perturbations, not VCDM's. The difference is expected to be small
  (mainly the large-scale ISW), but it is unchecked.

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

**The data measure the link stretch exponent at s ≈ 1 (0.9–1.0 ± 0.15).** Within this counting, links are stretching
with space today, not being replaced by new cells. Combined with the Stretch limit (E6), cell creation must have stopped
at some z_f between ~5 and 3.5×10⁷, after which the number of cells is fixed. This is conditional on the √N exponent
(β = 1/2) and on this counting picture.

## 3. Size: the counting fixes the shape of dark energy, not its amount (theory/magnitude.py)
With Planck-length links, N ≈ 4×10⁶⁰ and √N ≈ 2×10³⁰. The observed ρ_DE = 1.1×10⁻¹²³ Planck densities then needs a kick
of 6×10⁻¹⁵⁴ Planck densities per crossing (3×10⁻¹⁵⁰ for the largest allowed cells). Nothing in the model sets that
number. The cosmological-constant size problem remains as unsolved as in every other model.

## 4. Where the law must stop
ρ ∝ ȧ^(−1/2) diverges where ȧ = 0, at the Bounce. There the stretch time is infinite and the counting is meaningless, so
the law is a late-universe law that must be cut off near the Bounce. In the future ȧ grows (a ∝ t⁵), so ρ_DE falls and
there is no problem.

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
- The data then say, through the law, that the links are stretching today (s ≈ 1) and that cells are not being added
  (s = 0 excluded at ≥ 5σ). This is a new, testable statement about the grid.
- The amount of dark energy is not explained.
