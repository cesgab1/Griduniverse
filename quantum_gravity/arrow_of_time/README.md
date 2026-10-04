# Arrow of time: Penrose's gravitational entropy and the Weyl curvature hypothesis, tested on the grid (iterations 16–19)

**16. Does dark energy's turnover track gravitational-entropy production?** (iter16) **NO.**
Peak redshifts of the standard measures:

| measure | peak / turnover redshift |
|---|---|
| Weyl/Ricci ratio (linear) | 0.79 |
| galaxy halos | 1.39 |
| cluster halos | ~0 |
| black-hole entropy (from the measured star-formation history) | 1.19 |
| horizon entropy | 0.65 |

The measures disagree. The horizon value is circular: d(1/H²)/dt = 2(1 + q)/H depends only on the expansion history.
There is no independent link to Claim 1 (crossing at 0.68) or the toy (0.46). **Idea dropped.**

**17. What Weyl curvature means on the grid.** (iter17)
- Random cell irregularity (spread 0.20 per cell) averages away as 1/√n (0.009 over 512 cells), so it is not Weyl
  curvature. A random mosaic is compatible with Penrose's smooth start.
- Stretched cells carry a coherent strain linearly (shift = 1.25e).
- Cells REBUILT by the mosaic rule erase it almost completely (retained ≈ 0).
- **Lesson:** a grid that renews its cells cannot store tidal distortion in cell shapes; it must live in the link lengths
  (deficit angles). This matches the Regge result (quantum_gravity/regge_gauge.py): moving or renewing nodes is a gauge
  change, so the curvature survives renewal. Gravitational entropy = coherent link-length distortion.

**18. Generalised second law: horizon entropy (∝ 1/H²) must never decrease.** (iter18) **PASSES** for Claim 1
(despite its phantom phase at a < 0.59: matter keeps the total ρ + p ≥ 0.013), for the toy with f = 1 and f = 1.5, and
for ΛCDM.
- **Distinctive difference:** in ΛCDM, H → 0.83 H₀, so horizon entropy saturates (a finite "heat death" ceiling). In ours
  H keeps falling (0.28 / 0.20 / 0.03 H₀ by a = 100), so horizon entropy grows without bound.

**19. The jostled tension in a CONTRACTING universe.** (iter19)
- Hubble friction becomes anti-friction, and the tension energy grows as a⁻⁶ (stiff, w = +1). That is exactly the rate of
  anisotropic shear (Weyl-type energy, the BKL chaos), so it cannot smooth a collapse.
- **Lesson:** the dark-energy mechanism is intrinsically time-asymmetric. It calms the tension only while space expands.
  A collapse would end with high Weyl curvature, unlike our smooth start: Penrose's Bang/Crunch asymmetry, built in.

**Cross-applied.**
- 19 supports dropping the bounce: contraction amplifies the tension.
- 14 + 18: no recollapse, and horizon entropy grows forever.
- 17 + graviton tests: gravity and its entropy live in link lengths.
- 16 closes the timing idea.
- **Arrow-of-time picture:** smooth start (no coherent link distortion) → clumping builds coherent distortion
  (gravitational entropy) → expansion makes the tension's memory dissipative → horizon entropy rises without limit.
