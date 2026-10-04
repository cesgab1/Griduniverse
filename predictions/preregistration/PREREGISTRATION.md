# Grid Universe: pre-registered predictions (frozen 3 October 2026)

These predictions are written down **before** DESI DR3, Euclid, the full SPT-3G survey and Simons Observatory release
their results. The git commit that adds this file timestamps it. After the new data arrive, nothing below may be
re-tuned. Any later change of β, κ, the counting or the reading counts as a new model, not a rescue.

Numbers come from `make_predictions.py` → `predictions.json` / `predictions_table.txt`, with the figure in
`prediction_wz.png`. All parameters are from stored fits to DESI DR2 BAO + Planck 2018 distance priors + one supernova
set. Pantheon+, DES-Dovekie and Union3 agree to ±0.002 in w.

## P1. When dark energy crosses w = −1 (the main test)

| | Claim 1 (instant Tension-Rate Law, β = ½) | Memory version (κ = 3 from Hubble friction, quadratic) | Λ |
|---|---|---|---|
| crossing redshift | **0.68** (= acceleration onset) | **0.46** | never |
| w today | −0.911 | −0.897 | −1 |
| bin 0.1–0.4 | −0.947 | −0.954 | −1 |
| bin 0.4–0.6 | −0.980 | −1.007 | −1 |
| bin 0.6–0.8 | −1.001 | −1.042 | −1 |
| bin 0.8–1.1 | −1.022 | −1.075 | −1 |
| bin 1.1–1.6 | −1.045 | −1.111 | −1 |
| bin 1.6–2.1 | −1.061 | −1.136 | −1 |

- Both versions predict w > −1 at low redshift and w < −1 at high redshift: a single crossing, at the stated redshift.
- The two versions differ most in the 0.6–1.6 bins (by 0.04–0.07). Separating them needs σ(w_bin) ≈ 0.015–0.02.
- BAO distances: D_H/r_d differs from Λ by only ~0.3–0.5% at z < 1 (all values are in predictions.json). The binned w,
  or w0/wa, is the sharper test.

**Kill criteria**
- **Both versions die** if binned w(z) with σ ≲ 0.03 shows no crossing: w ≥ −1 everywhere, or w consistent with −1
  everywhere at > 3σ against the pattern above.
- **Claim 1 dies** if the crossing lies outside z = 0.5–0.9.
- **The memory version dies** if the crossing lies outside z = 0.3–0.6, or if the 0.8–1.1 bin is above −1.04 at > 3σ.

## P2. Claim 2: the electron was heavier when the CMB was released

- Prediction: m_e(recombination)/m_e(today) = **1.004–1.011**. The current measurement is 1.0099 ± 0.0049 (Planck +
  ACT DR6 + SPT-3G D1 + DESI DR2 + DES).
- Equivalent "slower early light" reading: α_rec/α_0 = 1.0031 ± 0.0013. This is partly the same signal, not independent.
- Switch to today's mass at z ≈ 100–200. 21-cm signature: a step of at most 0.7 mK at 7–14 MHz. Not reachable soon.
- **Kill criterion:** SPT-3G full survey or Simons Observatory finds m_e(rec)/m_e(0) consistent with 1.000 at
  σ ≤ 0.002, i.e. excluding 1.004 at > 2σ.

## P3. A little extra radiation (weak test)

- 0 < ΔN_eff < 0.107 (the dark radiation that jostles the layers).
- This is a weak test, because the amount is not fixed. Only a firm ΔN_eff ≲ 0.01 would hurt it, and that is beyond
  planned surveys.

## P4. Other stated consequences (for the record)

- Growth: σ8 about 0.6% below Λ with the same early universe.
- Hubble constant: 67.5 (Claim 1 alone, this pipeline), rising to ~69.5 with Claim 2. **The Hubble tension with SH0ES
  (73) is NOT solved.**
- Large-angle CMB: no excess required. The layer stack needs N_eff ≥ 2.9×10⁷ (real Planck low-ℓ likelihood).
- No energy-dependent photon delay down to the cell-size window (≤ 5.7×10⁻²⁸ m). Gravitational waves travel at c.
- Far future: no recollapse. With the mildly preferred negative baseline (f > 1), the universe coasts (aH → constant)
  instead of approaching de Sitter. Not testable soon.

## What does NOT count as success later

- Refitting β, κ, the jostled fraction f or the layer count after DR3.
- Adding a new ingredient without a new prediction (that is an accommodation; see the self-audit in IDEAS_LEDGER.md).
- Claiming the size of dark energy: in this model it is a free normalisation, not a prediction.
