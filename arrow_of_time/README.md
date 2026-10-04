
## Iteration 30: the grid's 'watermark' (primordial ripples). Pre-registered at d793c64, then compared
Data: n_s = 0.974 +/- 0.003 (ACT DR6 P-ACT-LB), Planck 2018 0.965 +/- 0.004; running dn_s/dlnk = 0.0062 +/- 0.0052 (P-ACT-LB).
- W1, cell jitter (local random nudges, thermal noise): n_s = 4. EXCLUDED overwhelmingly. The watermark cannot be noise
  made locally after the start; it needs correlations beyond the horizon.
- W2, Horava scaling with z = 3 exactly (BORROWED, Mukohyama 2009): perfectly scale invariant, n_s = 1. EXCLUDED at 8.7 sigma
  (ACT) / 8.8 sigma (Planck).
- W2 with z = 3 - eta: matching the tilt needs eta = 0.013 (ACT) to 0.017 (Planck), i.e. z = 2.983-2.987. One fitted number
  for one measured number = accommodation. Its extra prediction, no running (|alpha_s| < 1e-3), agrees with data (zero within
  1.2 sigma): a weak PASS. Amplitude not predicted (it needs a curvaton-type conversion with a free efficiency).
- Verdict: the grid's preferred time slicing CAN print a nearly scale-invariant watermark without inflation (BORROWED), but
  only with a small fitted departure from z = 3. Not a success; it is a consistent option. A derivation of eta from the grid
  (e.g. an anomalous scaling from the layers) would turn it into a prediction.

## Iteration 51: where could eta come from? (pre-registered: PREREG_51.md)  -> T1-T4 as expected
- Frozen-ripple spectra computed numerically (iter51_eta_sources.*). Formulas confirmed within 1.5-4%.
- Route 1, 'the z = 3 phase is ending' (crossover omega^2 = p^2 + p^6/M^4): the tilt comes with a large NEGATIVE running,
  alpha_s = (4/3)[eps/(1 - eps/3)](n_s - 1): radiation era -0.22 (43 sigma from ACT's +0.006 +/- 0.005); even a nearly
  inflating era (eps = 0.5) is 5 sigma off. EXCLUDED.
- Route 2, a constant anomalous exponent z = 3 - eta: n_s - 1 = -(eta/3) eps/(1 - eps/3), no running. SURVIVES.
  Radiation era: eta = 0.013 (ACT) - 0.018 (Planck); the frame field's small-scale spectral dimension must be 2.004-2.006,
  not exactly 2. A frame-free causal-set grid reduces to exactly 2 (BORROWED) -> it cannot be the source (consistent with the
  two-sector picture: the watermark belongs to the frame field).
- eta NOT derived (needs a loop calculation in the frame field's z = 3 theory). Narrowing, not a prediction. The running
  prediction (alpha_s ~ 0) is now a sharper test: a future |alpha_s| > ~0.01 would exclude route 2 as well.
