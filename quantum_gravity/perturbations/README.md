# Iteration 50: does the frame field change how structure grows? (Oct 2026)
Pre-registration: PREREG_50.md. Equations BORROWED from VCDM (De Felice, Mukohyama & Pookkillath, PLB 816, 136201, 2021):
LCDM perturbation equations except the momentum constraint, which picks up
    F(k) = [k^2 - 3 a^2 dH/dt] / [k^2 + (9/2) a^2 rho_m]     (F = 1 when dark energy has w = -1).

Results (iter50_horizon_growth.txt):
- P1 code check: F = 1 on LCDM to 2.5e-4 (finite differences) -- after one RECORDED revision: the first run left radiation in
  dH/dt but not in F's denominator (|F - 1| = 0.075); radiation was removed from this calculation (iter50_first_run_with_radiation.txt).
- P2 PASS: our model vs the same expansion history with ordinary gravity differ by < 1e-5 for k >= 0.01 h/Mpc (all galaxy
  surveys), 0.05% at k = 0.001, 1.2% at k = 3e-4, 15% only at k = 1e-4 (larger than the visible universe's horizon).
- P3 PASS (as a consistency check, NOT a new prediction): late-ISW potential decay on the largest scales changes by -3.6%
  (k = 2e-4) and -0.5% (k = 1e-3) -- far inside cosmic variance (the lowest CMB multipoles scatter by tens of percent).
- For comparison, the expansion-history effect alone (already in our fits) is larger on observed scales: growth 0.5% below
  LCDM at k = 0.01, potential decay 1.7% faster.
Verdict: the frame field's extra effect is confined to near-horizon scales and is unobservable; the all-background fits we ran
were justified. Claim 1's testable content stays in the expansion history (DESI 2027) and the growth it implies.
