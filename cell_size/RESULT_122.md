# RESULT 122: grid cell size from GRB photon timing (no Planck input)

Pre-registration: PREREG_122.md (committed before running). Code: iter122_cellsize.py. Numbers: RESULT_122_numbers.txt.

## Checks
- C1 reproduction: our formula gives 1.190 for Fermi row (a) (paper 1.19). PASS.
- C2 independent rebuild from the single 28 GeV photon: done (table below).

## Results (upper limits on cell size, metres; only hbar*c measured, H0 and distances used)
| Data | Linear graininess | Regular grid, diagonal (headline) | Regular grid, axis |
|---|---|---|---|
| GRB 090510, (a) any emission start | 1.4e-35 | 4.0e-26 | 2.3e-26 |
| GRB 090510, (d) >1 GeV start | 1.6e-36 | 1.4e-26 | 7.8e-27 |
| GRB 221009A (LHAASO, 95%) | **2.0e-36** | **1.7e-27** | 9.9e-28 |

- Cosmology spread on the headline: 1.65e-27 to 1.73e-27 m (negligible).
- In seconds: a regular grid with 1e-27 m cells would make 7 TeV photons from GRB 221009A arrive 4 s after 0.2 TeV
  photons. 1e-26 m would make it 404 s; the data allow at most ~12 s.

## Expectations vs outcome
- E1 (regular grid < ~1e-27 m from LHAASO; Fermi ~1e-26 to 1e-25): HIT (1.7e-27; Fermi 1.4e-26 to 4.0e-26).
- E2 (linear < ~1e-35 Fermi, ~2e-36 LHAASO): HIT (1.4e-35, 2.0e-36).

## What it means
- A regular "pixel" grid: cells must be smaller than ~2e-27 m, about a trillion times smaller than a proton. This
  is a measured limit with no Planck length used. It is still ~1e8 times larger than the Planck length, so a regular
  grid is NOT pushed to Planck size by these data.
- Linear graininess (light slowed in proportion to its energy): must be < 2e-36 m, below the Planck length (0.12 x).
  Any grid of that type with Planck-size cells is already excluded.
- The current Grid Universe grid is random and Lorentz-friendly (M3): it predicts zero average delay, so photon
  timing cannot size it. Sizing it needs a different signal (spread/blur of arrival times, or black-hole echoes).
- No lower limit on cell size exists from any measurement.

## Corrections logged
- lorentz/lorentz_check.py quoted LHAASO 1.47e20 / 1.2e12 GeV from memory; the paper's ML numbers fetched today are
  1.0e20 / 6.9e11 GeV. The grid bound there (5.7e-28 m, axis factor, older number) becomes 9.9e-28 m (axis) /
  1.7e-27 m (diagonal). Conclusions there unchanged (10-micron fluid channels still excluded by ~22 orders).

## Limits of this test
- Uses published delay limits and LHAASO's own fit, not raw photons (raw LHAASO photons are not public; the Fermi
  archive is not reachable from this workspace).
- Assumes the grid affects only propagation, not emission timing at the source (intrinsic lags could cancel a delay;
  two bursts at different distances agreeing makes that less likely).
