# PREREG 122: how big can a grid cell be? (photon timing, no Planck input)

Written and committed BEFORE any code is run.

## Question (Coalesce)
Measure the size of one bit (one grid cell) of space directly, in metres, without using the Planck length.

## Idea in one line
If space is a grid, short-wavelength (high-energy) light should feel the grid and travel slightly slower.
Gamma-ray bursts emit all energies at nearly the same instant from billions of light-years away.
So: no measured lag between high- and low-energy photons -> an upper limit on the cell size.

## Data (published numbers, quoted, not re-fitted)
A. GRB 090510, Fermi (Abdo et al. 2009, Nature 462, 331; Table 2):
   z = 0.900 (1-sigma low end of 0.903 +/- 0.003); cosmology [Omega_L, Omega_M, h] = [0.73, 0.27, 0.71];
   highest-energy photon 28.0 GeV (1-sigma low end of ~31 GeV) at T0 + 0.829 s.
   Delay limits by assumed emission start: (a) any <1 MeV emission: |dt| < 859 ms (paper: > 1.19 M_Planck)
   (b) main <1 MeV emission: < 299 ms; (c) main >0.1 GeV: < 181 ms; (d) >1 GeV: < 99 ms;
   (g) lag analysis of >1 GeV spikes: |dt/dE| < 30 ms/GeV.
B. GRB 221009A, LHAASO (Cao et al. 2024, PRL; arXiv:2402.06009): z = 0.151; H0 = 67.36, Omega_M = 0.315;
   64,000 photons 0.2-7 TeV. 95% limits (ML method, subluminal): E_QG,1 > 1.0e20 GeV, E_QG,2 > 6.9e11 GeV.
   (Note: an older repo note, lorentz/lorentz_check.py, quoted 1.47e20 / 1.2e12 from memory; this test uses the
   numbers fetched today, which are the more conservative ones.)

## Delay formula (Jacob-Piran, standard)
dt = (1+n)/2 * (E_h^n - E_l^n) / E_QG^n * (1/H0) * Integral_0^z (1+z')^n / sqrt(Om (1+z')^3 + OL) dz'

## Mapping from cell size l (metres) to E_QG -- uses only hbar*c (measured), no Planck length
M1 regular cubic grid (the natural "pixel" grid): group speed v = c [1 - (E l / hbar c)^2 * f / 8],
   f = 1 for light travelling along a grid axis, f = 1/3 along a body diagonal (weakest effect).
   -> E_QG,2 = sqrt(12/f) * hbar c / l, i.e. l = sqrt(12/f) * hbar c / E_QG,2.
M2 linear graininess (no simple grid gives this; kept for completeness): v = c (1 - E l / hbar c) -> l = hbar c / E_QG,1.
M3 random, Lorentz-friendly grid (the current Grid Universe choice, ledger iterations ~1409-1415 and PREREG_53 B3):
   predicts ZERO average delay -> this test cannot size it. Stated, not run.

## Free choices (counted): 4
1. Direction factor f in M1 (headline = diagonal f = 1/3, the most conservative; axis reported too).
2. Coefficient 1 in M2.
3. Cosmology inside the integral (headline = each paper's own; spread shown for H0 67-73, Om 0.27-0.32).
4. Which LHAASO number (headline = ML subluminal, the more conservative of those quoted).
No significance search -> no look-elsewhere correction needed; results are upper limits.

## Pre-registered checks and expectations
C1 (reproduction, must pass or STOP): our formula with row (a) inputs must give E_QG,1 = 1.19 x 1.22e19 GeV within 5%.
C2 (independent rebuild): from the single 28 GeV photon alone (rows a-d, g) compute l limits for M1 and M2 ourselves.
E1 (expectation, written before running): M1 (regular grid) l < ~1e-27 m, set by LHAASO; Fermi alone ~1e-26 to 1e-25 m.
E2: M2 (linear) l < ~1e-35 m (Fermi row a) and ~2e-36 m (LHAASO).
E3: give the delay in seconds for a few trial cell sizes, so the answer can be checked by eye.
Planck length appears ONLY as a comparison in the final printed line, never in a calculation.
