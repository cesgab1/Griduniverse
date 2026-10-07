# RESULT 123: "IT is space" -- the three ways out of the bit-count puzzle

Pre-registration: PREREG_123.md. Code: iter123_three_avenues.py, iter123_t2_sensitivity.py. Numbers: RESULT_123_*.txt.
Sources: SOURCES_123.md.

## T1 LOCKED BITS (bits fill space but mostly move together)
| Signal | Independent bits (b = 1/2) | Locked bits (b = 1/3) |
|---|---|---|
| GRB 090510 spikes <~10 ms (our calc) | l < 9.5e-14 m | l < 5.3e-4 m (half a millimetre!) |
| Holometer (Hogan's locked/shear model only) | -- | scale < 4.0e-36 m (2 sigma) |
| Image blurring (Perlman+2015) -- CONTESTED | ruled out | l < 2.2e-40 m if the method is right |
Also: a grid that randomly scatters light speed in proportion to energy: l < 5.8e-36 m (Fermi, 95%).
Expectation (locked bits hardest to see, weak timing limit): HIT.
Lesson: locking the bits HIDES them. Timing alone allows locked cells up to ~0.5 mm. Only specific locked models
(Hogan's) are sub-Planck excluded; the general idea survives unless the contested blurring method is right.

## T2 SURFACE BITS (a black hole's surface grows in tiles a = alpha*hbar*G/c^3)
Using GW250114 (spin 0.68 +/- 0.01, f220 +/- 2.4%) and GW150914; narrow-line assumption (Foit-Kleban):
| alpha | tile | verdict | robust to tolerance 2-5%? |
|---|---|---|---|
| 4 ln 2 = 2.77 (one bit per tile; Coalesce) | 7.2e-70 m^2 | allowed (untestable: lines denser than errors) | yes |
| 4 ln 3 = 4.39 (Hod) | | excluded at 2.0-2.4%, allowed at >= 2.9% | NO -> marginal |
| 8 ln 2 = 5.55 | | allowed | yes |
| 8 pi = 25.1 (Bekenstein/Maggiore) | | EXCLUDED | yes (2-5%); NOT under line-width criterion |
- About 90% of alpha in 0.5-120 excluded by GW250114 alone; everything below ~3.4-4.4 allowed.
- Alternative criterion (line only needs to fall within the ringing's own width, 15.6%): all alpha < 25.6 allowed,
  including 8 pi. The whole test rests on whether quantum lines are narrow (Foit-Kleban) or broadened (Agullo+).
- Published (pre-GW250114): Laghi+2021 uninformative; GWTC-3 no echoes.
Expectation ("small alpha untestable; only alpha >~ 10-30 testable"): first half HIT; second half MISS --
GW250114 tests down to alpha ~ 4 (better than expected) under narrow lines.
Lesson: Coalesce's one-bit-per-tile is still alive and is the hardest version to test; the large-tile proposal
(8 pi) fails if lines are narrow. Not support -- just survival.

## T3 THE AREA RULE BREAKS AT SMALL SCALES
- Correction on record: the area rule uses G at the HORIZON (km), which gravitational waves test; the small-scale
  extrapolation is inside Hawking's derivation (very short wavelengths near the horizon).
- Newton verified down to 3.9e-5 m; gap to a 1.7e-27 m cell: 22 powers of ten; to Planck size: 30.
- Area never decreases: GW150914 97%; GW250114 4.4 sigma -> classical area rule holds.
- Grainy-medium check: BEC analog black hole had Hawking energy = 0.12 of its grain energy and still gave the
  thermal spectrum with no free parameters; a real 62.7 Msun hole with 1.7e-27 m cells has ratio 7e-34 -- the
  analog was ~2e32 times more grainy (relative) and the law still held. Unruh 1995 found the same numerically.
Expectation (no evidence of breakdown; analog ~1e30 more grainy): HIT.
Lesson: avenue 3 loses support (graininess does not easily break the rule) but is not excluded (no test below 40 um).

## Overall
- Bit counting: none of the three is excluded. Locked bits are least constrained in general; surface bits with one
  bit per tile survive but are untestable now; small-scale breakdown has no evidence.
- Free choices: T1 2 (spike width, distance measure), T2 2 (tolerance, events) -- verdicts shown across choices.
- No detections claimed; no look-elsewhere needed.
