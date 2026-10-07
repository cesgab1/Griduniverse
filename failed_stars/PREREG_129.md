# PREREG 129: one planet-building recipe from stars down to failed stars?
Committed before data and code.
Data: per-object disk dust masses (from millimetre flux) and central masses in young regions, stars to brown dwarfs
(and planetary-mass objects if any mm data exist). Same-region comparisons preferred (e.g. Lupus, Cha I, Taurus,
Ophiuchus); Upper Sco (older) analysed separately.
Model: log M_dust = a + b log M_central. Fit STARS ONLY (M > 0.1 Msun); extrapolate to brown dwarfs (0.013-0.08).
Test: median residual of brown dwarfs from the star line (dex), with bootstrap uncertainty, non-detections handled by
(i) detections only and (ii) treating limits as values (bracketing); a BREAK is claimed only if |residual| > 0.3 dex
AND > 2 sigma in both treatments and in >= 2 regions (independent confirmation).
CONVENTION drill: dust mass from flux needs an assumed dust temperature. Two conventions: (A) fixed 20 K;
(B) T = 25 K (L/Lsun)^0.25 (cooler around faint objects). Report both; a break that appears in only one
convention is a convention artefact.
Expectation (written now): star slope b ~ 1.3-2.0; brown dwarfs near or slightly ABOVE the extrapolated line with
fixed 20 K and closer still with (B) -> one recipe; planetary-mass objects lack mm detections -> untested.
Free choices: 3 (temperature convention, limit treatment, region set).

## ADDENDUM (before running, after seeing only the data COVERAGE, not results)
Machine-readable full star catalogues were not reachable (CDS blocked; arXiv tables truncated). Available: L1688
(Testi 2022, 10 objects > 0.1 Msun spanning only 0.11-0.20 Msun; 14 below) and Lupus (Sanchis 2020, 6 stars
0.23-0.86 Msun with masses; 9 objects <= 0.10). Star-only slopes from these few rows are too poorly anchored, so:
headline star line = PUBLISHED star-only slope for that region (Lupus 1.73 +/- 0.25, Sanchis 2020; L1688 1.9 +/- 0.4
for M >= 0.15, Testi 2022), normalisation fitted to our star rows; our own free-slope fit reported as alternative.
Brown dwarfs here = M <= 0.08 Msun (0.1 boundary objects excluded from both groups). Convention B (luminosity-scaled
temperature) only possible for L1688 (Lupus table lacks luminosities). Upper Sco (brown dwarfs only) reported
descriptively. Free choices now 4 (adds: published vs own slope).
