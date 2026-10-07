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
