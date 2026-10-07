# RESULT 132: does the extra-pull threshold grow with the expansion rate?
Prereg: PREREG_132.md. Code: iter132_threshold.py. Numbers: RESULT_132_numbers.txt.
Models: C (compression) a_t = a0 * H(z)/H0; K (constant) a_t = a0 (a0 = 1.2 +/- 0.26 e-10, SPARC).

(1) Our own fit, Genzel+2017 six z = 0.85-2.38 galaxies (g_bar 1-9x a0): z<2 a = 4.3e-10 (+/-0.54 dex), z>2
    a = 2.3e-10 (+/-1.1 dex); chi^2 C 0.2, K 1.1 -> no preference. Too few galaxies, too uncertain. (A draft note
    misstated g_bar as 10-30x a0; corrected before results were written.)
(2) Direct published measurement (Ciocan+2026, MUSE-DARK III, A&A 709 L16; 79 galaxies, z 0.33-1.44):
    a(z~1) = 2.38 (+0.12/-0.10, 95%) e-10 m/s^2; linear fit a = (1.0 + 1.59 z) e-10.
    K (no change): excluded at 4.4 sigma (SPARC anchor) to 24 sigma (their own anchor).
    C (tracks H): matches with the SPARC local anchor (+0.6 sigma); with the survey's own extrapolated local value
    (1.0e-10) C is TOO SLOW (9.5 sigma) -- the authors: 'faster than H(z)'.
    Caveats (authors): baryons and total pull from the same model; molecular gas +-0.2 dex; reconciling with no
    evolution would need +0.2-0.45 dex more stellar mass.
(3) Baryonic Tully-Fisher zero point (indirect): Ubler+2017 -0.44 (z 0.9), -0.27 (z 2.3) dex -- fits neither C nor K
    (non-monotonic); Jeanneau+2026 (MUSE-DARK II, same team) 0.00 +/- 0.06 at z~1 -> fits K, rejects C at 4 sigma.
    The studies differ mainly in how gas is counted (convention) -> conflicting, no verdict.
(4) Milgrom 2017 on z~2 curves: ~4 a0 at z~2 'all but excluded'; C predicts 3.2 a0 at z = 2.2 (allowed, marginal).

Project rules: the evolving threshold is significant in ONE survey; not independently confirmed; the same team's
Tully-Fisher analysis finds no evolution; not yet a population across surveys -> PROMISING, NOT ESTABLISHED.
Expectation ('likely inconclusive / weak'): our own data HIT; MISS overall -- a direct, strong (if unconfirmed)
measurement of a rising threshold exists, which I did not anticipate.

LESSON 132: the single direct measurement says the threshold DOES rise with lookback time -- the direction
Coalesce's compression picture predicts, roughly H(z)-sized with the standard local anchor, but faster with the
survey's own anchor. Constant-threshold (standard MOND) is in trouble if it holds; dark-matter halos have no such
threshold to track. Deciders: independent RAR measurements at z~1-2 (other surveys) and a consistent gas convention.
