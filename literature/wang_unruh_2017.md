# Wang & Unruh 2017 (PRD 95, 103504; arXiv:1703.00543) and Carlip 2019 "Hiding the cosmological constant"

Claim: the huge vacuum energy does gravitate, but makes every point of space oscillate (expand/contract) out of phase;
the net, averaged expansion is tiny: H ~ Lam exp(-beta sqrt(G) Lam), Lam = high-energy cutoff.

Size check (wang_unruh_scale.py): matching today's H0 needs a cutoff of ~70-300 x the Planck energy, and a 1% change in
that cutoff changes dark energy by a factor ~18. The 10^-120 tuning becomes a tuning of a number inside an exponential.
The averaged expansion is de Sitter-like (behaves as w = -1, a constant).

Criticism (arXiv:1801.00138, PRD 97, 068301): the result relies on a cutoff that breaks Lorentz invariance; with
physically sensible regularisation the sign needed for the resonance flips and the mechanism fails. A comment on Carlip's
version was also published (PRL 125, 089001).

Verdict for Griduniverse: supports "a huge microscopic energy can be hidden by cancelling fluctuations" (same idea as
Blitz-Majid), but (1) does not explain the size without tuning, (2) predicts a constant (w = -1), whereas our fits prefer
the dynamic beta = 1/2 law, (3) is disputed. Not adopted. Our scale-setting candidate remains the sqrt(N) / everpresent
idea; note the raw everpresent-Lambda random walk failed BAO+SN badly (cosmology_fits/sorkin/sorkin_results.txt) while
its smoothed beta = 1/2 version beat LCDM (Delta chi2 = -6.7 background-only).
