# Pre-registration, iteration 52 (committed BEFORE running): the curvature-squared coefficients a and b
Averaged action (E15):  S = (1/16 pi G) INT sqrt(g) [ R + eps^2 (a R^2 + b Ric_ij Ric^ij) + ... ]  (3-D slice, our reader:
R_eps = (30/eps^2) x average over +/- three axes of the Ollivier curvature between uniform eps-balls at separation eps/2).
Method: on HOMOGENEOUS spaces derivative terms vanish, so R_eps = R + eps^2 (a R^2 + b Ric^2) + O(eps^3) exactly per point.
Two spaces give two equations (in 3-D, a and b span all curvature-squared terms since Riemann is fixed by Ricci):
   unit 3-sphere S^3:       R = 6, Ric = 2 g        -> R^2 = 36, Ric^2 = 12
   unit S^2 x line:         R = 2, Ric = diag(1,1,0) -> R^2 = 4,  Ric^2 = 2   (average over the three axes as in the reader)
Ollivier curvature computed by exact optimal transport between finely discretised balls (same point pattern carried by the
isometry, so flat space gives exactly 0). Expectations:
 C1 flat check: kappa = 0 to < 1e-6.
 C2 leading term: kappa/eps^2 -> Ric(v,v)/10 within 1% (Ollivier's theorem, already confirmed on the random grid at ~10%).
 C3 a and b come out as definite O(0.01-0.1) numbers with fit uncertainty < 20% (no prior value; sign not predicted).
 C4 consequence: real-universe corrections stay below |a|,|b| x 2e-61 (unmeasurable); the point is completeness of the equation.

## Outcome (after running; see iter52_first_runs.txt, iter52_curvature_squared_m*.txt, iter52b_convergence.txt)
- C1 PASS (flat kappa ~ 1e-16). C2 PASS after one recorded revision (hard-edged lattice balls had the wrong second moment;
  fractional-volume cells + a free overall factor in the fit): leading term converges to Ric/10 (factor 1.016 -> 1.006 as m = 5 -> 9).
- C3 FAIL: the eps-dependence that defines a and b drifts with resolution and changes sign (S3: -0.0136 ... +0.0016 between
  eps = 0.2 and 0.6 for m = 5 ... 13); finer runs exceed memory. a and b are NOT determined. Ollivier's theorem only bounds
  the next term as O(eps^3) for uniform balls, so even the FORM of the correction (curvature-squared at eps^2, or a
  non-analytic eps^1 |Ric|^(3/2) piece) is open. E15's 'a R^2 + b Ric^2' is therefore an assumption, now labelled as such.
- C4 unaffected: even an eps^1 correction would be ~ eps/L ~ 5e-31 at the most curved tested places; any form is unmeasurable.
