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
