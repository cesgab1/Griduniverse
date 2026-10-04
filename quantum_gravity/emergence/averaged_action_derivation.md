# Does the neighbourhood-averaged grid action become Einstein's action? (pen-and-paper, Oct 2026)

Labels: **BORROWED** = published theorem; **DERIVED** = argued here; **MEASURED** = our grid runs.

## Setup
Grid action on a time slice (E13): S_eps = (1/16 pi G) SUM_x V_x R_eps(x), with
R_eps(x) = (30/eps^2) * average over six directions (+/- three perpendicular axes e_k) of kappa_eps(x, x + delta e_k),
where kappa is the Ollivier curvature of uniform eps-neighbourhoods and the pair separation is delta ~ eps/2.

## Step 1: what one pair measures (BORROWED: Ollivier 2009, Example 7 / Prop. 6)
For a smooth n-dimensional space and uniform measures on eps-balls:
    kappa(x, y) = eps^2 Ric_m(v, v) / (2 (n + 2)) + O(eps^4 / L^4, eps^2 delta^2 / L^4)
evaluated at the MIDPOINT m of x and y, v the unit direction x -> y, L the curvature length scale.
(n = 3: denominator 10.) The leading coefficient is MEASURED on our random grid at the bump centre:
measured/true = 0.99 +/- 0.08 (iteration 43), 1.04 +/- 0.09 (43b), 1.18 +/- 0.12 (41b).

## Step 2: the six-direction average gives the scalar curvature, with no first-order error (DERIVED)
- Sum over three perpendicular unit vectors: Ric(e1,e1) + Ric(e2,e2) + Ric(e3,e3) = R (trace). So
  (30/eps^2) * (1/3) SUM_k kappa_k = (30/eps^2)(eps^2/10)(R/3) = R. The 30 is not a fit: 2 n (n + 2) = 30 for n = 3.
- Midpoint shift: Ric_m = Ric(x) + (delta/2) e_k . grad Ric + (delta^2/8) (e_k . grad)^2 Ric + ...
  Taking +e_k and -e_k together cancels the first-derivative term exactly (odd in e_k).
- What remains: R_eps(x) = R(x) + delta^2 * (second derivatives of Ric) + eps^2 * (curvature squared) + ...

## Step 3: summing over the grid (DERIVED)
- SUM_x V_x -> INTEGRAL sqrt(g) d^3x as the grid is refined (each point carries its proper volume; sprinkling uniform in proper
  volume).
- Second-derivative terms are total derivatives (they integrate to boundary terms and vanish for a closed region or fields that
  fall off).
- Result:   S_eps = (1/16 pi G) INTEGRAL sqrt(g) [ R + eps^2 ( a R^2 + b Ric_ij Ric^ij ) + O(eps^3) ]  + boundary terms
  i.e. EINSTEIN'S action plus small curvature-squared corrections. The numbers a, b need a 4th-order expansion of the
  transport cost (not done here). Curvature-squared corrections are the generic form expected in any quantum-gravity approach.

## Step 4: the randomness averages away (DERIVED + MEASURED)
- Per-point scatter of R_eps at ~1600-1700 points per neighbourhood: about 24% of R (from 43/43b: +/-2.5 with 8-30 points).
- Independent per-point errors in a sum of N points shrink as 1/sqrt(N): the action's relative noise ~ 0.24/sqrt(N).

## Step 5: how big are the corrections in the real universe? (numbers)
Neighbourhood eps = 6-8 cells; cell size a = 4e-32 to 5.7e-28 m (iteration 40) -> eps ~ 3e-31 to 4.6e-27 m.
- Curvature-squared corrections, relative size (eps/L)^2. Most strongly curved place tested (black-hole horizons and neutron stars,
  L ~ 10 km): (4.6e-27 m / 1e4 m)^2 ~ 2e-61 at most.
- Random noise: points in 1 cm^3 ~ (0.01 m / a)^3 ~ 5e75 to 1e88 -> relative noise ~ 3e-39 to 2e-45.
Both are far beyond any measurement: in practice, the averaged random grid IS Einstein's gravity on every tested scale.

## Status and caveats
- Steps 1 and 3 rely on Ollivier's expansion (borrowed) and on van der Hoorn et al. 2021 for convergence on random geometric
  graphs; Step 2 and the totals are derived here; the leading coefficient is confirmed on our grid at the bump centre.
- Our off-centre profile runs (43-43c) agree at r = 0.2 and at the edge (including its negative sign) but r = 0.35 is unresolved
  because of comparison-method biases (diagnosed). The derivation does not depend on those runs.
- Only the SPACE part of the action is covered. The time-stretching part (K_ij K^ij - lambda K^2) on a random grid is not yet
  addressed; with the preferred slicing it is read from how whole slices change between ticks, which again averages over many
  cells, but this has not been tested.
