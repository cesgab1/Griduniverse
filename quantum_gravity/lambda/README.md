# Why lambda = 1? (Oct 2026)

lambda weights the two stretching terms, K_ij K^ij - lambda K^2. Physically: G_cosmo / G_N = 1 / (1 + 3(lambda - 1)/2), i.e.
lambda != 1 makes gravity for the cosmic expansion differ from gravity for planets.

## What the data say (lambda_numbers.txt)
- Published: big-bang nucleosynthesis |G_cosmo/G_N - 1| <~ 1/8 -> 0 <~ lambda - 1 <~ 0.1; 0 <~ lambda - 1 <~ 0.01 when the
  other preferred-frame parameters vanish (Herrero-Valea, 'The status of Horava gravity', arXiv:2307.13039). Other parameters:
  graviton speed |beta'| <~ 1e-15 (GW170817), |alpha| <~ 1e-7 (pulsars/PPN).
- Our estimate from N_eff = 2.990 +/- 0.070 (treating a G change as extra radiation; approximate at the CMB):
  G_cosmo/G_N - 1 = -0.007 +/- 0.009 -> lambda - 1 in [-0.008, +0.017] (2 sigma).
- So lambda = 1 to about 1-2%.

## What the grid says
- CORRECTION: the earlier 'lambda ~ -0.5' estimate omitted contact terms and gives the same tensor structure for any cutoff, so it
  cannot discriminate grids; WITHDRAWN as a number. What stands: a frame-free (covariant) cutoff gives lambda = 1 exactly by
  symmetry; a preferred-frame cutoff generically gives lambda != 1 (value needs a full one-loop calculation, not done).
  This is the known general problem: Lorentz violation at the cutoff 'percolates' into low-energy physics through loops
  (Collins, Perez, Sudarsky, Urrutia & Vucetich, PRL 93, 191301, 2004). A grid with a built-in 'now' does not by itself give
  lambda = 1. This is a genuine problem for the model, recorded as such.

## Three routes to lambda = 1
1. **Grid random in SPACETIME, not just space** (BORROWED principle: a random sprinkling of points in spacetime is Lorentz
   invariant on average; Bombelli, Henson & Sorkin 2009). Then induced gravity (iteration 40) is covariant and lambda = 1
   exactly (Sakharov: a covariant cutoff induces R(4) = R + K_ij K^ij - K^2). The preferred 'now' would NOT be built into the
   grid; it would be carried by a field (the khronon / VCDM slicing) that enters only through the dark-energy function P(K).
   Cost: the Horava z = 3 small-scale rule (watermark, iteration 30/35; dimensional reduction, iteration 42) would no longer
   come from the grid's slicing; causal-set-style grids show their own small-scale dimensional reduction, untested here.
2. **Absorbed into the dark-energy function** (our VCDM home): any (lambda - 1) K^2 piece is just a K^2 term inside the
   function of K that already carries Claim 1. Claim 1's jostling produces K^(-1/2), not K^2, and our fits assumed
   G_cosmo = G_N; the data cap any K^2 piece at ~1-2%. This makes lambda a measurable number, not an independent mystery,
   but does not by itself explain why the K^2 piece is small.
3. **Flow towards lambda = 1 at large scales** (BORROWED, partial): in 2-D Horava gravity the renormalization flow runs towards
   lambda = 1 at long distances (but with loss of control as it gets there); not established in 3+1-D.

## Verdict
lambda = 1 is NOT explained by a grid with a built-in now (preferred-frame cutoff -> lambda != 1 generically). The cleanest fix is route 1: make the grid random in
spacetime (statistically Lorentz invariant) and let the preferred slicing live only in the dark-energy sector (route 2). That
keeps Claim 1 and gives lambda = 1, at the price of re-deriving the small-scale (z = 3) results. Proposed next test: compute
lambda induced by a field on a spacetime-random grid vs a space-random grid with a preferred time step.

## Trying all routes (Coalesce: 'time has to be its own dimension as well')
- Route 1 premise TESTED (iteration 46b): a grid random in space AND time has no preferred frame (diamond counts identical at all
  rapidities); a synchronised-tick grid carries a frame (weakly visible in counts). With time as its own random dimension, induced
  gravity is frame-free -> lambda = 1 by symmetry.
- Route 2 (VCDM): lambda - 1 is a K^2 piece of the dark-energy function; data: lambda - 1 in [-0.008, +0.017].
- Route 3 (flow to lambda = 1): borrowed 2-D evidence only; not computed here.
- Adopted as the working foundation: cells scattered randomly through space AND time (time a real dimension). The preferred 'now'
  (the book's pages) is not built into the grid; it is set by the universe's contents (the slicing of constant expansion rate,
  York time -- the same K that Claim 1 uses), entering only through the dark-energy function. Cost: the z = 3 small-scale results
  (watermark tilt, dimensional reduction) must be re-derived for a spacetime-random grid.
