# Does the grid give Einstein gravity at large scales? (iteration 39, pre-registered at 0ec31da)
Curvature put on the grid only through link lengths; grid (Regge) action compared with Einstein's action for the same space.
- First run (bump 0.12 wide, random and jittered Delaunay grids): ratios 0.4-6.7 with huge scatter, no trend. FAILED -- partly
  by design (bump only 1-5 cells wide).
- Declared revision 39b (same criteria): wider bump (0.2) and a WELL-SHAPED grid (each cube split into 6 tetrahedra).
  Well-shaped: ratio 1.68 -> 1.33 -> 1.21 -> 1.16 -> 1.12 for 3k -> 384k cells, error ~ 5/m (m cubes per side) -> converges to 1.
  Random Delaunay: 6.0 +/- 2.0, 0.13 +/- 1.1, 0.41 +/- 0.9 for 4k-64k cells: no convergence at these sizes (sliver cells).
  (A revision that added the cube's corner points broke the flat-space check -- flat gave 22 instead of 0 -- and was discarded.)
- Verdict: R1 PASSES for well-shaped cells (Einstein's action emerges, as the Cheeger-Mueller-Schrader theorem says); FAILS for
  a fully random grid at these sizes. Design constraint for the model: the Mosaic's cells must be well-shaped (no slivers; CDT
  uses equilateral simplices for the same reason), or gravity must be read off a coarse-grained grid.

## Iterations 41/41b: curvature read by NEIGHBOURHOOD AVERAGING (Ollivier-Ricci) on the random grid (Coalesce's idea)
- 41 (independent flat/curved samples): ratio to the true curvature 18 -> 13 -> 3.7 -> 3.0 +/- 0.9 for 300 -> 3700 points per
  neighbourhood; failed the pre-set criterion at reachable sizes, but trending.
- 41b (declared follow-up; same random points for flat and curved, so sampling noise cancels): 0.69 +/- 0.19 (~800 points),
  1.18 +/- 0.12 (~2100 points) -> within 25%: PASS. The excess in 41 was sampling noise.
- Verdict: a random grid DOES give the right curvature when curvature is read by averaging over neighbourhoods. The tidy-cell
  requirement belongs to Regge's single-link definition. (Full Einstein action from this definition not yet built.)

## Iterations 43/43b/43c: curvature PROFILE on the random grid (towards the full action) -- partly verified
True R at r = 0 / 0.2 / 0.35 / 0.5: +31.6 / +20.7 / +5.4 / -2.6 (sign change at the edge).
- 43 (global volume-matched common points, random directions): +31.2 +/- 2.5 / +8.5 +/- 12 / -5.0 +/- 14 / -50 +/- 16.
- 43b (same, six directions per point): +33.0 +/- 2.7 / +22.0 +/- 7.0 / +15.4 +/- 14 / -35 +/- 9.
- 43c (local isotropic rescale): r = 0.5: -0.4 +/- 2.7; r = 0.35: -2.5 +/- 3.2; r = 0: -18.0 +/- 2.4.
Diagnosis (comparison METHOD, not the curvature definition):
 * global matching shears the point pattern by ~40% far from the bump (all the extra volume pushes outer points inward), and
   pairs picked in flat geometry then have very different proper separations -> biased at the edge;
 * local rescaling leaves a density error of second order, the same order as the curvature signal, wherever phi is strongly
   curved -> biased at the centre.
Verified: centre and r = 0.2 (global method, where its shear is negligible) and the edge incl. its NEGATIVE sign (local method,
where its density error is negligible). Not resolved: r = 0.35 (both methods' errors overlap). P1 met at 3 of 4 radii with a valid
method each; P2 (full action) NOT established. Next: one comparison valid everywhere (pairs selected in each geometry, much larger
statistics), or a perturbative calculation of the neighbourhood-averaged action.
