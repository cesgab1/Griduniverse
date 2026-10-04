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
