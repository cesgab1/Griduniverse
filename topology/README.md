# Iteration 62: does space wrap around? Planck SMICA 2015 temperature map (user-supplied), largest scales
Pre-registration: PREREG_62.md. Key fact: the last-scattering sphere has diameter 27.74 Gpc; a cube larger than that makes NO
matched circles, so our predicted range (30-60 Gpc) can only show up as subtle reshaping of the largest patterns.
Method: 20-degree smoothing, HEALPix Nside 8, SMICA confidence mask (580 of 768 pixels), monopole/dipole removed, exact Gaussian
likelihood; cubic-torus covariance from the discrete modes with CAMB transfer functions (multipoles 2-12; isotropic above);
amplitude free; best of 3000 cube orientations per size. Pixel window not available offline (same filter on data and model).
Runs:
- iter62_run.py (first run, 48 orientations): METHOD CHECK FAILED -- skies simulated with a 20 Gpc cube were not recognised.
  Output kept in iter62_first_run_coarse_orientations.txt.
- iter62b_run.py (RECORDED REVISION): exact rotations of the harmonic covariance, 3000 orientations; 50 isotropic sims and 10 cube
  sims each at 20, 30, 35 Gpc for calibration. Result: iter62b_result.txt.
Results:
- Method check PASSES now: 20 Gpc cube skies give Delta chi2 -17 to -52 (null median -2).
- REAL SKY at 20 Gpc: -4.0, outside the whole range of 20 Gpc cube skies -> a 20 Gpc cube is ruled out by our own pipeline
  (agrees with Planck's published searches). T1 as expected.
- 30 and 35 Gpc: cube skies are INDISTINGUISHABLE from infinite-space skies with this test (medians -8.7 vs -8.6; -2.5 vs -2.8).
  The data sit in the middle of both (chance to look as cube-like 0.36 / 0.58). INCONCLUSIVE, as pre-registered (T2).
- 40-60 Gpc: data slightly LESS cube-like than typical infinite skies (p ~ 0.98); nothing preferred. Look-elsewhere over all
  sizes: p = 0.52. No hint (T3 not triggered).
Learned: the book's 'testable half' is not testable with large-scale temperature maps -- a cube of 30+ Gpc hides beyond the
last-scattering sphere. Remaining handles: CMB polarisation (E-modes from reionisation probe a slightly different sphere),
higher-resolution searches for subtle correlations, or future 21-cm surveys of a much larger volume. Status of the finite
'book' route: not tested, not excluded.
