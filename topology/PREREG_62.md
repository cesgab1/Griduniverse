# Pre-registration, iteration 62 (committed BEFORE touching the data): does space wrap around? Planck SMICA map test
Prediction under test (iterations 32, 38, 53): flat 3-torus (cube) with side L ~ 30-60 Gpc (allowed 27.5-95 Gpc).
KEY FACT written first: the last-scattering sphere has diameter D = 2 x 13.87 Gpc = 27.74 Gpc. A cube LARGER than D produces NO
matched circles (the sphere never meets its own copies). So a matched-circle search cannot see our predicted range at all;
the only handle is the subtle reshaping of the largest-scale patterns (correlations between a_lm that are zero in infinite
space). Test chosen accordingly:
 Method: Planck 2015 SMICA temperature (user-supplied file), smoothed to 20 deg, HEALPix Nside 8, its own confidence mask
 (TMASK), monopole + dipole projected out. Exact Gaussian pixel likelihood. Covariance = torus modes k = 2 pi n / L (cubic, all n),
 CAMB temperature transfer functions (includes ISW and Doppler), multipoles 2-12 from the torus sum, 13-30 isotropic (shared),
 amplitude free. Orientation of the cube: 48 random orientations per L (maximum taken). Compared with infinite flat space.
 L values: 20, 25, 27.74, 30, 35, 40, 50, 60 Gpc.
 Calibration (done before reading the result): 30 infinite-space simulations through the identical pipeline -> null spread of
 the best-orientation Delta chi2; 10 simulations of a cube with L = 20 Gpc -> must be detected (method check).
Expectations:
 T1 L = 20 Gpc is strongly disfavoured by the data (already excluded by Planck's own searches) and detected in its simulations.
 T2 L = 30-60 Gpc: Delta chi2 inside the null spread -> INCONCLUSIVE (expected; the effect is too weak beyond ~1.1 D).
 T3 If instead the data PREFER some L in 30-60 Gpc beyond the null 95% level, that is a hint (not a detection) -- recorded.
 Caveats stated now: Gaussian smoothing + coarse pixels, foreground residuals near the mask edge, cubic torus only.
