# Pre-registration, iteration 78: do stars / black holes set dark energy's clock? (Oct 4 2026, written BEFORE running)
Measured input: cosmic star-formation history psi(z) = 0.015 (1+z)^2.7 / (1 + ((1+z)/2.9)^5.6) (Madau & Dickinson 2014).
Only the SHAPE is used (no free shape parameters); today's amount is matched to the measured dark energy as for every model.
Shapes (rho_DE(z)/rho_DE(0)):
 S1 'stars feed it': proportional to the total mass ever turned into stars, integral_z^inf psi dt. (Black-hole growth follows the
    same shape to good accuracy -- Madau & Dickinson 2014 -- so this also covers 'black holes grow it'.)
 S2 'coupled black holes' (Farrah et al. 2023, k = 3: black-hole mass grows as a^3 after birth): integral_z^inf psi (1+z')^3 dt.
 Cosmic time uses the Lambda-CDM clock (convention; effect small).
Fit: DESI DR2 BAO + Planck CMB distance priors + supernovae (Pantheon+, DES-Dovekie, Union3), same pipeline as Claim 1 (TAB model,
3 parameters each, like Lambda and Claim 1). 2 shapes x 3 supernova sets.
Rules: SUPPORTED if it beats Lambda by delta chi2 <= -4 in all three sets; EXCLUDED if worse than Lambda by > +9 in all three.
Expectations: S1 grows fast late (stellar mass doubled since z ~ 1) -> phantom-like, against the data's preference -> EXCLUDED.
S2 mostly built by z ~ 2, nearly constant afterwards -> Lambda-like: delta chi2 between -2 and +6 vs Lambda; not better than
Claim 1 (-5 to -7). Note: whether black holes can supply the full AMOUNT (Farrah et al.'s claim) is not tested here, only the shape.
