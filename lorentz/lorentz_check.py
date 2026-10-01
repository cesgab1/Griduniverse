"""
Lorentz-invariance check for the grid.
(1) Light on a grid of spacing l: omega = (2c/l) sin(k l/2) -> v = c(1 - (k l)^2/8)  (quadratic, n = 2).
    Standard form v = c[1 - (3/2)(E/E_QG2)^2]  ->  E_QG2 = sqrt(12) hbar c / l.
    LHAASO GRB 221009A (95%): E_QG1 > 1.47e20 GeV (subluminal), E_QG2 > 1.2e12 GeV (subluminal).
    A random (Poisson) mosaic removes the systematic shift on average but leaves comparable random spread; same order bound.
(2) AeST preferred frame, NAIVE estimate: AeST's vector term -(K_B/2)F^2 is Einstein-aether with c1 = K_B, c3 = -K_B, c2 = c4 = 0.
    Foster-Jacobson: alpha1 = -8(c3^2 + c1 c4)/(2c1 - c1^2 + c3^2) = -4 K_B.  Bound |alpha1| < 1e-4 (Sagi 2009 compilation).
    Ignores AeST's scalar coupling (2-K_B) J.grad(phi), which may change this -> needs the real AeST calculation.
"""
import numpy as np
hbarc = 1.9733e-16   # GeV m
lP = 1.616e-35
EQG2 = 1.2e12
lmax = np.sqrt(12)*hbarc/EQG2
print(f"(1) quadratic photon-speed bound -> grid spacing seen by light must be < {lmax:.1e} m  (Planck length {lP:.1e} m: allowed by {np.log10(lmax/lP):.1f} orders)")
for lab, l in (("Planck length", lP), ("superfluid channel spacing (BFK fit, ~10 micron)", 1e-5)):
    E = np.sqrt(12)*hbarc/l
    print(f"    {lab:48s} -> E_QG2 = {E:.1e} GeV  ({'allowed' if E > EQG2 else f'EXCLUDED by {np.log10(EQG2/E):.0f} orders'})")
print("    -> the 10-micron fluid channels cannot be what light travels through: light must see a Planck-fine grid, the fluid must be dark.")
print("\n(2) naive AeST preferred-frame parameter alpha1 = -4 K_B  (bound |alpha1| < 1e-4):")
for KB in (0.5, 0.05, 1e-3, 2.5e-5):
    a1 = -4*KB
    print(f"    K_B = {KB:<8g} alpha1 = {a1:+.1e}  {'passes' if abs(a1) < 1e-4 else 'fails'}")
