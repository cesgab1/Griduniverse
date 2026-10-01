"""
Preferred-frame parameter alpha1 for AeST (Skordis & Zlosnik 2021), weak-field Newtonian regime.

1. AeST vector term  R - (K_B/2) F^2  is Einstein-aether (signature -+++, L = -M (grad u)(grad u),
   M = c1 g g + c2 dd + c3 dd - c4 u u g)  with  c1 = K_B, c3 = -K_B, c2 = 0  (so c13 = 0: gravitational waves at c).
2. AeST scalar terms  2(2-K_B) J.grad(phi) - (2-K_B)(1+lambda_s) Y   (Newtonian 'tracking' regime J -> lambda_s Y).
   For the slow-motion problem phi enters quadratically; eliminating it (grad phi = J_L/(1+lambda_s)) leaves
   + (2-K_B)/(1+lambda_s) |J_L|^2.  At first order in the source/aether velocity J_L^2 and J^2 give identical equations
   (the zeroth-order J is a pure gradient), so this is an Einstein-aether c4 = (2-K_B)/(1+lambda_s).
   Screening regime (J has Y^p, p >= 3/2): effective lambda_s -> infinity, c4 -> 0.
3. Check: Einstein-aether gives G_N = G/(1 - c14/2) = G (1+lambda_s) / ((1-K_B/2) lambda_s), which is exactly AeST's own
   Newtonian constant from its quasistatic action (eq. 6 of the PRL). Mapping confirmed.
4. Foster-Jacobson: alpha1 = -8(c3^2 + c1 c4)/(2c1 - c1^2 + c3^2)  ->  alpha1 = -4 c14 = -4[K_B + (2-K_B)/(1+lambda_s)].
   Bounds: |alpha1| < 1e-4 (lunar laser ranging); pulsars ~ 4e-5.  alpha2: Einstein-aether formula is singular here (c123 = 0),
   AeST's scalar carries that sector -> needs its own calculation (open).
"""
import numpy as np
a1 = lambda KB, ls: -4*(KB + (2-KB)/(1+ls))
GN = lambda KB, ls: (1+ls)/((1-KB/2)*ls)
print(f"{'K_B':>8s} {'lambda_s':>10s} {'alpha1':>11s}  verdict (|alpha1| < 1e-4)")
for KB in (0.5, 0.1, 1e-3, 2e-5):
    for ls in (1.0, 10.0, 1e5, np.inf):
        v = a1(KB, ls) if np.isfinite(ls) else -4*KB
        print(f"{KB:8.0e} {ls:10.0e} {v:+11.2e}  {'passes' if abs(v) < 1e-4 else 'fails'}")
print("\nAllowed corner: K_B < 2.5e-5 AND Solar System screened (effective lambda_s > ~8e4).")
print("Side effect (Skordis-Zlosnik 2022, eq.29): scalar sound speed^2 = (2-K_B)(1+K_B lambda_s/2)/(K2 K_B)")
for KB in (0.5, 2e-5):
    print(f"   K_B={KB}: with K2 = 7.5e3 -> c_s^2 = {(2-KB)/(7.5e3*KB):.3g} (lambda_s small); keeping c_s <= 1 needs K2 >= {(2-KB)/KB:.2g}")
