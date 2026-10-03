"""
Coalesce's question: do the vibrations of the jostling tell us the HEIGHT of dark energy?
1. Height: rho = (coupling)^2 x (energy scale) x Theta/(kappa H). Theta is free (0 < Delta N_eff < 0.107), the coupling is
   free: the vibrations fix how the height SCALES, not its value.
2. But the vibrations also make dark energy LUMPY, and the ratio lump/mean does not depend on the coupling at all.
   Spatial correlation of the jostled tension between two points a distance d apart (thermal part, memory rate gamma):
     C(d) ∝ ∫ dk sinc(k d) /(gamma^2 + k^2)   -> correlation length ~ c/gamma = c/(kappa H)
   With quadratic energy rho ∝ tau^2 (Gaussian tau), each correlated patch has delta rho / rho ~ sqrt(2).
3. Rough size of the gravitational potential such lumps make, and the large-angle CMB temperature they would imprint.
"""
import numpy as np
from scipy.integrate import quad
g = 1.0
C = lambda d: quad(lambda k: (np.sin(k*d)/(k*d) if d > 0 else 1.0)/(g*g + k*k), 0, 2000, limit=4000)[0]
c0 = C(0); out = ["Iteration 5 (question from Coalesce): what the vibrations of the jostling tell us", "",
                  "Spatial correlation of the tension, distance in units of c/(memory rate):"]
for d in (0.1, 0.3, 1, 3, 10):
    out.append(f"   d = {d:5.1f}:  C(d)/C(0) = {C(d)/c0:.3f}    (1 - e^-d)/d = {(1 - np.exp(-d))/d:.3f}")
kappa, OmDE = 3.0, 0.69
L = 1/kappa                                     # correlation length in Hubble lengths
Phi = 1.5*OmDE*np.sqrt(2)*L**2                  # Poisson: Phi ~ (3/2) Omega_DE (delta rho/rho) (L H/c)^2
out += ["", f"Correlation length: c/(kappa H0) = Hubble length / {kappa:.0f} ≈ {4400/kappa:.0f} Mpc today",
        f"Lumps: delta rho_DE / rho_DE ~ sqrt(2) on that scale (quadratic energy; coupling-free)",
        f"Potential of such lumps: Phi ~ (3/2) Omega_DE x sqrt(2) x (1/{kappa:.0f})^2 ≈ {Phi:.2f}",
        f"Large-angle CMB temperature imprint ~ Phi ~ {Phi:.1e}; observed large-angle anisotropy ~ 1e-5",
        f"-> the toy, as built, predicts dark-energy lumps ~{Phi/1e-5:.0e} times too strong. Escape routes:",
        "   (a) N independent jostled components: lumps shrink by 1/sqrt(N); need N >~ (Phi/1e-5)^2 ≈ " + f"{(Phi/1e-5)**2:.0e}",
        "   (b) the tension field responds only to the GLOBAL expansion (no local lumps), as in the VCDM embedding of Claim 1;",
        "       then the 'jostling' must act on the universe as a whole, not patch by patch",
        "   (c) the correlation length is much shorter than c/(kappa H) (needs noise not dominated by long waves)"]
txt = "\n".join(out); print(txt); open("iter5_vibrations.txt", "w").write(txt + "\n")
