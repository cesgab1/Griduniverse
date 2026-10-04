"""Light-Cone Model: can dark energy's size C be DERIVED from what the equation already contains (c, G, curvature radius a, H)?"""
import numpy as np
from scipy.optimize import brentq
Ok_meas, Ok_err, OL = 0.0023, 0.0011, 0.685
out = ["Light-Cone Model -- deriving C from the equation's own quantities (today's data: Omega_k = 0.0023 +/- 0.0011, Omega_DE = 0.685)", ""]
# D1: C equals the curvature energy density scale (3 c^2 / 8 pi G a^2)
d1 = lambda ok: ok*ok**0.25            # rho_DE0/rho_crit = Omega_k * (a0 H0/c)^(-1/2), a0 H0/c = Omega_k^(-1/2)
# D2: C = Planck density x (Planck length / curvature radius)^2 = c^2/(G a^2) = (8 pi/3) x curvature density
d2 = lambda ok: (8*np.pi/3)*ok**1.25
for name, f in (("D1 C = curvature energy density", d1), ("D2 C = Planck density x (Planck length/curvature radius)^2", d2)):
    pred = f(Ok_meas); need = brentq(lambda ok: f(ok) - OL, 1e-6, 10)
    out.append(f"{name}: predicts Omega_DE = {pred:.4f} (measured {OL}); matching would need Omega_k = {need:.3f} "
               f"(measured {Ok_meas} +/- {Ok_err}: {(need - Ok_meas)/Ok_err:.0f} sigma away)")
out += ["D3 C tied to the age since the first event (C ~ c^2/(G (c tau)^2)): dark energy would follow the total density at all times",
        "   -> already excluded by today's data (iteration 37: chi2 154 for 5 points).",
        "", "Conclusion: c, G, curvature and age cannot produce the size; a further ingredient is needed (quantum + something countable,",
        "or particle physics). If one is found, C stops being free and the equation predicts Omega_DE from matter and curvature."]
txt = "\n".join(out); print(txt); open("derive_C_options.txt", "w").write(txt + "\n")
