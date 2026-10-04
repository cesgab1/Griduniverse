"""Light-Cone Model: what sets C, the starting strength of fading dark energy? Numbers only."""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
Om, OL, Ok = 0.315, 0.685, 0.0023
out = ["What sets C (fading dark energy, H^2 = (8 pi G/3) rho_m + c^2/a^2 + C (aH/c)^(-1/2), a = curvature radius)", ""]
aHc = 1/np.sqrt(Ok)
out.append(f"1. Today aH/c = 1/sqrt(Omega_k) = {aHc:.1f} -> C = rho_DE0 x {np.sqrt(aHc):.2f} = {OL*np.sqrt(aHc):.2f} x critical density today.")
out.append("   Meaning: C is the dark-energy density the equation would have in the EMPTY light-cone state (aH = c, Special Relativity).")
out.append("2. The fading law itself has no built-in scale (a pure power): its SHAPE cannot fix C.")
# far-future attractor a ~ t^5: rho R_eh^2 = (3/8pi)(25/16) c^4/G ; quantum geometric-mean rule rho = hbar c/(2 l^2 L^2) with L = R_eh
k_att = 3/(8*np.pi)*25/16
l2 = 1/(2*k_att)
out.append(f"3. Far future (only fading dark energy): rho x horizon^2 = {k_att:.4f} c^4/G for ANY C. Matching the quantum rule")
out.append(f"   rho = hbar c/(2 l^2 L^2) there fixes the smallest length: l^2 = {l2:.2f} l_P^2 (l = {np.sqrt(l2):.2f} l_P) -- between Planck (1)")
out.append("   and the black-hole value (2) -- a derived candidate, independent of C.")
# today: how far from the attractor
H0R = quad(lambda x: 1/(np.exp(x)*np.sqrt(Om*np.exp(-3*x) + OL)), 0, 60, limit=400)[0]
today = OL*3/(8*np.pi)*H0R**2
out.append(f"4. Today rho x horizon^2 = {today:.4f} c^4/G = {today/k_att:.2f} of the far-future value: we are part-way to the attractor.")
out.append("   So C does not set the size scale; it sets WHEN dark energy takes over (acceleration began z ~ 0.68): C is the 'why now' number.")
txt = "\n".join(out); print(txt); open("what_sets_C.txt", "w").write(txt + "\n")
