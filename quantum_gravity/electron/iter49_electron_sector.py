"""
ITERATION 49 (pre-registered in quantum_gravity/dynamics/PREREG_47_49.md, E1-E2): does Claim 2 survive laboratory and
astronomical tests of the electron?
"""
import numpy as np
from scipy.optimize import brentq
out = ["ITERATION 49: electron-sector check for Claim 2 (expectations pre-registered)", ""]
# E1 Lorentz violation
out += ["E1 Lorentz violation: Claim 2 as a mass that depends on the tension field only as a scalar, m_e(phi) or m_e(averaged K).",
        "   A scalar mass term has no direction: no Lorentz-violating electron terms at tree level. Loop leakage suppressed as in",
        "   iteration 48 (c ~ 0). LEP synchrotron bound |kappa_tr - 4/3 c_00| < 5e-15 (Altschul 2009): PASS.",
        "   Requirement made explicit: the electron must NOT couple to the frame's direction u^mu (e.g. psi-bar gamma^mu u_mu psi);",
        "   that would be a direct Lorentz-violating term, bounded at 1e-15 or better.", ""]
# E2 time variation: switch used in every fit
delta = 0.0099                       # best fit (EQUATIONS.md E5)
mu_obs, mu_err = -1.0e-7, np.hypot(0.8e-7, 1.0e-7)        # Bagdonaite et al. 2013, z = 0.89, mu = m_p/m_e
mu_lim = abs(mu_obs) + 2*mu_err                          # 2 sigma two-sided envelope ~ 3.6e-7
clock = 8e-17 + 2*36e-18                                 # atomic clocks, mu_dot/mu = (-8 +/- 36)e-18 /yr (Lange et al. 2021, PRL 126, 011102), 2 sigma envelope ~ 1.5e-16 /yr
H0yr = 67.5e3/3.086e22*3.156e7
def frac(z, zs, dz): return 0.5*(1 + np.tanh((z - zs)/dz))
def Ez(z): return np.sqrt(0.315*(1+z)**3 + 0.685)
out.append(f"E2 switch m_e = 1 + delta/2 [1 + tanh((z - z_s)/dz)], delta = {delta}; bounds: |Delta mu/mu|(z=0.89) < {mu_lim:.1e},")
out.append(f"   |mu_dot/mu| today < {clock:.0e}/yr. (A heavier past electron makes mu = m_p/m_e smaller.)")
for zs in (100, 150, 200):
    for dz in (30, 15):
        dm = delta*frac(0.89, zs, dz); dfdz = delta*(1 - np.tanh((0 - zs)/dz)**2)/(2*dz)
        rate = dfdz*H0yr*Ez(0)                           # |d ln m_e/dt| today = |d delta_m/dz| * H (1+z), z = 0
        dzmax = brentq(lambda d: delta*frac(0.89, zs, d) - mu_lim, 1, 200)
        out.append(f"   z_s = {zs:3d}, dz = {dz:2d}: residual at z=0.89 = {dm:.1e} ({'FAIL' if dm > mu_lim else 'pass'}); "
                   f"drift today = {rate:.1e}/yr ({'FAIL' if rate > clock else 'pass'}); widest allowed dz = {dzmax:.1f}")
out += ["", "Verdict: the tanh tail used in the fits (dz = 30) is EXCLUDED by the z = 0.89 methanol bound for z_s = 100 and 150",
        "(allowed for z_s = 200). The CMB result does not depend on the width (switch_check), so the fix is a sharper switch:",
        "z_s/dz >~ 5. Width is a free choice (accommodation) unless the grid fixes it. Prediction made sharper: if the switch is",
        "near the allowed edge, a residual Delta mu/mu ~ 1e-7 at z ~ 1-3 would appear in next-generation (ELT/ALMA) measurements;",
        "a sharp (step-like) switch predicts exactly zero there."]
txt = "\n".join(out); print(txt); open("iter49_electron_sector.txt", "w").write(txt + "\n")
