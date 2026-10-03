"""
ITERATION 3: can CMB PHOTONS be the jostlers of iteration 2? If not, what can, and what would it predict?

A. Photons, linear coupling. Gauge invariance: the tension can only couple to E or B (not to the potential A).
   ∫ E dt = -ΔA is a boundary term, and B carries a factor k. Test: filtered variance of the domain-averaged E field in a
   thermal photon gas vs memory 1/gamma. A random walk needs Var ∝ 1/gamma.
B. Photons, quadratic coupling (energy density, F^2): noise ∝ fluctuations of the photon energy in the domain. Scaling
   only: zero-frequency spectrum S ~ <δrho^2>·(R/c) ∝ Theta^5/R^2 -> Var ∝ Theta^5/(R^2 gamma) -> rho_lin ∝ (Theta^5/H)^(1/2)
   -> d ln rho/d ln a = (q + 1 - 5)/2: link-stretch exponent s = 5. Compare: our data scan excludes s >= 1.75 by > 5 sigma.
C. A massless conformally coupled SCALAR relic (dark radiation) works (iteration 2). If it was ever in thermal contact with
   ordinary matter, it adds to the radiation density: Delta N_eff = (4/7) (g*(nu decoupling)/g*(its decoupling))^(4/3).
"""
import numpy as np
from scipy.integrate import quad
def var_E(gamma, Theta, R=1.0):
    # one Cartesian component of E, transverse photons (2 pol): <|E_k,i|^2> per mode ∝ (2/3) * k/2 * coth(k/2Theta); thermal part only
    f = lambda k: (4*np.pi*k*k/(2*np.pi)**3)*(2/3)*(k/2)*(1/np.tanh(k/(2*Theta)) - 1)*np.exp(-k*k*R*R)/(gamma**2 + k*k)
    return quad(f, 0, 40/R, limit=400, points=[gamma, 1/R])[0]
def var_phi(gamma, Theta, R=1.0):
    f = lambda k: (k/(4*np.pi**2))*(1/np.tanh(k/(2*Theta)) - 1)*np.exp(-k*k*R*R)/(gamma**2 + k*k)
    return quad(f, 0, 40/R, limit=400, points=[gamma, 1/R, 2*Theta])[0]
out = ["ITERATION 3: who jostles the grid tension?", "", "A. Thermal photons, tension driven by the electric field (Theta = 1, R = 1):"]
for g in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5):
    out.append(f"   memory {1/g:6.0e}:  photon E: Var = {var_E(g, 1.0):.4e}     scalar (iteration 2): Var = {var_phi(g, 1.0):.4e}")
out.append("   -> photon case SATURATES (no random walk: ∫E dt = -ΔA); the scalar grows ∝ memory. Same holds for B (extra factor k).")
out.append("")
out.append("B. Photon energy-density coupling: Var ∝ Theta^5/(R^2 gamma) -> stretch exponent s = 5; data (stretch_scan) allow")
out.append("   s = 0.72-1.06; s = 1.75 already costs Delta chi2 >= +20 vs the best. EXCLUDED.")
out.append("")
out.append("C. A massless conformally coupled scalar relic ('grid radiation'): Delta N_eff if it last shared a temperature with ordinary")
out.append("   matter when the number of particle types was g*:")
for lab, gs in (("above the electroweak scale (all SM particles, g* = 106.75)", 106.75), ("between electroweak and QCD (g* ~ 61.75)", 61.75),
                ("just below QCD (g* = 17.25)", 17.25), ("never: born with the grid, same temperature as SM at the Planck era, SM-only", 106.75)):
    dn = 4/7*(10.75/gs)**(4/3); out.append(f"   {lab:72s}: Delta N_eff = {dn:.3f}")
out.append("   Measured (BBN + Planck/ACT/SPT + DESI DR2, arXiv:2603.13226): N_eff = 2.990 ± 0.070, Delta N_eff < 0.107 (95%).")
out.append("   -> allowed only if the scalar decoupled ABOVE the QCD transition: prediction Delta N_eff ≈ 0.027-0.057")
out.append("      (smaller only if many unknown heavy particles existed). Its temperature today: (10.75/106.75)^(1/3) x 1.95 K = "
           f"{(10.75/106.75)**(1/3)*1.945:.2f} K for the 106.75 case.")
out.append("   Reach: Simons Observatory sigma(N_eff) ~ 0.05, CMB-S4-class ~ 0.03: the upper half of the range is testable.")
txt = "\n".join(out); print(txt); open("iter3_who_jostles.txt", "w").write(txt + "\n")
