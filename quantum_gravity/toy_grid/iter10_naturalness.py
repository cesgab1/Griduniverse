"""
ITERATION 10: is the coupling that sets the HEIGHT of dark energy small enough to keep the wobble FEATHERWEIGHT?
(Coalesce: layers vibrate on their own since the hot start; pinned but free in some directions.)

Minimal model per layer (hbar = c = 1, energies in eV):
  stretch field Phi (tension tau = dPhi/dt), wobble h along a FREE direction (massless when g = 0), coupling L_int = g Phi h
  (the only way to make the wobble push the tension directly:  tau' = -3H tau + g h ;  g has dimension mass^2)
  energy of the tension = (1/2) tau^2 per layer (ordinary motion energy, canonical field); N layers.
Height: <tau^2> = g^2 T_s / (4 pi * 3H)  (iteration 2 white-noise result, memory 1/(3H))  ->  rho_DE = N g^2 T_s / (24 pi H)
  extra radiation: N (4/7)(T_s/T_nu)^4 <= 0.107  ->  N T_s is largest at the largest allowed N (species window <= 1.2e15)
Stability & featherweight: mass matrix [[m_Phi^2, -g], [-g, m_h^2]] must be positive: g <= m_Phi m_h. The tension must be a
  Hubble-damped motion (m_Phi <~ 3H) and the wobble must be light enough to give zero-frequency noise (m_h <~ 3H):
  g_max = (3 H0)^2.
"""
import numpy as np
H0 = 67.4e3/3.0857e22*6.582e-16           # eV
rho = 0.69*3*(67.4e3/3.0857e22)**2/(8*np.pi*6.674e-11)*(2.998e8)**2   # J/m^3
rho_eV4 = rho/1.602e-19*(1.9733e-7)**3    # eV^4  (hbar c = 1.9733e-7 eV m)
Tnu = 1.945*8.617e-5
C = 0.107*7/4*Tnu**4
out = ["ITERATION 10: height vs featherweight (minimal linear-coupling model)", "",
       f"observed dark energy rho = {rho_eV4:.2e} eV^4 = ({rho_eV4**0.25*1e3:.2f} meV)^4;  H0 = {H0:.2e} eV", ""]
gmax = (3*H0)**2
out.append(f"stability + featherweight: g <= (3 H0)^2 = {gmax:.1e} eV^2")
for N in (1e7, 1e10, 1.2e15):
    NT = C**0.25*N**0.75                                  # max N*T_s allowed by Delta N_eff
    gmin = np.sqrt(24*np.pi*H0*rho_eV4/NT)
    out.append(f"N = {N:7.1e} layers: height needs g >= {gmin:.1e} eV^2   -> too strong by {gmin/gmax:.0e}x")
out += ["",
        "-> In this minimal model the coupling needed for the height is ~1e40x too strong to keep the wobble light and the",
        "   system stable: 'technically natural' FAILS here. The height and featherweight problems do NOT merge this way.",
        "",
        "Other routes checked in the same model:",
        "  * coupling to a PINNED-direction wobble (heavy, m_h >> 3H): no zero-frequency noise -> no random walk (dead).",
        "  * derivative coupling (g dPhi/dt h): the push becomes dh/dt, whose zero-frequency noise vanishes (iteration 3) (dead).",
        "  * dark energy stored as potential energy of the wandering field instead of motion: smaller, not larger (m_Phi <~ H).",
        "  * a large 'stiffness' multiplying the motion energy (Z tau^2/2): rescaling Phi -> sqrt(Z) Phi gives the same physics",
        "    with g -> g/sqrt(Z); the requirement is identical (no escape).",
        "  * instability if g > m_Phi m_h: growth time ~ 1/sqrt(g) ~ milliseconds for the g needed (catastrophic).",
        "What could still work (OPEN, each a new assumption): a NONLINEAR link between the wobbles and the tension energy that",
        "  is not a mass mixing; or the noise reaching the tension through gravity/the time stamp instead of a direct coupling."]
txt = "\n".join(out); print(txt); open("iter10_naturalness.txt", "w").write(txt + "\n")
