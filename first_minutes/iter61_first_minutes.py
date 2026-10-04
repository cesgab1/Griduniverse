"""ITERATION 61 (pre-registered in PREREG_61.md)."""
import numpy as np
from scipy.optimize import brentq
alpha = 1/137.035999084; hbar_c = 197.3269804e-15; me = 0.51099895
r_e = alpha*hbar_c/me; pop = hbar_c/me
out = ["ITERATION 61: first-minutes interactions -> timing and size rule (expectations pre-registered)", "",
       f"S1 classical radius r_e = {r_e:.4e} m, pop size = {pop:.4e} m, ratio = {r_e/pop:.9f} = alpha ({alpha:.9f}) exactly",
       f"   Thomson cross-section 8 pi r_e^2 / 3 = {8*np.pi*r_e**2/3:.4e} m^2 (textbook 6.6524e-29)", ""]
# S2: protons
mp = 938.272; out.append(f"S2 proton scattering strength vs electron: (m_e/m_p)^2 = {(me/mp)**2:.1e}; proton size/pop = 4 (iteration 60): no pair picture")
# positron survival: pairs per photon vs leftover electrons per photon (eta x electrons per baryon)
eta = 6.1e-10; ne_per_b = 0.877
def pairs_per_photon(T):   # non-relativistic e+ density / photon density (g=2 per species), Maxwell-Boltzmann
    return (2*(me*T/(2*np.pi))**1.5*np.exp(-me/T))/(2*1.2020569/np.pi**2*T**3)
T_x = brentq(lambda T: pairs_per_photon(T) - eta*ne_per_b, 0.005, 0.2)
t_x = 1.32*(1/T_x)**2           # s, after annihilation (g* = 3.36); rough
out += [f"   positrons fall below the leftover electrons at T = {1000*T_x:.0f} keV, t ~ {t_x/60:.0f} minutes -> after that, only electrons hold pops", ""]
# S3 zero-parameter prediction (iteration 54 machinery)
H0 = 67.4e3/3.0857e22; h = 0.674; Om, Or = 0.315, 9.1e-5; OL = 1 - Om - Or
def E1(zz):
    a = 1/(1 + zz); g = lambda le: np.exp(2*le) - Om*a**-3 - Or*a**-4 - OL*(a*np.exp(le))**-0.5
    return np.exp(brentq(g, -10, 40))
rho_DE0 = OL*1.0537e4*h*h; n_e0 = 0.0224*1.878e-29/1.6726e-24*(1 - 0.245/2)
zs = 124.0
pred = alpha*511e3*n_e0*(1+zs)**3/((E1(zs)/(1+zs))**-0.5)
out += [f"S3 zero-parameter prediction (delta = alpha, payday at Compton heating = expansion, z = {zs:.0f}):",
        f"   dark energy today = {pred:.0f} eV/cm^3 vs measured {rho_DE0:.0f} -> ratio {pred/rho_DE0:.2f}; 21-cm step {1420.406/(1+zs):.1f} MHz",
        "S4 nucleosynthesis with the heavier electron: 0.7-1.5 sigma (iteration 57)."]
txt = "\n".join(out); print(txt); open("iter61_first_minutes.txt", "w").write(txt + "\n")
