"""ITERATION 58 (pre-registered in PREREG_58.md): open-pop mechanism; proton test."""
import numpy as np
from scipy.optimize import brentq
e, eps0, hbar, c = 1.602176634e-19, 8.8541878128e-12, 1.054571817e-34, 2.99792458e8
me, mp = 9.1093837015e-31, 1.67262192e-27
alpha = e**2/(4*np.pi*eps0*hbar*c)
U = e**2/(4*np.pi*eps0*(hbar/(me*c)))
out = ["ITERATION 58: open-pop mechanism (expectations pre-registered)", "",
       f"M1 Coulomb energy of a pair at the electron's pop size / (m_e c^2) = {U/(me*c**2):.9f}; alpha = {alpha:.9f} (exact identity)",
       f"   single charge's field energy outside the pop size / (m_e c^2) = {U/2/(me*c**2):.6f} (alpha/2, iteration 57)",
       "M2 photon cloud: alpha/pi x log factors -- not exact (no single natural number): fails the 'exact' test", ""]
# Proton test: protons' own direct thermal contact with photons (Compton heating per proton) vs expansion
H0 = 67.4e3/3.0857e22; Om, Or = 0.315, 9.1e-5; OL = 1 - Om - Or
H = lambda z: H0*np.sqrt(Om*(1+z)**3 + Or*(1+z)**4 + OL)
sT, arad = 6.6524e-29, 7.5657e-16
GCp = lambda z: 8*sT*(me/mp)**2*arad*(2.7255*(1+z))**4/(3*mp*c)            # proton's own Compton heating rate (fully ionised)
zp = brentq(lambda z: GCp(z)/H(z) - 1, 1e2, 1e12)
# energy protons would hand over, compared with what dark energy may hold then (Claim 1) and with all matter then
h = 0.674; rho_DE0 = OL*1.0537e4*h*h                                         # eV/cm^3
n_p0 = 0.0224*1.878e-29/1.6726e-24*(1 - 0.245/2)                              # free protons ~ electrons (He nuclei aside)
def E1(zz):
    a = 1/(1 + zz); g = lambda le: np.exp(2*le) - Om*a**-3 - Or*a**-4 - OL*(a*np.exp(le))**-0.5
    return np.exp(brentq(g, -10, 60))
pay = alpha*938.272e6*n_p0*(1 + zp)**3
allowed = rho_DE0*((E1(zp))/(1 + zp))**-0.5
out += [f"Proton test: a proton loses its OWN thermal contact with light at z = {zp:.2e} (~160 years after the start)",
        f"   if it paid alpha m_p c^2 then: {pay:.2e} eV/cm^3 vs the dark energy Claim 1 allows then {allowed:.2e} -> {pay/allowed:.1e} x too much",
        "   -> protons must NOT pay: the lock must be electron-mediated (photons couple directly to electrons; protons only via Coulomb)."]
txt = "\n".join(out); print(txt); open("iter58_open_pop.txt", "w").write(txt + "\n")
