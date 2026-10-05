"""Iteration 80: dark energy size from stellar physics, no cosmological data. See PREREG_80.md."""
import numpy as np
G, c, hbar = 6.674e-11, 2.998e8, 1.0546e-34
Msun, Lsun = 1.989e30, 3.828e26; Gyr = 3.156e16
rDE = 5.85e-27
mp, mHe = 1.007825, 4.002602     # atomic masses (u): hydrogen-1, helium-4
eff = (4*mp - mHe)/(4*mp)
out = ["Iteration 80: dark energy's size from stellar physics (no cosmological data). PREREG_80.md", ""]
t_sun = 0.1*eff*Msun*c**2/Lsun
out.append(f"1. fusion efficiency {eff:.4f}; Sun's main-sequence lifetime = {t_sun/Gyr:.1f} Gyr")
# Eddington standard model with electron scattering only (Thomson opacity, X = 0.7), mu = 0.61
kappa = 0.02*(1 + 0.7)            # m^2/kg
mu = 0.61
f = lambda b: 1 - b - 0.00298*mu**4*b**4*(1.0)**2
from scipy.optimize import brentq
beta = brentq(f, 0.5, 1)
L_es = 4*np.pi*c*G*Msun*(1 - beta)/kappa
out.append(f"   Eddington Sun with electron-scattering opacity only: L = {L_es/Lsun:.1f} x real -> lifetime {t_sun/Gyr*Lsun/L_es:.2f} Gyr "
           "(explains the old 670x: the particle-constant formula ignores atomic absorption, which dims real stars)")
links = {"L1 own clock": 1.0, "L2 equals matter": (2/3*np.arcsinh(1))**2, "L3 acceleration begins": (2/3*np.arcsinh(1/np.sqrt(2)))**2}
out.append("2-3. predicted / measured dark energy:")
vals = []
for M in (0.8, 1.0, 1.2):
    t = t_sun*M**-2.5
    row = []
    for k, fct in links.items():
        r = fct*3/(8*np.pi*G*t**2)/rDE; vals.append(r); row.append(f"{k} {r:.2f}x")
    out.append(f"   star {M} M_sun (t_obs = {t/Gyr:.1f} Gyr): " + " | ".join(row))
lo, hi = min(vals), max(vals)
out.append("")
out.append(f"RESULT: predicted size = {lo:.2f}x to {hi:.2f}x measured (central Sun, link L2: {links['L2 equals matter']*3/(8*np.pi*G*t_sun**2)/rDE:.2f}x). "
           f"{'Consistent (1 inside range)' if lo < 1 < hi else 'NOT consistent'}. Observer argument, not a mechanism.")
txt = "\n".join(out); print(txt); open("iter80_stellar_size.txt", "w").write(txt + "\n")
