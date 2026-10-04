"""ITERATION 65 (pre-registered in PREREG_65.md): the light route for the excess payment."""
import numpy as np
from scipy.optimize import brentq
exec(open("iter64_size_no_threshold.py").read().split("Om, Or =")[0])          # z (high->low), heff, Tb, alpha, r
H0 = 67.4e3/3.0857e22; Om, Or = 0.315, 9.1e-5; OL = 1 - Om - Or
Hz = lambda zz: H0*np.sqrt(Om*(1+zz)**3 + Or*(1+zz)**4 + OL)
nH0 = 0.0224*1.878e-29/1.6726e-24*(1 - 0.245)*1e6                              # hydrogen atoms /m^3 today
ne_per_H = (1 - 0.245/2)/(1 - 0.245)
E_pop = alpha*511e3                                                             # eV
zm = 0.5*(z[1:] + z[:-1]); dz = -np.diff(z); dh = -np.diff(heff)                 # release per dz
dhdt = dh/dz*(1 + zm)*Hz(zm)                                                    # release rate per second (dz/dt = -(1+z)H)
# L1: absorption check for a 3.7 keV photon (hydrogen photo-ionisation ~ 6.3e-22 m^2 (E/13.6 eV)^-3)
sig = 6.3e-22*(E_pop/13.6)**-3
for zz in (380, 140, 40):
    out_abs = nH0*(1+zz)**3*sig*2.998e8/Hz(zz)
    print(f"z = {zz}: a {E_pop/1e3:.2f} keV photon is absorbed {out_abs:.0f}x per expansion time (energy redshift aside)")
alphaB = lambda T: 2.6e-19*(T/1e4)**-0.7                                        # m^3/s, case-B recombination
Tg = np.interp(zm[::-1], z[::-1], Tb[::-1])[::-1]
def extra_tau(f):
    inj = f*E_pop*ne_per_H*dhdt                                                 # eV per H atom per second
    ion_rate = inj/(3*13.6)                                                     # ionisations per atom per second (1/3 efficiency)
    nH = nH0*(1+zm)**3
    xe = np.clip(np.sqrt(ion_rate/(nH*alphaB(np.maximum(Tg, 10)))), 0, 1.0)      # ionisation-recombination balance
    dt = dz/((1 + zm)*Hz(zm))
    return np.sum(xe*nH*ne_per_H*6.6524e-29*2.998e8*dt), xe
out = ["ITERATION 65: can the excess go into light? (expectations pre-registered)", ""]
for f in (0.89, 0.96):
    t, xe = extra_tau(f)
    out.append(f"L1 X-ray route, light share {f:.2f}: ionised fraction at z = 140: {xe[np.argmin(abs(zm-140))]:.2f}; extra optical depth = {t:.2f} "
               f"(Planck allows ~0.015) -> {'EXCLUDED' if t > 0.015 else 'allowed'}")
fmax = brentq(lambda f: extra_tau(f)[0] - 0.015, 1e-12, 1.0)
out.append(f"   largest light share allowed: f_max = {fmax:.1e} (vs 0.89-0.96 needed)")
# L2: soft light into the CMB
rho_g = lambda zz: 0.2606*(1+zz)**4
n_e0 = 0.0224*1.878e-29/1.6726e-24*(1 - 0.245/2)
w = dh; frac = np.sum(w*E_pop*n_e0*(1+zm)**3/rho_g(zm))
out += ["", f"L2 soft light: total Delta rho/rho = {frac:.1e}; with 89-96% as light {0.89*frac:.1e}-{0.96*frac:.1e} vs FIRAS ~6e-5 -> allowed,"
        f" but needs ~{E_pop/(2.7*8.617e-5*(1+140)):.0e} microwave photons per electron (no mechanism)"]
txt = "\n".join(out); print(txt); open("iter65_light_route.txt", "w").write(txt + "\n")
