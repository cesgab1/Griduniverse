"""ITERATION 64 (pre-registered in PREREG_64.md): dark-energy size, threshold-free."""
import numpy as np, camb, ctypes
from camb.baseconfig import camblib
from scipy.optimize import brentq
alpha = 1/137.035999084
camblib.set_vconst.argtypes = [ctypes.c_double]*2; camblib.set_vswitch.argtypes = [ctypes.c_double]*2
camblib.set_vconst(1.0, 1 + alpha); camblib.set_vswitch(-1.0, 30.0)
p = camb.set_params(H0=67.4, ombh2=0.0224, omch2=0.1200, YHe=0.245, tau=0.054); p.Recomb.use_rosenbrock = False
r = camb.get_results(p)
z = np.logspace(np.log10(5.0), np.log10(1500.0), 20000)[::-1]                # high -> low z
Tb = r.get_background_redshift_evolution(z, ["T_b"], format="array")[:, 0]
s = np.gradient(np.log(Tb), np.log(1 + z)); h = np.clip(2 - s, 0, 1)
heff = np.minimum.accumulate(h)                                               # one-way: never restored
w = -np.diff(heff); zm = 0.5*(z[1:] + z[:-1])                                 # release weights (sum = h(start) - h(end))
Om, Or = 0.315, 9.1e-5; OL = 1 - Om - Or
def E1(zz, beta):
    a = 1/(1 + zz); g = lambda le: np.exp(2*le) - Om*a**-3 - Or*a**-4 - OL*(a*np.exp(le))**-beta
    return np.exp(brentq(g, -10, 40))
n_e0 = 0.0224*1.878e-29/1.6726e-24*(1 - 0.245/2)                              # /cm^3
rho_meas = OL*1.0537e4*0.674**2                                               # eV/cm^3
eV_cm3_to_kg_m3 = 1.602176634e-19/2.99792458e8**2*1e6
def today(beta):
    sel = w > 0
    return float(np.sum(w[sel]*alpha*511e3*n_e0*(1 + zm[sel])**3*np.array([(E1(x, beta)/(1 + x))**beta for x in zm[sel]])))
out = ["ITERATION 64: dark-energy size with no chosen threshold (expectations pre-registered)", ""]
mean_z = np.sum(w*zm)/np.sum(w); q = np.cumsum(w)/np.sum(w)
out.append(f"hold fraction from the gas temperature: released between z = {np.interp(0.1, q, zm):.0f} (10%) and z = {np.interp(0.9, q, zm):.0f} (90%),"
           f" median z = {np.interp(0.5, q, zm):.0f}; total released fraction {np.sum(w):.3f}")
for b in (0.5,):
    v = today(b); out.append(f"beta = 1/2 (model law): today's dark energy = {v*eV_cm3_to_kg_m3:.2e} kg/m^3 = {v/rho_meas:.2f} x measured")
vals = {b: today(b)/rho_meas for b in (0.41, 0.63, 0.85)}
out.append(f"beta measured 0.63 (+/-0.22): {vals[0.63]:.2f} x measured (range {vals[0.41]:.2f} - {vals[0.85]:.2f})")
out.append(f"measured value: {rho_meas*eV_cm3_to_kg_m3:.2e} kg/m^3")
txt = "\n".join(out); print(txt); open("iter64_size_no_threshold.txt", "w").write(txt + "\n")

# SECONDARY CHECK (added AFTER seeing the result, not pre-registered): a second threshold-free hold measure,
# h2 = Gamma/(Gamma + H) (energy-trading rate vs expansion rate), same one-way rule and same weighting.
xe = r.get_background_redshift_evolution(z, ["x_e"], format="array")[:, 0]
Hs = r.hubble_parameter(z)*1e3/3.0857e22; Tc = 2.7255*(1 + z); fHe = 0.245/(4*(1 - 0.245)); me = 1 + alpha
G = 8*6.6524e-29/me**2*7.5657e-16*Tc**4/(3*9.109e-31*me*2.998e8)*xe/(1 + fHe + xe)
h2 = np.minimum.accumulate(G/(G + Hs)); w = -np.diff(h2); q = np.cumsum(w)/np.sum(w)
v2 = today(0.5)
out2 = [f"SECONDARY (post-hoc) Gamma/(Gamma+H): released z = {np.interp(0.1, q, zm):.0f} (10%) .. {np.interp(0.9, q, zm):.0f} (90%), "
        f"median {np.interp(0.5, q, zm):.0f}; beta = 1/2: {v2/rho_meas:.2f} x measured"]
print("\n".join(out2)); open("iter64_size_no_threshold.txt", "a").write("\n".join(out2) + "\n")
