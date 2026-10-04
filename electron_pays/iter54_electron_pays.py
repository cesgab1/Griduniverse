"""
ITERATION 54 (pre-registered in PREREG_54.md): does the electron's ~1% mass handover pay for dark energy, with the switch set by
the gas losing thermal contact with the CMB?
"""
import numpy as np, camb
from scipy.optimize import brentq
h = 0.674; Om = 0.315; ombh2 = 0.0224; Yp = 0.245; OL = 1 - Om
p = camb.set_params(H0=100*h, ombh2=ombh2, omch2=Om*h*h - ombh2 - 0.00064, YHe=Yp, tau=0.054)
r = camb.get_results(p)
z = np.logspace(-1.5, 3.3, 6000)
ev = r.get_background_redshift_evolution(z, ["x_e", "T_b"], format="array"); xe, Tb = ev[:, 0], ev[:, 1]; Tc = 2.7255*(1 + z)
H = r.hubble_parameter(z)*1e3/3.0857e22
sT, mec, arad = 6.6524e-29, 9.109e-31*2.998e8, 7.5657e-16; fHe = Yp/(4*(1 - Yp))
CH = 8*sT*arad*Tc**4/(3*mec)*xe/(1 + fHe + xe)/H                      # Compton heating rate / expansion rate
slope = np.gradient(np.log(Tb), np.log(1 + z))
rho_DE0 = OL*1.0537e4*h*h                                               # eV/cm^3
n_e0 = ombh2*1.878e-29/1.6726e-24*(1 - Yp/2)                            # electrons /cm^3 today
def E1(zz):   # Claim 1 background, beta = 1/2
    a = 1/(1 + zz); g = lambda le: np.exp(2*le) - Om*a**-3 - 9.1e-5*a**-4 - OL*(a*np.exp(le))**-0.5
    return np.exp(brentq(g, -10, 30))
need_beta = lambda zs: rho_DE0*((E1(zs))/(1 + zs))**-0.5/(n_e0*(1 + zs)**3*511e3)
need_lam = lambda zs: rho_DE0/(n_e0*(1 + zs)**3*511e3)
def cross(y, t, zmin=20):
    m = z > zmin; zz, yy = z[m], y[m]; i = np.where(np.diff(np.sign(yy - t)))[0]
    return [float(np.interp(t, [yy[j], yy[j+1]], [zz[j], zz[j+1]])) for j in i]
d_obs, d_err = 0.0099, 0.0049
out = ["ITERATION 54: does the electron's handover pay for dark energy? (expectations pre-registered)", "",
       f"electrons today {n_e0:.3e}/cm^3; dark energy today {rho_DE0:.0f} eV/cm^3; measured delta = {d_obs} +/- {d_err}", "",
       " switch criterion                       z_s    delta needed (Claim 1)   (pure Lambda)   z-score   21-cm step"]
crit = {"Compton heating rate = expansion rate": cross(CH, 1.0), "gas 10% colder than CMB": cross(Tb/Tc, 0.9),
        "halfway to free cooling (slope 1.5)": cross(slope, 1.5), "gas 50% colder than CMB": cross(Tb/Tc, 0.5)}
for k, v in crit.items():
    zs = max(v); nb = need_beta(zs); nl = need_lam(zs)
    out.append(f" {k:38s} {zs:5.0f}    {100*nb:6.2f}%                 {100*nl:6.2f}%      {(d_obs - nb)/d_err:+5.1f}    {1420.4/(1+zs):5.1f} MHz")
# E2 reionisation
low = z < 15
out += ["", f"E2 after reionisation: max Compton/H for z < 15 = {CH[low].max():.3f} (at z = {z[low][np.argmax(CH[low])]:.1f}); "
        f"at z = 0.89: {np.interp(0.89, z, CH):.3f}; at z = 0: {np.interp(0.0316, z, CH):.3f}"]
ch089 = np.interp(0.89, z, CH)
out.append(f"   smooth reversible switch (extra mass ~ delta x min(1, Compton/H)): Delta mu/mu at z = 0.89 ~ {d_obs*min(1, ch089):.1e} "
           f"vs bound 3.6e-7 -> {'EXCLUDED' if d_obs*min(1, ch089) > 3.6e-7 else 'allowed'}")
out.append("   threshold (Compton/H must exceed 1) or one-way handover: no re-coupling, methanol bound safe.")
# E3
q = n_e0*0.0099*511e3*(1 + 124)**3
out += ["", f"E3 energy released at z = 124: {q:.0f} eV/cm^3 = {q/(0.2606*(1+124)**4):.1e} of the CMB energy (FIRAS 6e-5) -- into gas it would"
        " heat it to ~1e7 K and re-ionise everything (optical depth ~4 vs 0.054): it must go to the grid."]
txt = "\n".join(out); print(txt); open("iter54_electron_pays.txt", "w").write(txt + "\n")
