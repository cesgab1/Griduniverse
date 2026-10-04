"""
ITERATION 56 (pre-registered in PREREG_56.md): countable 'pop' quantities vs the electron's needed 0.60%.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
target, lo, hi = 0.0060, 0.0040, 0.0090
alpha = 1/137.035999
MP = 1.22091e19                                       # GeV
ferm = {"e": (0.000511, 1, -1), "mu": (0.10566, 1, -1), "tau": (1.77686, 1, -1),
        "u": (0.00216, 3, 2/3), "d": (0.00467, 3, -1/3), "s": (0.0934, 3, -1/3), "c": (1.27, 3, 2/3),
        "b": (4.18, 3, -1/3), "t": (172.7, 3, 2/3)}
res = []
den = sum(Nc*Q**2*np.log(MP/m) for m, Nc, Q in ferm.values())
res.append(("A  vacuum charged-pop share (polarisation)", np.log(MP/0.000511)/den))
inv3 = sum(1/m**3 for m, Nc, Q in ferm.values()); num3 = sum(m**3 for m, Nc, Q in ferm.values())
res += [("B  share of Compton volumes (1/m^3)", (1/0.000511**3)/inv3), ("B  share of pop number density (m^3)", 0.000511**3/num3)]
res += [("C  alpha", alpha), ("C  alpha/pi", alpha/np.pi), ("C  alpha/(2 pi) (Schwinger)", alpha/(2*np.pi)), ("C  3 alpha/(4 pi)", 3*alpha/(4*np.pi))]
# D: Thomson hits per electron between z = 1090 and z = 124
H0 = 67.4e3/3.0857e22; Om = 0.315; Or = 9.1e-5; OL = 1 - Om - Or
Hz = lambda z: H0*np.sqrt(Om*(1+z)**3 + Or*(1+z)**4 + OL)
rate = lambda z: 411e6*(1+z)**3*6.6524e-29*2.998e8           # photons/m^3 x sigma_T x c
hits = quad(lambda z: rate(z)/((1+z)*Hz(z)), 124, 1090)[0]
res += [("D  Thomson hits per electron (z 1090 -> 124)", hits), ("D  1 / hits", 1/hits), ("D  1 / ln(hits)", 1/np.log(hits))]
res += [("E  1/N, N = 163 (largest allowed)", 1/163), ("E  1/N, N = 6 (smallest)", 1/6)]
cell = 1.6e-35*1.0; lamC = 2.426e-12
res += [("F  cell / Compton wavelength", cell/lamC), ("F  (cell / Compton wavelength)^2", (cell/lamC)**2)]
out = ["ITERATION 56: countable pop quantities vs the electron's needed 0.60% (expectations pre-registered)", "",
       f" HIT window {100*lo:.2f}-{100*hi:.2f}%   (data 0.99 +/- 0.49%)", ""]
for name, v in res:
    hit = "HIT" if lo <= v <= hi else ""
    out.append(f" {name:46s} {v:12.4e}   = {100*v:10.4g} %   {hit}")
# consequences if delta = alpha
def E1(zz):
    a = 1/(1 + zz); g = lambda le: np.exp(2*le) - Om*a**-3 - Or*a**-4 - OL*(a*np.exp(le))**-0.5
    return np.exp(brentq(g, -10, 40))
h = 0.674; rho_DE0 = OL*1.0537e4*h*h; n_e0 = 0.0224*1.878e-29/1.6726e-24*(1 - 0.245/2)
need = lambda zs: rho_DE0*((E1(zs))/(1 + zs))**-0.5/(n_e0*(1 + zs)**3*511e3)
zs_alpha = brentq(lambda zs: need(zs) - alpha, 50, 400)
out += ["", f"IF delta = alpha: payday at z = 124 gives dark energy = {alpha/need(124):.2f} x measured;",
        f"   or, with dark energy as measured, payday at z = {zs_alpha:.0f} (21-cm step at {1420.406/(1+zs_alpha):.1f} MHz);",
        f"   the CMB value 0.99 +/- 0.49% vs alpha = {100*alpha:.3f}%: {(0.99 - 100*alpha)/0.49:+.1f} sigma."]
txt = "\n".join(out); print(txt); open("iter56_pops.txt", "w").write(txt + "\n")
