"""Iteration 74: dark energy's share at the Bang's own clock moments. See PREREG_74.md."""
import numpy as np
from scipy.optimize import brentq
ODE, Ok, h = 0.685, 0.0023, 0.674
Org = 2.469e-5*(1 + 0.2271*3.046)/h**2       # radiation today (photons + neutrinos)
Om = 1 - ODE - Ok - Org
T0K, kB_GeV = 2.7255, 8.617e-14              # K, GeV per K
def gstar(TGeV):
    if TGeV > 160: return 106.75, 106.75
    if TGeV > 0.15: return 61.75, 61.75
    if TGeV > 5e-4: return 10.75, 10.75
    return 3.36, 3.91
def a_of_T(TGeV):                           # entropy conservation
    gs = gstar(TGeV)[1]; return (T0K*kB_GeV/TGeV)*(3.91/gs)**(1/3)
def Orad(a):                                # radiation density / critical today, with g* steps
    T = T0K*kB_GeV/a                         # rough T for choosing the step (iterate once)
    for _ in range(3):
        g, gs = gstar(T); T = T0K*kB_GeV/a*(3.91/gs)**(1/3)
    g, gs = gstar(T)
    return Org*(g/3.36)*(3.91/gs)**(4/3)/a**4
def E_and_share(a):
    base = Om/a**3 + Orad(a) + Ok/a**2
    lnE = brentq(lambda l: np.exp(2*l) - base - ODE*(np.exp(l)*a)**-0.5, -100, 200)
    E = np.exp(lnE); rde = ODE*(E*a)**-0.5
    return E, rde/E**2, rde
S_today = 3.3e122
TP = 1.2209e19                               # Planck energy, GeV
aeq = brentq(lambda la: Om/np.exp(la)**3 - Orad(np.exp(la)), np.log(1e-6), np.log(1e-2)); aeq = np.exp(aeq)
moments = [("Planck time", a_of_T(TP)), ("electroweak 160 GeV", a_of_T(160)), ("QCD 0.15 GeV", a_of_T(0.15)),
           ("first 3 minutes 0.07 MeV", a_of_T(7e-5)), ("matter-radiation equality", aeq), ("transparency z=1090", 1/1091)]
rhoc0_planck = 3*(67.4e3/3.0857e22)**2/(8*np.pi*6.674e-11) / (2.998e8**5/(1.0546e-34*6.674e-11**2))
out = ["Iteration 74: dark energy's share at the Bang's clock moments (fading law extrapolated back). PREREG_74.md", ""]
hits = []
for name, a in moments:
    E, sh, rde = E_and_share(a)
    t1 = abs(np.log10(sh/0.5)) < np.log10(3); t2 = abs(np.log10(sh*S_today)) < np.log10(3)
    if t1: hits.append((name, "T1"))
    if t2: hits.append((name, "T2"))
    out.append(f"{name:27s} a = {a:.2e}  share = {sh:.2e}  rho_DE = {rde*rhoc0_planck:.2e} Planck density  "
               f"| T1 {'HIT' if t1 else 'miss'} | T2 {'HIT' if t2 else 'miss'} (share x 3.3e122 = {sh*S_today:.1e})")
out.append("")
out.append(f"Today: share 0.685, rho_DE = {ODE*rhoc0_planck:.2e} Planck density. Early growth law check: share ~ a^4.5 in radiation era.")
out.append("Hits: " + (", ".join(f"{n} ({t})" for n, t in hits) if hits else "none"))
txt = "\n".join(out); print(txt); open("iter74_bang_dial.txt", "w").write(txt + "\n")
