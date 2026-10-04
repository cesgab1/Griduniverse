"""
How much dark energy is there today, and how much is left later, under three readings fitted to the SAME data
(DESI DR2 BAO + Planck distance priors + supernovae): Lambda (constant), Claim 1 (our law), and DESI-style w0wa (a fitting formula,
NOT a law -- extrapolating it into the future is unreliable). Concentration only; nothing is 'spent'.
"""
import numpy as np
from scipy.optimize import brentq
G = 6.674e-11; Mpc = 3.0857e22; Gyr = 3.156e16; mH = 1.674e-27
models = {"Lambda (Pantheon+)": ("L", 0.3016, 0.6845, None),
          "Claim 1 (Pantheon+)": ("C", 0.31, 0.6749, None),
          "w0wa Pantheon+": ("W", 0.3107, 0.6765, (-0.8563, -0.4997)),
          "w0wa DES-Dovekie": ("W", 0.3127, 0.6749, (-0.8242, -0.6194)),
          "w0wa Union3": ("W", 0.3262, 0.6612, (-0.6908, -0.9673))}
def rde_ratio(kind, Om, OL, a, h_of_a=None, w=None):
    if kind == "L": return 1.0
    if kind == "W": w0, wa = w; return a**(-3*(1 + w0 + wa))*np.exp(-3*wa*(1 - a))
def H_of_a(kind, Om, OL, a, w=None):
    Or = 9.1e-5
    if kind == "C":
        f = lambda e: e*e - Om*a**-3 - Or*a**-4 - OL*(a*e)**-0.5
        e = brentq(f, 1e-9, 10*np.sqrt(Om*a**-3 + Or*a**-4) + 10); return e, OL*(a*e)**-0.5
    r = OL*rde_ratio(kind, Om, OL, a, w=w); return np.sqrt(Om*a**-3 + Or*a**-4 + r), r
out = ["How much dark energy: today, and later (concentration relative to today)", ""]
for name, (k, Om, h, w) in models.items():
    OL = 1 - Om - 9.1e-5; rc = 3*(h*100e3/Mpc)**2/(8*np.pi*G); today = OL*rc
    # time for a: integrate dt = da/(a H)
    agrid = np.geomspace(1e-6, 10, 40000); Hs = np.array([H_of_a(k, Om, OL, x, w)[0] for x in agrid])
    t = np.concatenate([[0], np.cumsum(np.diff(agrid)/(agrid[1:]*Hs[1:]))])*977.8/(100*h)
    t0 = np.interp(1, agrid, t)
    row = f"{name:22s} today {today:.2e} kg/m^3 (~{today/mH:.1f} H-atom masses/m^3); left at "
    for tt in (25, 40, 60):
        a = np.interp(tt, t, agrid)
        row += f"{tt} Gyr: {H_of_a(k, Om, OL, a, w)[1]/OL:.2f}  "
    out.append(row + f"(t0 = {t0:.1f} Gyr)")
out += ["", "Reading: today's concentration: 5.5-6.2e-27 kg/m^3 across readings (each reading individually ~2-3%).",
        "What is LEFT later depends on the reading: Lambda keeps all of it; Claim 1 fades gently; the w0wa formula, if trusted",
        "beyond the data, fades fast. Only the fitted past is measured; the future part is a projection."]
txt = "\n".join(out); print(txt); open("how_much_left.txt", "w").write(txt + "\n")
