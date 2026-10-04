"""
ITERATION 13 (lesson of iteration 11 applied): dark energy = constant BASELINE (size set elsewhere, e.g. at the bounce)
+ a JOSTLED part (memory kappa = 3, quadratic energy: the best zero-parameter toy of iteration 4). What fraction f can
the jostled part be?  rho_DE(a) = OL [ (1 - f) + f * rho_toy(a)/rho_toy(1) ].
Also the multiplicative alternative is just f = 1 with a free normalisation (that is the toy as fitted before).
Data: DESI DR2 BAO + Planck distance priors + each SN set (fit_law machinery).
"""
import numpy as np, os, sys, json
os.environ.setdefault("SNSET", "PANTHEON")
src = open("iter1_memory.py").read().split("# ---------- Part A")[0]
g = {"__name__": "x", "__file__": os.path.abspath("iter1_memory.py")}; exec(src, g)
fl = g["g"]; ZG, W_R = fl["ZG"], fl["W_R"]
F = float(os.environ.get("FRAC", "1"))
def E_mix(model, Om, h, extra):
    Or = W_R/h**2; OL = 1 - Om - Or; zp = 1 + ZG; zi = ZG[ZG <= 30]; a = 1/(1 + zi)
    r = g["rho_hist"](Om, h, 3.0, "quad", a); r = r/r[0]                 # shape normalised to today (a[0] = 1)
    rde = np.full_like(ZG, OL*(1 - F)); rde[:len(zi)] = OL*((1 - F) + F*r)
    return np.sqrt(Om*zp**3 + Or*zp**4 + rde)
L = fl["best"]("LCDM", [[0.31, 0.68, 0.0224]])           # LCDM with the ORIGINAL expansion code (fix: was patched first)
fl["E_of_z"] = E_mix
r = fl["best"]("MIX", [[0.31, 0.68, 0.0224]])
print(f"RESULT {os.environ['SNSET']} f={F}: dchi2 vs LCDM {r.fun - L.fun:+.2f}")
