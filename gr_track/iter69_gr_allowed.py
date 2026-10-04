"""ITERATION 69 (pre-registered in PREREG_69.md): GR-allowed (non-crossing) dark-energy histories vs today's data."""
import os, sys, numpy as np
from scipy.optimize import minimize
os.chdir("/home/claude/griduniverse/cosmology_fits")
src = open("fit_law.py").read().split("Ndata = len(bao)")[0]
exec(src)
def c2(x):   # x = Om, h, ombh2, w0, wend  (wend = w0 + wa)
    Om, h, ob, w0, we = x
    if w0 < -1 or we < -1: return 1e9
    return chi2("w0wa", [Om, h, ob, w0, we - w0])
best = None
for st in ([0.31, 0.68, 0.0224, -0.95, -0.95], [0.31, 0.68, 0.0224, -0.8, -1.0], [0.30, 0.69, 0.0224, -1.0, -0.8], [0.32, 0.66, 0.0224, -0.7, -0.9]):
    r = minimize(c2, st, method="Nelder-Mead", options=dict(maxiter=6000, xatol=1e-6, fatol=1e-6))
    if best is None or r.fun < best.fun: best = r
lc = min(minimize(lambda p: chi2("LCDM", p), s, method="Nelder-Mead", options=dict(maxiter=4000)).fun for s in ([0.31, 0.68, 0.0224], [0.30, 0.69, 0.0223]))
print(f"RESULT {SNSET}: LCDM {lc:.2f}  GR-allowed w>=-1: {best.fun:.2f} (delta {best.fun - lc:+.2f}) at w0 = {best.x[3]:.3f}, w(past) = {best.x[4]:.3f}")
