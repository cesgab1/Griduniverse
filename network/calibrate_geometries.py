"""
Fairness follow-up: each geometry carries flux differently per link, so the same link rule gives each its own effective a0.
Calibrate per geometry, then compare like with like:
  sigma   : effective conductivity from the static run (g = M / (sigma r^2))
  a0_eff  : fitted so that the adaptive profile best matches MOND (hard switch: g = gN above a0, sqrt(gN a0) below)
  fit     : rms of ln g about that MOND curve over r = 2-13, for two masses (shape quality; 0 = exact MOND)
  inner   : adaptive / static pull inside half the MOND radius (r = 2-4), mass chosen so r_M ~ 8 (Newton inside: 1.00).
            (A slope is not used: shell-averaging next to a point source on a discrete network distorts it in both runs.)
"""
import numpy as np
from scipy.optimize import minimize_scalar
from compare_geometries import run
def mond(r, M, sig, a):
    gN = M/(sig*r**2); return np.where(gN >= a, gN, np.sqrt(gN*a))
for geom in ("cubic", "mosaic", "web", "fungal"):
    mid, g, _, _ = run(geom, "static", 10.0); m = (mid >= 3) & (mid <= 12) & (g > 0)
    sig = np.exp(np.median(np.log(10.0/(g[m]*mid[m]**2))))
    prof = {M: run(geom, "adaptive", M)[:2] for M in (8.0, 64.0)}
    def cost(la):
        tot = 0
        for M, (mid, g) in prof.items():
            mm = (mid >= 2) & (mid <= 13) & (g > 0); tot += np.mean((np.log(g[mm]) - np.log(mond(mid[mm], M, sig, np.exp(la))))**2)
        return tot/len(prof)
    b = minimize_scalar(cost, bounds=(-12, 4), method="bounded"); a = np.exp(b.x); rms = np.sqrt(b.fun)
    Mbig = 64*a*sig                                                       # r_M = sqrt(M/(sig a)) = 8
    mid, g, _, _ = run(geom, "adaptive", Mbig); _, gs, _, _ = run(geom, "static", Mbig); mm = (mid >= 2) & (mid <= 4)
    inner = np.mean(g[mm]/gs[mm])
    mo = (mid >= 10) & (mid <= 13.5) & (g > 0); outer = np.polyfit(np.log(mid[mo]), np.log(g[mo]), 1)[0]
    print(f"{geom:7s} sigma {sig:6.3f} | a0_eff {a:8.2e} | MOND shape rms {rms:5.2f} | with r_M = 8: inner adaptive/static {inner:5.2f}, outer slope {outer:5.2f}")
