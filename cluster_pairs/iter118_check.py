"""Iteration 118 check (post-hoc, flagged): same-size scatter within one team's sample only (Sun+2009 groups),
and how big unmodelled star-mass variation would have to be to explain it."""
import os, numpy as np, pandas as pd
from scipy.optimize import minimize
d = pd.read_csv("clusters.csv", comment="#"); d = d[d["sample"] == "S09"]
fs = 0.025*d.M500.values**-0.37; fv = d.fg.values + fs; B = 1/fv
sB = ((d.fg_lo + d.fg_hi)/2).values/fv/np.log(10); x = np.log10(d.r500.values/700); y = np.log10(B)
def nll(p):
    a, b, ls = p; v = sB**2 + np.exp(2*ls); return np.sum((y - a - b*x)**2/v + np.log(v))
best = minimize(nll, [1, -0.5, np.log(0.05)], method="Nelder-Mead", options=dict(xatol=1e-8, fatol=1e-8, maxiter=20000))
zero = minimize(lambda p: nll([p[0], p[1], -20]), best.x[:2], method="Nelder-Mead").fun
s = np.exp(best.x[2])
out = [f"Sun+2009 groups only (N={len(d)}): intrinsic scatter {s:.3f} dex, vs zero delta(-2lnL) {zero-best.fun:.1f} -> {np.sqrt(max(zero-best.fun,0)):.1f} sigma",
       f"star-mass variation needed to produce it alone: visible-fraction scatter {np.log(10)*s*100:.0f}% -> stars (~{np.median(fs/fv)*100:.0f}% of visible) would have to vary by ~{np.log(10)*s/np.median(fs/fv)*100:.0f}% group to group"]
print("\n".join(out)); open("iter118_check.txt", "w").write("\n".join(out) + "\n")
