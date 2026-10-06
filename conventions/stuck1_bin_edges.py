"""Stuck item S1: are the free-bin dark-energy features robust to the choice of bin edges?
For each edge set: best-fit w per bin, and delta-chi2 when that one bin is forced to w = -1 (all else refitted)."""
import os, numpy as np
from scipy.optimize import minimize
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "drill3_binned_w.py")).read().split("res = {}")[0]
EDGES_ENV = [float(x) for x in os.environ["EDGES"].split(",")] + [1e9]
src = src.replace("EDGES = [0.0, 0.4, 0.8, 1.5, 1e9]", "EDGES = EDGES_ENV")
exec(src)
best = fit("BINW", [[0.31, 0.68, 0.0224, -0.9, -1.0, -1.0, -1.0], [0.31, 0.68, 0.0224, -1, -1, -1, -1]])
lines = [f"SNSET={SNSET} edges {EDGES[:-1]}: chi2 {best.fun:.2f}; w = {np.round(best.x[3:], 2)}"]
for i in range(4):
    f = lambda q: c2("BINW", np.r_[q[:3+i], -1.0, q[3+i:]])
    r = minimize(f, np.delete(best.x, 3+i), method="Nelder-Mead", options=dict(maxiter=20000, xatol=1e-7, fatol=1e-6))
    r = minimize(f, r.x, method="Nelder-Mead", options=dict(maxiter=20000, xatol=1e-7, fatol=1e-6))
    d = r.fun - best.fun
    lines.append(f"   bin {EDGES[i]:.1f}-{EDGES[i+1] if EDGES[i+1] < 100 else 'inf'}: w {best.x[3+i]:+.2f}; forcing -1 costs {d:.1f} -> {np.sqrt(max(d,0)):.1f} sigma")
print("\n".join(lines), flush=True)
