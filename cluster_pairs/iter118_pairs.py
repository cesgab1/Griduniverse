"""Iteration 118: same-size clusters/groups vs 'visible matter alone sets gravity' (PREREG_118.md)."""
import os, numpy as np, pandas as pd
from scipy.optimize import minimize_scalar, minimize
here = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(os.path.join(here, "clusters.csv"), comment="#")
G, MS, KPC = 6.674e-11, 1.989e30, 3.086e19
nu = lambda x: 1/(1 - np.exp(-np.sqrt(x)))
out = ["Iteration 118 -- 31 clusters and groups (Vikhlinin+2006, Sun+2009)"]
for label, fstar in [("stars f* = 0.025 (M/1e14)^-0.37", lambda M: 0.025*(M)**-0.37), ("stars f* = 0", lambda M: 0*M), ("stars f* = 0.02", lambda M: 0.02 + 0*M)]:
    fs = fstar(d.M500.values); fv = d.fg.values + fs
    B = 1/fv; efg = (d.fg_lo + d.fg_hi).values/2; sB = efg/fv/np.log(10)               # error on log10 B
    r = d.r500.values*KPC; gobs = G*d.M500.values*1e14*MS/r**2; gbar = gobs*fv
    sM = ((d.M_lo + d.M_hi)/2/d.M500).values/np.log(10); sg = np.hypot(sM, 0)           # error on log g_obs (mass)
    # T1: best acceleration scale
    def c2(la, sint=0.0):
        p = np.log10(gbar*nu(gbar/10**la)); return np.sum((np.log10(gobs) - p)**2/(sg**2 + sint**2))
    la = minimize_scalar(c2, bounds=(-12, -7), method="bounded").x
    grid = np.linspace(-12, -7, 2001); cc = np.array([c2(x) for x in grid]); ok = grid[cc <= cc.min() + 1]
    off = np.median(gobs/(gbar*nu(gbar/1.2e-10)))
    # T2: intrinsic scatter of log B at fixed size: fit log B = p0 + p1 log r500 + intrinsic scatter (max likelihood)
    x = np.log10(d.r500.values/700); y = np.log10(B)
    def nll(p):
        a, b, ls = p; v = sB**2 + np.exp(2*ls); return np.sum((y - a - b*x)**2/v + np.log(v))
    best = minimize(nll, [1, -0.5, np.log(0.05)], method="Nelder-Mead", options=dict(xatol=1e-8, fatol=1e-8, maxiter=20000))
    zero = minimize(lambda p: nll([p[0], p[1], -20]), best.x[:2], method="Nelder-Mead").fun
    sig_int = np.exp(best.x[2]); dchi = zero - best.fun
    # T3: slope of log f_vis vs log r500 = -slope of log B
    out += ["", f"[{label}]",
            f"  boost B = total / visible: median x{np.median(B):.1f}, range x{B.min():.1f}-x{B.max():.1f}",
            f"  T1 best acceleration scale a = {10**la:.2e} m/s^2 (1-sigma {10**ok.min():.1e}-{10**ok.max():.1e}) vs galaxies 1.2e-10 -> x{10**la/1.2e-10:.0f};"
            f" with galaxies' a the rule under-predicts total pull by x{off:.2f} (median)",
            f"  T2 same size: intrinsic scatter of log B beyond errors = {sig_int:.3f} dex (x{10**sig_int:.2f}); vs zero: delta(-2lnL) = {dchi:.1f}"
            f" -> {np.sqrt(max(dchi,0)):.1f} sigma",
            f"  T3 visible fraction vs size: slope {-best.x[1]:.2f} (rule predicts ~0.7-1)"]
    if label.startswith("stars f* = 0.025"):
        out.append("  bins of size (r500): " + "; ".join(
            f"{lo}-{hi} kpc: N {s.sum()}, B x{np.median(B[s]):.1f} ({B[s].min():.1f}-{B[s].max():.1f})"
            for lo, hi in [(350, 500), (500, 650), (650, 800), (800, 1200), (1200, 1500)] for s in [(d.r500.values >= lo) & (d.r500.values < hi)] if s.sum()))
txt = "\n".join(out); print(txt); open(os.path.join(here, "iter118_pairs.txt"), "w").write(txt + "\n")
