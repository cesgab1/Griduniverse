"""
ITERATION 46 (route 1 premise; Coalesce: 'time has to be its own dimension'): does a grid with no built-in 'now' look the same to
every moving observer, while a grid laid down in synchronised ticks does not?
Two 1+1-D grids with the same density: (S) points scattered randomly in space AND time (time a dimension like space);
(T) points scattered randomly in space but laid down on synchronised time ticks (built-in 'now').
Test: boost the grid by speed v (rapidity 1, 2) and measure, in a fixed window, the distribution of the time gaps between each point
and its nearest point in the future light cone (the grid's own 'tick' as seen by that observer). Compare with the unboosted
distribution (two-sample Kolmogorov-Smirnov).
Expectations written BEFORE running: (S) statistically identical after boosts (KS p > 0.05); (T) clearly different (p < 1e-6).
"""
import numpy as np
from scipy.stats import ks_2samp
rng = np.random.default_rng(46)
def sprinkle(kind, rho=400.0, T=30.0, X=30.0):
    n = rng.poisson(rho*T*X); t = rng.uniform(-T/2, T/2, n); x = rng.uniform(-X/2, X/2, n)
    if kind == "T": t = np.round(t*np.sqrt(rho))/np.sqrt(rho)          # synchronised ticks, same density
    return t, x
def boost(t, x, eta): c, s = np.cosh(eta), np.sinh(eta); return c*t - s*x, c*x - s*t
def gaps(t, x, W=3.0):
    sel = np.where((np.abs(t) < W) & (np.abs(x) < W))[0]; o = np.argsort(t); ts, xs = t[o], x[o]; out = []
    for i in sel[:3000]:
        j0 = np.searchsorted(ts, t[i], side="right"); dt = ts[j0:j0 + 4000] - t[i]; dx = np.abs(xs[j0:j0 + 4000] - x[i])
        ok = (dt > dx) & (dt > 0)
        if ok.any(): out.append(np.sqrt((dt[ok]**2 - dx[ok]**2)).min())    # proper time to nearest future-cone point
    return np.array(out)
lines = ["ITERATION 46: does the grid have a preferred 'now'? (expectations committed before running)", ""]
for kind, name in (("S", "random in space AND time"), ("T", "random in space, synchronised ticks")):
    t, x = sprinkle(kind); g0 = gaps(t, x)
    for eta in (1.0, 2.0):
        tb, xb = boost(t, x, eta); g1 = gaps(tb, xb); p = ks_2samp(g0, g1).pvalue
        lines.append(f"{name:38s} boost rapidity {eta}: median gap {np.median(g0):.4f} -> {np.median(g1):.4f}, KS p = {p:.2e}")
txt = "\n".join(lines); print(txt); open("iter46_boost_test.txt", "w").write(txt + "\n")
