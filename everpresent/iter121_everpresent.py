"""Iteration 121: simplest everpresent Lambda (Lambda = W(V)/V, random walk over past 4-volume) vs data (PREREG_121)."""
import numpy as np
from scipy.integrate import quad
Om, Or = 0.31, 9.0e-5; OL = 1 - Om - Or
E = lambda z: np.sqrt(Om*(1+z)**3 + Or*(1+z)**4 + OL)
t = lambda z: quad(lambda x: 1/((1+x)*E(x)), z, np.inf, limit=500)[0]       # age in 1/H0
zs = [1100, 1.5, 0.8, 0.4, 0.0]
V = np.array([t(z)**4 for z in zs])                                        # past 4-volume ~ t^4 (arbitrary units)
rng = np.random.default_rng(121); N = 200000
dV = np.diff(np.r_[0, V]); W = np.cumsum(rng.standard_normal((N, len(V)))*np.sqrt(dV), axis=1)
L = W/V; r = L/L[:, -1:]                                                   # dark-energy density relative to today
a = L[:, -1] > 0
b = np.all((r[:, 1:4] > 0.5) & (r[:, 1:4] < 2), axis=1)
c = np.abs(r[:, 0])*(OL/Om)/1101**3 < 0.05
out = ["Iteration 121 -- simplest everpresent Lambda (Sorkin) vs data; 200000 realisations",
       f"age ratios t(z)/t0: " + ", ".join(f"z={z}: {t(z)/t(0):.3g}" for z in zs[:-1]),
       f"typical |dark energy(z)| / today (median of |ratio|): " + ", ".join(f"z={z}: {np.median(np.abs(r[:,i])):.3g}" for i, z in enumerate(zs[:-1])),
       f"(a) positive today: {a.mean():.3f}",
       f"(b) 0.5-2x today at z = 0.4, 0.8, 1.5: {b.mean():.4f}  (given (a): {b[a].mean():.4f})",
       f"(c) < 5% of matter at z = 1100: {c.mean():.4f}",
       f"(a)+(b)+(c) together: {(a & b & c).mean():.5f}   (constant Lambda: 1)",
       f"needed damping of past fluctuations for (c) in half of realisations: early jitter must be smaller by x{np.median(np.abs(r[:,0])*(OL/Om)/1101**3)/0.05:.0f}"]
txt = "\n".join(out); print(txt); open("iter121_everpresent.txt", "w").write(txt + "\n")
