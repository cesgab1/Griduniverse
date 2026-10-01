"""
Two extra same-fashion checks:
 1. Systematic direction preference: fit the shell potential at r = 10 (adaptive, M = 8) with the cubic-symmetry pattern
    K4 = (x^4 + y^4 + z^4)/r^4 - 3/5. Random scatter averages away on large scales; a systematic pattern does not.
    Reported as the K4 amplitude in units of the potential drop per unit radius.
 2. A second mass that is not at the network's centre (source at (6,0,0)): outer slope measured around the source
    (r = 3-8 from it). Only the web has a built-in centre; the others should not care.
"""
import numpy as np
from compare_geometries import prepare, solve
from geometries import BUILDERS
def adapt(B, Lng, free, n, src, M, it=60):
    s = np.ones(len(Lng)); x = None
    for _ in range(it):
        Lap_cond = s/Lng; p, x = solve_src(B, Lap_cond, free, n, M, src, x)
        Q = np.abs(Lap_cond*(B @ p)); new = np.clip(np.sqrt(Q), 1e-4, 1); ch = np.max(abs(new - s)); s = 0.5*s + 0.5*new
        if ch < 1e-4: break
    return p
import scipy.sparse as sp, scipy.sparse.linalg as spl
def solve_src(B, cond, free, n, M, src, x0):
    Lap = (B.T @ sp.diags(cond) @ B).tocsr(); A = Lap[free][:, free]; s = np.zeros(n); s[src] = 4*np.pi*M
    x, _ = spl.cg(A, s[free], x0=x0, rtol=1e-9, maxiter=40000, M=sp.diags(1/A.diagonal())); p = np.zeros(n); p[free] = x; return p, x
for geom in ("cubic", "mosaic", "web", "fungal"):
    P, E, Lng, B, r, bnd = prepare(*BUILDERS[geom](0)); n = len(P); free = ~bnd
    p = adapt(B, Lng, free, n, 0, 8.0)
    sh = (r > 9.5) & (r < 10.5); u = P[sh]/r[sh, None]; K4 = (u**4).sum(1) - 0.6
    sh2 = (r > 10.5) & (r < 11.5); gdrop = p[sh].mean() - p[sh2].mean()
    A = np.vstack([np.ones(sh.sum()), K4]).T; coef = np.linalg.lstsq(A, p[sh], rcond=None)[0][1]
    src = np.argmin(np.linalg.norm(P - np.array([6.0, 0, 0]), axis=1)); free2 = free.copy(); free2[src] = True
    p2 = adapt(B, Lng, free2, n, src, 2.0); d = np.linalg.norm(P - P[src], axis=1)
    bins = np.arange(2.5, 8.6, 1.0); mid = 0.5*(bins[1:] + bins[:-1])
    pm = np.array([p2[(d >= a) & (d < b)].mean() for a, b in zip(bins[:-1], bins[1:])]); g = -np.gradient(pm, mid)
    ok = g > 0; sl = np.polyfit(np.log(mid[ok]), np.log(g[ok]), 1)[0] if ok.sum() > 1 else np.nan
    print(f"{geom:7s} cubic-pattern amplitude {coef/gdrop:+.3f} (radial units) | off-centre mass: outer slope {sl:5.2f}")
