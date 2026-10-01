"""
Mosaic network (random Delaunay mesh), flow-thickening tubes. Question: is the fluid held IN the thickened tubes the MOND
halo mass that Khronon's one-fluid bookkeeping needs?
Fluid per tube = (length) x (cross-section); conductance grows with cross-section as A^n:
  n = 2 (pipe flow, Poiseuille), n = 1 (conductance ~ area), n = 1/2.  Thickness factor s = conductance / max -> A ~ s^(1/n).
Phantom (MOND halo) mass inside r, from the measured pull: M_ph(<r) = sigma (g_adaptive - g_static) r^2  (network units).
Test: one universal fluid density C per tube (same for all masses) with M_tubes(<r) = C sum len s^(1/n) must match M_ph(<r)
for every mass (M = 4, 8, 16, 32) over r = 3-13.
Analytic expectation: far out s ~ g/a0 ~ sqrt(M)/r, so tube fluid density ~ (sqrt(M)/r)^(1/n); the MOND halo density is
~ sqrt(M)/r^2. Matching the r-dependence needs n = 1/2, matching the M-dependence needs n = 1 -> no single n works.
"""
import numpy as np
from scipy.optimize import minimize_scalar
from compare_geometries import prepare, solve, profile
from geometries import BUILDERS
P, E, Lng, B, r, bnd = prepare(*BUILDERS["mosaic"](0)); n = len(P); free = ~bnd
rmid = np.linalg.norm(0.5*(P[E[:, 0]] + P[E[:, 1]]), axis=1)
sigma = 3.138                                                               # from calibrate_geometries (static conductivity)
res = {}
for M in (4.0, 8.0, 16.0, 32.0):
    ps, _ = solve(B, 1/Lng, free, n, M); mid, gs = profile(ps, r)
    s = np.ones(len(E)); x = None
    for _ in range(80):
        p, x = solve(B, s/Lng, free, n, M, x); Q = np.abs((s/Lng)*(B @ p)); new = np.clip(np.sqrt(Q), 1e-4, 1)
        ch = np.max(abs(new - s)); s = 0.5*s + 0.5*new
        if ch < 1e-4: break
    _, g = profile(p, r)
    Mph = sigma*(g - gs)*mid**2
    res[M] = (mid, Mph, s)
rs = np.arange(3, 13.5, 1.0)
print("MOND halo mass inside r (network units):")
for M, (mid, Mph, s) in res.items():
    print(f"  M = {M:4.0f}: " + " ".join(f"{np.interp(R_, mid, Mph):7.1f}" for R_ in rs[::2]) + f"   (r = {rs[::2].astype(int).tolist()})")
for nexp in (2.0, 1.0, 0.5):
    curves = {}
    for M, (mid, Mph, s) in res.items():
        cont = np.array([np.sum(Lng[rmid < R_]*s[rmid < R_]**(1/nexp)) - np.sum(Lng[rmid < R_]*1e-4**(1/nexp)) for R_ in rs])
        curves[M] = (cont, np.interp(rs, mid, Mph))
    def cost(lc):
        return np.mean([np.mean((np.log(np.exp(lc)*c) - np.log(np.maximum(ph, 1e-6)))**2) for c, ph in curves.values()])
    b = minimize_scalar(cost, bounds=(-20, 10), method="bounded")
    perM = {M: np.exp(np.median(np.log(np.maximum(ph, 1e-6)/c))) for M, (c, ph) in curves.items()}
    shape = {M: np.polyfit(np.log(rs), np.log(c), 1)[0] - np.polyfit(np.log(rs), np.log(np.maximum(ph, 1e-6)), 1)[0]
             for M, (c, ph) in curves.items()}
    print(f"n = {nexp}: universal-C mismatch rms {np.sqrt(b.fun):.2f} (ln units) | best C per mass " +
          ", ".join(f"{M:.0f}: {v:.3g}" for M, v in perM.items()) +
          " | slope error vs r: " + ", ".join(f"{v:+.2f}" for v in shape.values()))
