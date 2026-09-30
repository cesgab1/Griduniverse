"""
Neighbours stiffen the grid (external field effect) — internal test on SPARC.
Tension rule: slack depends on TOTAL strain |g_internal + g_neighbours|. A galaxy sitting in a neighbour's slope (external pull e*a0)
gets less boost in its outskirts, so its rotation curve SAGS at large radius instead of staying flat.
Model (Chae et al. 2020 analytic form, simple interpolation):  g = g_N [ 1/2 - A/y + sqrt((1/2 - A/y)^2 + B/y) ],  y = g_N/a0,
   A = e(1+e/2)/(1+e), B = 1+e.   e = 0 is the isolated case (the rule already tested, a0 = 1.15e-10).
Per galaxy: fit stellar mass-to-light (log-normal prior 0.5, 0.1 dex; bulge 0.7) with e = 0, then with e free (0-1).
Measured: how many galaxies prefer e > 0, the summed improvement, the typical e, and whether e rises where the curve is measured far out.
"""
import numpy as np, glob, os, json
from scipy.optimize import minimize
a0 = 1.15e-10; kpc = 3.0857e19
def model(r, Vg, Vd, Vb, Yd, e):
    gN = (Vg*np.abs(Vg) + Yd*Vd**2 + 1.4*Yd*Vb**2)*1e6/(r*kpc)
    gN = np.clip(gN, 1e-16, None); y = gN/a0
    if e <= 0: nu = 0.5 + np.sqrt(0.25 + 1/y)
    else:
        A = e*(1 + e/2)/(1 + e); B = 1 + e
        nu = 0.5 - A/y + np.sqrt((0.5 - A/y)**2 + B/y)
    return np.sqrt(gN*nu*r*kpc)/1e3
res = []
for f in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(f, comments="#"); d = d[None] if d.ndim == 1 else d
    r, V, eV, Vg, Vd, Vb = d.T[:6]
    ok = (V > 0) & (eV > 0)
    if ok.sum() < 8: continue
    r, V, eV, Vg, Vd, Vb = r[ok], V[ok], np.hypot(eV[ok], 0.03*V[ok]), Vg[ok], Vd[ok], Vb[ok]
    def chi(p, e):
        Yd = 10**p[0]; return np.sum(((V - model(r, Vg, Vd, Vb, Yd, e))/eV)**2) + ((p[0] - np.log10(0.5))/0.1)**2
    f0 = minimize(lambda p: chi(p, 0.0), [np.log10(0.5)], method="Nelder-Mead")
    best = None
    for e0 in (0.01, 0.05, 0.15, 0.4):
        fe = minimize(lambda p: chi(p[:1], abs(p[1])) + (1e6 if abs(p[1]) > 1 else 0), [f0.x[0], e0], method="Nelder-Mead")
        if best is None or fe.fun < best.fun: best = fe
    gN_out = (Vg[-1]*abs(Vg[-1]) + 0.5*Vd[-1]**2 + 0.7*Vb[-1]**2)*1e6/(r[-1]*kpc)/a0
    res.append(dict(name=os.path.basename(f)[:-11], n=int(ok.sum()), chi0=float(f0.fun), chie=float(best.fun), e=float(abs(best.x[1])), yout=float(gN_out)))
dchi = np.array([x["chi0"] - x["chie"] for x in res]); e = np.array([x["e"] for x in res]); yout = np.array([x["yout"] for x in res])
N = len(res)
print(f"{N} galaxies with >= 8 points")
print(f"prefer a neighbour pull (Delta chi2 > 1): {np.mean(dchi > 1)*100:.0f}%;  strongly (> 4): {np.mean(dchi > 4)*100:.0f}%;  (> 9): {np.mean(dchi > 9)*100:.0f}%")
print(f"summed improvement {dchi.sum():.0f} for {N} extra parameters;  AIC-type net {dchi.sum() - 2*N:.0f}")
print(f"fitted neighbour pull e = g_ext/a0: median {np.median(e):.3f}, 16-84% {np.percentile(e,16):.3f}-{np.percentile(e,84):.3f}")
far = yout < 0.1; near = yout > 0.3
print(f"galaxies measured deep into the weak-pull region (outer g_N < 0.1 a0): {far.sum()}, prefer e>0 (dchi2>4): {np.mean(dchi[far] > 4)*100:.0f}%, median e {np.median(e[far]):.3f}")
print(f"galaxies whose curve stops early (outer g_N > 0.3 a0): {near.sum()}, prefer e>0: {np.mean(dchi[near] > 4)*100:.0f}%, median e {np.median(e[near]):.3f}")
print("expected for comparison (Chae et al. 2020/2021): e ~ 0.02-0.05 typical, larger for galaxies in crowded environments")
top = sorted(res, key=lambda x: -(x["chi0"] - x["chie"]))[:8]
print("strongest preferences:", ", ".join(f"{x['name']} (e={x['e']:.2f}, dchi2={x['chi0']-x['chie']:.0f})" for x in top))
json.dump(res, open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/neighbour_test.json", "w"), indent=1)
