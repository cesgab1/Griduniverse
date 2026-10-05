import os, json, numpy as np
from scipy.optimize import minimize
os.chdir("/home/claude/griduniverse/cosmology_fits")
src = open("fit_law.py").read().split("Ndata = len(bao)")[0]
src = src.replace("mu = 5 * np.log10((1 + zhel) * o[\"DM\"](z_sn)) + 25", "mu = 5 * np.log10((1 + zhel) * o[\"DM\"](z_sn) * (1 + z_sn)**EPS) + 25")
assert "EPS" in src
EPS = 0.0
exec(src)
MODEL = os.environ.get("MODEL", "LCDM")
out = {}
for eps in (None, 0.0):
    def f(p):
        global EPS
        if eps is None: EPS = p[3]; q = p[:3]
        else: EPS = 0.0; q = p
        return chi2(MODEL, q)
    starts = [[0.31, 0.68, 0.0224, 0.0], [0.30, 0.69, 0.0224, 0.05], [0.32, 0.67, 0.0224, -0.05]] if eps is None else [[0.31, 0.68, 0.0224]]
    best = min((minimize(f, s, method="Nelder-Mead", options=dict(maxiter=8000, xatol=1e-7, fatol=1e-6)) for s in starts), key=lambda r: r.fun)
    out["free" if eps is None else "zero"] = dict(chi2=float(best.fun), x=[float(v) for v in best.x])
# width of eps from a profile
bx = out["free"]["x"]; grid = np.linspace(bx[3] - 0.08, bx[3] + 0.08, 9); prof = []
for e in grid:
    EPS = e
    r = minimize(lambda q: chi2(MODEL, q), bx[:3], method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-7, fatol=1e-6))
    prof.append(r.fun)
a2 = np.polyfit(grid, prof, 2)[0]; out["eps_sigma"] = float(1/np.sqrt(a2))
print("RESULT", SNSET, MODEL, json.dumps(out))
