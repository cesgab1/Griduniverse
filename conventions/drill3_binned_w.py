"""Conventions drill D3 (exploratory): dark energy WITHOUT the w0-wa (CPL) straight-line convention.
w(z) free in 4 redshift bins (0-0.4, 0.4-0.8, 0.8-1.5, >1.5); fit BAO (DESI DR2) + SN (SNSET) + CMB priors.
Compare: LCDM, w0-wa, our law (beta 1/2), and the law's own effective w in the same bins."""
import os, numpy as np
from scipy.optimize import minimize
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../cosmology_fits/fit_law.py")).read()
exec(src.split("Ndata =")[0])
EDGES = [0.0, 0.4, 0.8, 1.5, 1e9]
_E = E_of_z
def lnrho_binned(ws):
    lz = np.log1p(ZG); out = np.zeros_like(ZG)
    for w, a, b in zip(ws, EDGES[:-1], EDGES[1:]):
        out += 3*(1 + w)*(np.clip(lz, np.log1p(a), np.log1p(b)) - np.log1p(a))
    return out
def E_of_z(model, Om, h, extra):
    if model != "BINW": return _E(model, Om, h, extra)
    Or = W_R/h**2; OL = 1 - Om - Or; zp = 1 + ZG
    return np.sqrt(Om*zp**3 + Or*zp**4 + OL*np.exp(lnrho_binned(extra)))
def c2(model, p):
    if model in ("BINW", "w0wa") and np.any(np.abs(np.array(p[3:])) > 4): return 1e12
    return chi2(model, p)
def fit(model, starts):
    b = None
    for x0 in starts:
        r = minimize(lambda p: c2(model, p), x0, method="Nelder-Mead", options=dict(maxiter=20000, maxfev=40000, xatol=1e-6, fatol=1e-5))
        r = minimize(lambda p: c2(model, p), r.x, method="Nelder-Mead", options=dict(maxiter=20000, maxfev=40000, xatol=1e-7, fatol=1e-6))
        if b is None or r.fun < b.fun: b = r
    return b
res = {}
res["LCDM"] = fit("LCDM", [[0.31, 0.68, 0.0224]])
res["LAW"] = fit("LAW", [[0.31, 0.68, 0.0224]])
res["w0wa"] = fit("w0wa", [[0.31, 0.68, 0.0224, -0.8, -0.6], [0.30, 0.69, 0.0224, -0.95, -0.1]])
res["BINW"] = fit("BINW", [[0.31, 0.68, 0.0224, -0.9, -1.0, -1.0, -1.0], [0.31, 0.68, 0.0224, -1, -1, -1, -1]])
out = [f"Drill 3 -- dark energy shape without the w0-wa convention (SNSET={SNSET})"]
for m, r in res.items():
    k = len(r.x); out.append(f"  {m:5s} chi2 {r.fun:9.2f}  (vs LCDM {r.fun - res['LCDM'].fun:+7.2f}; extra parameters {k-3})  params {np.round(r.x, 3)}")
# 1-sigma errors on binned w: profile each bin
wb = res["BINW"].x
errs = []
for i in range(4):
    def prof(wi):
        f = lambda q: c2("BINW", np.r_[q[:3+i], wi, q[3+i:]])
        x0 = np.delete(wb, 3+i); return minimize(f, x0, method="Nelder-Mead", options=dict(maxiter=8000, xatol=1e-6, fatol=1e-5)).fun
    grid = wb[3+i] + (np.linspace(-0.25, 0.25, 51) if i < 2 else np.linspace(-0.8, 0.8, 65)); pr = np.array([prof(g) for g in grid])
    assert abs(prof(wb[3+i]) - res['BINW'].fun) < 0.05, 'profile not converged'
    ok = grid[pr <= res["BINW"].fun + 1]; errs.append((ok.min(), ok.max()))
# law's effective w in the same bins
Om, h = res["LAW"].x[:2]; E = _E("LAW", Om, h, []); Or = W_R/h**2; zp = 1 + ZG
rde = E**2 - Om*zp**3 - Or*zp**4; sel = ZG < 3
weff = -1 + np.gradient(np.log(np.clip(rde[sel], 1e-30, None)), np.log(zp[sel]))/3
out.append("  binned w (1-sigma, profiled)  |  our law's own effective w  |  LCDM -1")
for i, (a, b) in enumerate(zip(EDGES[:-1], EDGES[1:])):
    s = (ZG[sel] >= a) & (ZG[sel] < min(b, 3))
    out.append(f"   z {a:.1f}-{min(b,3):.1f}: w = {wb[3+i]:+.2f} [{errs[i][0]:+.2f}, {errs[i][1]:+.2f}]   law {np.mean(weff[s]):+.2f}   "
               f"-1 is {'inside' if errs[i][0] <= -1 <= errs[i][1] else 'OUTSIDE'} 1 sigma")
txt = "\n".join(out); print(txt); open(f"drill3_binned_w_{SNSET}.txt", "w").write(txt + "\n")
