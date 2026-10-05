"""Iteration 82: dark energy slice by slice. See PREREG_82.md."""
import os, json, numpy as np
from scipy.optimize import minimize, brentq
os.chdir("/home/claude/griduniverse/cosmology_fits")
exec(open("fit_law.py").read().split("Ndata = len(bao)")[0])
NODES = np.array([0.0, 0.4, 0.8, 1.4, 2.5])
_E_old = E_of_z
def E_of_z(model, Om, h, extra):
    if model != "BIN": return _E_old(model, Om, h, extra)
    Or = W_R/h**2; OL = 1 - Om - Or; zp = 1 + ZG
    vals = np.r_[1.0, extra]
    sh = np.interp(ZG, NODES, vals)          # beyond 2.5 held at last value
    if np.any(vals <= 0): return np.full_like(ZG, np.nan)
    return np.sqrt(Om*zp**3 + Or*zp**4 + OL*sh)
def c2(p): return chi2("BIN", p)
best = None
for st in ([0.31, 0.68, 0.0224, 1, 1, 1, 1], [0.31, 0.68, 0.0224, 1.05, 1.05, 1.0, 1.0], [0.30, 0.69, 0.0224, 0.95, 0.95, 1.0, 1.0]):
    r = minimize(c2, st, method="Nelder-Mead", options=dict(maxiter=40000, maxfev=40000, xatol=1e-7, fatol=1e-6, adaptive=True))
    r = minimize(c2, r.x, method="Nelder-Mead", options=dict(maxiter=40000, maxfev=40000, xatol=1e-8, fatol=1e-7, adaptive=True))
    if best is None or r.fun < best.fun: best = r
x0 = best.x; n = len(x0); step = np.array([2e-3, 2e-3, 5e-5, 1e-2, 1e-2, 2e-2, 4e-2])
Hs = np.zeros((n, n)); f0 = c2(x0)
for i in range(n):
    for j in range(i, n):
        ei = np.zeros(n); ej = np.zeros(n); ei[i] = step[i]; ej[j] = step[j]
        Hs[i, j] = Hs[j, i] = (c2(x0+ei+ej) - c2(x0+ei-ej) - c2(x0-ei+ej) + c2(x0-ei-ej))/(4*step[i]*step[j])
cov = 2*np.linalg.inv(Hs)
r = x0[3:]; C = cov[3:, 3:]; err = np.sqrt(np.diag(C))
# Claim 1 pattern at nodes (fading law, same Om)
Om = x0[0]; ODE = 1 - Om
Ef = lambda z: np.exp(brentq(lambda l: np.exp(2*l) - Om*(1+z)**3 - ODE*(np.exp(l)/(1+z))**-0.5, -5, 10))
claim = np.array([(Ef(z)/(1+z))**-0.5 for z in NODES[1:]])
Ci = np.linalg.inv(C)
chi_const = (r-1) @ Ci @ (r-1); chi_claim = (r-claim) @ Ci @ (r-claim)
print("RESULT", SNSET, json.dumps(dict(chi2=float(best.fun), params=[float(v) for v in x0], values=[float(v) for v in r],
      errors=[float(v) for v in err], claim=[float(v) for v in claim], chi_const=float(chi_const), chi_claim=float(chi_claim),
      corr=np.round(C/np.outer(err, err), 2).tolist())))
