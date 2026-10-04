"""
ITERATION 23 (Coalesce: the size is tied to our horizon): holographic dark energy with the FUTURE event horizon
(Li 2004): rho_DE = 3 c^2 M_P^2 / R_h^2  ->  w = -1/3 - 2 sqrt(Omega_DE)/(3c),  d ln rho_DE / d ln a = -2 + 2 sqrt(Omega_DE)/c.
One extra parameter c. Same data and code as Claim 1 (DESI DR2 BAO + Planck 2018 distance priors + one SN set).
Compare: LCDM (0 extra), Claim 1 (0 extra), holographic (1 extra), w0wa (2 extra).
"""
import numpy as np, os, json
from scipy.integrate import solve_ivp
os.environ["ONLY"] = "NONE"; FIT = "/home/claude/griduniverse/cosmology_fits/fit_law.py"
g = {"__name__": "x", "__file__": FIT}
try: exec(open(FIT).read(), g)
except SystemExit: pass
ZG, W_R = g["ZG"], g["W_R"]
def E_hde(model, Om, h, extra):
    c = extra[0]; Or = W_R/h**2; OL = 1 - Om - Or; zp = 1 + ZG; zi = ZG          # FULL history (fix: was z <= 30; HDE is not negligible early)
    def rhs(lna, y):
        a = np.exp(lna); rde = np.exp(y[0]); Ode = rde/(rde + Om*a**-3 + Or*a**-4)
        return [-2 + 2*np.sqrt(Ode)/c]
    s = solve_ivp(rhs, [0, -np.log(1 + zi[-1])], [np.log(OL)], t_eval=-np.log(1 + zi), rtol=1e-9, atol=1e-11, method='LSODA')
    rde = np.zeros_like(ZG); rde[:len(zi)] = np.exp(s.y[0])
    return np.sqrt(Om*zp**3 + Or*zp**4 + rde)
SN = os.environ.get("SNSET", "PANTHEON")
L = g["best"]("LCDM", [[0.31, 0.68, 0.0224]])
E0 = g["E_of_z"]; g["E_of_z"] = lambda m, Om, h, ex: E_hde(m, Om, h, ex) if m == "HDE" else E0(m, Om, h, ex)
orig_chi2 = g["chi2"]
def chi2_hde(model, p, parts=False):
    if model == "HDE" and not (0.2 < p[3] < 3.0): return 1e12
    return orig_chi2(model, p, parts)
g["chi2"] = chi2_hde
H = g["best"]("HDE", [[0.31, 0.68, 0.0224, 0.7], [0.30, 0.68, 0.0224, 0.9], [0.32, 0.67, 0.0224, 0.55]])
Om, h, wb, c = H.x; Ode0 = 1 - Om
w0 = -1/3 - 2*np.sqrt(Ode0)/(3*c)
out = {"SN": SN, "LCDM": float(L.fun), "HDE": float(H.fun), "params": [float(x) for x in H.x], "w0": float(w0)}
print("RESULT", json.dumps(out))
tot, cb, cs, cc = g["chi2"]("HDE", H.x, parts=True); tl, lb, ls, lc = orig_chi2("LCDM", L.x, parts=True) if False else (None,)*4
g["E_of_z"] = E0; tl, lb, ls, lc = orig_chi2("LCDM", L.x, parts=True)
print(f"PARTS {SN}: HDE BAO {cb:.1f} SN {cs:.1f} CMB {cc:.1f} | LCDM BAO {lb:.1f} SN {ls:.1f} CMB {lc:.1f}")
