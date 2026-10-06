"""Iteration 111: GRB (A118, Amati relation) Hubble diagram -- LAW vs LCDM, plus EdS and BARE-OPEN (PREREG_111.md)."""
import os, numpy as np, pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize
here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, "../cosmology_fits/fit_law.py")).read(); exec(src.split("def chi2")[0])
d = pd.read_csv(os.path.join(here, "A118.csv"), comment="#", dtype={"name": str})
z, Ep, eEp, S, eS = (d[c].values for c in ["z", "Ep", "Ep_err", "Sbolo", "Sbolo_err"]); S, eS = S*1e-5, eS*1e-5
y, sy = np.log10(Ep), eEp/(Ep*np.log(10)); sS = eS/(S*np.log(10))
MPC = 3.0857e24; h = 0.70
def DM_of(model, Om):
    zp = 1 + ZG
    if model in ("LCDM", "LAW"): E = E_of_z(model, Om, h, []); Ok = 0.0
    elif model == "EdS": E = np.sqrt(zp**3); Ok = 0.0
    else: Ok = 1 - Om; E = np.sqrt(Om*zp**3 + Ok*zp**2)
    chi = cumulative_trapezoid(1/E, ZG, initial=0)*C_KMS/(100*h)
    DMt = chi if Ok < 1e-8 else C_KMS/(100*h)/np.sqrt(Ok)*np.sinh(np.sqrt(Ok)*chi*100*h/C_KMS)
    return np.interp(z, ZG, DMt)
def m2lnL(model, Om):
    dL = (1 + z)*DM_of(model, Om)*MPC
    x = np.log10(4*np.pi*dL**2*S/(1 + z)/1e52)
    def f(p):
        a, b, ls = p; v = np.exp(2*ls) + sy**2 + b*b*sS**2; r = y - a - b*x
        return np.sum(r*r/v + np.log(2*np.pi*v))
    best = min((minimize(f, p0, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-7, maxiter=20000)) for p0 in ([2.0, 1.2, np.log(0.4)], [2.3, 0.8, np.log(0.3)])), key=lambda r: r.fun)
    return best.fun, best.x
out = ["Iteration 111 -- A118 gamma-ray bursts, Amati relation fitted jointly with each cosmology (H0 absorbed)", f"N = {len(z)}, z {z.min():.2f}-{z.max():.2f}"]
res = {}
for model in ("LCDM", "LAW"):
    oms = np.linspace(0.06, 0.99, 94); L = []
    for o in oms:
        try: L.append(m2lnL(model, o)[0])
        except Exception: L.append(np.nan)
    L = np.array(L); i = np.nanargmin(L); ok = oms[L <= L[i] + 1]
    f, p = m2lnL(model, oms[i]); res[model] = f
    L03 = m2lnL(model, 0.31)[0]
    out.append(f"  {model:9s} best Om {oms[i]:.2f} (1-sigma range {ok.min():.2f}-{ok.max():.2f}); -2lnL {f:.3f}; Amati a {p[0]:.3f} b {p[1]:.3f} s_int {np.exp(p[2]):.3f} dex;"
               f" at Om = 0.31: -2lnL {L03:.3f} (delta {L03-f:+.2f})")
for model, Om in (("EdS", 1.0), ("BARE-OPEN", 0.05)):
    f, p = m2lnL(model, Om); res[model] = f
    out.append(f"  {model:9s} (Om fixed {Om}): -2lnL {f:.3f}; Amati b {p[1]:.3f} s_int {np.exp(p[2]):.3f}")
out.append("")
for m in ("LAW", "EdS", "BARE-OPEN"):
    dl = res[m] - res["LCDM"]
    out.append(f"  {m:9s} minus LCDM: delta(-2lnL) = {dl:+.3f}" + ("  (LCDM has 1 more free parameter here)" if m != "LAW" else ""))
txt = "\n".join(out); print(txt); open(os.path.join(here, "iter111_grb.txt"), "w").write(txt + "\n")
