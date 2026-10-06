"""Stuck item S2: supernova standardisation -- does our law's edge over constant dark energy depend on how the
host-galaxy 'mass step' is handled? Pantheon+ (m_b_corr already includes the published standardisation).
Variants: published; + free extra step (split log M* = 10); + free step that grows with redshift (step x (1 + k z)),
the worry being that changing host galaxies over cosmic time could imitate evolving dark energy."""
import os, numpy as np
from scipy.optimize import minimize
os.environ["SNSET"] = "PANTHEON"
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../cosmology_fits/fit_law.py")).read()
exec(src.split("Ndata =")[0])
hm = sn["HOST_LOGMASS"].values; s = np.where(hm >= 10, 0.5, -0.5); ok = (hm > 0) & (hm < 13)
s = np.where(ok, s, 0.0)
mb0 = mb.copy()
def c2(model, p, extra):
    global mb
    d, k = (list(extra) + [0, 0])[:2]
    mb = mb0 + (d + k*z_sn)*s; v = chi2(model, p); mb = mb0; return v
def fit(model, nx, step_free):
    best = None
    for x0 in ([0.31, 0.68, 0.0224], [0.30, 0.69, 0.0223]):
        f = lambda q: c2(model, q[:3], q[3:]) if nx else c2(model, q, [])
        st = x0 + [0.0]*nx
        r = minimize(f, st, method="Nelder-Mead", options=dict(maxiter=20000, xatol=1e-7, fatol=1e-6))
        r = minimize(f, r.x, method="Nelder-Mead", options=dict(maxiter=20000, xatol=1e-7, fatol=1e-6))
        if best is None or r.fun < best.fun: best = r
    return best
out = [f"S2 -- host-galaxy step, Pantheon+ ({ok.sum()} of {len(s)} supernovae have host masses)"]
for lab, nx in [("published standardisation", 0), ("+ free extra step", 1), ("+ free step changing linearly with redshift (d0 + d1 z)", 2)]:
    L = fit("LCDM", nx, True); W = fit("LAW", nx, True)
    ex = "" if nx == 0 else f"; fitted step {L.x[3]:+.3f} mag" + (f", change per unit z {L.x[4]:+.3f} mag" if nx == 2 else "") + " (LCDM fit)"
    out.append(f"  {lab}: law minus LCDM chi2 = {W.fun - L.fun:+.2f}{ex}")
txt = "\n".join(out); print(txt); open("stuck2_host_step.txt", "w").write(txt + "\n")
