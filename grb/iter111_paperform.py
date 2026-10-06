"""Iteration 111 check: the paper's own form log Eiso = beta + gamma log Ep (scatter on Eiso) -- transcription check + same comparison."""
import os, numpy as np
here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, "iter111_grb.py")).read().split('out = ["Iteration')[0]; exec(src)
from scipy.optimize import minimize
def m2lnL_p(model, Om):
    dL = (1 + z)*DM_of(model, Om)*MPC; X = np.log10(4*np.pi*dL**2*S/(1 + z))
    def f(p):
        bt, g, ls = p; v = np.exp(2*ls) + sS**2 + g*g*sy**2; r = X - bt - g*y
        return np.sum(r*r/v + np.log(2*np.pi*v))
    r = min((minimize(f, p0, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-7, maxiter=20000)) for p0 in ([49.5, 1.2, np.log(0.4)], [50, 1.0, np.log(0.5)])), key=lambda r: r.fun)
    return r.fun, r.x
out = ["Iteration 111 -- paper's form: log Eiso = beta + gamma log Ep"]; res = {}
for model in ("LCDM", "LAW"):
    oms = np.linspace(0.06, 0.99, 94); L = np.array([m2lnL_p(model, o)[0] for o in oms]); i = np.argmin(L); ok = oms[L <= L[i] + 1]
    f, p = m2lnL_p(model, oms[i]); res[model] = f; L03 = m2lnL_p(model, 0.31)[0]
    out.append(f"  {model:9s} best Om {oms[i]:.2f} (1-sigma {ok.min():.2f}-{ok.max():.2f}); -2lnL {f:.3f}; beta {p[0]:.2f} gamma {p[1]:.3f} s_int {np.exp(p[2]):.3f} dex; Om=0.31 delta {L03-f:+.2f}")
for model, Om in (("EdS", 1.0), ("BARE-OPEN", 0.05)):
    f, p = m2lnL_p(model, Om); res[model] = f; out.append(f"  {model:9s} -2lnL {f:.3f}; gamma {p[1]:.3f} s_int {np.exp(p[2]):.3f}")
for m in ("LAW", "EdS", "BARE-OPEN"): out.append(f"  {m:9s} minus LCDM: delta(-2lnL) = {res[m]-res['LCDM']:+.3f}")
txt = "\n".join(out); print(txt); open(os.path.join(here, "iter111_paperform.txt"), "w").write(txt + "\n")
