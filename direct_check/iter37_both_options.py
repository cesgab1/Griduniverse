"""
ITERATION 37: test BOTH size options against the direct leftover (DESI DR2 expansion rate - matter), allowing dark energy to
rise OR fall, and checking whether a possible measurement 'event' (the galaxy-type switch, iteration-33 wobble check) changes it.
Models for the leftover omega_DE(z):
  A  tracking   : dark energy a fixed fraction of the total -> follows matter, ~ (1+z)^3   (amplitude free)
  B  stamp      : constant                                                                 (amplitude free)
  P  free shape : k (1+z)^p, p free (p < 0: rose toward today; p > 0: fell toward today)   (2 free)
  C1 Claim 1    : fixed prediction, nothing free
Each fitted to (i) all six DESI bins, (ii) without the two luminous-red-galaxy bins (z 0.51, 0.71).
Expectations written BEFORE running:
  X1 A (tracking) is strongly disfavoured at low z (it needs dark energy ~3-4x larger at z ~ 1 than today; the leftover is flat
     or falling into the past).
  X2 B (constant) and C1 both fit acceptably in (i) and (ii).
  X3 With all bins, the free shape prefers p < 0 (higher recently) at ~2 sigma, driven by the red-galaxy bins; without them p is
     consistent with 0. If so, the 'rise' is not robust and no revision of the model is warranted.
"""
import numpy as np
from scipy.optimize import minimize
src = open("leftover.py").read().split("# Claim 1 prediction")[0]
g = {}; exec(src, g); X, z = g["X"], g["z"]
from scipy.optimize import brentq
wr = g["wr"]; Om, h = 0.31, 0.6749; Or = wr/h**2; OL = 1 - Om - Or
def Ec(zz):
    a = 1/(1 + zz); f = lambda e: e*e - Om*a**-3 - Or*a**-4 - OL*(a*e)**-0.5
    return brentq(f, 1e-6, 1e4)
C1 = np.array([OL*h*h*(Ec(zz)/(1 + zz))**-0.5 for zz in z])
def analyse(mask, label):
    m = X[:, mask].mean(0); Cv = np.cov(X[:, mask].T); W = np.linalg.inv(Cv); zz = z[mask]
    chi = lambda r: r @ W @ r
    res = {}
    for name, shape in (("A tracking", (1 + zz)**3), ("B constant", np.ones(mask.sum()))):
        k = (shape @ W @ m)/(shape @ W @ shape); res[name] = (chi(m - k*shape), 1, k)
    fP = lambda p: chi(m - p[0]*(1 + zz)**p[1])
    r = minimize(fP, [0.33, 0.0], method="Nelder-Mead"); res["P free shape"] = (r.fun, 2, r.x)
    # 1-sigma on p by profile
    ps = np.linspace(-6, 6, 1201); prof = [minimize(lambda k: fP([k[0], p]), [0.33], method="Nelder-Mead").fun for p in ps]
    prof = np.array(prof); ok = ps[prof <= r.fun + 1]
    res["C1 Claim 1"] = (chi(m - C1[mask]), 0, None)
    lines = [f"--- {label} ({mask.sum()} bins) ---"]
    for k_, (c2, npar, par) in res.items():
        dof = mask.sum() - npar
        extra = f"  amplitude {par:.3f}" if npar == 1 else (f"  p = {par[1]:+.2f} (68%: {ok.min():+.2f} to {ok.max():+.2f})" if npar == 2 else "")
        lines.append(f"  {k_:14s} chi2 = {c2:6.2f} for {dof} dof{extra}")
    return lines
allm = np.ones(len(z), bool); nolrg = ~np.isin(np.round(z, 2), [0.51, 0.71])
out = ["ITERATION 37: both size options + free shape vs the direct leftover (expectations committed before running)", ""]
out += analyse(allm, "all DESI bins") + [""] + analyse(nolrg, "without the red-galaxy bins")
txt = "\n".join(out); print(txt); open("iter37_both_options.txt", "w").write(txt + "\n")
