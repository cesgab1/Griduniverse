"""Iteration 99: timescape (uneven expansion) vs LCDM vs LAW on late-universe data (PREREG_99.md)."""
import os, numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize_scalar, minimize
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../cosmology_fits/fit_law.py")).read()
exec(src.split("def chi2")[0])
ZMAX = 3.0
zgrid = np.linspace(0, ZMAX, 3001)

def timescape(f):
    """returns H0*D_M/c (dressed H0 units) on zgrid"""
    B = 2*(1 - f)*(2 + f)/(9*f); b = B**(1/3); y0 = ((2 + f)/3)**(1/3)
    y = np.linspace(y0, 0.02, 20000)
    opz = 2**(4/3)*y*(y**3 + B)/(f**(1/3)*y**3*(2*y**3 + 3*B)**(4/3))
    F = lambda y: 2*y + b/6*np.log((y + b)**2/(y**2 - b*y + b**2)) + b/np.sqrt(3)*np.arctan((2*y - b)/(np.sqrt(3)*b))
    dL_bare = opz**2*y**2*(F(y0) - F(y))                  # Hbar0 d_L / c
    ratio = (4*f*f + f + 4)/(2*(2 + f))                    # H0 / Hbar0
    DM = dL_bare/opz*ratio                                 # H0 D_M / c
    z = opz - 1
    return np.interp(zgrid, z, DM)

def flat(model, Om):
    E = E_of_z(model, Om, 0.68, [])
    chi = cumulative_trapezoid(1/E, ZG, initial=0)
    return np.interp(zgrid, ZG, chi)                       # H0 D_M / c

def sn_chi2(DM):
    dl = (1 + zhel)*np.interp(z_sn, zgrid, DM)
    r = mb - 5*np.log10(dl); B = (sn_icov @ r).sum()
    return r @ sn_icov @ r - B**2/sn_icov_sum

def bao_chi2(DM):
    DH = np.gradient(DM, zgrid)
    def c2(s):                                             # s = r_d H0 / c
        pred = []
        for z, q in zip(bao.z, bao.q):
            dm, dh = np.interp(z, zgrid, DM), np.interp(z, zgrid, DH)
            pred.append(dm/s if q == "DM_over_rs" else dh/s if q == "DH_over_rs" else (z*dm*dm*dh)**(1/3)/s)
        d = bao.v.values - np.array(pred); return d @ bao_icov @ d
    r = minimize_scalar(c2, bounds=(0.01, 0.1), method="bounded"); return r.fun

# self-check: low-z limit
for f in (0.6, 0.8):
    DM = timescape(f); i = np.searchsorted(zgrid, 0.01)
    print(f"self-check f_v0={f}: H0 D_M/c at z={zgrid[i]:.3f} = {DM[i]:.5f} (expect ~ z - small) ratio {DM[i]/zgrid[i]:.4f}")

out = [f"Iteration 99 -- SNSET={SNSET}"]
def best(fun, lo, hi):
    r = minimize_scalar(fun, bounds=(lo, hi), method="bounded", options=dict(xatol=1e-5)); return r.x, r.fun
models = {"timescape": (lambda p: timescape(p), 0.30, 0.99), "LCDM": (lambda p: flat("LCDM", p), 0.1, 0.6),
          "LAW": (lambda p: flat("LAW", p), 0.1, 0.6)}
res = {}
for name, (mk, lo, hi) in models.items():
    ps, cs = best(lambda p: sn_chi2(mk(p)), lo, hi)
    pb, cb = best(lambda p: bao_chi2(mk(p)), lo, hi)
    pc, cc = best(lambda p: sn_chi2(mk(p)) + bao_chi2(mk(p)), lo, hi)
    res[name] = (cs, cb, cc)
    out.append(f"  {name:9s} SN only: chi2 {cs:9.2f} (param {ps:.3f}) | BAO shape: chi2 {cb:6.2f} (param {pb:.3f}) | SN+BAO: chi2 {cc:9.2f} (param {pc:.3f})")
L = res["LCDM"]
for name in ("timescape", "LAW"):
    out.append(f"  {name:9s} minus LCDM: SN {res[name][0]-L[0]:+7.2f}   BAO {res[name][1]-L[1]:+7.2f}   SN+BAO {res[name][2]-L[2]:+7.2f}")
txt = "\n".join(out); print(txt)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"iter99_{SNSET}.txt"), "w").write(txt + "\n")
