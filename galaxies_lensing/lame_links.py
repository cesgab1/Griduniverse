"""
'Lame links': far from a galaxy the grid is so slack it can no longer carry tension, so the pull fades out.
Grid mechanics (tension_grid.py): with tension T and local anchoring k, the pull reaches only lambda = sqrt(T/k): force x (1 + r/lambda) e^(-r/lambda).
Test: KiDS-1000 lensing around isolated galaxies in 4 stellar-mass bins, 35 kpc - 2.6 Mpc (Brouwer et al. 2021, Fig. 3 data, full covariance).
Models (the galaxy's total ordinary mass M is left free per bin, so star-mass uncertainty and hot gas do not bias the SHAPE test):
   A  slack rule, pull reaches forever (a0 = 1.15e-10)
   B  slack rule with lame-link fade-out at lambda (one lambda for all bins)
   C  Newton, visible matter only (for reference)
Predicted excess surface density computed by full projection of the model's pull (no isothermal shortcut).
"""
import numpy as np
from scipy.optimize import minimize
Gk = 4.30091e-6                     # kpc (km/s)^2 / Msun
a0 = 1.15e-10*3.0857e19/1e6         # (km/s)^2/kpc
D = "/home/claude/kids/"
bins = []
for k in range(1, 5):
    d = np.loadtxt(D + f"Fig-3_Lensing-rotation-curves_Massbin-{k}.txt")
    bins.append((d[:, 0]*1e3, d[:, 1]/d[:, 4]))        # R in kpc, ESD Msun/pc^2
c = np.loadtxt(D + "Fig-3_Lensing-rotation-curves_Massbins_covmatrix.txt")
mv = np.unique(c[:, 0]); Rv = [b[0]/1e3 for b in bins]
def ix(m, R):
    k = int(np.where(np.isclose(mv, m))[0][0]); j = int(np.argmin(abs(Rv[k] - R))); return k*len(Rv[0]) + j
n = sum(len(b[0]) for b in bins); C = np.zeros((n, n))
for row in c: C[ix(row[0], row[2]), ix(row[1], row[3])] = row[4]/row[6]
C = 0.5*(C + C.T); print("covariance min eigenvalue", np.linalg.eigvalsh(C).min()); Ci = np.linalg.inv(C)
yobs = np.concatenate([b[1] for b in bins])
r = np.geomspace(1, 1e5, 3000)      # kpc
def esd_profile(M, R, lam=None, newton=False):
    gN = Gk*M/r**2
    g = gN if newton else gN/(-np.expm1(-np.sqrt(gN/a0)))
    if lam: g = g*(1 + r/lam)*np.exp(-r/lam)
    Menc = g*r**2/Gk                                   # effective enclosed mass (kpc, Msun)
    dM = np.gradient(Menc, r)                          # = 4 pi r^2 rho
    out = []
    for RR in R:
        m = r > RR
        rr = r[m]; w = dM[m]
        M2D = np.interp(RR, r, Menc) + np.trapezoid(w*(1 - np.sqrt(1 - (RR/rr)**2)), rr)
        # Sigma(R) = int dM/(2 pi r) * 1/sqrt(r^2-R^2) ... = int 4pi r^2 rho * 2/(4 pi r sqrt(r^2-R^2)) dr
        u = np.geomspace(1e-4, 1, 400)                 # substitute r = R/sqrt(1-t^2) handled via rr grid with softening
        Sig = np.trapezoid(w/(2*np.pi*rr*np.sqrt(rr**2 - RR**2 + (0.002*RR)**2)), rr)
        out.append((M2D/(np.pi*RR**2) - Sig)/1e6)
    return np.array(out)
def predict(p, lam=None, newton=False):
    out = []
    for (R, _), lm in zip(bins, p):
        out.append(esd_profile(10**lm, R, lam, newton))
    return np.concatenate(out)
def chi2(y): r_ = yobs - y; return r_ @ Ci @ r_
res = {}
for lab, lam_free, newton in (("A  slack rule, reaches forever", False, False), ("C  Newton, visible matter only", False, True), ("B  slack rule + lame-link fade-out", True, False)):
    p0 = [10.4, 10.7, 10.9, 11.1] + ([np.log10(1000)] if lam_free else [])
    f = lambda p: chi2(predict(p[:4], 10**p[4] if lam_free else None, newton))
    b = minimize(f, p0, method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-3, fatol=1e-2))
    res[lab] = b
    extra = f", fade-out distance lambda = {10**b.x[4]/1e3:.2f} Mpc" if lam_free else ""
    print(f"{lab:40s} chi2 = {b.fun:6.1f} for {len(yobs)} points; fitted ordinary mass per bin log M = {', '.join(f'{x:.2f}' for x in b.x[:4])}{extra}")
# profile of lambda
print("\nchi2 vs fade-out distance (masses refitted):")
for lam in (300, 500, 1000, 2000, 4000, 1e5):
    f = lambda p: chi2(predict(p, lam))
    b = minimize(f, res["B  slack rule + lame-link fade-out"].x[:4], method="Nelder-Mead", options=dict(maxiter=2000, xatol=1e-3, fatol=1e-2))
    print(f"   lambda = {lam/1e3:6.1f} Mpc: chi2 = {b.fun:.1f}")
print("\nstellar-mass bin limits (log M*): 8.5-10.3, 10.3-10.6, 10.6-10.8, 10.8-11.0")
