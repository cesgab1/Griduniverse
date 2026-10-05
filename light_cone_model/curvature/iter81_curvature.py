"""Iteration 81: curvature free, constant vs Claim 1. See PREREG_81.md."""
import os, numpy as np, json
from scipy.optimize import minimize
from scipy.integrate import solve_ivp, cumulative_trapezoid
os.chdir("/home/claude/griduniverse/cosmology_fits")
exec(open("fit_law.py").read().split("Ndata = len(bao)")[0])
def E_k(model, Om, h, Ok):
    Or = W_R/h**2; OL = 1 - Om - Or - Ok; zp = 1 + ZG
    if model == "LCDM": return np.sqrt(Om*zp**3 + Or*zp**4 + Ok*zp**2 + OL)
    def rhs(lna, y):
        a = np.exp(lna); rde = np.exp(y[0]); rm, rr, rk = Om*a**-3, Or*a**-4, Ok*a**-2
        return [0.5*(rm/2 + rr - rde)/(rm + rr + rk + rde)]
    zi = ZG[ZG <= 30]
    sol = solve_ivp(rhs, [0, -np.log(31)], [np.log(OL)], t_eval=-np.log(1 + zi), rtol=1e-8, atol=1e-10)
    rde = np.zeros_like(ZG); rde[:len(zi)] = np.exp(sol.y[0])
    return np.sqrt(Om*zp**3 + Or*zp**4 + Ok*zp**2 + rde)
def chi2k(model, p):
    Om, h, wb, Ok = p
    if not (0.1 < Om < 0.6 and 0.5 < h < 0.9 and 0.018 < wb < 0.027 and -0.05 < Ok < 0.05): return 1e12
    E = E_k(model, Om, h, Ok)
    if not np.all(np.isfinite(E)) or np.any(E <= 0): return 1e12
    H0 = 100*h; chi = cumulative_trapezoid(1/E, ZG, initial=0)*C_KMS/H0
    if abs(Ok) > 1e-8:
        s = np.sqrt(abs(Ok))*H0/C_KMS
        dm = np.sinh(s*chi)/s if Ok > 0 else np.sin(s*chi)/s
    else: dm = chi
    DM = lambda z: np.interp(z, ZG, dm); DH = lambda z: C_KMS/(H0*np.interp(z, ZG, E))
    wm = Om*h**2; wcb = wm - W_NU
    rd = 55.154*np.exp(-72.3*(W_NU + 0.0006)**2)/(wcb**0.25351*wb**0.12807)
    g1 = 0.0783*wb**-0.238/(1 + 39.5*wb**0.763); g2 = 0.560/(1 + 21.1*wb**1.81)
    zs = 1048*(1 + 0.00124*wb**-0.738)*(1 + g1*wm**g2)
    Rb = 3*wb/(4*W_GAMMA)/(1 + ZG); cs = C_KMS/np.sqrt(3*(1 + Rb))/(H0*E); sel = ZG >= zs
    rs = np.trapezoid(np.r_[np.interp(zs, ZG, cs), cs[sel]], np.r_[zs, ZG[sel]]) + cs[sel][-1]*(1 + ZG[-1])
    R = np.sqrt(Om)*H0*DM(zs)/C_KMS; lA = np.pi*DM(zs)/rs*LA_CAL
    pred = []
    for z, q in zip(bao.z, bao.q):
        if q == "DM_over_rs": pred.append(DM(z)/rd)
        elif q == "DH_over_rs": pred.append(DH(z)/rd)
        else: pred.append((z*DM(z)**2*DH(z))**(1/3)/rd)
    d = bao.v.values - np.array(pred); cb = d @ bao_icov @ d
    mu = 5*np.log10((1 + zhel)*DM(z_sn)) + 25; r = mb - mu
    A = r @ sn_icov @ r; B = (sn_icov @ r).sum(); csn = A - B**2/sn_icov_sum
    dc = np.array([R, lA, wb]) - cmb_mean
    return cb + csn + dc @ cmb_icov @ dc
def fit(model, fixk=None):
    best = None
    for st in ([0.31, 0.68, 0.0224, 0.002], [0.30, 0.69, 0.0223, -0.002], [0.32, 0.67, 0.0225, 0.0]):
        if fixk is not None:
            f = lambda p: chi2k(model, [*p, fixk]); x0 = st[:3]
        else:
            f = lambda p: chi2k(model, p); x0 = st
        r = minimize(f, x0, method="Nelder-Mead", options=dict(maxiter=8000, xatol=1e-7, fatol=1e-5))
        r = minimize(f, r.x, method="Nelder-Mead", options=dict(maxiter=8000, xatol=1e-8, fatol=1e-6))
        if best is None or r.fun < best.fun: best = r
    return best
res = {}
for m in ("LCDM", "LAW"):
    free = fit(m); flat = fit(m, fixk=0.0)
    ok = free.x[3]
    # 1-sigma width from the chi2 profile: scan Ok
    grid = np.linspace(ok - 0.006, ok + 0.006, 25)
    prof = [fit(m, fixk=g).fun for g in grid]
    a2, a1, a0 = np.polyfit(grid, prof, 2); sig = 1/np.sqrt(a2)
    res[m] = dict(chi2_free=float(free.fun), chi2_flat=float(flat.fun), Ok=float(ok), sigma=float(sig), params=[float(v) for v in free.x])
print("RESULT", SNSET, json.dumps(res))
