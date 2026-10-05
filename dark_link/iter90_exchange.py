import os, json, numpy as np
from scipy.optimize import minimize
from scipy.integrate import solve_ivp, cumulative_trapezoid
os.chdir("/home/claude/griduniverse/cosmology_fits")
exec(open("fit_law.py").read().split("Ndata = len(bao)")[0])
BETA0 = float(os.environ.get("BASE_BETA", 0.0))
def background(Om, h, wb, xi):
    Or = W_R/h**2; Ob = wb/h**2; Oc0 = Om - Ob; OL = 1 - Om - Or
    def rhs(x, y):
        rc, rde = y; a = np.exp(x); rm = rc + Ob/a**3; rr = Or/a**4
        q = (rm/2 + rr - rde)/(rm + rr + rde)
        return [-3*rc + xi*rde, BETA0*q*rde - xi*rde]
    xs = -np.log(1 + ZG[ZG <= 1e4])
    sol = solve_ivp(rhs, [0, xs[-1]], [Oc0, OL], t_eval=xs, rtol=1e-9, atol=1e-14, method="LSODA")
    if not sol.success: return None, None
    rc = np.zeros_like(ZG); rde = np.zeros_like(ZG)
    n = len(xs); rc[:n], rde[:n] = sol.y
    C = rc[n-1]*(1 + ZG[n-1])**-3                     # early a^-3 coefficient
    rc[n:] = C*(1 + ZG[n:])**3
    zp = 1 + ZG
    E2 = rc + Ob*zp**3 + Or*zp**4 + rde
    return np.sqrt(E2), (C + Ob)*h**2
def chi2x(p, xi):
    Om, h, wb = p
    if not (0.1 < Om < 0.6 and 0.5 < h < 0.9 and 0.018 < wb < 0.027): return 1e12
    E, wm = background(Om, h, wb, xi)
    if E is None or not np.all(np.isfinite(E)) or np.any(E <= 0) or wm <= W_NU: return 1e12
    H0 = 100*h; chi = cumulative_trapezoid(1/E, ZG, initial=0)*C_KMS/H0
    DM = lambda z: np.interp(z, ZG, chi); DH = lambda z: C_KMS/(H0*np.interp(z, ZG, E))
    wcb = wm - W_NU
    rd = 55.154*np.exp(-72.3*(W_NU + 0.0006)**2)/(wcb**0.25351*wb**0.12807)
    g1 = 0.0783*wb**-0.238/(1 + 39.5*wb**0.763); g2 = 0.560/(1 + 21.1*wb**1.81)
    zs = 1048*(1 + 0.00124*wb**-0.738)*(1 + g1*wm**g2)
    Rb = 3*wb/(4*W_GAMMA)/(1 + ZG); cs = C_KMS/np.sqrt(3*(1 + Rb))/(H0*E); sel = ZG >= zs
    rs = np.trapezoid(np.r_[np.interp(zs, ZG, cs), cs[sel]], np.r_[zs, ZG[sel]]) + cs[sel][-1]*(1 + ZG[-1])
    R = np.sqrt(wm)/h*H0*DM(zs)/C_KMS; lA = np.pi*DM(zs)/rs*LA_CAL
    pred = []
    for z, qq in zip(bao.z, bao.q):
        if qq == "DM_over_rs": pred.append(DM(z)/rd)
        elif qq == "DH_over_rs": pred.append(DH(z)/rd)
        else: pred.append((z*DM(z)**2*DH(z))**(1/3)/rd)
    d = bao.v.values - np.array(pred); cb = d @ bao_icov @ d
    mu = 5*np.log10((1 + zhel)*DM(z_sn)) + 25; r = mb - mu
    A = r @ sn_icov @ r; B = (sn_icov @ r).sum(); csn = A - B**2/sn_icov_sum
    dc = np.array([R, lA, wb]) - cmb_mean
    return cb + csn + dc @ cmb_icov @ dc
prof = []
x0 = [0.31, 0.68, 0.0224]
for xi in [-0.4, -0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3, 0.4]:
    best = None
    for st in (x0, [0.30, 0.69, 0.0224], [0.32, 0.67, 0.0224]):
        r = minimize(lambda p: chi2x(p, xi), st, method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-7, fatol=1e-5))
        if best is None or r.fun < best.fun: best = r
    prof.append((xi, float(best.fun), [float(v) for v in best.x])); x0 = list(best.x)
    print("xi", xi, best.fun, best.x, flush=True)
print("RESULT", SNSET, BETA0, json.dumps(prof))
