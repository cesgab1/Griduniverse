"""
Direct fit to public data: DESI DR2 BAO + Pantheon+ supernovae + Planck 2018 CMB distance priors.
Models
  LCDM : constant dark energy                                   params  Om, h, wb
  w0wa : DESI's two-parameter dark energy (CPL)                 params  Om, h, wb, w0, wa
  LAW  : d ln rho_DE / d ln a = BETA q   (w = -1 - BETA q/3), rho_DE ~ adot^-BETA, BETA = 1/2 (env BETA); no free DE parameter   params  Om, h, wb
Data
  DESI DR2 BAO: CobayaSampler/bao_data desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb (13 points + covariance)
  Pantheon+: PantheonPlusSH0ES/DataRelease, zHD > 0.01, STAT+SYS covariance, absolute magnitude marginalised
  CMB: Chen, Huang & Wang 2019 Planck 2018 distance priors (wCDM case): R, l_A, omega_b (n_s marginalised)
"""
import numpy as np, pandas as pd, sys, json
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import minimize
C_KMS = 299792.458
W_GAMMA = 2.469e-5
W_R = W_GAMMA * (1 + 0.2271 * 3.046)
import os
W_NU = float(os.environ.get("MNU", 0.06)) / 93.14
BETA = float(os.environ.get("BETA", 0.5))   # rho_DE ~ adot^-BETA  <=>  d ln rho_DE/d ln a = BETA q  (the law: BETA = 1/2)
STRETCH = float(os.environ.get("STRETCH", 1.0))   # link stretch exponent s: physical link length ∝ a^s; law d ln rho/d ln a = BETA (q + 1 - s); s = 1 is Claim 1
LA_CAL = float(os.environ.get("LA_CAL", 1.0))   # set below: CAMB l_A / this code's l_A at the Planck 2018 best-fit LCDM (recalibrated Oct 2026 after r_s fix)

# ---------------- data ----------------
BAO = "/home/claude/cobayasampler/bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_"
bao = pd.read_csv(BAO + "mean.txt", sep=r"\s+", comment="#", header=None, names=["z", "v", "q"])
bao_icov = np.linalg.inv(np.loadtxt(BAO + "cov.txt"))
import os
SNSET = os.environ.get("SNSET", "PANTHEON")
if SNSET == "PANTHEON":
    PP = "/home/claude/pantheonplussh0es/datarelease/Pantheon+_Data/4_DISTANCES_AND_COVAR/"
    sn = pd.read_csv(PP + "Pantheon+SH0ES.dat", sep=r"\s+")
    raw = np.loadtxt(PP + "Pantheon+SH0ES_STAT+SYS.cov")
    N = int(raw[0]); cov = raw[1:].reshape(N, N)
    m = (sn["zHD"] > 0.01).values
    sn = sn[m].reset_index(drop=True); cov = cov[np.ix_(m, m)]
    sn_icov = np.linalg.inv(cov)
    z_sn, zhel, mb = sn["zHD"].values, sn["zHEL"].values, sn["m_b_corr"].values
elif SNSET == "UNION3":   # Union3 + UNITY1.5 compressed distances (Rubin et al. 2023): row 0 = z, col 0 = mu, rest = inverse covariance
    from astropy.io import fits
    d = fits.open("/home/claude/union3_release/mu_mat_union3_cosmo=2_mu.fits")[0].data
    z_sn = d[0, 1:].astype(float); zhel = z_sn.copy(); mb = d[1:, 0].astype(float); sn_icov = d[1:, 1:].astype(float)
else:   # DES-Dovekie (recalibrated DES 5-year sample); file stores the INVERSE covariance, upper triangle
    DD = "/home/claude/des-science/des-sn5yr/4_DISTANCES_COVMAT/"
    rows = [l.split() for l in open(DD + "DES-Dovekie_HD.csv") if l.startswith("SN:")]
    z_sn = np.array([float(r[3]) for r in rows]); zhel = np.array([float(r[4]) for r in rows]); mb = np.array([float(r[5]) for r in rows])
    d = np.load(DD + "STAT+SYS.npz"); n = int(d["nsn"][0])
    sn_icov = np.zeros((n, n)); sn_icov[np.triu_indices(n)] = d["cov"]
    il = np.tril_indices(n, -1); sn_icov[il] = sn_icov.T[il]
sn_icov_sum = sn_icov.sum()
cmb_mean = np.array([1.7493, 301.462, 0.02239])
cmb_sig = np.array([0.00465, 0.0895, 0.00015])
cmb_corr = np.array([[1.0, 0.47, -0.66], [0.47, 1.0, -0.34], [-0.66, -0.34, 1.0]])
cmb_icov = np.linalg.inv(cmb_corr * np.outer(cmb_sig, cmb_sig))
print(f"data: {len(bao)} BAO points, {len(z_sn)} supernovae, 3 CMB priors")

# ---------------- expansion histories ----------------
ZG = np.concatenate([np.linspace(0, 3, 1500)[:-1], np.geomspace(3, 1e8, 5000)])
def E_of_z(model, Om, h, extra):
    Or = W_R / h**2; OL = 1 - Om - Or; zp = 1 + ZG
    if model == "LCDM":
        return np.sqrt(Om * zp**3 + Or * zp**4 + OL)
    if model == "w0wa":
        w0, wa = extra
        return np.sqrt(Om * zp**3 + Or * zp**4 + OL * zp**(3 * (1 + w0 + wa)) * np.exp(-3 * wa * ZG / zp))
    if model == "TAB":   # fixed, pre-computed dark-energy shape rho_DE(z)/rho_DE(0) from file TABFILE (no DE parameter)
        tz, ts = np.loadtxt(os.environ["TABFILE"], unpack=True)
        sh = np.interp(ZG, tz, ts, right=ts[-1])
        return np.sqrt(Om * zp**3 + Or * zp**4 + OL * sh)
    # LAW: integrate d ln rho / d ln a = BETA (q + 1 - STRETCH) backwards; dark energy is negligible above z ~ 20
    def rhs(lna, y):
        a = np.exp(lna); rde = np.exp(y[0]); rm, rr = Om * a**-3, Or * a**-4
        if os.environ.get('EXACTQ'): return [BETA * ((rm / 2 + rr - rde - BETA * (1 - STRETCH) * rde / 2) / (rm + rr + rde + BETA * rde / 2) + 1 - STRETCH)]   # self-consistent q (w of DE included)
        return [BETA * ((rm / 2 + rr - rde) / (rm + rr + rde) + 1 - STRETCH)]
    zi = ZG[ZG <= 30]
    sol = solve_ivp(rhs, [0, -np.log(31)], [np.log(OL)], t_eval=-np.log(1 + zi), rtol=1e-8, atol=1e-10)
    rde = np.zeros_like(ZG); rde[:len(zi)] = np.exp(sol.y[0])
    return np.sqrt(Om * zp**3 + Or * zp**4 + rde)

def observables(model, p):
    Om, h, wb = p[:3]; extra = p[3:]
    E = E_of_z(model, Om, h, extra)
    if not np.all(np.isfinite(E)) or np.any(E <= 0): return None
    GD, GA = float(os.environ.get("GDELTA", 0)), float(os.environ.get("GAT", 1e-4))
    if GD: E = E * np.sqrt(1 + GD / (1 + ((1 / (1 + ZG)) / GA)**4))    # gravity stronger by GD before a ~ GA (grid stretching phase), today's value after
    H0 = 100 * h
    chi = cumulative_trapezoid(1 / E, ZG, initial=0) * C_KMS / H0           # comoving distance, Mpc
    DM = lambda z: np.interp(z, ZG, chi)
    DH = lambda z: C_KMS / (H0 * np.interp(z, ZG, E))
    wm = Om * h**2; wcb = wm - W_NU
    rd = 55.154 * np.exp(-72.3 * (W_NU + 0.0006)**2) / (wcb**0.25351 * wb**0.12807)
    # CMB: recombination redshift, sound horizon, shift parameter, acoustic scale
    g1 = 0.0783 * wb**-0.238 / (1 + 39.5 * wb**0.763); g2 = 0.560 / (1 + 21.1 * wb**1.81)
    zs = 1048 * (1 + 0.00124 * wb**-0.738) * (1 + g1 * wm**g2)
    Rb = 3 * wb / (4 * W_GAMMA) / (1 + ZG)
    cs_over_H = C_KMS / np.sqrt(3 * (1 + Rb)) / (H0 * E)
    sel = ZG >= zs
    zz = np.r_[zs, ZG[sel]]; ff = np.r_[np.interp(zs, ZG, cs_over_H), cs_over_H[sel]]    # start exactly at z* (fixed Oct 2026: was first grid point above z*)
    rs = np.trapezoid(ff, zz) + cs_over_H[sel][-1] * (1 + ZG[-1])       # tail: radiation era, c_s/H ~ (1+z)^-2
    R = np.sqrt(Om) * H0 * DM(zs) / C_KMS; lA = np.pi * DM(zs) / rs * LA_CAL   # calibrated once to Planck 2018 LCDM (same for all models)
    return dict(DM=DM, DH=DH, rd=rd, R=R, lA=lA, wb=wb)

if "LA_CAL" not in os.environ:   # calibrate l_A once: CAMB (Planck 2018 best fit H0 67.36, ombh2 0.02237, omch2 0.1200, mnu 0.06) gives 301.7286
    _Om = (0.02237 + 0.1200 + W_NU * 0.6736**2) / 0.6736**2
    LA_CAL = 301.7286 / observables("LCDM", [_Om, 0.6736, 0.02237])["lA"]
print(f"l_A calibration factor {LA_CAL:.5f}")

def chi2(model, p, parts=False):
    Om, h, wb = p[:3]
    if not (0.1 < Om < 0.6 and 0.5 < h < 0.9 and 0.018 < wb < 0.027): return 1e12
    o = observables(model, p)
    if o is None: return 1e12
    # BAO
    pred = []
    for z, q in zip(bao.z, bao.q):
        if q == "DM_over_rs": pred.append(o["DM"](z) / o["rd"])
        elif q == "DH_over_rs": pred.append(o["DH"](z) / o["rd"])
        else: pred.append((z * o["DM"](z)**2 * o["DH"](z))**(1 / 3) / o["rd"])
    d = bao.v.values - np.array(pred); c_bao = d @ bao_icov @ d
    # supernovae (absolute magnitude marginalised analytically)
    mu = 5 * np.log10((1 + zhel) * o["DM"](z_sn)) + 25
    r = mb - mu; A = r @ sn_icov @ r; B = (sn_icov @ r).sum()
    c_sn = A - B**2 / sn_icov_sum
    # CMB
    dc = np.array([o["R"], o["lA"], o["wb"]]) - cmb_mean; c_cmb = dc @ cmb_icov @ dc
    tot = c_bao + c_sn + c_cmb
    return (tot, c_bao, c_sn, c_cmb) if parts else tot

def best(model, starts):
    bestres = None
    for x0 in starts:
        r = minimize(lambda p: chi2(model, p), x0, method="Nelder-Mead",
                     options=dict(maxiter=6000, maxfev=12000, xatol=1e-6, fatol=1e-4))
        r = minimize(lambda p: chi2(model, p), r.x, method="Nelder-Mead", options=dict(maxiter=6000, xatol=1e-7, fatol=1e-5))
        if bestres is None or r.fun < bestres.fun: bestres = r
    return bestres

Ndata = len(bao) + len(z_sn) + 3
results = {}
ONLY = os.environ.get("ONLY")
for model, starts in [("LCDM", [[0.31, 0.68, 0.0224], [0.30, 0.69, 0.0223]]),
                      ("TAB", [[0.31, 0.68, 0.0224], [0.30, 0.67, 0.0223]]),
                      ("LAW", [[0.31, 0.68, 0.0224], [0.30, 0.67, 0.0223]]),
                      ("w0wa", [[0.31, 0.68, 0.0224, -0.8, -0.6], [0.32, 0.66, 0.0224, -0.7, -1.0], [0.30, 0.69, 0.0224, -0.95, -0.1]])]:
    if ONLY and model not in ONLY.split(","): continue
    r = best(model, starts)
    tot, cb, cs, cc = chi2(model, r.x, parts=True); k = len(r.x) + 1   # +1 for the SN absolute magnitude
    results[model] = dict(params=[float(x) for x in r.x], chi2=float(tot), bao=float(cb), sn=float(cs), cmb=float(cc),
                          k=k, AIC=float(tot + 2 * k), BIC=float(tot + k * np.log(Ndata)))
    print(f"{model:5s}  chi2 = {tot:9.2f}  (BAO {cb:6.2f}, SN {cs:8.2f}, CMB {cc:5.2f})  params {np.round(r.x, 4)}")
if ONLY:
    print("ONLYRESULT", json.dumps(results)); sys.exit()
L = results["LCDM"]
for m_ in ["LAW", "w0wa"]:
    R = results[m_]
    print(f"{m_:5s} vs LCDM:  delta chi2 = {R['chi2'] - L['chi2']:7.2f}   delta AIC = {R['AIC'] - L['AIC']:7.2f}   delta BIC = {R['BIC'] - L['BIC']:7.2f}")
R, W = results["LAW"], results["w0wa"]
print(f"LAW vs w0wa: delta chi2 = {R['chi2'] - W['chi2']:.2f} with 2 fewer parameters;  delta AIC = {R['AIC'] - W['AIC']:.2f},  delta BIC = {R['BIC'] - W['BIC']:.2f}")
print(f"BETA = {BETA}")
