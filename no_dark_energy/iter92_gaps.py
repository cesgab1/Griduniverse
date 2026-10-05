"""
Iteration 92: gaps of the dark-energy-free universe (see PREREG_92.md).
Reuses data + E_of_z + l_A calibration from cosmology_fits/fit_law.py (exec of its top part),
adds curvature-aware distances so the OPEN (matter + curvature) universe can be tested.
"""
import os, numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.optimize import minimize
src = open(os.path.join(os.path.dirname(__file__), "../cosmology_fits/fit_law.py")).read()
exec(src.split("def chi2")[0])          # data, ZG, E_of_z, observables, LA_CAL

GYR = 977.792                            # 1/(H0 in km/s/Mpc) in Gyr

def Ez(model, Om, h):
    Or = W_R / h**2; zp = 1 + ZG
    if model == "EdS":  return np.sqrt((1 - Or) * zp**3 + Or * zp**4), 0.0
    if model == "OPEN": Ok = 1 - Om - Or; return np.sqrt(Om * zp**3 + Or * zp**4 + Ok * zp**2), Ok
    return E_of_z(model, Om, h, []), 0.0

def obs(model, p):
    if model == "EdS": h, wb = p; Om = 1 - W_R / h**2
    else: Om, h, wb = p
    E, Ok = Ez(model, Om, h)
    if not np.all(np.isfinite(E)) or np.any(E <= 0): return None
    H0 = 100 * h; DHub = C_KMS / H0
    chi = cumulative_trapezoid(1 / E, ZG, initial=0) * DHub
    DMt = DHub / np.sqrt(Ok) * np.sinh(np.sqrt(Ok) * chi / DHub) if Ok > 1e-8 else chi
    DM = lambda z: np.interp(z, ZG, DMt); DH = lambda z: DHub / np.interp(z, ZG, E)
    wm = Om * h**2; wcb = wm - W_NU
    rd = 55.154 * np.exp(-72.3 * (W_NU + 0.0006)**2) / (wcb**0.25351 * wb**0.12807)
    g1 = 0.0783 * wb**-0.238 / (1 + 39.5 * wb**0.763); g2 = 0.560 / (1 + 21.1 * wb**1.81)
    zs = 1048 * (1 + 0.00124 * wb**-0.738) * (1 + g1 * wm**g2)
    Rb = 3 * wb / (4 * W_GAMMA) / (1 + ZG); csH = C_KMS / np.sqrt(3 * (1 + Rb)) / (H0 * E)
    sel = ZG >= zs; zz = np.r_[zs, ZG[sel]]; ff = np.r_[np.interp(zs, ZG, csH), csH[sel]]
    rs = np.trapezoid(ff, zz) + csH[sel][-1] * (1 + ZG[-1])
    R = np.sqrt(Om) * H0 * DM(zs) / C_KMS; lA = np.pi * DM(zs) / rs * LA_CAL
    # age: t0 = int dz / ((1+z) H)
    age = np.trapezoid(1 / ((1 + ZG) * E), ZG) * GYR / H0
    return dict(DM=DM, DH=DH, rd=rd, R=R, lA=lA, wb=wb, Om=Om, h=h, E=E, age=age)

def parts(model, p, use=("bao", "sn", "cmb")):
    o = obs(model, p)
    if o is None: return None
    out = {}
    pred = []
    for z, q in zip(bao.z, bao.q):
        if q == "DM_over_rs": pred.append(o["DM"](z) / o["rd"])
        elif q == "DH_over_rs": pred.append(o["DH"](z) / o["rd"])
        else: pred.append((z * o["DM"](z)**2 * o["DH"](z))**(1 / 3) / o["rd"])
    d = bao.v.values - np.array(pred); out["bao"] = d @ bao_icov @ d
    mu = 5 * np.log10((1 + zhel) * o["DM"](z_sn)) + 25
    r = mb - mu; A = r @ sn_icov @ r; B = (sn_icov @ r).sum(); out["sn"] = A - B**2 / sn_icov_sum
    dc = np.array([o["R"], o["lA"], o["wb"]]) - cmb_mean; out["cmb"] = dc @ cmb_icov @ dc
    out["tot"] = sum(out[k] for k in use)
    return out, o

def fit(model, use, starts):
    lo = dict(EdS=[(0.2, 1.0), (0.005, 0.05)], OPEN=[(0.05, 1.0), (0.2, 1.0), (0.005, 0.05)])
    bounds = lo.get(model, [(0.1, 0.6), (0.5, 0.9), (0.018, 0.027)])
    def f(p):
        if any(not (a < x < b) for x, (a, b) in zip(p, bounds)): return 1e12
        r = parts(model, p, use); return 1e12 if r is None else r[0]["tot"]
    best = None
    for x0 in starts:
        r = minimize(f, x0, method="Nelder-Mead", options=dict(maxiter=8000, xatol=1e-7, fatol=1e-5))
        r = minimize(f, r.x, method="Nelder-Mead", options=dict(maxiter=8000, xatol=1e-8, fatol=1e-6))
        if best is None or r.fun < best.fun: best = r
    return best.x, parts(model, best.x, use)

def growth(o):
    """linear growth D(a), normalised D = a deep in matter era; return D(1)/1 (suppression vs pure matter growth)."""
    Om, h, E = o["Om"], o["h"], o["E"]; Or = W_R / h**2
    lnE = np.log(E); lna_g = -np.log(1 + ZG)
    def H(lna): return np.exp(np.interp(-lna, -lna_g, lnE))   # lna_g decreasing -> flip
    def rhs(lna, y):
        a = np.exp(lna); e = H(lna); eps = 1e-4
        dlnH = (np.log(H(lna + eps)) - np.log(H(lna - eps))) / (2 * eps)
        Oma = Om * a**-3 / e**2
        return [y[1], -(2 + dlnH) * y[1] + 1.5 * Oma * y[0]]
    a0 = 1e-3
    s = solve_ivp(rhs, [np.log(a0), 0], [a0, a0], rtol=1e-8, atol=1e-12)
    return s.y[0, -1]

starts3 = [[0.31, 0.68, 0.0224], [0.30, 0.69, 0.0223]]
startsE = [[0.5, 0.0224], [0.7, 0.0224], [0.35, 0.0224]]
startsO = [[0.3, 0.68, 0.0224], [0.6, 0.5, 0.0224], [0.9, 0.4, 0.0224]]
models = [("LCDM", starts3), ("LAW", starts3), ("EdS", startsE), ("OPEN", startsO)]

lines = [f"Iteration 92 -- gaps of the dark-energy-free universe  (SNSET={SNSET})", ""]
# --- (A) supernova shape alone ---
lines.append("(A) Supernova distance SHAPE alone (absolute magnitude marginalised; omega_b irrelevant)")
snres = {}
for m, st in models:
    x, (pp, o) = fit(m, ("sn",), st); snres[m] = pp["sn"]
    lines.append(f"  {m:5s} chi2_SN = {pp['sn']:8.2f}   Om = {o['Om']:.3f}")
for m in ["LAW", "EdS", "OPEN"]:
    lines.append(f"  {m:5s} minus LCDM: delta chi2_SN = {snres[m] - snres['LCDM']:+8.2f}")
lines.append("")
# --- (B) all three data sets combined ---
lines.append("(B) BAO + SN + CMB combined")
comb = {}
for m, st in models:
    x, (pp, o) = fit(m, ("bao", "sn", "cmb"), st); comb[m] = (pp, o)
    D = growth(o); comb[m] = (pp, o, D)
    lines.append(f"  {m:5s} chi2 = {pp['tot']:9.2f} (BAO {pp['bao']:8.2f}, SN {pp['sn']:8.2f}, CMB {pp['cmb']:8.2f})"
                 f"  Om={o['Om']:.3f} h={o['h']:.3f} age={o['age']:.2f} Gyr  growth D0={D:.3f}")
L = comb["LCDM"]
lines.append("")
lines.append("(C) Gap table (each model at its own best combined fit)")
lines.append(f"  {'model':5s} {'d chi2 tot':>10s} {'d BAO':>9s} {'d SN':>9s} {'d CMB':>9s} {'age gap (sigma)':>16s} {'sigma8 pred':>11s}")
for m in ["LCDM", "LAW", "EdS", "OPEN"]:
    pp, o, D = comb[m]
    agesig = (13.3 - o["age"]) / 0.5
    s8 = 0.811 * D / L[2]
    lines.append(f"  {m:5s} {pp['tot']-L[0]['tot']:10.2f} {pp['bao']-L[0]['bao']:9.2f} {pp['sn']-L[0]['sn']:9.2f} {pp['cmb']-L[0]['cmb']:9.2f}"
                 f" {agesig:+16.2f} {s8:11.3f}")
lines.append("  age gap = (13.3 Gyr oldest stars - model age)/0.5 Gyr; positive = universe younger than its stars")
lines.append("  sigma8 pred = 0.811 x growth(model)/growth(LCDM), same early amplitude; transfer-shape change ignored")
lines.append("")
# --- (D) EdS with H0 forced to local value: the classic age problem in its simplest form ---
for h in (0.674, 0.73):
    lines.append(f"(D) EdS age at h={h}: {2/3*GYR/(100*h):.2f} Gyr  (oldest stars 13.3 +/- 0.5)")
txt = "\n".join(lines); print(txt)
open(os.path.join(os.path.dirname(__file__), f"iter92_gaps_{SNSET}.txt"), "w").write(txt + "\n")
