"""Iteration 93: the bare universe -- atoms, light, neutrinos; no dark energy, no dark matter (PREREG_93.md)."""
import os, numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.optimize import minimize_scalar
src = open(os.path.join(os.path.dirname(__file__), "../cosmology_fits/fit_law.py")).read()
exec(src.split("def chi2")[0])
GYR = 977.792; WB = 0.0224

def build(model, h):
    Or = W_R / h**2; zp = 1 + ZG
    if model == "LCDM":   Om, Ok = (0.02237 + 0.1200 + W_NU * h**2) / h**2, 0.0; OL = 1 - Om - Or; wb = 0.02237
    elif model == "BARE-OPEN": Om = (WB + W_NU) / h**2; Ok = 1 - Om - Or; OL = 0.0; wb = WB
    else:                 Om = 1 - Or; Ok = 0.0; OL = 0.0; wb = Om * h**2 - W_NU   # BARE-FLAT: atoms fill to critical
    E = np.sqrt(Om * zp**3 + Or * zp**4 + Ok * zp**2 + OL)
    H0 = 100 * h; DHub = C_KMS / H0
    chi = cumulative_trapezoid(1 / E, ZG, initial=0) * DHub
    DMt = DHub / np.sqrt(Ok) * np.sinh(np.sqrt(Ok) * chi / DHub) if Ok > 1e-8 else chi
    wm = Om * h**2
    Rb = 3 * wb / (4 * W_GAMMA) / zp; csH = C_KMS / np.sqrt(3 * (1 + Rb)) / (H0 * E)
    def horizon(zstart):
        sel = ZG >= zstart; zz = np.r_[zstart, ZG[sel]]; ff = np.r_[np.interp(zstart, ZG, csH), csH[sel]]
        return np.trapezoid(ff, zz) + csH[sel][-1] * (1 + ZG[-1])
    g1 = 0.0783 * wb**-0.238 / (1 + 39.5 * wb**0.763); g2 = 0.560 / (1 + 21.1 * wb**1.81)
    zs = 1048 * (1 + 0.00124 * wb**-0.738) * (1 + g1 * wm**g2)                       # Hu & Sugiyama
    b1 = 0.313 * wm**-0.419 * (1 + 0.607 * wm**0.674); b2 = 0.238 * wm**0.223
    zd = 1291 * wm**0.251 / (1 + 0.659 * wm**0.828) * (1 + b1 * wb**b2)              # Eisenstein & Hu drag epoch
    rd = horizon(zd); lA = np.pi * np.interp(zs, ZG, DMt) / horizon(zs) * LA_CAL
    age = np.trapezoid(1 / (zp * E), ZG) * GYR / H0
    return dict(E=E, DMt=DMt, DHub=DHub, rd=rd, lA=lA, age=age, Om=Om, Ok=Ok, wb=wb, h=h, zs=zs, zd=zd)

RD_CAL = 147.09 / build("LCDM", 0.6736)["rd"]       # one calibration of the integrated sound horizon, applied to all

def chis(o):
    DM = lambda z: np.interp(z, ZG, o["DMt"]); DH = lambda z: o["DHub"] / np.interp(z, ZG, o["E"]); rd = o["rd"] * RD_CAL
    pred = [DM(z) / rd if q == "DM_over_rs" else DH(z) / rd if q == "DH_over_rs" else (z * DM(z)**2 * DH(z))**(1/3) / rd
            for z, q in zip(bao.z, bao.q)]
    d = bao.v.values - np.array(pred); cb = d @ bao_icov @ d
    mu = 5 * np.log10((1 + zhel) * DM(z_sn)) + 25; r = mb - mu
    B = (sn_icov @ r).sum(); cs = r @ sn_icov @ r - B**2 / sn_icov_sum
    return cb, cs, (o["lA"] - cmb_mean[1]) / cmb_sig[1]

def growth_since(o, z_from):
    Om, h, E = o["Om"], o["h"], o["E"]; lna_g = -np.log(1 + ZG)
    H = lambda x: np.exp(np.interp(-x, -lna_g, np.log(E)))
    def rhs(x, y):
        a = np.exp(x); e = 1e-4; dl = (np.log(H(x + e)) - np.log(H(x - e))) / (2 * e)
        return [y[1], -(2 + dl) * y[1] + 1.5 * Om * a**-3 / H(x)**2 * y[0]]
    a0 = 1 / (1 + z_from); s = solve_ivp(rhs, [np.log(a0), 0], [a0, a0], rtol=1e-8, atol=1e-14)
    return s.y[0, -1] / a0

out = [f"Iteration 93 -- bare universe (atoms + light + neutrinos; nothing dark)   SNSET={SNSET}", ""]
ref = build("LCDM", 0.6736); cbL, csL, lAL = chis(ref)
out.append(f"LCDM reference (Planck 2018 params, not refitted): BAO chi2 {cbL:.1f}, SN chi2 {csL:.1f}, l_A pull {lAL:+.1f} sigma, age {ref['age']:.2f} Gyr")
for model in ["BARE-OPEN", "BARE-FLAT"]:
    f = lambda h: sum(chis(build(model, h))[:2]) + chis(build(model, h))[2]**2
    h = minimize_scalar(f, bounds=(0.3, 1.0), method="bounded").x
    o = build(model, h); cb, cs, lp = chis(o)
    hs = minimize_scalar(lambda h: chis(build(model, h))[1], bounds=(0.3, 1.0), method="bounded").x
    csn = chis(build(model, hs))[1]
    out += ["", f"== {model} ==",
            f"  best combined h = {h:.3f}: Om = {o['Om']:.4f}, Omega_k = {o['Ok']:.4f}, omega_b = {o['wb']:.4f}",
            f"  BAO chi2 {cb:10.1f}  (LCDM {cbL:.1f})  -> delta {cb-cbL:+.1f}",
            f"  SN  chi2 {csn:10.1f}  at SN-only best h={hs:.3f}  -> delta vs LCDM(Planck) {csn-csL:+.1f}",
            f"  CMB acoustic angle l_A = {o['lA']:.1f} vs 301.46 +/- 0.09  -> {lp:+.0f} sigma",
            f"  sound horizon r_d = {o['rd']*RD_CAL:.1f} Mpc (LCDM 147.1)"]
    for hh in (0.674, 0.73):
        oo = build(model, hh)
        out.append(f"  at h={hh}: age {oo['age']:.2f} Gyr (stars 13.3+/-0.5 -> {(13.3-oo['age'])/0.5:+.1f} sigma), Omega_k {oo['Ok']:.3f}, omega_b {oo['wb']:.3f}")
    oo = build(model, 0.674)
    out.append(f"  flatness: Omega_k {oo['Ok']:.3f} vs measured 0.002 +/- 0.0013 -> {(oo['Ok']-0.002)/0.0013:.0f} sigma")
    out.append(f"  atoms: omega_b {oo['wb']:.4f} vs BBN 0.0222 +/- 0.0005 -> {(oo['wb']-0.0222)/0.0005:.0f} sigma")
    G = growth_since(oo, oo["zs"])
    out.append(f"  growth of lumps since recombination (z={oo['zs']:.0f}): x{G:.0f};  seeds ~1e-5..1e-4 -> today {1e-5*G:.1e}..{1e-4*G:.1e} (needed ~1)")
G_L = growth_since(ref, ref["zs"]) ; out.append(f"\n(LCDM: dark-matter lumps grow x{G_L:.0f} after recombination, and had a head start before it)")
txt = "\n".join(out); print(txt)
open(os.path.join(os.path.dirname(__file__), f"iter93_bare_{SNSET}.txt"), "w").write(txt + "\n")
