"""
ITERATION 24, part 1 (PREDICTIONS ONLY, no cosmological data touched): the 'trampoline' idea.
Idea (Synthesis + Coalesce, Oct 2026): the Ocean (cold dark matter, in the gaps) is never drained (the sky says DM is
constant to a few % since the CMB), so dark energy cannot be paid out of it. What CAN change is how the Ocean is ARRANGED:
when it gathers into Pools (halos), the Pools sag the sheets, and the sag adds tension. Then 'why now' = 'after the Pools
formed'.  This script turns that into a dark-energy history BEFORE any fit.

Free choices (counted, per the self-audit rules):
  C1 what measures 'sag'      : T1 collapsed mass fraction f_coll(>M_min)   [tension ~ how much Ocean sits in Pools]
                                T2 sag energy density rho_m * <|Phi|>/c^2   [tension ~ total depth of all dents]
  C2 halo mass function       : Sheth-Tormen (standard)
  C3 smallest Pool that counts: M_min = 1e10 Msun/h (galaxy-hosting); sensitivity 1e8 and 1e12 shown
  C4 halo density             : 200 x mean matter density
  C5 growth                   : LCDM linear growth (the sag history barely feeds back; checked by the size, part C)
Amplitude (the SIZE) is NOT predicted unless an amplifier is derived; part C reports the bare size.
Outputs: w(z), crossing redshifts, bin-averaged w in the pre-registration bins, and the bare size ratio.
"""
import numpy as np, json, camb
from scipy.integrate import solve_ivp, quad

h = 0.6736; Om0 = 0.3153; ombh2, omch2 = 0.02237, 0.1200
pars = camb.set_params(H0=100*h, ombh2=ombh2, omch2=omch2, mnu=0.06, As=2.1e-9, ns=0.9649)
pars.set_matter_power(redshifts=[0.0], kmax=200.0)
res = camb.get_results(pars)
kh, _, pk = res.get_matter_power_spectrum(minkh=1e-4, maxkh=200, npoints=2000)
pk = pk[0]; s8 = res.get_sigma8_0()
lnk = np.log(kh)

def sigma_R(R):   # R in Mpc/h, z = 0, top hat
    x = kh*R; W = 3*(np.sin(x) - x*np.cos(x))/x**3
    return np.sqrt(np.trapezoid(kh**3*pk*W**2/(2*np.pi**2), lnk))

rho_m0 = 2.775e11*Om0                    # (Msun/h) / (Mpc/h)^3
lnM = np.linspace(np.log(1e6), np.log(1e16), 400); M = np.exp(lnM)
R = (3*M/(4*np.pi*rho_m0))**(1/3)
sig0 = np.array([sigma_R(r) for r in R])
dlns_dlnM = np.gradient(np.log(sig0), lnM)

# LCDM linear growth D(a), D(1) = 1
OL = 1 - Om0
def Ea(a): return np.sqrt(Om0*a**-3 + OL)
def grow(lna, y):
    a = np.exp(lna); E2 = Om0*a**-3 + OL; dlnE = -1.5*Om0*a**-3/E2
    return [y[1], -(2 + dlnE)*y[1] + 1.5*Om0*a**-3/E2*y[0]]
lna_g = np.linspace(np.log(1e-3), 0, 3000)
sol = solve_ivp(grow, [lna_g[0], 0], [1e-3, 1e-3], t_eval=lna_g, rtol=1e-9)
Dg = sol.y[0]/sol.y[0][-1]
D = lambda a: np.interp(np.log(a), lna_g, Dg)

dc = 1.686; A_, a_, p_ = 0.3222, 0.707, 0.3
def st_fnu(nu):  # Sheth-Tormen mass fraction per ln(nu)... returns dF/dlnM given nu and slope
    nup = np.sqrt(a_)*nu
    return A_*np.sqrt(2/np.pi)*nup*(1 + nup**(-2*p_))*np.exp(-nup**2/2)

G = 4.30091e-9    # Mpc (km/s)^2 / Msun
c_kms = 299792.458
def measures(a, Mmin):
    sig = sig0*D(a); nu = dc/sig
    dF = st_fnu(nu)*np.abs(dlns_dlnM)             # mass fraction per ln M
    sel = M >= Mmin
    fcoll = np.trapezoid(dF[sel], lnM[sel])
    # potential depth of a halo of mass M at 200 x mean density, physical units (M in Msun: divide /h)
    rho_phys = rho_m0*h**2*a**-3                  # Msun / Mpc^3 (physical)
    Mphys = M/h; Rv = (3*Mphys/(4*np.pi*200*rho_phys))**(1/3)
    phi = G*Mphys/Rv/c_kms**2                      # |Phi|/c^2
    phibar = np.trapezoid((dF*phi)[sel], lnM[sel])  # mass-weighted over ALL matter (uncollapsed matter contributes 0)
    return fcoll, phibar

out = ["ITERATION 24 part 1: trampoline (Pools sag the sheets) -- predicted dark-energy history BEFORE fitting", "",
       f"Linear power from CAMB (Planck 2018 LCDM), sigma8 = {s8:.4f}; Sheth-Tormen halos; growth = LCDM", ""]
z = np.linspace(0, 6, 1201); a = 1/(1 + z)
bins = [(0.1, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.1), (1.1, 1.6), (1.6, 2.1)]
preds = {}
for Mmin in (1e8, 1e10, 1e12):
    F = np.array([measures(ai, Mmin) for ai in a])
    fc, pb = F[:, 0], F[:, 1]
    for name, S in (("T1_fcoll", fc), ("T2_sag", Om0*a**-3*pb)):
        lnS = np.log(S); w = -1 - np.gradient(lnS, np.log(a))/3
        cross = [float(0.5*(z[i] + z[i+1])) for i in range(len(z)-1) if (w[i] + 1)*(w[i+1] + 1) < 0]
        wb = {f"{lo}-{hi}": float(np.mean(w[(z >= lo) & (z < hi)])) for lo, hi in bins}
        key = f"{name}_Mmin{Mmin:.0e}"
        preds[key] = dict(w0=float(w[0]), crossing=cross, w_bins=wb, S_over_S0={f"{zz}": float(np.interp(zz, z, S/S[0])) for zz in (0.5, 1, 2, 3, 5)})
        out.append(f"{key:22s} w0 {w[0]:+.3f}  crossing {np.round(cross, 2).tolist()}  bins " + " ".join(f"{k}:{v:+.3f}" for k, v in wb.items()))
        out.append(f"{'':22s} rho_DE(z)/rho_DE(0): " + "  ".join(f"z={k}: {v:.3g}" for k, v in preds[key]['S_over_S0'].items()))
        if Mmin == 1e10:
            np.savetxt(f"iter24_shape_{name}.txt", np.c_[z, S/S[0]], header="z  rho_DE/rho_DE0 (pre-registered shape, Mmin=1e10)")
    if Mmin == 1e10:
        size_ratio = Om0*pb[0]/OL
        out += ["", f"C. bare SIZE today (T2, Mmin 1e10): rho_m <|Phi|>/c^2 = {Om0*pb[0]:.2e} of critical; dark energy needs {OL:.2f}",
                f"   -> the sag is {size_ratio:.1e} of what is needed; an amplifier of ~{1/size_ratio:.0e} is required (not derived)", ""]
out += ["Reference (pre-registration sheet): Claim 1 crossing 0.68, w0 -0.911, bins -0.947 -0.980 -1.001 -1.022 -1.045 -1.061",
        "Data direction (DESI DR2 era): w > -1 today, w < -1 at z ~ 1 (dark energy PEAKED in the past and is now fading)."]
txt = "\n".join(out); print(txt)
open("iter24_part1_predictions.txt", "w").write(txt + "\n")
json.dump(preds, open("iter24_predictions.json", "w"), indent=1)
