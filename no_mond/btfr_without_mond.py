"""
R2 (truth T2): baryonic Tully-Fisher relation, v_flat^4 ~ M_baryon, slope ~3.85, small scatter. Can 'ocean = cold dark matter'
(Newtonian gravity, NFW halos placed by abundance matching, no MOND) produce the SAME slope and scatter from the SAME galaxies?
Observed: SPARC galaxies with a flat outer curve; v_flat = mean of the outer 3 points; M_b = stars (M/L 0.5 disk, 0.7 bulge) +
1.33 x gas, from the enclosed mass at the last point.
Model: v at the same outer radii from baryons + NFW halo, halo mass from stellar mass (Moster+2013) with its scatter
(0.15 dex in M* at fixed halo -> converted to halo scatter with the local slope), concentration scatter 0.11 dex
(Dutton & Maccio 2014). 300 random realisations. Also: MOND at fixed M/L for reference.
"""
import numpy as np, glob
from scipy.optimize import brentq
rng = np.random.default_rng(3)
import os
FLAT = float(os.environ.get("FLAT", 0.10)); ERR = float(os.environ.get("ERR", 1.0)); print(f"cuts: flat within {FLAT:.0%}, mean outer error < {ERR:.0%}")
G = 4.30e-6; rhoc = 136.0; a0 = 1.2e-10*3.086e19/1e6
nu = lambda y: 1/(1 - np.exp(-np.sqrt(y)))
def shmr(M): M1, N, b, g = 10**11.59, 0.0351, 1.376, 0.608; return 2*N*M/((M/M1)**-b + (M/M1)**g)
def Mh_of(Ms): return np.exp(brentq(lambda lm: np.log(shmr(np.exp(lm))) - np.log(Ms), np.log(1e8), np.log(1e16)))
def conc(M): return 10**(0.905 - 0.101*np.log10(M*0.7/1e12))
def nfw_v2(r, M, dc=0.0):
    R200 = (3*M/(4*np.pi*200*rhoc))**(1/3); c = conc(M)*10**dc; x = r*c/R200; m = lambda x: np.log(1 + x) - x/(1 + x)
    return G*M*m(x)/m(c)/r
gals = []
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(fn, comments="#", ndmin=2); r, V, eV, Vg, Vd, Vb, _, _ = d.T
    ok = (r > 0) & (V > 0); r, V, eV, Vg, Vd, Vb = r[ok], V[ok], eV[ok], Vg[ok], Vd[ok], Vb[ok]
    if len(r) < 6: continue
    o = slice(-3, None)
    if abs(V[-1] - V[-3])/V[-3] > FLAT or np.mean(eV[o]/V[o]) > ERR: continue   # outer curve not flat / poorly measured
    Ms = (0.5*Vd[-1]**2 + 0.7*Vb[-1]**2)*r[-1]/G; Mg = 1.33*Vg[-1]**2*r[-1]/G
    if Ms < 1e6: continue
    vb2 = Vg*abs(Vg) + 0.5*Vd**2 + 0.7*Vb**2
    gals.append(dict(r=r[o], V=V[o].mean(), Ms=Ms, Mb=Ms + Mg, vb2=vb2[o], Mh=Mh_of(Ms), ev=np.mean(eV[o]/V[o])/np.sqrt(3)/np.log(10)))
N = len(gals); lMb = np.log10([g["Mb"] for g in gals])
def fit(lv, lm):
    s, b = np.polyfit(lv, lm, 1); si = 1/np.polyfit(lm, lv, 1)[0]; sb = np.tan((np.arctan(s) + np.arctan(si))/2)   # forward, inverse, bisector
    res = lv - (lm - np.mean(lm - sb*lv))/sb
    return sb, np.std(res)                                    # slope of log Mb vs log v; scatter in log v (dex)
lVobs = np.log10([g["V"] for g in gals])
so, eo = fit(lVobs, lMb)
lVmond = np.log10([np.mean(np.sqrt(nu(np.maximum(g["vb2"], 1e-9)/g["r"]/a0)*np.maximum(g["vb2"], 1e-9))) for g in gals])
sm, em = fit(lVmond, lMb)
def slope_ms(Mh): e = 1e-3; return (np.log(shmr(Mh*(1 + e))) - np.log(shmr(Mh)))/np.log(1 + e)
S, E, S0, E0 = [], [], None, None
for k in range(301):
    lv = []
    for g in gals:
        if k == 0: dh, dc = 0.0, 0.0
        else: dh, dc = rng.normal(0, 0.15)/slope_ms(g["Mh"]), rng.normal(0, 0.11)
        lv.append(np.log10(np.mean(np.sqrt(g["vb2"] + nfw_v2(g["r"], g["Mh"]*10**dh, dc)))))
    s, e = fit(np.array(lv), lMb)
    if k == 0: S0, E0 = s, e
    else: S.append(s); E.append(e)
S, E = np.array(S), np.array(E)
print(f"R2 baryonic Tully-Fisher, {N} SPARC galaxies with flat outer curves (same galaxies, same baryons for every model; bisector slope of log M_b vs log v)")
print(f"  observed                          slope {so:5.2f}   scatter {eo:.3f} dex in v (includes measurement errors)")
print(f"  MOND, fixed M/L (0 free)          slope {sm:5.2f}   scatter {em:.3f}")
print(f"  CDM halos, abundance matched, no scatter        slope {S0:5.2f}   scatter {E0:.3f}")
print(f"  CDM halos + SHMR & concentration scatter        slope {S.mean():5.2f} +/- {S.std():.2f}   scatter {E.mean():.3f} +/- {E.std():.3f}")
ev = np.array([g["ev"] for g in gals]); print(f"  typical measurement error in log v_flat: {np.median(ev):.3f} dex (v only; distance and M/L errors add to M_b)")
print(f"  literature: Lelli+2016/2019 slope 3.85 +/- 0.09, intrinsic scatter ~0.025 dex in v (0.10 dex in M)")
bs = []
for _ in range(2000):
    i = rng.integers(0, N, N); bs.append(fit(lVobs[i], lMb[i])[0])
print(f"  observed slope bootstrap error: +/- {np.std(bs):.2f}  -> CDM-halo slope is {abs(so - S.mean())/np.hypot(np.std(bs), S.std()):.1f} sigma below")
