"""
Heated-halo equilibrium (follow-up to valve_vibration.py).
The shaken fluid is a fluid: once heated to a speed spread sigma_h it settles into hydrostatic (isothermal) equilibrium,
    rho_f(r) = rho_bg * exp(DeltaPhi(r) / sigma_h^2),
where rho_bg is the heated fluid density in the surroundings and DeltaPhi is the well depth relative to the surroundings.
This replaces the crude 'P(v < V_esc at R500)' rule: it automatically counts the deeper binding near halo centres.
Potential: MOND from observed baryons (nu(y) = 1/(1-exp(-sqrt y)), a0 = 1.2e-10) + Newtonian self-gravity of the fluid
(solved iteratively). MOND's logarithmic potential is cut off by the external field of large-scale structure, g_e ~ 0.02-0.03 a0
(literature estimates for SPARC), beyond which the well flattens.  Units: kpc, Msun, km/s.
Outputs: extra fluid mass within R500 as a fraction of the cosmic share (target values from gas_retention_test),
and for a Milky-Way-like galaxy the fluid within 30 kpc (SPARC) and 300 kpc (KiDS lensing) relative to its baryons.
"""
import numpy as np
from scipy.optimize import brentq
G = 4.30e-6; a0 = 1.2e-10*3.086e19/1e6; rhoc = 128.0; Of = 0.261; fcos = 5.4
r = np.logspace(-1, 4.3, 1500)                                     # 0.1 kpc .. 20 Mpc
nu = lambda y: 1/(1 - np.exp(-np.sqrt(np.maximum(y, 1e-12))))
target = {12.0: 0.0, 13.0: 0.02, 13.5: 0.54, 14.0: 0.71, 14.5: 0.57, 15.0: 0.38}
fb_obs = {12.0: 0.03, 13.0: 0.097, 13.5: 0.083, 14.0: 0.089, 14.5: 0.116, 15.0: 0.167}

def hernq(M, a): return M*r**2/(r + a)**2
def baryons(lM):
    M500 = 10**lM; R500 = (3*M500/(4*np.pi*500*rhoc))**(1/3)
    if lM < 12.5:                                                   # Milky-Way-like: stars + hot gas (Brouwer-style, rho~r^-2 to 143 kpc)
        Ms = 5e10; Mb = hernq(Ms, 3.0) + Ms*np.minimum(r, 143)/143
    else:
        Mb500 = fb_obs[lM]*M500; fst = np.interp(lM, [13, 15], [0.35, 0.10]); rc = 0.1*R500
        gas = lambda x: x - rc*np.arctan(x/rc)                     # beta = 2/3 profile, M(<r)
        Mg = (1 - fst)*Mb500*gas(np.minimum(r, 3*R500))/gas(R500)
        Mb = Mg + hernq(fst*Mb500, 0.015*R500)
    return M500, R500, Mb

def well(gtot, ge):
    re = r[np.argmax(gtot < ge)] if np.any(gtot < ge) else r[-1]   # external field takes over where internal g < g_e
    g = np.where(r < re, gtot, 0.0)
    dphi = np.cumsum((g*np.gradient(r))[::-1])[::-1]               # integral from r to re
    return dphi + np.interp(re, r, gtot)*re                        # Newtonian-like tail beyond re

def solve(lM, sig, ge_frac=0.025, env=1.0, it=60):
    M500, R500, Mb = baryons(lM)
    gN = G*Mb/r**2; gb = nu(gN/a0)*gN; rho_bg = env*Of*rhoc
    Mf = np.zeros_like(r)
    for _ in range(it):
        gt = gb + G*Mf/r**2
        x = np.minimum(well(gt, ge_frac*a0)/sig**2, 60)
        rho = rho_bg*(np.exp(x) - 1)                                # excess over the surroundings
        Mnew = np.cumsum(4*np.pi*r**2*rho*np.gradient(r))
        if np.any(~np.isfinite(Mnew)) or Mnew[-1] > 1e17: return None # no equilibrium (fluid collapses: isothermal catastrophe)
        Mf = 0.5*Mf + 0.5*Mnew
    return M500, R500, Mb, Mf

def frac(lM, sig, **kw):
    s = solve(lM, sig, **kw)
    if s is None: return np.inf
    M500, R500, Mb, Mf = s
    return np.interp(R500, r, Mf)/(fcos*np.interp(R500, r, Mb))

if __name__ == "__main__":
    sigs = [250, 300, 350, 400, 500, 600, 800]
    for ge in (0.015, 0.025, 0.04):
      for env in (1.0, 3.0):
        print(f"\nexternal field g_e = {ge} a0, surroundings = {env:.0f} x mean fluid density")
        print("  sigma_h | fraction of cosmic share within R500 at log M = " + " ".join(f"{k:>5}" for k in target) + " | MW fluid/baryons <30kpc <300kpc | chi2")
        best = None
        for s in sigs:
            fr = [frac(lM, s, ge_frac=ge, env=env) for lM in target]
            sol = solve(12.0, s, ge_frac=ge, env=env)
            if sol is None: mw = (np.inf, np.inf)
            else: _, _, Mb, Mf = sol; mw = (np.interp(30, r, Mf)/np.interp(30, r, Mb), np.interp(300, r, Mf)/np.interp(300, r, Mb))
            chi = sum(((min(f, 3) - t)/0.15)**2 for f, t in zip(fr, target.values()))
            print(f"  {s:6d}  | {'':48s}" + " ".join(f"{min(f,99):5.2f}" for f in fr) + f" | {min(mw[0],99):6.2f} {min(mw[1],99):7.2f} | {chi:7.1f}")
    print("\nTargets (as measured):        " + " ".join(f"{t:5.2f}" for t in target.values()))
    print("Targets (20% hydrostatic bias): 0.00  0.50  1.09  1.23  0.97  0.65")
    print("MW limits: SPARC needs fluid/baryons << 1 inside 30 kpc; KiDS (MOND + hot gas fits, chi2 103) tolerates little at 300 kpc.")

# ---- Supply-limited version. 'No equilibrium' (99 above) means the heated fluid stays bound as a self-gravitating halo:
# the system simply keeps what it had, i.e. its full cosmic share (fraction 1). It cannot gather more than its supply.
# Real systems sit in different surroundings: average over external fields g_e = 0.01-0.05 a0 and surroundings 1-3x mean.
def pop_frac(lM, sig):
    vals = [min(frac(lM, sig, ge_frac=ge, env=en), 1.0) for ge in np.geomspace(0.01, 0.05, 5) for en in (1.0, 2.0, 3.0)]
    return np.mean(vals)
if __name__ == "__main__":
    tb = {12.0: 0.0, 13.0: 0.50, 13.5: 1.0, 14.0: 1.0, 14.5: 0.97, 15.0: 0.65}
    print("\nSupply-limited, averaged over surroundings:")
    print("  sigma_h | log M = " + " ".join(f"{k:>5}" for k in target) + " | chi2 (as measured) | chi2 (20% bias)")
    for s in (250, 300, 350, 400, 450, 500, 550, 600, 700):
        fr = [pop_frac(lM, s) for lM in target]
        c1 = sum(((f - t)/0.15)**2 for f, t in zip(fr, target.values())); c2 = sum(((f - t)/0.15)**2 for f, t in zip(fr, tb.values()))
        print(f"  {s:6d}  |        " + " ".join(f"{f:5.2f}" for f in fr) + f" | {c1:8.1f}           | {c2:6.1f}")
