"""
Slack regions of the grid vs 175 SPARC galaxy rotation curves.
Tension picture: G_eff = G * T0/T_local. Where the grid is slacker, the same mass pulls harder (and, with gamma = 1, bends light equally).
Observed acceleration g_obs = V_obs^2/r; Newtonian pull of the visible matter g_bar = (Vgas|Vgas| + Yd Vdisk^2 + Yb Vbul^2)/r.
If slack explains the rotation curves: g_obs = g_bar * S, where S = T0/T_local is the slack factor.
Three rules for where the grid is slack:
  A. Random patches, unrelated to matter           -> S scatters freely from galaxy to galaxy / radius to radius
  B. Slack where visible matter is dense            -> S rises with the local surface density of stars
  C. Slack where the grid is barely strained (weak pull): taut near mass, loose far away -> S = nu(g_bar/a0), one universal a0
Stellar mass-to-light 0.5 (disk) / 0.7 (bulge) at 3.6 micron; points with errV/V > 10% dropped (no quality table available).
"""
import numpy as np, glob, os, json
from scipy.optimize import minimize_scalar, minimize
Yd, Yb = 0.5, 0.7
kpc = 3.0857e19
rows = []
for f in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(f, comments="#")
    if d.ndim == 1: d = d[None]
    r, V, eV, Vg, Vd, Vb, SBd, SBb = d.T
    gbar = (Vg*np.abs(Vg) + Yd*Vd**2 + Yb*Vb**2) * 1e6 / (r*kpc)
    gobs = V**2 * 1e6 / (r*kpc)
    ok = (V > 0) & (eV/V < 0.1) & (gbar > 0) & (r > 0)
    name = os.path.basename(f).replace("_rotmod.dat", "")
    for i in np.where(ok)[0]:
        rows.append((name, r[i], gobs[i], gbar[i], 2*eV[i]/V[i]/np.log(10), Yd*SBd[i] + Yb*SBb[i]))
names = np.array([x[0] for x in rows]); r = np.array([x[1] for x in rows]); go = np.array([x[2] for x in rows])
gb = np.array([x[3] for x in rows]); elog = np.array([x[4] for x in rows]); sig = np.array([x[5] for x in rows])
lgo, lgb = np.log10(go), np.log10(gb); gal = np.unique(names)
print(f"{len(gal)} galaxies, {len(go)} points")
res0 = lgo - lgb
print(f"\nNo slack (Newton with visible matter only): observed pull exceeds it by a median factor {10**np.median(res0):.2f}; scatter {np.std(res0):.2f} dex")
print(f"   at the weakest pulls (g_bar < 1e-11 m/s^2): factor {10**np.median(res0[gb<1e-11]):.1f};   strongest (> 1e-9): {10**np.median(res0[gb>1e-9]):.2f}")
# ---- Rule A: random patches. If S were random per galaxy, the required S would not depend on g_bar; test correlation
c_A = np.corrcoef(lgb, res0)[0, 1]
per_gal = np.array([np.median(res0[names == g]) for g in gal])
print(f"\nA. Random slack patches: required slack vs strength of pull, correlation r = {c_A:.2f} (random patches predict r ~ 0)")
# ---- Rule B: slack tied to stellar density
ms = sig > 0
c_B = np.corrcoef(np.log10(sig[ms]), res0[ms])[0, 1]
print(f"B. Slack where stars are dense: required slack vs stellar surface density, correlation r = {c_B:.2f} (rule B predicts r > 0)")
# ---- Rule C: slack where the pull is weak. Two forms of nu:
nus = {"McGaugh 2016 form  1/(1-exp(-sqrt(y)))": lambda y: 1/(1 - np.exp(-np.sqrt(y))),
       "'simple' form  1/2+sqrt(1/4+1/y)": lambda y: 0.5 + np.sqrt(0.25 + 1/y)}
out = {}
for lab, nu in nus.items():
    def chi(la0): 
        pred = lgb + np.log10(nu(gb/10**la0)); return np.sum(((lgo - pred)/np.hypot(elog, 0.1))**2)
    b = minimize_scalar(chi, bounds=(-11, -9), method="bounded")
    a0 = 10**b.x; resC = lgo - (lgb + np.log10(nu(gb/a0)))
    # intrinsic scatter: total minus observational
    obs = np.median(elog)
    intr = np.sqrt(max(np.std(resC)**2 - np.mean(elog**2), 0))
    trend = np.polyfit(lgb, resC, 1)[0]
    print(f"\nC. Slack where the pull is weak, {lab}:\n   best a0 = {a0:.2e} m/s^2 (one number for all {len(gal)} galaxies)")
    print(f"   residual scatter {np.std(resC):.3f} dex, of which measurement errors explain {np.sqrt(np.mean(elog**2)):.3f}; intrinsic {intr:.3f} dex; leftover trend with pull {trend:+.3f}")
    print(f"   galaxies whose median residual is within 0.1 dex: {np.mean([abs(np.median(resC[names==g]))<0.1 for g in gal])*100:.0f}%")
    out[lab] = dict(a0=a0, scatter=float(np.std(resC)), intrinsic=float(intr))
H0 = 67.5e3/3.0857e22; c = 2.998e8
print(f"\n   compare a0 with c*H0/(2 pi) = {c*H0/(2*np.pi):.2e} m/s^2 (a known coincidence noted by Milgrom)")
# Rule C with a per-galaxy free a0: does a universal a0 suffice?
nu = list(nus.values())[0]; a0u = out[list(nus)[0]]["a0"]
la = []
for g in gal:
    m = names == g
    if m.sum() < 5: continue
    b = minimize_scalar(lambda x: np.sum(((lgo[m] - lgb[m] - np.log10(nu(gb[m]/10**x)))/np.hypot(elog[m], 0.05))**2), bounds=(-12, -8), method="bounded")
    la.append(b.x)
la = np.array(la)
print(f"   per-galaxy best a0: median {10**np.median(la):.2e}, spread {np.std(la):.2f} dex (16-84%: {10**np.percentile(la,16):.1e} to {10**np.percentile(la,84):.1e})")
# binned relation for the page
bins = np.arange(-12.2, -8.4, 0.3); prof = []
for lo in bins:
    m = (lgb >= lo) & (lgb < lo + 0.3)
    if m.sum() > 10: prof.append((lo + 0.15, float(np.median(lgo[m])), float(np.std(lgo[m])), int(m.sum())))
json.dump(dict(out=out, prof=prof, cA=c_A, cB=c_B, n_gal=len(gal), n_pts=len(go), a0_spread=float(np.std(la))),
          open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/slack_sparc.json", "w"), indent=1)
print("\nbinned: log g_bar -> median log g_obs (scatter, N):")
for p in prof: print(f"   {p[0]:6.2f} -> {p[1]:6.2f} ({p[2]:.2f}, {p[3]})")
