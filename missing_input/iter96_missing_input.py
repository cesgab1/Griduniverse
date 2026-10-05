"""Iteration 96: missing measured input? a0_eff = a0 (X/Xref)^n, fitted inside SPARC, extrapolated to clusters (PREREG_96.md)."""
import numpy as np, glob
from scipy.optimize import minimize_scalar
G = 4.30e-6; conv = 3.086e19/1e6
nu = lambda y: 1/(1 - np.exp(-np.sqrt(y)))
def fit_a0(gb, go, w):
    f = lambda la: np.sum(w*(np.log10(nu(gb/10**la)*gb) - np.log10(go))**2)
    return minimize_scalar(f, bounds=(-12.5, -8.5), method="bounded").x
rows = []
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(fn, comments="#", ndmin=2); r, V, eV, Vg, Vd, Vb, _, _ = d.T
    ok = (r > 0) & (V > 0); r, V, eV, Vg, Vd, Vb = r[ok], V[ok], eV[ok], Vg[ok], Vd[ok], Vb[ok]
    vb2 = Vg*abs(Vg) + 0.5*Vd**2 + 0.7*Vb**2; ok = vb2 > 0
    if ok.sum() < 6: continue
    r, V, eV, vb2, Vg, Vd, Vb = r[ok], V[ok], eV[ok], vb2[ok], Vg[ok], Vd[ok], Vb[ok]
    Ms = (0.5*Vd[-1]**2 + 0.7*Vb[-1]**2)*r[-1]/G; Mg = 1.33*max(Vg[-1], 0)*abs(Vg[-1])*r[-1]/G
    if Ms < 1e6: continue
    el = 2*np.sqrt(eV**2 + (0.05*V)**2)/V/np.log(10)
    la = fit_a0(vb2/r/conv, V**2/r/conv, 1/el**2)
    Mb = Ms + Mg; R = r[-1]
    rows.append(dict(la=la, M=Mb, v2=np.max(V)**2, R=R, fg=Mg/Mb, S=Mb/(np.pi*(R*1e3)**2)))
la = np.array([x["la"] for x in rows]); la0 = np.median(la)
Mcl, vcl, Rcl = 3e14, 1400.0, 2000.0
cl = dict(M=Mcl, v2=vcl**2, R=Rcl, fg=0.83, S=Mcl/(np.pi*(Rcl*1e3)**2))
names = dict(M="baryonic mass", v2="well depth v^2", R="size", fg="gas fraction", S="surface density")
out = [f"Iteration 96: {len(la)} SPARC galaxies, median a0 {10**la0:.2e} m/s^2, spread {np.std(la):.2f} dex", "",
       f"{'input X':16s} {'SPARC range (dex)':>17s} {'cluster vs median (dex)':>24s} {'SPARC slope n':>15s} {'n needed x5.3':>14s} {'x17':>7s} {'tension x5.3':>13s} {'x17':>7s}"]
for k in names:
    x = np.log10([r_[k] for r_ in rows]); xc = np.log10(cl[k])
    A = np.vstack([x - x.mean(), np.ones_like(x)]).T
    coef, res, *_ = np.linalg.lstsq(A, la, rcond=None); s2 = res[0]/(len(x) - 2)
    n, en = coef[0], np.sqrt(s2/np.sum((x - x.mean())**2))
    dx = xc - np.median(x); rng_ = np.percentile(x, 95) - np.percentile(x, 5)
    nr = [np.log10(t)/dx if abs(dx) > 1e-3 else np.inf for t in (5.3, 17)]
    ts = [(q - n)/en for q in nr]
    out.append(f"{names[k]:16s} {rng_:17.2f} {dx:+24.2f} {n:+8.3f}+/-{en:.3f} {nr[0]:+14.2f} {nr[1]:+7.2f} {ts[0]:+12.1f}s {ts[1]:+6.1f}s")
    pct = np.mean(x >= xc)*100
    out.append(f"{'':16s} (clusters sit at the {100-pct:.0f}th percentile of galaxies for this input)")
txt = "\n".join(out); print(txt); open("iter96_missing_input.txt", "w").write(txt + "\n")
