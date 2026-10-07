"""PREREG 132: does the extra-pull threshold scale with H(z)?  C: a_t = a0*E(z);  K: a_t = a0."""
import numpy as np, pandas as pd
from scipy.optimize import brentq
G=6.674e-11; Msun=1.989e30; kpc=3.086e19
E = lambda z, Om=0.3: np.sqrt(Om*(1+z)**3 + 1-Om)
a0 = 1.2e-10; a0_err = 0.26e-10
out=[]; P=lambda s="": (print(s), out.append(s))
def fit_a(gbar, gobs):
    if gobs <= gbar*1.0001: return 0.0
    f = lambda la: gbar/(1-np.exp(-np.sqrt(gbar/10**la))) - gobs
    return 10**brentq(f, -14, -6)
P("(1) Our own: Genzel+2017 six galaxies (enclosed baryons = half of total at R1/2; free choice)")
g = pd.read_csv("data/genzel2017_nature_table1.csv", comment="#")
rows=[]
for _, r in g.iterrows():
    for lab, Mb in (("prior", r.Mbar_prior_1e11), ("fit", r.Mbar_fit_1e11)):
        R = r.R12_H_kpc*kpc; gb = G*0.5*Mb*1e11*Msun/R**2; go = (r.vc_R12_kms*1e3)**2/R
        rows.append((r.galaxy, r.z, lab, gb, go, fit_a(gb, go)))
d = pd.DataFrame(rows, columns="gal z mass gbar gobs a".split())
for _, r in d[d.mass=="prior"].iterrows():
    P(f"   {r.gal:11s} z={r.z:.2f}: g_bar={r.gbar:.1e}, g_obs={r.gobs:.1e} (ratio {r.gobs/r.gbar:.2f}) -> a = {r.a:.1e}"
      f"   [C predicts {a0*E(r.z):.1e}, K {a0:.1e}]")
P("   note: g_bar = 0.9e-10 to 1.1e-9 m/s^2 = 1-9x a0 (corrected; an earlier draft said 10-30x)")
for lab in ("prior","fit"):
    s = d[d.mass==lab]
    P(f"   [{lab} masses] median g_obs/g_bar = {np.median(s.gobs/s.gbar):.2f}; galaxies with g_obs<=g_bar: {(s.gobs<=s.gbar*1.0001).sum()}/6")
rng=np.random.default_rng(132)
s0=d[d.mass=="prior"].reset_index(drop=True); gg=g.reset_index(drop=True)
P("   Monte Carlo (baryonic mass +- quoted error, speed +-10%, assumed), bins z<2 and z>2:")
chi={"C":0.0,"K":0.0}
for lab, sel in (("z<2", gg.z<2), ("z>2", gg.z>2)):
    meds=[]
    for _ in range(3000):
        aa=[]
        for i in np.where(sel)[0]:
            r=gg.iloc[i]; Mb=max(0.05, rng.normal(r.Mbar_prior_1e11, r.Mbar_prior_err)); v=r.vc_R12_kms*rng.normal(1,0.1)
            R=r.R12_H_kpc*kpc; gb=G*0.5*Mb*1e11*Msun/R**2; go=(v*1e3)**2/R; aa.append(max(fit_a(gb,go),1e-12))
        meds.append(np.median(np.log10(aa)))
    m, e = np.median(meds), np.std(meds); zc = gg.z[sel].median()
    for k, pred in (("C", np.log10(a0*E(zc))), ("K", np.log10(a0))):
        chi[k] += ((m-pred)/e)**2
    P(f"     {lab}: median log a = {m:.2f} +/- {e:.2f} (a = {10**m:.1e}); C {a0*E(zc):.1e}, K {a0:.1e}")
P(f"   chi^2: C = {chi['C']:.1f}, K = {chi['K']:.1f} (2 bins) -> delta = {chi['K']-chi['C']:+.1f} "
  f"({'prefer C' if chi['K']-chi['C']>4 else 'prefer K' if chi['C']-chi['K']>4 else 'no preference (|delta|<4)'})")

P("\n(2) Published direct measurement: Ciocan+2026 (MUSE-DARK III, 79 galaxies, 0.33<z<1.44)")
for z in (0.5, 1.0, 1.44):
    lin = 1.0e-10 + 1.59e-10*z
    P(f"   z={z:.2f}: measured (linear fit) {lin:.2e};  C with SPARC anchor {a0*E(z):.2e};  C with their own a(0)=1.0e-10 "
      f"{1.0e-10*E(z):.2e};  K {a0:.1e}")
meas, merr = 2.38e-10, 0.055e-10   # 95% CI ~ +0.12/-0.10 -> 1-sigma ~0.055
for lab, pred, perr in (("C, SPARC anchor", a0*E(1.0), a0_err*E(1.0)), ("C, own anchor", 1.0e-10*E(1.0), 0.02e-10*E(1.0)),
                        ("K, SPARC anchor", a0, a0_err), ("K, own anchor", 1.0e-10, 0.02e-10)):
    P(f"   z~1: {lab:16s} predicts {pred:.2e} -> deviation {(meas-pred)/np.hypot(merr,perr):+.1f} sigma")

P("\n(3) Published baryonic Tully-Fisher zero point vs local (dex in baryonic mass at fixed speed; C predicts -log E(z))")
for lab, z, dB, err in (("Ubler+2017", 0.9, -0.44, 0.06), ("Ubler+2017", 2.3, -0.27, 0.07), ("Jeanneau+2026 (MUSE-DARK II)", 1.0, 0.00, 0.06)):
    pC = -np.log10(E(z))
    P(f"   {lab:28s} z={z}: measured {dB:+.2f} +/- {err:.2f};  C {pC:+.2f} ({(dB-pC)/err:+.1f} sigma);  K +0.00 ({dB/err:+.1f} sigma)")
P("   (Ubler errors approximated from quoted zero-point errors; gas treatment differs between studies)")

P("\n(4) Milgrom 2017 (Genzel z~2 curves): 'all but exclude ~4 a0 at z~2'; C predicts {:.1f} a0 at z=2.2; Ciocan linear"
  " extrapolation {:.1f} a0".format(E(2.2), (1.0+1.59*2.2)/1.2))
open("RESULT_132_numbers.txt","w").write("\n".join(out)+"\n")
