"""
ITERATION 12: joint large-angle test with the REAL Planck 2018 low-l TT likelihood (l = 2-29, Gibbs, native python).
Two effects on the same multipoles:
  (a) bounce before inflation: primordial power suppressed on the largest scales,
      P(k) -> P(k) [1 - exp(-(k/k_c)^3.35)]   (standard cutoff form used for pre-inflation bounces; 1 parameter k_c)
  (b) leftover layer lumps (iterations 8-9, with matter response): extra D_l = D_l^(one layer) / N_eff  (1 parameter N_eff)
High-l parameters fixed at the Planck 2018 best fit (only the largest angles are tested here).
Report chi2 differences vs plain LCDM, counting parameters honestly.
"""
import numpy as np, camb, os
from cobaya.likelihoods.planck_2018_lowl.TT import TT
like = TT({"packages_path": "/home/claude/cobaya_packages"})
T0 = 2.7255e6
lump = np.load("lump_cl_one_layer_with_matter.npy")
def dl_with_cut(kc):
    pars = camb.set_params(H0=67.36, ombh2=0.02237, omch2=0.1200, mnu=0.06, tau=0.0544, lmax=60)
    ks = np.geomspace(1e-6, 10, 2000); P = 2.1e-9*(ks/0.05)**(0.9649 - 1)
    if kc > 0: P = P*(1 - np.exp(-(ks/kc)**3.35))
    pars.set_initial_power_table(ks, P, effective_ns_for_nonlinear=0.9649)
    r = camb.get_results(pars); return r.get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=False)["total"][:, 0]
def chi2(dl): return -2*like.log_likelihood(dl)
base = dl_with_cut(0); c0 = chi2(base)
out = ["ITERATION 12: bounce suppression vs layer lumps, real Planck low-l TT likelihood (l = 2-29)", "",
       f"plain LCDM (Planck best fit): -2 lnL = {c0:.2f}", ""]
kcs = [1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 7e-4]
dls = {kc: dl_with_cut(kc) for kc in kcs}
out.append("(a) bounce alone (1 parameter):")
best_b = min(kcs, key=lambda k: chi2(dls[k]))
for kc in kcs: out.append(f"   k_c = {kc:.0e} /Mpc:  delta chi2 = {chi2(dls[kc]) - c0:+.2f}   D_2 = {dls[kc][2]:.0f} muK^2")
out.append("")
out.append("(b) layers alone (1 parameter; they only ADD power):")
def add_lumps(dl, N):
    d = dl.copy(); l = lump[:, 0].astype(int); d[l] += l*(l + 1)*lump[:, 1]/N/2/np.pi*T0**2; return d
Ns = [1e6, 3e6, 1e7, 3e7, 1e8, 1e9]
for N in Ns: out.append(f"   N_eff = {N:.0e}:  delta chi2 = {chi2(add_lumps(base, N)) - c0:+.2f}")
def Nbound(dl):
    cref = chi2(dl); g = np.geomspace(1e5, 1e10, 200); c = np.array([chi2(add_lumps(dl, N)) - cref for N in g])
    ok = g[c < 4.0]; return ok.min() if len(ok) else np.inf
nb0 = Nbound(base); out.append(f"   2-sigma lower bound (delta chi2 = 4 above LCDM): N_eff >= {nb0:.1e}")
out.append("")
out.append("(c) joint (2 parameters): best over the grid")
best = None
for kc in kcs:
    for N in Ns + [1e12]:
        c = chi2(add_lumps(dls[kc], N)) - c0
        if best is None or c < best[0]: best = (c, kc, N)
out.append(f"   best: k_c = {best[1]:.0e}, N_eff = {best[2]:.0e}: delta chi2 = {best[0]:+.2f}  (2 parameters; AIC change = {best[0] + 4:+.2f})")
nbb = Nbound(dls[best_b]); out.append(f"   with the best bounce suppression in place, the layer bound becomes N_eff >= {nbb:.1e}")
out += ["", "Reading: if the bounce alone improves chi2 by less than ~2 per parameter, the low quadrupole is not evidence for it;",
        "the layer bound with the real likelihood replaces the cosmic-variance estimate of iterations 8-9."]
txt = "\n".join(out); print(txt); open("iter12_bounce_vs_layers.txt", "w").write(txt + "\n")
