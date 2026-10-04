"""
ITERATION 51 (pre-registered in PREREG_51.md): spectra of the frame field's frozen ripples for two ways of tilting them.
Units M = 1. Era: H = a^(-epsilon) (epsilon = 2 radiation). Mode k freezes at a_* where omega(k/a_*) = H(a_*); P = p_*^3/omega_*.
"""
import numpy as np
from scipy.optimize import brentq
def spectrum(omega, lnk, eps):
    P = []
    for lk in lnk:
        g = lambda la: omega(lk - la) + eps*la                         # ln omega(p) - ln H, H = a^-eps (omega takes ln p)
        grid = np.linspace(-1000, 1000, 40001); gv = g(grid); i = np.argmax(gv < 0)     # FIRST freeze (omega drops below H)
        if gv[i] >= 0: P.append(np.nan); continue                                         # never freezes
        la = brentq(g, grid[i-1], grid[i]); lp = lk - la; P.append(3*lp - omega(lp))
    P = np.array(P); ns = 1 + np.gradient(P, lnk); al = np.gradient(ns, lnk); ok = np.isfinite(al); return ns[ok], al[ok]
out = ["ITERATION 51: sources of the watermark tilt (expectations pre-registered)", ""]
lnk = np.linspace(-30, 30, 6001)
# constant exponent z = 3 - eta
for eps in (2.0, 1.5, 0.5):
    for eta in (0.015,):
        z = 3 - eta; ns, al = spectrum(lambda lp: z*lp, lnk, eps)
        pred = -(eta/3)*eps/(1 - eps/3)
        out.append(f"constant z = 3 - {eta}: era epsilon = {eps}: n_s - 1 = {np.median(ns) - 1:+.5f} (formula {pred:+.5f}); "
                   f"running = {np.median(al):+.1e}")
out.append("")
# crossover omega^2 = p^2 + p^6
alpha_obs, alpha_err = 0.0062, 0.0052
for eps in (2.0, 1.5, 0.5):
    ns, al = spectrum(lambda lp: 0.5*np.logaddexp(2*lp, 6*lp), lnk, eps)
    for target in (0.974, 0.965):
        i = np.argmin(np.abs(ns - target)); pred = (4/3)*eps/(1 - eps/3)*(ns[i] - 1)
        out.append(f"crossover, epsilon = {eps}: where n_s = {ns[i]:.4f}, running alpha_s = {al[i]:+.4f} (formula {pred:+.4f}); "
                   f"vs ACT {alpha_obs} +/- {alpha_err}: {abs(al[i] - alpha_obs)/alpha_err:.0f} sigma")
out.append("")
# spectral dimension needed
for ns_obs, lab in ((0.974, "ACT"), (0.965, "Planck")):
    eta = (1 - ns_obs)/2; out.append(f"{lab}: radiation era needs eta = {eta:.4f}, z = {3-eta:.4f}, frame-field d_s = 1 + 3/z = {1 + 3/(3-eta):.4f}")
txt = "\n".join(out); print(txt); open("iter51_eta_sources.txt", "w").write(txt + "\n")
