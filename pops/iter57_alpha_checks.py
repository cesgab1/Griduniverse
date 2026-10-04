"""ITERATION 57 (pre-registered in PREREG_57.md): delta = alpha against known facts."""
import numpy as np
alpha = 1/137.035999; d = alpha
out = ["ITERATION 57: delta = alpha (0.730%) against things we already know (expectations pre-registered)", ""]
# K1 white dwarfs
for name, v, s in (("GD133", -2.7e-5, np.hypot(4.7e-5, 0.2e-5)), ("G29-38", -5.8e-5, np.hypot(3.8e-5, 0.3e-5))):
    out.append(f"K1 {name}: local version predicts Delta mu/mu = {-d:.1e}; measured {v:.1e} +/- {s:.1e} -> off by {abs(-d - v)/s:.0f} sigma")
out.append("   -> LOCAL version EXCLUDED; mechanism must be global (cosmic photon bath / cosmic frame) and one-way.")
# K2 CMB fits (predictions/README.md item 3)
me = np.array([1.000, 1.004, 1.008, 1.012, 1.016]); dchi = np.array([0, -2.0, -3.8, -3.8, -2.2]); H0 = np.array([68.2, 68.6, 69.3, 69.8, 70.4])
c = np.polyfit(me, dchi, 2); x = 1 + d
out += ["", f"K2 CMB+BAO+SN at m_e = {x:.4f}: Delta chi2 vs LCDM = {np.polyval(c, x):.1f} (best possible {np.polyval(c, -c[1]/(2*c[0])):.1f}); "
        f"H0 = {np.interp(x, me, H0):.1f} km/s/Mpc (SH0ES 73: tension eased, not solved)"]
# K3 BBN
for name, m, s in (("PRIMAT", 0.510, 0.007), ("NACRE II", 0.504, 0.007)):
    out.append(f"K3 BBN ({name}): predicted m_e = {0.51099895*(1+d):.4f} MeV vs {m} +/- {s} -> {(0.51099895*(1+d) - m)/s:+.1f} sigma")
# K4 mechanisms
out += ["", f"K4a field energy outside the pop size: alpha/2 = {100*d/2:.3f}% -> outside the 0.40-0.90% window: NOT the mechanism as stated",
        "K4b frame lock: right structure (global, one-way, ends when the cosmic photon fluid lets go), no numerical factor: OPEN"]
txt = "\n".join(out); print(txt); open("iter57_alpha_checks.txt", "w").write(txt + "\n")
