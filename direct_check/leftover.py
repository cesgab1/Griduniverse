"""
Coalesce's direct check: speed (expansion rate) + mass (matter) at each moment -> leftover = dark energy at that moment.
Friedmann: (H/100)^2 = omega_m (1+z)^3 + omega_r (1+z)^4 + omega_DE(z)   (omega = density x h^2, in units of the critical density
for H = 100 km/s/Mpc). Inputs: DESI DR2 D_H/r_d (expansion rate at each redshift, real data + covariance), sound horizon r_d and
matter omega_m from Planck 2018 (standard early universe). If dark energy is constant, omega_DE(z) is the same at every z.
"""
import numpy as np, pandas as pd
c = 299792.458
B = "/home/claude/cobayasampler/bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_"
d = pd.read_csv(B + "mean.txt", sep=r"\s+", comment="#", header=None, names=["z", "v", "q"]); cov = np.loadtxt(B + "cov.txt")
sel = (d.q == "DH_over_rs").values; z = d.z.values[sel]; v = d.v.values[sel]; C = cov[np.ix_(sel, sel)]
rd, srd = 147.09, 0.26; wm, swm = 0.1430, 0.0011; wr = 2.469e-5*(1 + 0.2271*3.046)
rng = np.random.default_rng(1); N = 200000
V = rng.multivariate_normal(v, C, N); RD = rng.normal(rd, srd, N); WM = rng.normal(wm, swm, N)
H = c/(V*RD[:, None])                                   # km/s/Mpc
X = (H/100)**2 - WM[:, None]*(1 + z)**3 - wr*(1 + z)**4  # leftover = dark energy density x h^2
m, s = X.mean(0), X.std(0)
# Claim 1 prediction for the same quantity (Pantheon+ best fit Om 0.31, h 0.6749)
from scipy.optimize import brentq
Om, h = 0.31, 0.6749; Or = wr/h**2; OL = 1 - Om - Or
def E(zz):
    a = 1/(1 + zz); f = lambda e: e*e - Om*a**-3 - Or*a**-4 - OL*(a*e)**-0.5
    return brentq(f, 1e-6, 1e4)
pred = [OL*h*h*(E(zz)/(1+zz))**-0.5 for zz in z]
out = ["Direct check: dark energy left over at each moment = (expansion rate)^2 - matter - radiation  [omega units]", "",
       "  z      leftover (data)        Claim 1 prediction     constant (Lambda, Planck)"]
for i in range(len(z)):
    out.append(f" {z[i]:4.2f}   {m[i]:+.3f} +/- {s[i]:.3f}       {pred[i]:.3f}                 {0.6736**2 - 0.1430:.3f}")
# test of constancy: fit one constant to all points with full MC covariance
Cx = np.cov(X.T); w = np.linalg.inv(Cx); one = np.ones(len(z))
k = (one @ w @ m)/(one @ w @ one); chi_const = (m - k) @ w @ (m - k)
chi_c1 = (m - np.array(pred)) @ w @ (m - np.array(pred))
out += ["", f"best single constant: {k:.3f} +/- {1/np.sqrt(one @ w @ one):.3f};  chi2 for 'constant' = {chi_const:.2f} for {len(z)-1} dof",
        f"chi2 for the Claim 1 prediction (no fitting) = {chi_c1:.2f} for {len(z)} points",
        "Caveat: assumes the standard early universe for r_d and omega_m (Claim 2's heavier electron would shift r_d by ~1%)."]
txt = "\n".join(out); print(txt); open("leftover.txt", "w").write(txt + "\n")
