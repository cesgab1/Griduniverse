"""
Test 2: expansion history out to the first billion years with quasars (Lusso et al. 2020 sample, 2,014 objects, distance moduli
from the X-ray/UV relation with fixed calibration, full covariance; Benetti et al. 2025 release). Absolute scale marginalised
(like supernovae), so only the SHAPE of distance vs redshift is tested.
Models: standard model (Planck-like Om = 0.315), the dark-energy law fit (Om 0.311), both with no free parameters in shape.
"""
import numpy as np
from scipy.integrate import quad, solve_ivp
D = "/home/claude/qso/MCMC/Cobaya/qso_data/"
d = np.genfromtxt(D + "lcparam_full_long_zhel.txt", usecols=(1, 2, 4, 5))
z, zh, mu_obs, dmu = d.T
C = np.loadtxt(D + "QSO_covmat.txt", skiprows=1).reshape(len(z), len(z)) + np.diag(dmu**2)
Ci = np.linalg.inv(C); one = np.ones(len(z))
print(f"{len(z)} quasars, z = {z.min():.3f} to {z.max():.2f}; beyond z = 5.7 (first billion years): {np.sum(z > 5.7)}; z > 3: {np.sum(z > 3)}")
def dl_lcdm(Om):
    E = lambda x: np.sqrt(Om*(1+x)**3 + (1-Om) + 9e-5*(1+x)**4)
    return np.array([(1+zhi)*quad(lambda x: 1/E(x), 0, zi)[0] for zi, zhi in zip(z, zh)])
def dl_law(Om, beta=0.5):
    orad = 9e-5
    zz = np.linspace(0, z.max()+0.1, 3000); lna = -np.log1p(zz)
    f = lambda x, y: [beta*(Om*np.exp(-3*x)/2 + orad*np.exp(-4*x) - np.exp(y[0]))/(Om*np.exp(-3*x) + orad*np.exp(-4*x) + np.exp(y[0]))]
    s = solve_ivp(f, [0, lna[-1]], [np.log(1-Om-orad)], t_eval=lna, rtol=1e-9); rde = np.exp(s.y[0])
    E = np.sqrt(Om*(1+zz)**3 + orad*(1+zz)**4 + rde)
    chi = np.concatenate([[0], np.cumsum(0.5*(1/E[1:] + 1/E[:-1])*np.diff(zz))])
    return (1+zh)*np.interp(z, zz, chi)
def score(dl, lab):
    r = mu_obs - 5*np.log10(dl); off = (one @ Ci @ r)/(one @ Ci @ one); r = r - off
    chi = r @ Ci @ r
    print(f"\n{lab}: chi2 = {chi:.1f} for {len(z)} quasars")
    edges = [0, 0.5, 1, 1.5, 2, 3, 4, 5.7, 8]
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (z >= lo) & (z < hi)
        if m.sum():
            w = 1/dmu[m]**2; mean = np.sum(w*r[m])/w.sum(); err = 1/np.sqrt(w.sum())
            print(f"   z {lo:3.1f}-{hi:3.1f} ({m.sum():4d}): distance offset {mean:+.3f} +- {err:.3f} mag  ({mean/err:+.1f} sigma)")
    return chi
c1 = score(dl_lcdm(0.315), "standard model (Om 0.315)")
c2 = score(dl_law(0.311), "dark-energy law (beta 1/2, Om 0.311)")
# what the quasars alone prefer
best = min(((score_q := None) or 0, om) for om in [0])  # placeholder
chis = []
for om in np.arange(0.2, 1.01, 0.1):
    dl = dl_lcdm(om); r = mu_obs - 5*np.log10(dl); r -= (one @ Ci @ r)/(one @ Ci @ one); chis.append((r @ Ci @ r, om))
print("\nquasars alone, flat standard model: chi2 vs Om:", ", ".join(f"{om:.1f}:{c:.1f}" for c, om in chis))
