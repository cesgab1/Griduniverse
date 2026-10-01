"""
Viscoelastic (Maxwell) dark fluid: elastic (stiff) for strain faster than tau, liquid for slower.
Halo collapse time t_ff = sqrt(3 pi / (32 G rho_vir)), rho_vir = 200 rho_crit(z_f) -> t_ff = 0.111 / H(z_f): depends only on formation epoch.
Formation epoch z_f (half of the mass assembled) from extended Press-Schechter (Lacey & Cole 1993), LCDM sigma(M) from CLASS.
Fluid kept OUT of a halo if t_ff(z_f) < tau (stiff); flows IN if t_ff(z_f) > tau (soft).
Requirements: galaxies (1e11-1e12.5 Msun) almost never capture (SPARC RAR scatter is ~0.1 dex, so <~5%), clusters (1e14-1e15) almost always.
Also: linear growth (timescale ~1/H) must be liquid at all z that matter (z < ~10), and the early-universe elastic phase (H tau > 1) needs
a vanishing modulus.
"""
import numpy as np
from classy import Class
from scipy.integrate import quad
h = 0.68; Om = 0.31; H0 = 100*h; H0inv_Gyr = 977.8/H0
c = Class(); c.set({"H0":68,"omega_b":0.0224,"omega_cdm":0.12,"n_s":0.968,"A_s":2.1e-9,"output":"mPk","P_k_max_1/Mpc":50,"N_ur":3.046})
c.compute()
rho_m = Om*2.775e11*h**2                       # Msun/Mpc^3
def sigma(M):
    R = (3*M/(4*np.pi*rho_m))**(1/3); return c.sigma(R, 0.0)
Hz = lambda z: np.sqrt(Om*(1+z)**3 + 1 - Om)
def D(z): return c.scale_independent_growth_factor(z)
dc = lambda z: 1.686/D(z)
def P_form_after(M, z):          # P(z_f < z) = 1 - P(half-mass progenitor already in place at z)  [Lacey & Cole 1993]
    S0, Sh = sigma(M)**2, sigma(M/2)**2; w = dc(z) - dc(0)
    if w <= 0: return 0.0
    # substitution y = w/sqrt(S-S0): integrand (M/M1) sqrt(2/pi) exp(-y^2/2), y from w/sqrt(Sh-S0) to infinity
    ylo = w/np.sqrt(Sh - S0)
    f = lambda y: (M/M1(S0 + (w/y)**2))*np.sqrt(2/np.pi)*np.exp(-y*y/2)
    return 1 - quad(f, ylo, ylo + 12, limit=200)[0]
Ms = np.logspace(9, 16, 300); Ss = np.array([sigma(m)**2 for m in Ms])
M1 = lambda S: np.exp(np.interp(S, Ss[::-1], np.log(Ms[::-1])))
zgrid = np.linspace(0, 8, 161)
def zf_samples(M):
    cdf = np.array([P_form_after(M, z) for z in zgrid]); cdf = np.maximum.accumulate(np.clip(cdf, 0, 1))
    return cdf
tff = lambda z: 0.111/Hz(z)*H0inv_Gyr      # Gyr
masses = [(1e11, "dwarf-ish galaxy"), (1e12, "Milky-Way-mass galaxy"), (3e12, "massive galaxy"), (1e13, "group"), (1e14, "cluster"), (1e15, "massive cluster")]
cdfs = {M: zf_samples(M) for M, _ in masses}
print("median formation redshift and collapse time:")
for M, lab in masses:
    zmed = np.interp(0.5, cdfs[M], zgrid); print(f"   {lab:22s} M={M:.0e}: z_f median {zmed:.2f}, t_ff {tff(zmed):.2f} Gyr")
print("\nfraction CAPTURING the fluid (t_ff > tau, i.e. formed later than z where t_ff = tau):")
print(f"{'tau (Gyr)':>9s} " + " ".join(f"{lab[:14]:>14s}" for _, lab in masses) + "   z where 1/H = tau (elastic before)")
for tau in (0.4, 0.6, 0.8, 1.0, 1.2, 1.5):
    zs = np.linspace(0, 8, 2001); zcut = zs[np.argmin(np.abs(tff(zs) - tau))] if tff(0) > tau else 0.0
    fr = [np.interp(zcut, zgrid, cdfs[M]) for M, _ in masses]          # P(z_f < zcut)
    zel = zs[np.argmin(np.abs(H0inv_Gyr/Hz(zs) - tau))]
    print(f"{tau:9.1f} " + " ".join(f"{f:14.2f}" for f in fr) + f"   z > {zel:.1f}")
