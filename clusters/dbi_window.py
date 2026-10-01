"""
Khronon no-go check. For any K(Q) with w << 1:  c_ad^2 = K'/(Q K''),  static mass term mu_eff^2 = K''/2  ->  c_ad^2 mu_eff^2 = K'/(2Q) ~ 4 pi G rho.
So the fluid's cosmological crossover (Jeans-like) length today = 2 pi / mu_eff: the same number that sets the extra pull around galaxies.
Scan DBI lambda_D (mu^-1 = 22.3 Mpc fixed): sigma8, P(k) at z=3 (Lyman-alpha), and isolated-galaxy lensing boost from the mass term.
"""
import numpy as np, warnings; warnings.filterwarnings("ignore")
exec(open("dbi_cosmology.py").read().split("ks = ")[0])
src = open("khronon_mass_term.py").read()
exec(src.split('print("Cluster')[0])
ks = np.array([1, 2, 5])
c = Class(); c.set({**cosmo, "omega_cdm":0.12}); c.compute(); p3 = np.array([c.pk(k,3) for k in ks]); s80 = c.sigma8(); c.struct_cleanup(); c.empty()
I0 = 3*0.12/0.68**2*(68/2.998e5)**2; mu = 1/22.3
print(f"{'lambda_D':>9s} {'1/mu_eff':>9s} {'crossover':>9s} {'sigma8':>7s}  P/P_LCDM z=3 (k=1,2,5)   galaxy lensing boost at 0.3 / 1 Mpc (R=3 Mpc)")
for lam in (1e8, 1e9, 1e10, 1e11, 1e12):
    s0 = I0/(2*mu**2); r = 1/np.sqrt(1+lam*s0**2); mu_eff_inv = (1/mu)*r**1.5
    c = Class(); c.set({**cosmo, **kh, "aest_dbi_mu": mu, "aest_dbi_lam": lam}); c.compute()
    s8 = c.sigma8(); q = np.array([c.pk(k,3) for k in ks])/p3; c.struct_cleanup(); c.empty()
    rr = np.geomspace(1*kpc, 3*Mpc, 4000)
    g0, _, _ = solve(galaxy_Mb, None, 3*Mpc, rr); g1, _, _ = solve(galaxy_Mb, mu_eff_inv, 3*Mpc, rr)
    b = [np.interp(x*Mpc, rr, g1/g0) for x in (0.3, 1.0)]
    print(f"{lam:9.0e} {mu_eff_inv:7.2f}Mpc {2*np.pi*mu_eff_inv:7.1f}Mpc {s8:7.3f}  " + " ".join(f"{v:.3f}" for v in q) + "            " + " / ".join(fmt(v) for v in b), flush=True)
print(f"LCDM sigma8 {s80:.3f}; measured S8 ~ 0.78-0.83. Lyman-alpha tolerates ~2% loss at k=5. KiDS galaxy lensing tolerates ~10% boost.")
