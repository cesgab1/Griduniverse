# DBI Khronon with large lambda_D (fluid stays in its pressureless 'early' branch until today): LCDM recovered?
# Also: effective static mass term today, mu_eff^2 = K''(Q_today)/2 = mu^2/(1 - lam x0^2)^(3/2)  (what galaxies would feel)
import numpy as np
exec(open("dbi_cosmology.py").read().split("ks = ")[0])
ks = np.array([0.01, 0.05, 0.1, 0.5, 1, 5])
c = Class(); c.set({**cosmo, "omega_cdm":0.12}); c.compute(); tt0 = c.lensed_cl(2500)["tt"][2:]; p0 = np.array([c.pk(k,0) for k in ks]); s80 = c.sigma8(); c.struct_cleanup(); c.empty()
I0 = 3*0.12/0.68**2*(68/2.998e5)**2
for lam in (1e8, 1e12, 1e16):
    mu = 1/22.3
    c = Class(); c.set({**cosmo, **kh, "aest_dbi_mu": mu, "aest_dbi_lam": lam})
    try:
        c.compute(); tt = c.lensed_cl(2500)["tt"][2:]; p = np.array([c.pk(k,0) for k in ks])
        s0 = I0/(2*mu**2); x0 = s0/np.sqrt(1+lam*s0**2); r = np.sqrt(1-lam*x0**2)
        print(f"lam={lam:.0e}: TT vs LCDM max {np.max(abs(tt/tt0-1))*100:.2f}% | sigma8 {c.sigma8():.3f} (LCDM {s80:.3f}) | P/P_LCDM "
              + " ".join(f"{v:.3f}" for v in p/p0) + f" | static mass term today: 1/mu_eff = {22.3*r**1.5*1e3:.3g} kpc", flush=True)
    except Exception as e: print(f"lam={lam:.0e}: failed {str(e)[:150]}")
    c.struct_cleanup(); c.empty()
