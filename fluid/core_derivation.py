"""
Derive the fluid core instead of picking it.
A fluid with P = kappa rho^2 (Gamma = 2; Khronon's simple quadratic K, dense channel gas) in hydrostatic equilibrium is an n = 1
polytrope: radius R = pi sqrt(kappa / (2 pi G)), independent of mass. For Khronon's K: kappa = 2 pi G / mu^2  ->  R = pi / mu.
SPARC needs core >~ 150 kpc  ->  1/mu ~ 48 kpc.  DBI version: dust (pressureless) above a transition density, so dense cluster
cores collapse like cold dark matter; early universe is dust.
The same mu sets the fluid's cosmological crossover length 2 pi / mu_eff. Check sigma8 and Lyman-alpha (P(k=5/Mpc, z=3)) in CLASS
for 1/mu = 48 kpc and the dust transition at z_t (eps = lam s0^2 = 1/(1+z_t)^6).
"""
import numpy as np, warnings; warnings.filterwarnings("ignore")
exec(open("/home/claude/griduniverse/clusters/dbi_cosmology.py").read().split("ks = ")[0])
c_ = Class(); c_.set({**cosmo, "omega_cdm":0.12}); c_.compute(); s80 = c_.sigma8(); p3 = np.array([c_.pk(k, 3) for k in (1, 2, 5)]); c_.struct_cleanup(); c_.empty()
I0 = 3*0.12/0.68**2*(68/2.998e5)**2
print(f"n=1 polytrope core radius pi/mu: 1/mu = 48 kpc -> core {np.pi*48:.0f} kpc (any mass)")
for mu_inv_kpc in (48.0, 20.0):
    mu = 1/(mu_inv_kpc*1e-3)
    for zt in (14.0, 30.0, 100.0):
        eps = (1+zt)**-6; s0 = I0/(2*mu**2); lam = eps/s0**2
        rho_tr = 1/np.sqrt(eps)        # transition density / today's mean fluid density
        c_ = Class(); c_.set({**cosmo, **kh, "aest_dbi_mu": mu, "aest_dbi_lam": lam})
        try:
            c_.compute(); s8 = c_.sigma8(); q = np.array([c_.pk(k, 3) for k in (1, 2, 5)])/p3
            print(f"1/mu {mu_inv_kpc:4.0f} kpc (core {np.pi*mu_inv_kpc:4.0f} kpc), dust before z={zt:5.0f} (transition at {rho_tr:.1e} x mean density): "
                  f"sigma8 {s8:.3f} (LCDM {s80:.3f}) | Lya P/P_LCDM z=3 k=1,2,5: " + " ".join(f"{v:.3f}" for v in q), flush=True)
        except Exception as e: print("failed", str(e)[:100])
        c_.struct_cleanup(); c_.empty()
