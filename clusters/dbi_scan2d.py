"""2D scan of DBI Khronon (mu, lambda_D): need sigma8 in [0.78,0.83], Lyman-alpha P(k=5,z=3)/LCDM >= 0.98,
isolated-galaxy lensing boost at 1 Mpc <= 1.10, and critical depth v_crit in ~[700, 1100] km/s (galaxies stay fluid, groups/clusters collapse).
Parameterised by mu and eps = lam s0^2 (s0 = I0/2mu^2): x0 = s0/sqrt(1+eps), x_lim = s0/sqrt(eps), 1/mu_eff = (1/mu)(1+eps)^(-3/4)."""
import numpy as np, warnings, json; warnings.filterwarnings("ignore")
exec(open("dbi_cosmology.py").read().split("ks = ")[0])
exec(open("khronon_mass_term.py").read().split('print("Cluster')[0])
c_ = 2.998e5; I0 = 3*0.12/0.68**2*(68/2.998e5)**2
cl = Class(); cl.set({**cosmo, "omega_cdm":0.12}); cl.compute(); p5 = cl.pk(5, 3); cl.struct_cleanup(); cl.empty()
rr = np.geomspace(1*kpc, 3*Mpc, 3000); g0, _, _ = solve(galaxy_Mb, None, 3*Mpc, rr)
import os
MUS = json.loads(os.environ.get("MUS", "[5.0, 22.3, 100.0]")); EPS = json.loads(os.environ.get("EPS", "[0.1, 0.3, 1.0, 3.0]"))
rows = []
for mu_inv in MUS:
    mu = 1/mu_inv; s0 = I0/(2*mu**2)
    for eps in EPS:
        lam = eps/s0**2; vc = c_*np.sqrt(s0*(1/np.sqrt(eps) - 1/np.sqrt(1+eps))); mue = mu_inv*(1+eps)**-0.75
        cl = Class(); cl.set({**cosmo, **kh, "aest_dbi_mu": mu, "aest_dbi_lam": lam})
        try:
            cl.compute(); s8 = cl.sigma8(); ly = cl.pk(5, 3)/p5
        except Exception as e:
            s8 = ly = np.nan
        cl.struct_cleanup(); cl.empty()
        g1, _, _ = solve(galaxy_Mb, mue, 3*Mpc, rr); lb = np.interp(Mpc, rr, g1/g0)
        ok = (0.78 <= s8 <= 0.83) and ly >= 0.98 and lb <= 1.10 and 700 <= vc <= 1100
        rows.append(dict(mu_inv=mu_inv, eps=eps, lam=lam, vcrit=vc, mu_eff_inv=mue, sigma8=s8, lya=ly, lens1Mpc=lb, all_ok=bool(ok)))
        print(f"1/mu={mu_inv:6.1f} eps={eps:4.1f} lam={lam:8.1e}: v_crit {vc:6.0f} km/s | 1/mu_eff {mue:6.2f} Mpc | sigma8 {s8:.3f} | Lya {ly:.3f} | lens@1Mpc {fmt(lb)} | {'ALL PASS' if ok else ''}", flush=True)
json.dump(rows, open(os.environ.get("OUT", "dbi_scan2d.json"), "w"), indent=1)
