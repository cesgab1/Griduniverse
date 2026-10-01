"""Khronon with the paper's DBI K(Q) (mu^-1 = 22.3 Mpc, lambda_D = 1): does its cosmology look like LCDM, as claimed?"""
import numpy as np
from classy import Class
cosmo = {"H0":68.0,"omega_b":0.0224,"n_s":0.968,"A_s":2.1e-9,"tau_reio":0.058,"N_ur":2.0328,"N_ncdm":1,"m_ncdm":0.06,
         "output":"tCl,pCl,lCl,mPk","lensing":"yes","l_max_scalars":2500,"P_k_max_1/Mpc":12.0,"z_max_pk":3.5,"gauge":"newtonian"}
kh = {"omega_cdm":0.0,"fluid_equation_of_state":"AEST","use_ppf":"no","aest_KB":0.5,"aest_K2":7.5e3,"aest_Q0":0.1,"aest_Z0":1e-9,
      "aest_khronon":"yes","aest_dbi":"yes","Omega_fld":0.12/0.68**2,"aest_de_beta":0.0}
ks = np.array([0.01, 0.05, 0.1, 0.5, 1, 5])
c = Class(); c.set({**cosmo, "omega_cdm":0.12}); c.compute()
tt0 = c.lensed_cl(2500)["tt"][2:]; p0 = np.array([c.pk(k, 0) for k in ks]); s80 = c.sigma8(); c.struct_cleanup(); c.empty()
for mu_inv in (2.0, 0.2, 0.02, 0.002):
    c = Class(); c.set({**cosmo, **kh, "aest_dbi_mu": 1/mu_inv, "aest_dbi_lam": 1.0})
    try:
        c.compute(); bg = c.get_background(); a = 1/(1+bg["z"]); w = bg["(.)w_fld"] if "(.)w_fld" in bg else (print([k for k in bg if "fld" in k]) or 0*bg["z"])
        tt = c.lensed_cl(2500)["tt"][2:]; p = np.array([c.pk(k, 0) for k in ks])
        wtxt = " ".join(f"a={aa:g}:{np.interp(np.log(aa), np.log(a[::-1]), w[::-1]):.1e}" for aa in (1e-4, 1e-3, 0.01, 0.1, 1))
        print(f"1/mu={mu_inv:6} Mpc: w(a) {wtxt} | TT vs LCDM max {np.max(abs(tt/tt0-1))*100:.2f}% | sigma8 {c.sigma8():.3f} (LCDM {s80:.3f}) | P/P_LCDM "
              + " ".join(f"{v:.3f}" for v in p/p0), flush=True)
    except Exception as e:
        print(f"1/mu={mu_inv}: failed {str(e)[:200]}", flush=True)
    c.struct_cleanup(); c.empty()
