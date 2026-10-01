# Sanity check: Khronon (GDM) branch vs plain LCDM with the same parameters; then Khronon + grid dark energy (beta=1/2, exchange)
import numpy as np
from classy import Class
cosmo = {"H0":68.0,"omega_b":0.0224,"n_s":0.968,"A_s":2.1e-9,"tau_reio":0.058,"N_ur":2.0328,"N_ncdm":1,"m_ncdm":0.06,
         "output":"tCl,pCl,lCl,mPk","lensing":"yes","l_max_scalars":2500,"P_k_max_1/Mpc":2.0}
kh = {"gauge":"newtonian","omega_cdm":0.0,"fluid_equation_of_state":"AEST","use_ppf":"no","aest_KB":0.5,"aest_K2":7.5e3,
      "aest_Q0":0.1,"aest_Z0":1e-9,"aest_khronon":"yes","Omega_fld":0.12/0.68**2}
def run(p):
    c = Class(); c.set(p); c.compute(); tt = c.lensed_cl(2500)["tt"][2:]; s8 = c.sigma8(); h = c.Hubble(0)
    c.struct_cleanup(); c.empty(); return tt, s8
ref, s8r = run({**cosmo, "omega_cdm":0.12, "gauge":"newtonian"})
for lab, extra in (("Khronon + Lambda", {"aest_de_beta":0.0}), ("Khronon + grid DE beta=1/2, exchange", {"aest_de_beta":0.5,"aest_exchange":"yes"}),
                   ("Khronon, K_B = 1e-5 (irrelevant: no vector)", {"aest_de_beta":0.0,"aest_KB":1e-5})):
    tt, s8 = run({**cosmo, **kh, **extra}); d = tt/ref-1
    print(f"{lab:46s}: TT vs LCDM max {np.max(abs(d))*100:.3f}%  (l=200 {d[198]*100:+.3f}%, l=1500 {d[1498]*100:+.3f}%)  sigma8 {s8:.4f} (LCDM {s8r:.4f})")
