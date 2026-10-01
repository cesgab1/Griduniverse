import numpy as np, sys
from classy import Class
base = {"gauge":"newtonian","N_ur":2.0328,"N_ncdm":1,"m_ncdm":0.06,"omega_cdm":0.0,"fluid_equation_of_state":"AEST","use_ppf":"no",
        "aest_K2":7.5e3,"aest_Q0":0.1,"aest_Z0":1e-9,"aest_de_beta":0.5,"aest_exchange":"yes",
        "H0":68.502588,"omega_b":0.022548605,"Omega_fld":0.11788114/0.68502588**2,"n_s":0.97434898,"A_s":1e-10*np.exp(3.0495784),
        "tau_reio":0.059290191,"output":"tCl,pCl,lCl","lensing":"yes","l_max_scalars":2500}
ref = None
for kb in [0.5, 0.2, 0.1, 0.05, 0.02]:
    c = Class(); c.set({**base, "aest_KB": kb})
    try:
        c.compute(); cl = c.lensed_cl(2500); tt = cl["tt"][2:2501]
        if ref is None: ref = tt
        d = tt/ref - 1
        print(f"K_B={kb}: ok, TT change vs K_B=0.5: max {np.nanmax(abs(d))*100:.2f}%  at l=200 {d[198]*100:+.2f}%, l=1000 {d[998]*100:+.2f}%  nan={np.isnan(tt).any()}", flush=True)
    except Exception as e:
        print(f"K_B={kb}: FAILED {str(e)[:300]}", flush=True)
    c.struct_cleanup(); c.empty()
