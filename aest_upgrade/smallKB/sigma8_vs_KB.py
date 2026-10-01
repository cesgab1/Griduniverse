# Does the early fast mode (rate ~ sqrt(3 Omega/K_B) aH) reach observables? sigma8 and P(k) vs K_B in CLASS (high precision)
import numpy as np
from classy import Class
base = {"gauge":"newtonian","N_ur":2.0328,"N_ncdm":1,"m_ncdm":0.06,"omega_cdm":0.0,"fluid_equation_of_state":"AEST","use_ppf":"no",
        "aest_K2":7.5e3,"aest_Q0":0.1,"aest_Z0":1e-9,"aest_de_beta":0.5,"aest_exchange":"yes",
        "H0":68.502588,"omega_b":0.022548605,"Omega_fld":0.11788114/0.68502588**2,"n_s":0.97434898,"A_s":1e-10*np.exp(3.0495784),
        "tau_reio":0.059290191,"output":"mPk,tCl","P_k_max_1/Mpc":2.0,
        "tol_perturbations_integration":1e-7,"perturbations_sampling_stepsize":0.02}
ks = np.array([0.001, 0.01, 0.05, 0.1, 0.5, 1.0]); ref = None
for kb in (0.5, 0.3, 0.2, 0.15, 0.1):
    c = Class(); c.set({**base, "aest_KB": kb})
    try:
        c.compute(); pk = np.array([c.pk(k, 0.) for k in ks]); s8 = c.sigma8()
        ref = pk if ref is None else ref
        print(f"K_B={kb:5}: sigma8={s8:.4f}  P(k)/P(k;K_B=0.5) at k={ks.tolist()}: " + " ".join(f"{r:.4f}" for r in pk/ref), flush=True)
    except Exception as e:
        print(f"K_B={kb}: failed {str(e)[:150]}", flush=True)
    c.struct_cleanup(); c.empty()
