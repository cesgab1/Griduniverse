"""
Where does Khronon's fluid clump?  GDM sound speed (Blanchet & Skordis eq 4.31): c_s^2 = c_ad^2 / (1 + c_ad^2 k^2 / (4 pi G a^2 rho (1+w))).
For c_ad^2 k^2 >> 4 pi G a^2 rho: c_s^2 k^2 -> 4 pi G a^2 rho (1+w): pressure cancels self-gravity -> the fluid stops self-clumping
below the crossover length lambda_* = 2 pi c_ad / sqrt(4 pi G rho) and only follows the baryons.
c_ad^2 = K'/(Q K'') = Z0/Q for the Cosh K (paper Z0 = 1e-9, Q0 = 0.1).
Measure: matter power spectrum vs LCDM at z = 0 and z = 3 (Lyman-alpha forest scales k ~ 0.5-5 /Mpc) for several Z0.
"""
import numpy as np
from classy import Class
cosmo = {"H0":68.0,"omega_b":0.0224,"n_s":0.968,"A_s":2.1e-9,"tau_reio":0.058,"N_ur":2.0328,"N_ncdm":1,"m_ncdm":0.06,
         "output":"mPk","P_k_max_1/Mpc":12.0,"z_max_pk":3.5,"gauge":"newtonian"}
kh = {"omega_cdm":0.0,"fluid_equation_of_state":"AEST","use_ppf":"no","aest_KB":0.5,"aest_K2":7.5e3,"aest_Q0":0.1,
      "aest_khronon":"yes","Omega_fld":0.12/0.68**2,"aest_de_beta":0.0}
ks = np.array([0.1, 0.5, 1, 2, 5, 10])
c = Class(); c.set({**cosmo, "omega_cdm":0.12}); c.compute()
ref = {z: np.array([c.pk(k, z) for k in ks]) for z in (0, 3)}; c.struct_cleanup(); c.empty()
H0 = 68/2.998e5; OmK = 0.12/0.68**2
for Z0 in (1e-12, 1e-10, 1e-9, 1e-8):
    cad = np.sqrt(Z0/0.1); lam = 2*np.pi*cad/np.sqrt(1.5*OmK*H0**2)
    c = Class(); c.set({**cosmo, **kh, "aest_Z0": Z0})
    try:
        c.compute()
        r0 = np.array([c.pk(k, 0) for k in ks])/ref[0]; r3 = np.array([c.pk(k, 3) for k in ks])/ref[3]
        print(f"Z0={Z0:.0e}: c_ad^2={Z0/0.1:.0e}, crossover length today ~{lam:.2g} Mpc | P/P_LCDM at k={ks.tolist()}: z=0 " +
              " ".join(f"{v:.3f}" for v in r0) + " | z=3 " + " ".join(f"{v:.3f}" for v in r3), flush=True)
    except Exception as e:
        print(f"Z0={Z0:.0e}: failed {str(e)[:120]}")
    c.struct_cleanup(); c.empty()
