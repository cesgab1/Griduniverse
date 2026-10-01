"""
Could ended black holes (stellar clusters of them, or supermassive ones) explain voids / the Hubble tension?
1. Have any ended?  bounce (M^2) / Hawking (M^3) lifetimes vs age (planck_star_times.txt).
2. Energy budget: even if every black hole had ended, how much energy is that vs what the Hubble tension needs?
3. Light speed: does a cold, empty region slow light, and would that change measured H0?
"""
import numpy as np
h = 0.675; rho_crit = 2.775e11*h**2          # Msun / Mpc^3
rho_smbh = 4.4e5                              # Msun/Mpc^3 (Shankar et al. 2009 range 3-5.5e5)
Om_smbh = rho_smbh/rho_crit
Om_allbh = 1e-4                               # order of magnitude incl. stellar-mass holes (literature ~1e-4 to 3e-4)
need = (73.0/67.5)**2 - 1                     # fractional change in H^2 the distance ladder implies
print(f"SMBH mass density: Omega = {Om_smbh:.1e};  all black holes: Omega ~ {Om_allbh:.0e}")
print(f"Hubble tension needs a ~{need*100:.0f}% change in H^2 (or in the distance calibration)")
print(f"   all SMBHs converted to energy: {Om_smbh/need:.1e} of what is needed;  all black holes: {Om_allbh/need:.1e}")
# light in a cold void: vacuum speed is exactly c; ionised gas slows radio waves (dispersion) ~ n_e
ne = 2e-7*1e6   # electrons per m^3 in a void (~0.2 per m^3... use 2e-7 cm^-3 x 1e6)
for nu_hz, lab in ((1.4e9, "radio 1.4 GHz"), (5e14, "visible light")):
    nu_p = 8.98*np.sqrt(ne)                   # plasma frequency, Hz
    print(f"   {lab}: speed deficit in a void (1 - v/c) = {0.5*(nu_p/nu_hz)**2:.1e}")
print("Distances in the Hubble diagram come from brightness (standard candles) and redshift, not from light travel time:")
print("   slowing light delays arrival but changes neither brightness nor redshift -> no effect on H0.")
