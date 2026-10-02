"""
Is there a grid-level channel that removes a halo's random motions in ~0.3 Gyr, only in the MOND regime (section 9's requirement)?
1. Gate: the fluid's motions couple to the Xi field with strength ~ [k^2/(k^2 + mu_eff^2)]^2. In the DBI dust era mu_eff is huge
   (crossover << halo) and the coupling vanishes; in the MOND era it is ~1. So a channel tied to Xi would be on only in the MOND
   regime, which is what is needed. This part comes for free from the DBI law.
2. Channel: radiation of Xi (khronon) waves by the halo's time-varying mass distribution, quadrupole-like:
      P ~ G M^2 sigma^6 / (c_k^5 R^2),  needed P_req = M sigma^2 / (2 tau).
   Solve for the wave speed c_k that gives tau = 0.3 Gyr, then evaluate P/P_req for c_k = c.
3. Constraint: vacuum gravitational-Cherenkov limits (Elliott, Moore & Stoica 2005, hep-ph/0505211): ultra-high-energy cosmic
   rays would radiate into any gravitationally coupled aether/khronon mode slower than themselves, so all such modes must
   travel at ~c. Slow (galactic-speed) modes are excluded.
"""
import numpy as np
G = 4.30e-6*1.0227**2; kms = 1.0227                                         # kpc^3/(Msun Gyr^2); km/s -> kpc/Gyr
M, sig, R, tau = 2e12, 200*kms, 200.0, 0.3
ck = (2*G*M*sig**4*tau/R**2)**0.2
print(f"wave speed needed for tau = {tau} Gyr: c_k ~ {ck/kms:.0f} km/s")
for c in (ck, 3e3*kms, 3e5*kms):
    P_ratio = (G*M**2*sig**6/(c**5*R**2))/(M*sig**2/(2*tau))
    print(f"   c_k = {c/kms:9.0f} km/s: radiated / needed = {P_ratio:.1e}  (dissipation time {tau/P_ratio:.1e} Gyr)")
print("Cherenkov limit: gravitationally coupled khronon/aether modes must travel at ~c, so the channel gives a dissipation time ~3e15 Gyr.")
