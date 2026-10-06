"""A violent, vibrating start (Coalesce): what would survive, and could it make dark matter (Planck-mass leftovers)?
Calculation, not a data fit. Limits marked (mem) are from memory."""
import numpy as np
from scipy.special import erfcinv
G, c, hb = 6.674e-11, 2.998e8, 1.0546e-34
mP = np.sqrt(hb*c/G); tP = np.sqrt(hb*G/c**5)
# 1 gravitational waves from the shaking: allowed energy vs light, from N_eff = 2.99 +/- 0.17 -> dN_eff < ~0.34 (2 sigma)
frac = 7/8*(4/11)**(4/3)*0.34
print(f"1. shaking left as gravitational waves: at most {frac:.0%} of the energy in light (from N_eff); the rest must have become heat")
# 2 Planck-mass leftovers (black holes that stop at the grid's density cap) as ALL dark matter
rho_dm = 0.26*8.5e-27; teq = 1.6e12
n = rho_dm/mP
beta = np.sqrt(tP/teq)                         # fraction of the universe collapsing at Planck time (grows ~a until equality)
dc = 0.45; sigma = dc/(np.sqrt(2)*erfcinv(beta))
print(f"2. Planck-mass leftover = {mP*1e9:.0f} micrograms; to be all dark matter today: {n:.1e} per m^3 (one per cube {n**(-1/3)/1e3:.0f} km on a side)")
print(f"   fraction of Planck-sized regions that may collapse at the start: {beta:.1e}; ripple size needed at the Planck scale: {sigma:.3f}")
print(f"   (measured ripples on galaxy scales: ~2e-5 -> the start must be ~{sigma/2e-5:.0f}x rougher at the Planck scale than on large scales)")
print(f"   if the start were fully chaotic (ripples ~0.3-1): leftovers would outweigh what is seen by x{1/beta*0.1:.0e} or more -> excluded unless they evaporate completely")
print(f"   through Earth: {n*2.5e5*np.pi*6.371e6**2:.1f} per second -- interacting only by gravity (quantum-sensor arrays proposed to catch them, mem)")
