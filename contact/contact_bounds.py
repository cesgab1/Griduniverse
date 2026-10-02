"""
Family 7 (new): the fluid is in CONTACT with visible matter (energy exchange), so halos settle to match the baryons.
Literature version: Famaey, Khoury & Penco 2018 (arXiv:1712.01316): DM-baryon energy exchange reproduces the RAR if the
cross-section obeys the 'master relation'  sigma = C a0 /(n_DM v^2),  C ~ 1/16  ->  sigma/m = C a0 /(rho_DM v^2).
Grid reading: each bit of visible matter trades energy with the fluid at a rate fixed per baryon (set by a0), however much
fluid is there. Quick checks: (1) the needed sigma/m in galaxies; (2) the same law in the early universe vs the Planck bound on
DM-baryon scattering (Gluscevic & Boddy 2018, velocity-independent, heavy DM: sigma/m < 8.8e-27 cm^2/GeV = 4.9e-3 cm^2/g);
(3) the exchange rate in the Lyman-alpha-era intergalactic gas: Gamma = rho_b (sigma/m) v = C a0 (rho_b/rho_DM)/v.
cgs units.
"""
import numpy as np
a0 = 1.2e-8; Msun_pc3 = 6.77e-23; H0 = 67.7e5/3.086e24; rho_dm0 = 0.26*3*H0**2/(8*np.pi*6.674e-8); fb = 0.19
bound = 4.9e-3
print("(1) galaxies (where halos must settle): needed sigma/m")
for C in (1/16, 1/10):
    for rho, v, lab in ((0.01, 100e5, "disc outskirts rho 0.01 Msun/pc^3, v 100 km/s"), (0.1, 200e5, "inner disc rho 0.1, v 200")):
        print(f"   C={C:.3f} {lab:46s} sigma/m = {C*a0/(rho*Msun_pc3*v**2):7.2f} cm^2/g")
print("(2) same law at recombination (z = 1100), relative DM-baryon speed 8.6 km/s (thermal) to 30 km/s (streaming):")
for C in (1/16, 1/10):
    for v in (8.6e5, 30e5):
        s = C*a0/(rho_dm0*1101**3*v**2); print(f"   C={C:.3f} v={v/1e5:4.1f} km/s: sigma/m = {s:.2e} cm^2/g = {s/bound:6.1f} x the Planck bound")
print("(3) exchange rate in intergalactic gas vs the expansion rate (Gamma/H):")
for z, v in ((3, 20e5), (100, 2e5), (1100, 8.6e5)):
    H = H0*np.sqrt(0.31*(1 + z)**3 + 0.69); G = (1/16)*a0*fb/v
    print(f"   z={z:5d}, v {v/1e5:4.1f} km/s: Gamma/H = {G/H:9.3g}")
