"""
Family 7 follow-up: the contact law sigma/m = C a0/(rho_DM v^2) switched on after recombination (z_s ~ 100-200).
For heavy dark matter (m >> proton, the case Famaey-Khoury-Penco need):
 - drag on the gas per unit time:      R_b   = rho_DM (sigma/m) v = C a0 / v       (independent of the particle mass)
 - gas heating/cooling rate:           Gam_b = 2 (m_p/m) C a0 / v                 (suppressed by m_p/m)
(A) 21-cm: the cold fluid can cool the gas at most by the fraction n_DM/(n_b+n_DM) of its heat (heat capacity), and at rate Gam_b.
(B) Lyman-alpha forest: gas pressure smooths the gas below the 'filtering scale' (measured ~80 kpc comoving, Rorai+2017, z~2).
    If the drag time 1/R_b is much shorter than the sound-crossing time of that scale, the gas is locked to the dark fluid and the
    smoothing disappears. v = relative speed or fluid velocity spread; the IGM fluid spread is < the gas sound speed, so using the
    gas sound speed is the most generous choice.
"""
import numpy as np
C = 1/16; a0 = 1.2e-8; mp = 0.938; H0 = 67.7e5/3.086e24; kpc = 3.086e21
print("(A) 21-cm (z ~ 17): maximum fractional cooling of the gas by the cold fluid, and rate vs expansion at z = 100")
for m in (1, 10, 100, 1000):
    nfrac = (5.3*mp/m)/(1 + 5.3*mp/m); H = H0*np.sqrt(0.31*101**3); Gb = 2*(mp/m)*C*a0/2e5
    print(f"   m = {m:5d} GeV: max gas cooling {nfrac:6.1%};  gas heating-rate/H at z=100: {Gb/H:6.3f}")
print("(B) Lyman-alpha filtering (most generous v = gas sound speed):")
for z, cs in ((2, 12e5), (3, 13e5), (5, 11e5)):
    R = C*a0/cs; lam = 80*kpc/(1 + z); tsound = lam/cs; H = H0*np.sqrt(0.31*(1 + z)**3 + 0.69)
    print(f"   z = {z}: drag time {1/R/3.15e13:6.1f} Myr, sound crossing of the filtering scale {tsound/3.15e13:6.0f} Myr -> gas locked "
          f"{tsound*R:5.1f}x faster than pressure can act; drag/expansion {R/H:5.0f}")
print("(C) inside a galaxy: a scattering contact also drags.  Rotating disc gas (v_rot) moving through a non-rotating halo (spread s):")
print("    drag rate on the gas = rho_DM (sigma/m) v_rel = C a0 v_rot / s^2")
for vr, s in ((100e5, 80e5), (200e5, 150e5), (250e5, 180e5)):
    R = C*a0*vr/s**2; print(f"   v_rot {vr/1e5:.0f} km/s, halo spread {s/1e5:.0f} km/s: gas spin-down time {1/R/3.15e16:5.2f} Gyr (discs are ~10 Gyr old)")
