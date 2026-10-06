"""Quark nuggets as dark matter in the verified window 1e7 - 1e15 kg: how they behave and how they could be caught."""
import numpy as np
rho_nuc, Msun = 4e17, 1.989e30
rho_loc = 0.4*1.783e-27/1e-6; v = 2.5e5; R_E = 6.371e6; rho_E = 5500.0; yr = 3.156e7; kt = 4.184e12
print(f"{'mass (kg)':>10s} {'size':>8s} {'collision area/mass (cm^2/g)':>29s} {'nuggets in Milky Way halo':>26s} {'Earth hit every':>16s} {'impact energy':>16s} {'deposited crossing Earth':>26s}")
for M in (1e7, 1e9, 1e11, 1e13, 1e15):
    R = (3*M/(4*np.pi*rho_nuc))**(1/3); A = np.pi*R**2
    sig = A/M*1e4/1e3                     # cm^2/g
    n_halo = 1e12*Msun/M
    every = 1/(rho_loc*v/M*np.pi*R_E**2*yr)
    KE = 0.5*M*v**2; dep = A*rho_E*v**2*2*R_E
    size = f"{R*1e3:.1f} mm" if R < 0.01 else f"{R*100:.1f} cm"
    print(f"{M:10.0e} {size:>8s} {sig:29.1e} {n_halo:26.1e} {every:>12.0f} yr {KE/kt:>11.0f} kt {dep/kt:>21.2f} kt")
print("""
Bullet-Cluster limit on collisions: < ~1 cm^2/g -> nuggets are ~1e13-1e16 times below it: they pass through like cold dark matter.
A galaxy holds 1e27-1e35 of them: smooth, so on cosmic and galaxy scales they act exactly like cold dark matter (same successes,
same galaxy-regularity puzzle). Distinguishing signal is only local: a nugget crossing Earth is a hypersonic (250 km/s) line source
depositing ~100 kt (1e7 kg) to ~2e7 kt (1e15 kg) along a 12,700 km track (entry and exit points), once per ~14 yr (1e7 kg) to ~1e9 yr (1e15 kg).
""")
