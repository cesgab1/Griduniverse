"""
Can the dark fluid's pressure keep it out of galaxy halos but let it into clusters?
Fluid equation of state P = kappa rho^Gamma (Gamma = 2: dense channel gas / Khronon quadratic K;  Gamma = 3: dilute Tonks channel gas,
the MOND-superfluid law). Sound speed c_s^2 = Gamma kappa rho^(Gamma-1) (units of c^2).
Allowed normalisation (largest kappa permitted by):
  (i)  early universe: w = kappa rho^(Gamma-1) <= 0.016 at a = 10^-4.5 (GDM bound quoted by Blanchet & Skordis 2024)
  (ii) Lyman-alpha: crossover length 2 pi c_s / sqrt(4 pi G rho) at z = 3, mean density, <= 0.3 Mpc comoving (condensate_jeans.txt)
Jeans mass at halo density (virial overdensity 200 x mean fluid density at formation z): M_J = (pi^(5/2)/6) c_s^3 / (G^(3/2) rho^(1/2)).
Fluid enters halos with M > M_J, stays out of lighter ones. Want: 1e12 < M_J < 1e14.5 Msun.
"""
import numpy as np
G, c, Msun, Mpc = 6.674e-11, 2.998e8, 1.989e30, 3.0857e22
h = 0.68; H0 = 100*h*1e3/Mpc; rhoc = 3*H0**2/(8*np.pi*G); OmK = 0.26
rho = lambda z: OmK*rhoc*(1+z)**3                                   # kg/m^3
def kappa_max(Gam):
    k1 = 0.016/rho(10**4.5 - 1)**(Gam-1)                               # w bound (kappa in units with c=1: P/(rho c^2))
    z = 3; lam = 0.3*Mpc/(1+z)                                         # physical
    cs2_max = (lam*np.sqrt(4*np.pi*G*rho(z))/(2*np.pi))**2/c**2
    k2 = cs2_max/(Gam*rho(z)**(Gam-1))
    return min(k1, k2), k1, k2
def MJ(Gam, kap, z, Delta=200):
    r = Delta*rho(z); cs = np.sqrt(Gam*kap*r**(Gam-1))*c
    return np.pi**2.5/6*cs**3/(G**1.5*np.sqrt(r))/Msun
for Gam in (2, 3):
    k, k1, k2 = kappa_max(Gam)
    print(f"Gamma = {Gam}: max kappa from early-w bound {k1:.2e}, from Lyman-alpha {k2:.2e} -> binding: {'early w' if k1 < k2 else 'Lyman-alpha'}")
    for zf, lab in ((2.0, "galaxy halos forming z~2"), (0.5, "clusters forming z~0.5")):
        print(f"   {lab:26s}: Jeans mass at virial density = {MJ(Gam, k, zf):.1e} Msun (largest allowed)")
print("\nWindow needed: fluid kept out of galaxies (M_J > 1e12 at z~2) yet held by clusters (M_J < 3e14 at z~0.5).")
