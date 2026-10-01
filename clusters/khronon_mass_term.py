"""
Clusters on the Khronon base. Quasi-static field equation (Blanchet & Skordis 2024, eq 3.21):
   div[(1+J_Y) grad phi] + mu^2 phi = 4 pi G rho_b
Spherical, phi < 0 in a well: f(g) g = G [M_b(<r) + M_mu(<r)] / r^2,  M_mu = (mu^2/G) Int |phi| r^2 dr   (extra 'condensate' mass)
phi(r) = -Int_r^R g dr', zero at the matching radius R where the object meets the cosmic background (environment: unknown -> scanned).
MOND function: the empirical McGaugh nu (same as our SPARC work). Iterated to convergence.
Data: clusters follow log g_obs = 0.52 log g_bar - 4.19 (CLASH BCG+cluster RAR, Tian et al. 2024; used in cluster_mix.py);
isolated-galaxy lensing (KiDS-1000) is MOND-like out to ~1 Mpc (our light-only ratio 0.91-0.94).
"""
import numpy as np, warnings
warnings.filterwarnings("ignore")
fmt = lambda v: f"{v:6.2f}" if np.isfinite(v) and v < 1e3 else "  runaway"
G, Msun, kpc, Mpc, a0 = 6.674e-11, 1.989e30, 3.0857e19, 3.0857e22, 1.2e-10
nu = lambda y: 1/(-np.expm1(-np.sqrt(np.maximum(y, 1e-30))))
def solve(Mb_of_r, mu_inv, R, r):
    mu2 = (1/(mu_inv*Mpc))**2 if mu_inv else 0.0
    Mmu = np.zeros_like(r)
    for _ in range(200):
        gN = G*(Mb_of_r(r) + Mmu)/r**2; g = nu(gN/a0)*gN
        phi = -np.concatenate([np.cumsum((g[1:]*np.diff(r))[::-1])[::-1], [0]])   # zero at r[-1] = R
        integrand = mu2/G*np.abs(phi)*r**2
        Mnew = np.concatenate([[0], np.cumsum(0.5*(integrand[1:]+integrand[:-1])*np.diff(r))])
        if np.max(abs(Mnew-Mmu)) < 1e-6*np.max(Mb_of_r(r)): break
        Mmu = 0.5*Mmu + 0.5*Mnew
    return g, gN, Mmu
def cluster_Mb(r):   # gas beta-model (r_c 150 kpc, beta 2/3, 1e14 Msun inside 1 Mpc) + BCG (Hernquist 1e12, a = 20 kpc)
    rc = 150*kpc; x = r/rc; m = x - np.arctan(x); m1 = 1e3*kpc/rc; m1 = m1 - np.arctan(m1)
    return 1e14*Msun*m/m1 + 1e12*Msun*r**2/(r + 20*kpc)**2
def galaxy_Mb(r):    # isolated L* galaxy: 1e11 Msun Hernquist a = 3 kpc (lensing regime r > 30 kpc)
    return 1e11*Msun*r**2/(r + 3*kpc)**2
print("Cluster: g_obs / (measured cluster RAR) at r = 0.1, 0.3, 0.6, 1.0 Mpc   [1.0 = matches data]")
for R_Mpc in (3.0, 10.0):
    r = np.geomspace(1*kpc, R_Mpc*Mpc, 4000)
    for mu_inv in (None, 22.0, 2.0, 1.0, 0.5, 0.25):
        g, gN, Mmu = solve(cluster_Mb, mu_inv, R_Mpc*Mpc, r)
        gb = G*cluster_Mb(r)/r**2; data = 10**(0.52*np.log10(gb) - 4.19)
        vals = [np.interp(x*Mpc, r, g/data) for x in (0.1, 0.3, 0.6, 1.0)]
        lab = "MOND only" if mu_inv is None else f"1/mu = {mu_inv:>4} Mpc"
        print(f"   R = {R_Mpc:>4} Mpc, {lab:16s}: " + " ".join(fmt(v) for v in vals) + "   (condensate/baryon mass at 1 Mpc " + fmt(np.interp(Mpc, r, Mmu)/cluster_Mb(Mpc)) + ")")
print("\nIsolated galaxy (lensing): g_obs / MOND-only at r = 0.1, 0.3, 0.6, 1.0 Mpc   [KiDS wants ~1]")
for R_Mpc in (1.5, 3.0):
    r = np.geomspace(1*kpc, R_Mpc*Mpc, 4000)
    g0, _, _ = solve(galaxy_Mb, None, R_Mpc*Mpc, r)
    for mu_inv in (22.0, 2.0, 1.0, 0.5, 0.25):
        g, gN, Mmu = solve(galaxy_Mb, mu_inv, R_Mpc*Mpc, r)
        print(f"   R = {R_Mpc:>4} Mpc, 1/mu = {mu_inv:>4} Mpc: " + " ".join(fmt(np.interp(x*Mpc, r, g/g0)) for x in (0.1, 0.3, 0.6, 1.0)))
