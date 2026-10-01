"""
Test of the 're-melting' prediction with KiDS-1000 isolated-galaxy lensing (Brouwer et al. 2021), 4 stellar-mass bins,
lens redshifts ~0.1-0.5 (KiDS-bright; median ~0.25).
Prediction under the DBI phase change: a galaxy keeps the halo it captured at z ~ 1-3 while its well is deeper than v_crit(z).
Models for g_obs(g_bar), same baryons (stars + cold gas), full covariance:
 A  MOND only (halos gone / never formed)
 B  MOND + retained fluid halo: NFW with abundance-matched mass (Moster+2013), MOND acting on the total field (what Khronon gives)
 C  LCDM reference: Newtonian baryons + the same NFW halo
"""
import numpy as np
G_pc = 4.52e-30; pcm = 3.086e16; D = "/home/claude/kids/"
G, Msun, kpc, Mpc, c = 6.674e-11, 1.989e30, 3.0857e19, 3.0857e22, 2.998e8
nu = lambda y: 1/(-np.expm1(-np.sqrt(np.maximum(y, 1e-30)))); a0 = 1.2e-10
cov = np.loadtxt(D + "Fig-9_RAR-KiDS-isolated_Massbins_covmatrix.txt")
logMs = [10.0, 10.45, 10.7, 10.9]; fgas = [0.4, 0.15, 0.1, 0.07]       # median stellar mass per bin; cold-gas fraction (rough)
def moster_halo(Ms):      # invert Moster+2013 (z=0) stellar-to-halo relation
    Mh = np.logspace(10, 14, 4000); N, M1, b, g = 0.0351, 10**11.59, 1.376, 0.608
    ms = Mh*2*N/((Mh/M1)**-b + (Mh/M1)**g); return np.interp(Ms, ms, Mh)
def nfw(Mh, z):
    H = 68e3/Mpc*np.sqrt(0.31*(1+z)**3+0.69); rho_c = 3*H**2/(8*np.pi*G)
    r200 = (3*Mh*Msun/(4*np.pi*200*rho_c))**(1/3); cc = 9.0*(Mh/1e12)**-0.1/(1+z)**0.5; rs = r200/cc
    f = lambda x: np.log(1+x) - x/(1+x)
    Menc = lambda r: Mh*Msun*f(r/rs)/f(cc)
    phi0 = G*Mh*Msun/(rs*f(cc))
    return Menc, np.sqrt(phi0)/1e3
# v_crit at the lens redshift (best DBI point)
mu = 1/(300*Mpc); lam = 9.2e6; rho_bar = 0.26*3*(68e3/Mpc)**2/(8*np.pi*G)
def vcrit(z):
    s = 8*np.pi*G*rho_bar*(1+z)**3/(2*mu**2*c**2); x = s/np.sqrt(1+lam*s*s); return c*np.sqrt(max(1/np.sqrt(lam) - x, 0))/1e3
zl = 0.25
print(f"lens redshift ~{zl}: v_crit = {vcrit(zl):.0f} km/s (z=0.1: {vcrit(0.1):.0f}, z=0.4: {vcrit(0.4):.0f})\n")
tot = {"A": 0, "B": 0, "C": 0}
for i in range(4):
    d = np.loadtxt(D + f"Fig-9_RAR-KiDS-isolated_Massbin-{i+1}.txt"); gb = d[:, 0]; go = 4*G_pc*d[:, 1]/d[:, 4]*pcm
    n = len(gb); C = cov[i*n*4 + 0:0]  # placeholder
    sel = (cov[:, 0] == cov[:, 0][0])
    # pick the block for bin i (covariance rows ordered by bin, then radius)
    mins = np.unique(cov[:, 0]); blk = cov[(cov[:, 0] == mins[i]) & (cov[:, 1] == mins[i])]
    C = (blk[:, 4]/blk[:, 6]).reshape(n, n)*(4*G_pc*pcm)**2; Ci = np.linalg.inv(C)
    Ms = 10**logMs[i]; Mb = Ms*(1+fgas[i]); Mh = moster_halo(Ms); Menc, vdepth = nfw(Mh, zl)
    r = np.sqrt(G*Mb*Msun/gb)                                         # radius where the (point-like) baryons give g_bar
    gh = G*Menc(r)/r**2
    models = {"A": nu(gb/a0)*gb, "B": nu((gb+gh)/a0)*(gb+gh), "C": gb + gh}
    chis = {k: float((go-m) @ Ci @ (go-m)) for k, m in models.items()}
    for k in tot: tot[k] += chis[k]
    kept = "kept (well deeper than v_crit)" if vdepth > vcrit(zl) else "re-melted"
    print(f"bin {i+1}: log M* {logMs[i]}, halo M200 {Mh:.1e} Msun, halo well depth {vdepth:.0f} km/s -> predicted {kept}")
    print(f"        chi2 ({n} pts): A MOND-only {chis['A']:7.1f} | B MOND+halo {chis['B']:8.1f} | C LCDM {chis['C']:7.1f}"
          f" | median g_obs/model: A {np.median(go/models['A']):.2f}  B {np.median(go/models['B']):.2f}  C {np.median(go/models['C']):.2f}")
print(f"\nTotal chi2: A {tot['A']:.1f} | B {tot['B']:.1f} | C {tot['C']:.1f}")
