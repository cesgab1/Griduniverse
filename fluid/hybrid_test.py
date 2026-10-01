"""
Hybrid 'ocean' model: superfluid halo whose MOND push acts through the geometry (Khronon-like), sourced by baryons only;
the fluid halo itself pulls with plain Newtonian gravity (no MOND boost), as in Berezhiani-Khoury.
   g_obs = nu(g_bar/a0) g_bar + g_halo
Halo: abundance-matched mass (Moster+2013) x f (f = 1: the full cosmic share), NFW, or cored (superfluid core ~ polytrope).
Tests: KiDS-1000 isolated-galaxy lensing (4 mass bins, full covariance) -> allowed f; clusters (CLASH/BCG RAR, Tian+2024) -> needed f.
"""
import numpy as np
from scipy.optimize import minimize_scalar
exec(open("/home/claude/griduniverse/clusters/remelt_kids_test.py").read().split("tot = {")[0].split("# v_crit at the lens")[0])
def halo_menc(Mh, z, core_kpc=0.0):
    Menc, _ = nfw(Mh, z)
    if core_kpc <= 0: return Menc
    rc = core_kpc*kpc
    return lambda r: Menc(r)*(r**3/(r**3 + rc**3))      # cored inside rc (superfluid-like core), NFW outside
zl = 0.25
data = []
mins = np.unique(cov[:, 0])
for i in range(4):
    d = np.loadtxt(D + f"Fig-9_RAR-KiDS-isolated_Massbin-{i+1}.txt"); gb = d[:, 0]; go = 4*G_pc*d[:, 1]/d[:, 4]*pcm; n = len(gb)
    blk = cov[(cov[:, 0] == mins[i]) & (cov[:, 1] == mins[i])]; Ci = np.linalg.inv((blk[:, 4]/blk[:, 6]).reshape(n, n)*(4*G_pc*pcm)**2)
    Ms = 10**logMs[i]; Mb = Ms*(1+fgas[i]); Mh = moster_halo(Ms); r = np.sqrt(G*Mb*Msun/gb)
    data.append((gb, go, Ci, Mh, r))
def chi2(f, core):
    tot = 0
    for gb, go, Ci, Mh, r in data:
        m = nu(gb/a0)*gb + f*G*halo_menc(Mh, zl, core)(r)/r**2; res = go - m; tot += res @ Ci @ res
    return tot
print("KiDS lensing, MOND(baryons) + f x Newtonian halo (60 points):")
for core in (0.0, 30.0, 100.0):
    b = minimize_scalar(lambda f: chi2(f, core), bounds=(0, 3), method="bounded")
    fs = np.linspace(0, 3, 601); c2 = np.array([chi2(f, core) for f in fs]); ok = fs[c2 < b.fun + 4]
    print(f"   halo {'NFW' if core == 0 else f'cored {core:.0f} kpc':14s}: chi2(f=0, MOND only) {chi2(0, core):6.1f} | best f = {b.x:.2f} chi2 {b.fun:6.1f}"
          f" | f allowed (2 sigma) {ok.min():.2f}-{ok.max():.2f} | chi2(f=1) {chi2(1, core):7.1f}")
# clusters: same rule
exec(open("/home/claude/griduniverse/clusters/khronon_mass_term.py").read().split('print("Cluster')[0].split("def solve")[0])
def cluster_Mb(r):
    rc = 150*kpc; x = r/rc; m = x - np.arctan(x); m1 = 1e3*kpc/rc; m1 = m1 - np.arctan(m1)
    return 1e14*Msun*m/m1 + 1e12*Msun*r**2/(r + 20*kpc)**2
rr = np.array([0.1, 0.3, 0.6, 1.0])*Mpc
gbc = G*cluster_Mb(rr)/rr**2; dat = 10**(0.52*np.log10(gbc) - 4.19)
Mb_tot = cluster_Mb(1.5*Mpc)/Msun; Mh_cl = Mb_tot*(0.31/0.048)          # cosmic ratio halo for the cluster
print(f"\nCluster (baryons {Mb_tot:.1e} Msun, cosmic-share halo {Mh_cl:.1e}): g_obs / measured cluster RAR at r = 0.1, 0.3, 0.6, 1.0 Mpc")
for f in (0.0, 0.5, 1.0):
    for core in (0.0, 100.0):
        g = nu(gbc/a0)*gbc + f*G*halo_menc(Mh_cl, 0.3, core)(rr)/rr**2
        print(f"   f = {f:.1f}, halo {'NFW' if core == 0 else 'cored 100 kpc':14s}: " + "  ".join(f"{v:.2f}" for v in g/dat))
