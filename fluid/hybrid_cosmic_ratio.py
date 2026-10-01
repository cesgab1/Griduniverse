"""
One rule for every system: fluid halo mass = R x (visible baryonic mass), R = Omega_fluid/Omega_b = 5.4 (cosmic ratio, not fitted),
pulling with plain Newtonian gravity; MOND acts on the baryons (through the geometry, so lensing sees it).
Halo profile: NFW (concentration of a halo of total mass (1+R) M_b) or cored (superfluid-like core).
Tests: KiDS lensing (R free, then R = 5.4 fixed), clusters (R = 5.4), SPARC rotation curves (R = 5.4: does the halo spoil MOND fits?).
"""
import numpy as np, glob
from scipy.optimize import minimize_scalar
exec(open("hybrid_test.py").read().split("def chi2")[0])
Rcos = 0.26/0.048
def chi2R(R, core):
    tot = 0
    for (gb, go, Ci, Mh, r), lM, fg in zip(data, logMs, fgas):
        Mb = 10**lM*(1+fg); Mf = R*Mb
        Menc = halo_menc(Mf + Mb, zl, core)                         # profile of the combined halo, scaled to fluid mass
        m = nu(gb/a0)*gb + G*Mf/(Mf+Mb)*Menc(r)/r**2; res = go - m; tot += res @ Ci @ res
    return tot
print("KiDS lensing (60 pts): MOND(baryons) + R x M_b fluid halo, Newtonian")
for core in (0.0, 30.0, 100.0):
    b = minimize_scalar(lambda R: chi2R(R, core), bounds=(0, 20), method="bounded")
    Rs = np.linspace(0, 20, 801); c2 = np.array([chi2R(R, core) for R in Rs]); ok = Rs[c2 < b.fun + 4]
    print(f"   {'NFW' if core == 0 else f'cored {core:.0f} kpc':14s}: MOND only {chi2R(0, core):6.1f} | best R = {b.x:.1f} (2 sigma {ok.min():.1f}-{ok.max():.1f}), chi2 {b.fun:6.1f} | R = 5.4 (cosmic): chi2 {chi2R(Rcos, core):6.1f}")
# clusters with R = 5.4
def cluster_Mb(r):
    rc = 150*kpc; x = r/rc; m = x - np.arctan(x); m1 = 1e3*kpc/rc; m1 = m1 - np.arctan(m1)
    return 1e14*Msun*m/m1 + 1e12*Msun*r**2/(r + 20*kpc)**2
rr = np.array([0.1, 0.3, 0.6, 1.0])*Mpc; gbc = G*cluster_Mb(rr)/rr**2; dat = 10**(0.52*np.log10(gbc) - 4.19)
Mb_cl = cluster_Mb(1.5*Mpc)/Msun
print(f"\nCluster, R = 5.4: g_obs / measured at 0.1, 0.3, 0.6, 1.0 Mpc")
for core in (0.0, 30.0, 100.0):
    Mf = Rcos*Mb_cl; Menc = halo_menc(Mf + Mb_cl, 0.3, core)
    g = nu(gbc/a0)*gbc + G*Mf/(Mf+Mb_cl)*Menc(rr)/rr**2
    print(f"   {'NFW' if core == 0 else f'cored {core:.0f} kpc':14s}: " + "  ".join(f"{v:.2f}" for v in g/dat))
# SPARC with R = 5.4
kpc_ = 3.0857e19
GO, GB, EL, PH = [], [], [], []
tot0 = tot1 = tot2 = 0; npt = 0
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(fn, comments="#", ndmin=2); r, V, eV, Vg, Vd, Vb, _, _ = d.T
    gb = (Vg*abs(Vg) + 0.5*Vd**2 + 0.7*Vb**2)*1e6/(r*kpc_); go = V**2*1e6/(r*kpc_); ok = (V > 0)&(eV/V < 0.1)&(gb > 0)&(r > 0)
    if ok.sum() < 3: continue
    Mb = (V[ok][-1]*0 + (Vg[-1]*abs(Vg[-1]) + 0.5*Vd[-1]**2 + 0.7*Vb[-1]**2)*1e6*r[-1]*kpc_/G)/Msun      # baryonic mass inside last point
    if Mb <= 0: continue
    el = np.hypot(2*eV[ok]/V[ok]/np.log(10), 0.1)
    mond = nu(gb[ok]/a0)*gb[ok]
    for core, acc in ((30.0, 1), (100.0, 2)):
        Mf = Rcos*Mb; Menc = halo_menc(Mf + Mb, 0.0, core); rm = r[ok]*kpc_
        mod = mond + G*Mf/(Mf+Mb)*Menc(rm)/rm**2
        c = np.sum(((np.log10(go[ok]) - np.log10(mod))/el)**2)
        if acc == 1: tot1 += c
        else: tot2 += c
    tot0 += np.sum(((np.log10(go[ok]) - np.log10(mond))/el)**2); npt += ok.sum()
print(f"\nSPARC rotation curves ({npt} pts, M/L 0.5): chi2 MOND only {tot0:.0f} | + 5.4 M_b halo cored 30 kpc {tot1:.0f} | cored 100 kpc {tot2:.0f}")
