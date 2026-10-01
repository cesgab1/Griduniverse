"""
'Salt' + 'current' rule, no free numbers:  dark fluid = 5.4 x the baryons a system retained, spread like its hot gas
(fluid shares the gas's fate: heated and expelled with it from shallow wells, kept in deep ones), Newtonian pull;
MOND acts on the baryons (through the geometry, so lensing sees it).
   g_obs = nu(g_b/a0) g_b + 5.4 x g_b,distributed-like-gas
A. Groups & clusters at R500: observed gas fraction f_gas,500 = 2.23e-7 M500^0.39 (Popesso et al. 2024, arXiv:2411.16555),
   stellar fraction (stars + ICL) f_* = 0.025 (M500/1e14)^-0.45 (approximate, Gonzalez et al. 2013). Prediction vs observed M500.
B. KiDS-1000 isolated galaxies (Brouwer+21): stars + cold gas (their eq: log f_cold = -0.69 log M* + 6.63) + hot gas = M*
   with rho ~ r^-2 truncated at 143 kpc (their nominal model); fluid = 5.4 x all of it, spread like the hot gas.
   chi2 with full covariance, plus Brouwer's 0.1 dex error floor, plus the reliable region only (g_bar > 1e-13).
C. SPARC: same rule (hot gas = stellar mass, r^-2 to 143 kpc): effect on rotation curves inside the measured radii.
"""
import numpy as np, glob
G, Msun, kpc, Mpc = 6.674e-11, 1.989e30, 3.0857e19, 3.0857e22
a0 = 1.2e-10; nu = lambda y: 1/(-np.expm1(-np.sqrt(np.maximum(y, 1e-30)))); R = 0.26/0.048
H0 = 70e3/Mpc; rhoc = 3*H0**2/(8*np.pi*G)
print("A. Groups and clusters at R500: predicted dynamical mass / observed M500")
for lM in (13.0, 13.5, 14.0, 14.5, 15.0):
    M5 = 10**lM; fg = 2.23e-7*M5**0.39; fs = 0.025*(M5/1e14)**-0.45; Mb = (fg+fs)*M5
    R5 = (3*M5*Msun/(4*np.pi*500*rhoc))**(1/3); gb = G*Mb*Msun/R5**2
    pred = (nu(gb/a0) + R)*Mb
    print(f"   M500 = 1e{lM}: f_gas {fg:.3f}, f_* {fs:.3f}, baryon fraction {fg+fs:.3f} (cosmic 0.157) | MOND boost {nu(gb/a0):.2f}"
          f" | predicted/observed: MOND alone {nu(gb/a0)*Mb/M5:.2f}, MOND + 5.4x retained baryons {pred/M5:.2f}")
# B. KiDS
G_pc = 4.52e-30; pcm = 3.086e16; D = "/home/claude/kids/"
cov = np.loadtxt(D + "Fig-9_RAR-KiDS-isolated_Massbins_covmatrix.txt"); mins = np.unique(cov[:, 0])
logMs = [10.0, 10.45, 10.7, 10.9]; Racc = 143*kpc
tot = {}
print("\nB. KiDS isolated galaxies (chi2: full | with 0.1 dex floor | reliable g_bar>1e-13 only, with floor)")
for i in range(4):
    d = np.loadtxt(D + f"Fig-9_RAR-KiDS-isolated_Massbin-{i+1}.txt"); x = d[:, 0]; go = 4*G_pc*d[:, 1]/d[:, 4]*pcm; n = len(x)
    blk = cov[(cov[:, 0] == mins[i]) & (cov[:, 1] == mins[i])]; C = (blk[:, 4]/blk[:, 6]).reshape(n, n)*(4*G_pc*pcm)**2
    Ms = 10**logMs[i]; fcold = 10**(-0.69*logMs[i] + 6.63); Mgal = Ms*(1+fcold); Mhot = Ms
    r = np.sqrt(G*Mgal*Msun/x)                                             # data x-axis = stars + cold gas, point mass
    Mhot_r = Mhot*np.minimum(r/Racc, 1.0); Mb_r = Mgal + Mhot_r
    gb = G*Mb_r*Msun/r**2
    Mfl_r = R*(Mgal + Mhot)*np.minimum(r/Racc, 1.0)
    models = {"MOND (stars+cold)": nu(x/a0)*x, "MOND (+hot gas)": nu(gb/a0)*gb, "MOND(+hot) + fluid 5.4x": nu(gb/a0)*gb + G*Mfl_r*Msun/r**2}
    floor = np.diag((0.1*np.log(10)*go)**2); rel = x > 1e-13
    for k, m in models.items():
        res = go - m
        c1 = res @ np.linalg.inv(C) @ res; c2 = res @ np.linalg.inv(C + floor) @ res
        Cr = (C + floor)[np.ix_(rel, rel)]; c3 = res[rel] @ np.linalg.inv(Cr) @ res[rel]
        t = tot.setdefault(k, np.zeros(3)); t += [c1, c2, c3]
for k, t in tot.items(): print(f"   {k:28s}: {t[0]:8.1f} | {t[1]:7.1f} | {t[2]:6.1f}   (60 / 60 / {sum((np.loadtxt(D+f'Fig-9_RAR-KiDS-isolated_Massbin-{i+1}.txt')[:,0] > 1e-13).sum() for i in range(4))} pts)")
# C. SPARC
kpc_ = kpc; c0 = c1 = 0; npt = 0; fr30 = []
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(fn, comments="#", ndmin=2); r, V, eV, Vg, Vd, Vb, _, _ = d.T
    gb = (Vg*abs(Vg) + 0.5*Vd**2 + 0.7*Vb**2)*1e6/(r*kpc_); go = V**2*1e6/(r*kpc_); ok = (V > 0)&(eV/V < 0.1)&(gb > 0)&(r > 0)
    if ok.sum() < 3: continue
    Mst = (0.5*Vd[-1]**2 + 0.7*Vb[-1]**2)*1e6*r[-1]*kpc_/G/Msun; Mgas = max(Vg[-1]*abs(Vg[-1]), 0)*1e6*r[-1]*kpc_/G/Msun
    Mhot = Mst; rr = r[ok]*kpc_
    gtot = gb[ok] + G*Mhot*Msun*np.minimum(rr/Racc, 1)/rr**2
    mond = nu(gtot/a0)*gtot; fl = G*R*(Mst + Mgas + Mhot)*Msun*np.minimum(rr/Racc, 1)/rr**2
    el = np.hypot(2*eV[ok]/V[ok]/np.log(10), 0.1)
    c0 += np.sum(((np.log10(go[ok]) - np.log10(nu(gb[ok]/a0)*gb[ok]))/el)**2); c1 += np.sum(((np.log10(go[ok]) - np.log10(mond + fl))/el)**2); npt += ok.sum()
    fr30.append(np.median(fl/(mond)))
print(f"\nC. SPARC ({npt} pts): chi2 MOND alone {c0:.0f} | MOND(+hot gas) + 5.4x fluid spread like hot gas {c1:.0f} | median fluid/MOND pull {np.median(fr30):.2f}")
