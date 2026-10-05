"""Iteration 95: cluster environment test (PREREG_95.md)."""
import numpy as np, glob, re
G, MSUN, A0, MPC = 6.674e-11, 1.989e30, 1.2e-10, 3.086e22
Gk = 4.30e-6   # kpc (km/s)^2 / Msun
kap = lambda s, M: (s*1e3)**4 / (G*M*MSUN*A0)
L = ["A. Virgo cluster (sigma 638 km/s)"]
for M in (1.65e13, 1.0e14):
    L.append(f"   M_b {M:.2e} Msun -> kappa {kap(638, M):.3f}  (= {kap(638, M)/kap(11.7, 3e7):.1f}x Fornax dwarf 0.039; Coma 0.21-0.42)")
gV = G*7e14*MSUN/(16.5*MPC)**2
L.append(f"   Virgo's pull at the Milky Way: {gV:.1e} m/s^2 = {gV/A0:.4f} a0")
# B. SPARC
rows = []
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    D = float(re.search(r"Distance = ([\d.]+)", open(fn).readline()).group(1))
    d = np.loadtxt(fn, comments="#", ndmin=2); r, V, eV, Vg, Vd, Vb = d.T[:6]
    ok = (r > 0) & (V > 0); r, V, eV, Vg, Vd, Vb = r[ok], V[ok], eV[ok], Vg[ok], Vd[ok], Vb[ok]
    if len(r) < 6 or abs(V[-1]-V[-3])/V[-3] > 0.10: continue
    Mb = (0.5*Vd[-1]**2 + 0.7*Vb[-1]**2 + 1.33*Vg[-1]*abs(Vg[-1]))*r[-1]/Gk
    if Mb <= 1e6: continue
    vf = V[-3:].mean()
    res = np.log10((vf*1e3)**4/(G*Mb*MSUN*A0))
    rows.append((fn.split("/")[-1][:-11], D == 18.0, res))
u = np.array([x[2] for x in rows if x[1]]); f = np.array([x[2] for x in rows if not x[1]])
diff = u.mean() - f.mean(); err = np.hypot(u.std(ddof=1)/np.sqrt(len(u)), f.std(ddof=1)/np.sqrt(len(f)))
L += ["", "B. SPARC: residual r = log10(v_flat^4 / G M_b a0)  (0 = exactly on the galaxy law)",
      f"   Ursa Major cluster members: N={len(u)}  mean {u.mean():+.3f}  scatter {u.std(ddof=1):.3f} dex",
      f"   field / other galaxies:     N={len(f)}  mean {f.mean():+.3f}  scatter {f.std(ddof=1):.3f} dex",
      f"   difference {diff:+.3f} +/- {err:.3f} dex  ({diff/err:+.1f} sigma);  in speed: {100*(10**(diff/4)-1):+.1f}%"]
# C. Milky Way
L += ["", "C. Milky Way (field galaxy)"]
for M in (5e10, 6e10, 7e10):
    L.append(f"   M_b {M:.0e} -> v_flat predicted {(G*M*MSUN*A0)**0.25/1e3:.0f} km/s   (observed ~200 at 25 kpc, Eilers+2019)")
txt = "\n".join(L); print(txt); open("iter95_environment.txt", "w").write(txt + "\n")
