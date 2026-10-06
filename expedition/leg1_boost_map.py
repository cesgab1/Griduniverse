"""Expedition leg 1 (exploratory, descriptive -- not a pre-registered test): if there were NO unseen matter, how much
stronger than Newton would gravity have to be, system by system? Boost = (pull measured) / (pull of visible matter)."""
import glob, numpy as np
KPC = 3.086e19
gb, go, rr = [], [], []
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(fn, comments="#", ndmin=2); r, V, eV, Vg, Vd, Vb = d.T[:6]
    vb2 = Vg*np.abs(Vg) + 0.5*Vd**2 + 0.7*Vb**2              # McGaugh+2016 mass-to-light 0.5 disk / 0.7 bulge
    ok = (vb2 > 0) & (V > 0) & (eV < 0.1*V)
    gb += list(vb2[ok]*1e6/(r[ok]*KPC)); go += list(V[ok]**2*1e6/(r[ok]*KPC)); rr += list(r[ok])
gb, go, rr = map(np.array, (gb, go, rr)); boost = go/gb
out = [f"SPARC: {len(gb)} points in 175 galaxies (velocity error < 10%)", "",
       "A. Boost sorted by the pull of visible matter (acceleration):"]
for lo in np.arange(-12.5, -8.4, 0.5):
    s = (np.log10(gb) >= lo) & (np.log10(gb) < lo + 0.5)
    if s.sum() > 10: out.append(f"   visible pull 1e{lo:+.1f}..1e{lo+0.5:+.1f} m/s^2: median boost x{np.median(boost[s]):.2f}  (scatter 16-84%: {np.percentile(boost[s],16):.2f}-{np.percentile(boost[s],84):.2f}, N={s.sum()})")
out.append("B. Same points sorted by DISTANCE from the galaxy centre:")
for lo, hi in [(0, 1), (1, 3), (3, 10), (10, 30), (30, 100)]:
    s = (rr >= lo) & (rr < hi)
    out.append(f"   {lo:3d}-{hi:3d} kpc: median boost x{np.median(boost[s]):.2f}  (16-84%: {np.percentile(boost[s],16):.2f}-{np.percentile(boost[s],84):.2f}, N={s.sum()})")
# spread comparison: how well does each variable organise the boost?
def spread(x):
    bins = np.quantile(x, np.linspace(0, 1, 9)); res = []
    for a, b in zip(bins[:-1], bins[1:]):
        s = (x >= a) & (x <= b); res.append(np.std(np.log10(boost[s])))
    return np.mean(res)
out.append(f"   scatter of log(boost) within bins: by acceleration {spread(np.log10(gb)):.3f} dex, by distance {spread(rr):.3f} dex")
out += ["", "C. Other systems (verified in this project unless marked):",
        "   lab (G measured, masses cm-m apart): x1.00000 +/- 0.00002",
        "   solar system (planet orbits, ranging): x1 to better than ~1e-9 (memory)",
        "   Sun's neighbourhood, vertical (iteration 117): x1.3-1.55",
        f"   Coma cluster at 2.5 Mpc (iteration 113, gas-measured / visible): x8.3; visible pull there ~{6.674e-11*2.03e14*1.989e30/(2500*KPC)**2:.1e} m/s^2",
        "   Bullet Cluster main galaxies, lensing (iteration 102): x5.3 -- AND located on the galaxies, not the gas",
        "   Milky Way + Andromeda timing (iteration 98): x28",
        "   whole universe, early (CMB, Planck): matter / atoms = x6.4",
        "   Segue 1 dwarf: x~1000+ (memory, Simon+2011)"]
txt = "\n".join(out); print(txt); open("leg1_boost_map.txt", "w").write(txt + "\n")
