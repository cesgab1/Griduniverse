"""Conventions drill 1 (exploratory robustness check): does the galaxy pattern (boost set by weakness of pull) survive the
borrowed mass-to-light convention? Disk 0.2-0.8 (bulge = 1.4 x disk), covering Chabrier/Kroupa/Salpeter-type choices."""
import glob, numpy as np
KPC = 3.086e19
D = []
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(fn, comments="#", ndmin=2); r, V, eV, Vg, Vd, Vb = d.T[:6]; ok = (V > 0) & (eV < 0.1*V)
    D.append((r[ok], V[ok], Vg[ok], Vd[ok], Vb[ok]))
out = ["Drill 1 -- galaxy pattern vs the mass-to-light convention (SPARC)"]
for Y in (0.2, 0.35, 0.5, 0.65, 0.8):
    gb, go = [], []
    for r, V, Vg, Vd, Vb in D:
        vb2 = Vg*np.abs(Vg) + Y*Vd**2 + 1.4*Y*Vb**2; k = vb2 > 0
        gb += list(vb2[k]*1e6/(r[k]*KPC)); go += list(V[k]**2*1e6/(r[k]*KPC))
    gb, go = np.array(gb), np.array(go); B = go/gb
    def spread(x):
        q = np.quantile(x, np.linspace(0, 1, 9)); return np.mean([np.std(np.log10(B[(x >= a) & (x <= b)])) for a, b in zip(q[:-1], q[1:])])
    bins = "; ".join(f"1e{lo:+.0f}: x{np.median(B[(np.log10(gb) >= lo-0.25) & (np.log10(gb) < lo+0.25)]):.1f}" for lo in (-12, -11, -10, -9))
    out.append(f"  M/L disk {Y:.2f}: boost by visible pull {bins}; scatter {spread(np.log10(gb)):.3f} dex; boost < 1 (impossible, too much star mass) in {np.mean(B < 0.8)*100:.0f}% of points")
txt = "\n".join(out); print(txt); open("drill1_mass_to_light.txt", "w").write(txt + "\n")
