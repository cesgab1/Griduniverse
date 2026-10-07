"""Iteration 120: galaxy pattern vs viewing angle (PREREG_120.md)."""
import glob, os, numpy as np, pandas as pd
KPC = 3.086e19
inc = pd.read_csv("sparc_inclinations.csv", comment="#").set_index("name")
P = []
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    n = os.path.basename(fn).replace("_rotmod.dat", "")
    d = np.loadtxt(fn, comments="#", ndmin=2); r, V, eV, Vg, Vd, Vb = d.T[:6]; ok = (V > 0) & (eV < 0.1*V)
    for ri, Vi, g, dk, b in zip(r[ok], V[ok], Vg[ok], Vd[ok], Vb[ok]):
        vb2 = g*abs(g) + 0.5*dk**2 + 0.7*b**2
        if vb2 > 0: P.append((n, inc.loc[n, "inc"], vb2*1e6/(ri*KPC), Vi**2*1e6/(ri*KPC), g*abs(g) > 2*0.8*(dk**2 + 1.4*b**2)))
P = pd.DataFrame(P, columns=["name", "inc", "gbar", "gobs", "gasdom"])
P["lb"] = np.log10(P.gobs/P.gbar); P["lg"] = np.log10(P.gbar)
edges = np.arange(-12.4, -8.4, 0.2); cen = []; med = []
for a, b in zip(edges[:-1], edges[1:]):
    s = (P.lg >= a) & (P.lg < b)
    if s.sum() >= 10: cen.append((a+b)/2); med.append(P.lb[s].median())
P["res"] = P.lb - np.interp(P.lg, cen, med)                      # residual from the all-galaxy curve
out = ["Iteration 120 -- galaxy pattern by viewing angle (SPARC, mass-to-light 0.5)"]
for lab, lo, hi in [("low tilt 30-50", 30, 50), ("middle 50-70", 50, 70), ("high tilt 70-90", 70, 91)]:
    for sub, m in [("all points", np.ones(len(P), bool)), ("gas-dominated", P.gasdom.values)]:
        s = (P.inc >= lo) & (P.inc < hi) & m
        boots = [np.median(np.random.default_rng(i).choice(P.res[s], s.sum())) for i in range(1000)]
        bins = "; ".join(f"1e{c:+.0f}: x{10**P.lb[s & (abs(P.lg-c) < 0.25)].median():.2f}" for c in (-12, -11, -10, -9) if (s & (abs(P.lg-c) < 0.25)).sum() >= 5)
        out.append(f"  {lab:16s} {sub:13s}: N {s.sum():4d} pts / {P.name[s].nunique():3d} galaxies; offset from all-galaxy curve "
                   f"{P.res[s].median():+.3f} dex (+/-{np.std(boots):.3f}); scatter {P.res[s].std():.3f} dex | {bins}")
txt = "\n".join(out); print(txt); open("iter120_angle.txt", "w").write(txt + "\n")
