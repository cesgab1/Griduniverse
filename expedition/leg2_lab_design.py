"""Expedition leg 2 -- lab design: does the galaxy gap appear between small lab masses whose OWN pull is weak?
H_total (gap set by total pull; Earth's 9.8 m/s^2 switches it off): signal = Newton everywhere.
H_self  (gap set by the source's own pull): signal = Newton x galaxy gap curve at that pull (measured SPARC curve).
Ratio method: same masses, two separations; G and calibration cancel in the ratio (systematics partly cancel)."""
import glob, numpy as np
KPC = 3.086e19; G = 6.674e-11
lg, lb = [], []
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(fn, comments="#", ndmin=2); r, V, eV, Vg, Vd, Vb = d.T[:6]
    vb2 = Vg*np.abs(Vg) + 0.5*Vd**2 + 0.7*Vb**2; ok = (V > 0) & (eV < 0.1*V) & (vb2 > 0)
    lg += list(np.log10(vb2[ok]*1e6/(r[ok]*KPC))); lb += list(np.log10(V[ok]**2/vb2[ok]))
lg, lb = np.array(lg), np.array(lb); c, m = [], []
for a in np.arange(-12.6, -8.2, 0.2):
    s = (lg >= a) & (lg < a + 0.2)
    if s.sum() >= 10: c.append(a + 0.1); m.append(np.median(lb[s]))
gap = lambda g: 10**np.interp(np.log10(g), c, m, left=m[0], right=0.0)   # measured curve; x1 above its range
SYS, STAT = 4e-11, 4e-12                                                    # Westphal+2021 demonstrated
out = ["Expedition leg 2 -- lab design (predictions put on record before any such measurement)",
       f"galaxy gap curve used: measured SPARC medians, range 1e{min(c):.1f}..1e{max(c):.1f} m/s^2", ""]
for lab, M, r1, r2 in [("90 mg gold (Westphal-class), 2.4 mm vs 10 mm", 9e-5, 2.4e-3, 10e-3),
                       ("90 mg gold, 2.4 mm vs 20 mm", 9e-5, 2.4e-3, 20e-3),
                       ("1 g tungsten, 5 mm vs 50 mm", 1e-3, 5e-3, 50e-3)]:
    g1, g2 = G*M/r1**2, G*M/r2**2
    s1, s2 = g1*gap(g1), g2*gap(g2)
    ratio_N, ratio_S = g2/g1, s2/s1
    extra = s2 - g2
    out.append(f"{lab}: own pull {g1:.1e} -> {g2:.1e} m/s^2")
    out.append(f"   far signal: Newton {g2:.1e}, H_self {s2:.1e} (x{s2/g2:.2f}); extra {extra:.1e} = {extra/STAT:.0f}x statistical, {extra/SYS:.1f}x systematic precision")
    out.append(f"   far/near ratio: Newton {ratio_N:.4f}, H_self {ratio_S:.4f} (x{ratio_S/ratio_N:.2f}) -- calibration and G cancel")
out += ["", "Quantum version (gravity-entanglement proposals, 1e-14 kg at 0.2 mm, memory): own pull ~1.7e-17 m/s^2 -- below the",
        "galaxy curve's measured range; extrapolating the deep-regime slope gives x~2600 faster entanglement (H_self) vs x1 (H_total)."]
txt = "\n".join(out); print(txt); open("leg2_lab_design.txt", "w").write(txt + "\n")
