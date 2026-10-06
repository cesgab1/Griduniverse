"""Conventions drills D1 + D2 (exploratory re-analyses).
D1: galaxy pattern from gas-dominated points only (gas from 21-cm radio, needs no mass-to-light convention).
D2: replace the assumed interpolation function nu with the galaxies' own measured boost curve, then apply it to the
    83 boundary-free cluster/group points of iteration 119."""
import glob, os, numpy as np
KPC = 3.086e19
D = []
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(fn, comments="#", ndmin=2); r, V, eV, Vg, Vd, Vb = d.T[:6]; ok = (V > 0) & (eV < 0.1*V)
    D.append((r[ok], V[ok], Vg[ok], Vd[ok], Vb[ok]))
def collect(Y, gas_dom=False):
    gb, go = [], []
    for r, V, Vg, Vd, Vb in D:
        stars = Y*Vd**2 + 1.4*Y*Vb**2; vb2 = Vg*np.abs(Vg) + stars
        k = (vb2 > 0) & ((Vg*np.abs(Vg) > 2*0.8*(Vd**2 + 1.4*Vb**2)) if gas_dom else True)   # gas > 2x stars even at M/L 0.8
        gb += list(vb2[k]*1e6/(r[k]*KPC)); go += list(V[k]**2*1e6/(r[k]*KPC))
    return np.array(gb), np.array(go)
out = ["D1 -- gas-dominated points only (gas > 2x stars even at the highest mass-to-light): convention-light pattern"]
for Y in (0.2, 0.5, 0.8):
    gb, go = collect(Y, True); B = go/gb
    out.append(f"  M/L {Y}: N {len(gb)}; boost " + "; ".join(
        f"1e{lo:+.1f}: x{np.median(B[(np.log10(gb) >= lo-0.25) & (np.log10(gb) < lo+0.25)]):.1f}" for lo in (-12.0, -11.5, -11.0, -10.5)
        if ((np.log10(gb) >= lo-0.25) & (np.log10(gb) < lo+0.25)).sum() > 5))
# D2: empirical galaxy curve (all points, M/L 0.5), median log boost in 0.2-dex bins
gb, go = collect(0.5); lg = np.log10(gb); lb = np.log10(go/gb)
edges = np.arange(-12.6, -8.4, 0.2); cen, med = [], []
for a, b in zip(edges[:-1], edges[1:]):
    s = (lg >= a) & (lg < b)
    if s.sum() >= 10: cen.append((a + b)/2); med.append(np.median(lb[s]))
cen, med = np.array(cen), np.array(med)
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../cluster_points/iter119_points.py")).read().split('out = [f"Iteration 119')[0]
src = src.replace("here = os.path.dirname(os.path.abspath(__file__))", "here = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../cluster_points')"); exec(src)
fv = P.fg + 0.025*(np.where(P["where"] == "r500", P.M, 2*P.M)/(1e14*MS))**-0.37
gobs = G*P.M/P.r**2; gbar = fv*gobs; lgb = np.log10(gbar)
inside = (lgb >= cen.min()) & (lgb <= cen.max())
pred = gbar*10**np.interp(lgb, cen, med); off = gobs/pred
out += ["", f"D2 -- galaxies' OWN measured boost curve (no assumed formula), range 1e{cen.min():.1f}..1e{cen.max():.1f} m/s^2",
        f"  cluster points inside that range: {inside.sum()} of {len(P)}",
        f"  inner points: measured / galaxy-curve x{np.median(off[inside & (P['where']=='r2500')]):.2f};  outer points: x{np.median(off[inside & (P['where']=='r500')]):.2f}",
        f"  (with the assumed formula, iteration 119: inner x2.44, outer x1.40)"]
txt = "\n".join(out); print(txt); open("drill2_conventionfree.txt", "w").write(txt + "\n")
