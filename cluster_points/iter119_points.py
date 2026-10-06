"""Iteration 119: boundary-free cluster/group points vs the galaxy rule (PREREG_119.md)."""
import os, numpy as np, pandas as pd
here = os.path.dirname(os.path.abspath(__file__))
G, MS, KPC, MPC = 6.674e-11, 1.989e30, 3.086e19, 3.086e22
nu = lambda x: 1/(1 - np.exp(-np.sqrt(x)))
def rhoc(z, h): return 3*(100*h*1e3/MPC)**2*(0.3*(1 + z)**3 + 0.7)/(8*np.pi*G)
pts = []   # (name, sample, label, r_m, M_kg, fg, efg)
# Vikhlinin+2006 r2500 (non-duplicate systems; M2500 in 1e14, fg2500) and z (Table 1)
V = [("A133", .0569, 1.24, .065, .002), ("A383", .1883, 1.68, .090, .005), ("A478", .0881, 4.23, .096, .004), ("A907", .1603, 2.30, .090, .003),
     ("A1413", .1429, 3.08, .092, .003), ("A1795", .0622, 2.66, .088, .003), ("A2029", .0779, 4.41, .090, .003), ("A2390", .2302, 3.50, .127, .005),
     ("USGCS152", .0153, 0.07, .043, .001)]
for n, z, M, fg, e in V:
    Mk = M*1e14*MS; r = (3*Mk/(4*np.pi*2500*rhoc(z, 0.72)))**(1/3); pts.append((n, "V06", "r2500", r, Mk, fg, e))
c = pd.read_csv(os.path.join(here, "../cluster_pairs/clusters.csv"), comment="#")
for _, x in c.iterrows():
    pts.append((x["name"], x["sample"], "r500", x.r500*KPC, x.M500*1e14*MS, x.fg, (x.fg_lo + x.fg_hi)/2))
s = pd.read_csv(os.path.join(here, "sun09_r2500.csv"), comment="#")
for _, x in s.iterrows():
    r = x.r*KPC; Mk = 2500*rhoc(x.z, 0.73)*4/3*np.pi*r**3; pts.append((x["name"], "S09", "r2500", r, Mk, x.fg, (x.fg_lo + x.fg_hi)/2))
P = pd.DataFrame(pts, columns=["name", "sample", "where", "r", "M", "fg", "efg"])
out = [f"Iteration 119 -- {len(P)} boundary-free points ({(P['where']=='r2500').sum()} inner, {(P['where']=='r500').sum()} outer)"]
for lab, fst in [("stars 0.025(M5/1e14)^-0.37", lambda M5: 0.025*M5**-0.37), ("stars 0", lambda M5: 0*M5), ("stars 0.02", lambda M5: 0.02 + 0*M5)]:
    M5 = np.where(P["where"] == "r500", P.M, 2*P.M)/(1e14*MS)          # stellar-fraction proxy uses ~M500
    fv = P.fg + fst(M5); gobs = G*P.M/P.r**2; gbar = fv*gobs
    off = gobs/(gbar*nu(gbar/1.2e-10))
    out += ["", f"[{lab}]"]
    for w in ("r2500", "r500"):
        m = P["where"] == w
        out.append(f"  {w:5s}: N {m.sum():2d}, visible pull {np.median(gbar[m]):.1e} m/s^2, total/visible x{np.median(1/fv[m]):.1f};"
                   f" measured / galaxy-rule prediction: median x{np.median(off[m]):.2f} (16-84%: {np.percentile(off[m],16):.2f}-{np.percentile(off[m],84):.2f})")
    for lo in (-12.5, -12.0, -11.5, -11.0, -10.5, -10.0):
        m = (np.log10(gbar) >= lo) & (np.log10(gbar) < lo + 0.5)
        if m.sum(): out.append(f"     visible pull 1e{lo:+.1f}..1e{lo+0.5:+.1f}: N {m.sum():2d}, offset x{np.median(off[m]):.2f}   [galaxies at same pull: offset x1 by definition]")
txt = "\n".join(out); print(txt); open(os.path.join(here, "iter119_points.txt"), "w").write(txt + "\n")
