"""Iteration 116: fuzzy (wave) dark matter vs all MW dwarfs with measured sigma (PREREG_116.md)."""
import os, numpy as np, pandas as pd
here = os.path.dirname(os.path.abspath(__file__))
hbar, c, eV, PC = 1.0546e-34, 2.998e8, 1.602e-19, 3.086e16
d = pd.read_csv(os.path.join(here, "simon2019_dwarfs.csv"), comment="#")
k = np.sqrt(3.93/3)
def mmin(sig, R): return k*hbar/(sig*1e3*4/3*R*PC)*c**2/eV
d["m_min"] = mmin(d.sig, d.R12); d["m_min_cons"] = mmin(d.sig + d.s_hi, d.R12 + d.R_hi)
d["sig_r"] = d.sig*4/3*d.R12
d = d.sort_values("m_min", ascending=False)
out = ["Iteration 116 -- minimum particle mass if each dwarf's dark core is no more compact than one quantum wave allows",
       f"N = {len(d)} dwarfs (Simon 2019, measured sigma)", "",
       "  dwarf              R1/2 pc  sigma   m_min (eV)   conservative (sigma+1, R+1 err)"]
for _, r in d.iterrows():
    out.append(f"  {r['name']:17s} {r.R12:6.0f} {r.sig:6.1f}   {r.m_min:.2e}     {r.m_min_cons:.2e}")
out += ["", "  dwarfs above:  " + ", ".join(f"{x:.0e} eV: {np.sum(d.m_min > x)} (conservative {np.sum(d.m_min_cons > x)})" for x in (1e-22, 1e-21, 2e-21, 1e-20, 1e-19)),
        f"  strongest three: " + ", ".join(f"{n} {m:.1e}" for n, m in zip(d.name[:3], d.m_min_cons[:3])) + " (conservative)",
        f"  sigma x r1/2 spread: {d.sig_r.min():.0f} to {d.sig_r.max():.0f} km/s pc = x{d.sig_r.max()/d.sig_r.min():.0f}"
        "  (one soliton mass for all would need this to be constant)",
        f"  mass that makes Fornax's core a soliton: {mmin(11.7, 792):.1e} eV -> excluded by {np.sum(d.m_min_cons > 3*mmin(11.7, 792))} dwarfs at >x3 (conservative)"]
txt = "\n".join(out); print(txt); open(os.path.join(here, "iter116_fuzzy.txt"), "w").write(txt + "\n")
