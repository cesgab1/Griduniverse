"""Iteration 117: local unseen mass at the Sun -- 'around' predicts zero (PREREG_117.md)."""
import os, numpy as np
tot, etot = 68.0, 4.0
out = ["Iteration 117 -- is the extra mass IN the Milky Way at the Sun, or only around it?"]
for lab, v, e in [("visible, McKee+2015 within 1.1 kpc", 43.8, 3.4*43.8/47.1),
                  ("visible, Bovy&Rix stars 38+/-4 + McKee gas 13.7+/-1.6", 38 + 13.7, np.hypot(4, 1.6))]:
    ex, ee = tot - v, np.hypot(etot, e)
    out.append(f"  {lab}: {v:.1f} +/- {e:.1f};  measured {tot:.0f} +/- {etot:.0f} -> unseen IN the disk column {ex:.1f} +/- {ee:.1f} Msun/pc^2 = {ex/ee:.1f} sigma above 0")
out.append(f"  cross-check: local dark density 0.013 +/- 0.003 Msun/pc^3 = {0.013/0.003:.1f} sigma; over 2.2 kpc of height ~ {0.013*2200:.0f} Msun/pc^2")
out.append("  'Around, not in' predicts 0 for every line above.")
txt = "\n".join(out); print(txt); open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "iter117_local.txt"), "w").write(txt + "\n")
