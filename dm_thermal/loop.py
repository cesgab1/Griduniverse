"""Dark matter on a loop (pieces keep breaking and re-forming): what the data allow (Coalesce)."""
import numpy as np
wb_bbn, e_bbn = 0.0222, 0.0005        # deuterium (mem, also used in iteration 93)
wb_cmb, e_cmb = 0.02239, 0.00015      # Planck compressed prior used in fit_law.py
wc = 0.120
d = wb_cmb - wb_bbn; e = np.hypot(e_bbn, e_cmb); up = d + 2*e
print(f"ordinary matter: element-making (3 min) {wb_bbn} vs early-universe glow (380,000 yr) {wb_cmb}: change {d:+.5f} +/- {e:.5f}")
print(f"-> at most {up:.5f} added (2 sigma) = {up/wb_bbn:.1%} of ordinary matter = {up/wc:.1%} of dark matter could have leaked into")
print("   ordinary matter between 3 minutes and 380,000 years")
print("a loop that never touches ordinary matter or light: unconstrained by this (a 'dark sector' with its own internal cycle)")
