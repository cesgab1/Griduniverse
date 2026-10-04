"""ITERATION 67 (pre-registered in PREREG_67.md)."""
import numpy as np
from scipy.optimize import brentq
rho, c, hbar_c, eV = 5.84e-27, 2.998e8, 3.1615e-26, 1.602176634e-19
m = (rho*c**2*hbar_c**3)**0.25/eV                                    # eV
dm21, dm31 = 7.42e-5, 2.51e-3
def masses(ml, order):
    if order == "normal": return ml, np.sqrt(ml**2 + dm21), np.sqrt(ml**2 + dm31)
    m2 = np.sqrt(ml**2 + dm31); return np.sqrt(m2**2 - dm21), m2, ml   # inverted: m3 lightest
out = ["ITERATION 67: dark energy sets the lightest neutrino mass? (expectations pre-registered)", "",
       f"rule: m c^2 = (rho_DE (hbar c)^3)^(1/4) = {m*1e3:.2f} meV   (h-version would be {m*(2*np.pi)**0.75*1e3:.1f} meV)"]
bound = 0.0642
for order in ("normal", "inverted"):
    ms = masses(m, order); s = sum(ms); mins = sum(masses(0.0, order))
    out.append(f"{order:8s}: masses {', '.join(f'{x*1e3:.1f}' for x in ms)} meV; sum {s*1e3:.1f} meV (oscillation minimum {mins*1e3:.1f}); "
               f"DESI+CMB bound {bound*1e3:.1f} -> {'PASS' if s < bound else 'EXCLUDED'}")
mmax = brentq(lambda ml: sum(masses(ml, 'normal')) - bound, 0, 0.1)
out.append(f"today's data allow (normal ordering) any lightest mass below {mmax*1e3:.1f} meV")
s_h = sum(masses(m*(2*np.pi)**0.75, "normal"))
out.append(f"h-version (not registered): sum {s_h*1e3:.1f} meV -> {'PASS' if s_h < bound else 'EXCLUDED'}")
lo, hi = 1e-4, 3e-2
out.append(f"look-elsewhere: of a log-uniform prior on the lightest mass from 0.1 to 30 meV, {np.log(mmax/lo)/np.log(hi/lo):.0%} passes "
           "today's bound -> passing is weak evidence")
txt = "\n".join(out); print(txt); open("iter67_neutrinos.txt", "w").write(txt + "\n")
