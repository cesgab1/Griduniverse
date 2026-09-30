"""
(1) Where does confinement happen, measured in grid spacings?  Run the measured strong coupling alpha_s(M_Z) = 0.1180 up to the grid
    scale (1/l_d, l_d = 0.60 l_P) with the Standard Model's two-loop beta function (quark thresholds at m_b, m_t),
    then down to the confinement scale Lambda_QCD. Grid coupling beta = 6/g^2 = 6/(4 pi alpha_s) for SU(3).
(2) Weak force on the grid: chiral (left-handed-only) fermions with a gauge field are consistent on a grid only if the gauge anomalies
    cancel (Luscher's construction). Check the Standard Model's charges generation by generation.
"""
import numpy as np
from scipy.integrate import solve_ivp
from fractions import Fraction as F
MZ, mb, mt = 91.1876, 4.18, 172.6
EP = 1.220890e19; Egrid = EP / 0.603                     # GeV
def beta(t, a, nf):
    b0 = 11 - 2 * nf / 3; b1 = 102 - 38 * nf / 3
    return [-(b0 / (2 * np.pi)) * a[0]**2 - (b1 / (8 * np.pi**2)) * a[0]**3]
def run(a0, mu0, mu1, nf):
    return solve_ivp(beta, [np.log(mu0), np.log(mu1)], [a0], args=(nf,), rtol=1e-10).y[0, -1]
a_t = run(0.1180, MZ, mt, 5); a_grid = run(a_t, mt, Egrid, 6)
print(f"alpha_s at the grid scale ({Egrid:.2e} GeV): {a_grid:.4f}  ->  SU(3) grid coupling beta = 6/(4 pi alpha) = {6/(4*np.pi*a_grid):.1f}")
a_b = run(0.1180, MZ, mb, 5)
# Lambda where alpha_s (nf=3 below m_b ~ approx) blows up: integrate down until alpha = 1
mu, a = mb, a_b
while a < 1.0 and mu > 0.05:
    mu2 = mu * 0.98; a = run(a, mu, mu2, 4 if mu > 1.27 else 3); mu = mu2
print(f"strong coupling reaches 1 at {mu*1e3:.0f} MeV (confinement scale; proton mass 938 MeV)")
print(f"confinement length / grid spacing = {Egrid/mu:.1e}   (= e^{np.log(Egrid/mu):.1f}: set by the logarithmic running, no fine-tuning)")
# sensitivity: what grid-scale coupling changes would move the proton mass by x2?
a_hi = run(a_t, mt, Egrid, 6)
for f in (0.97, 1.03):
    a_back = run(run(a_grid * f, Egrid, mt, 6), mt, MZ, 5)
    print(f"  grid coupling x{f}: alpha_s(M_Z) -> {a_back:.4f}")

print("\nGauge-anomaly cancellation, one generation (left-handed Weyl fields; hypercharge Y, Q = T3 + Y):")
fields = [("quark doublet Q", 3, 2, F(1, 6)), ("up-type u^c", 3, 1, F(-2, 3)), ("down-type d^c", 3, 1, F(1, 3)),
          ("lepton doublet L", 1, 2, F(-1, 2)), ("electron e^c", 1, 1, F(1, 1))]
Y3 = sum(c * w * y**3 for _, c, w, y in fields); Ygrav = sum(c * w * y for _, c, w, y in fields)
SU2Y = sum(c * y for _, c, w, y in fields if w == 2); SU3Y = sum(w * y for _, c, w, y in fields if c == 3)
nd = sum(c for _, c, w, _ in fields if w == 2)
print(f"  [U(1)_Y]^3 = {Y3}   [grav]^2 U(1)_Y = {Ygrav}   [SU(2)]^2 U(1)_Y = {SU2Y}   [SU(3)]^2 U(1)_Y = {SU3Y}   SU(2) doublets = {nd} (even: no Witten anomaly)")
print("  -> all cancel, so the Standard Model is exactly the kind of chiral theory a grid can hold; drop any one field and it cannot.")
fields2 = fields[:-1]
print(f"  e.g. without the electron: [U(1)_Y]^3 = {sum(c*w*y**3 for _, c, w, y in fields2)}  (anomalous)")
