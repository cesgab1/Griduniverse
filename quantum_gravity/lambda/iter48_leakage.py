"""
ITERATION 48 (pre-registered in quantum_gravity/dynamics/PREREG_47_49.md, L1-L3): the leakage factor c.
The one frame-dependent field enters only through Claim 1's function of the slice expansion K (plus the slicing time).
"""
import numpy as np
from scipy.optimize import brentq
Om, Or = 0.315, 9.1e-5; Ode = 1 - Om - Or
out = ["ITERATION 48: leakage of the frame-dependent field into gravity (expectations pre-registered)", ""]
# L1: tree level (BORROWED, VCDM / type-II minimally modified gravity): two graviton polarisations, c_GW = 1 exactly,
#     G_eff = G for sub-horizon modes, static stars = GR (TOV). Check that Claim 1 is of that form: its action density is a function
#     f(K, t) only; f carries no graviton (traceless) piece at linear order because K is a trace.
out += ["L1 tree level: Claim 1 = function of K (a trace) and the slicing time -> VCDM class (BORROWED): gravitational-wave speed",
        "   exactly 1, G_eff = G below the horizon, static stars as in GR. Tree-level leakage c = 0.", ""]
# L2: second derivative of f with respect to K shifts the K^2 weight: f = -(2/3) C (a K/3)^(-1/2)
#     d2f/dK2 = -rho_DE/(18 H^2) = -(M_P^2) Omega_DE/6  ->  Delta lambda = Omega_DE/6
def E(a):
    g = lambda le: np.exp(2*le) - Om*a**-3 - Or*a**-4 - Ode*(a*np.exp(le))**-0.5
    return np.exp(brentq(g, -20, 60, xtol=1e-14, maxiter=500))
def ode(z):
    a = 1/(1+z); e = E(a); return Ode*(a*e)**-0.5/e**2
K = 3.0; a = 1.0; C = 3*Ode; h = 1e-4                                          # units H0 = 1, 8 pi G = 1 (M_P^2 = 1)
f = lambda K: -(2/3)*C*(a*K/3)**-0.5
d2 = (f(K+h) - 2*f(K) + f(K-h))/h**2
out += [f"L2 numerical d2f/dK2 today = {d2:.5f}; formula -Omega_DE/6 = {-Ode/6:.5f}  -> Delta lambda = Omega_DE/6",
        "   Delta lambda along the history (horizon-scale perturbations only; no extra degree of freedom):"]
for z in (0, 0.5, 1, 3, 1100, 4e8):
    out.append(f"     z = {z:>9g}: Omega_DE = {ode(z):.3e}, Delta lambda = {ode(z)/6:.3e}")
out += ["   At nucleosynthesis it is ~1e-35: the BBN bound (|lambda-1| < 0.017) is untouched. Today 0.11, but only on horizon-sized",
        "   modes (momentum constraint), not in G for galaxies or the Solar System. To do: VCDM perturbations in a Boltzmann code.", ""]
# L3: where must K be read? cell-scale vacuum ripples: h ~ l_P / R, K ~ h c / R over a region of size R
lP, H0, c = 1.616e-35, 67.5e3/3.086e22, 2.998e8
K0 = 3*H0/c
for N, cell in ((2.9e7, 0.49*np.sqrt(2.9e7)*lP), (1e15, 0.49*np.sqrt(1e15)*lP)):
    out.append(f"L3 N = {N:.1e}: cell {cell:.1e} m; ripple K at the cell scale ~ l_P/cell^2 = {lP/cell**2:.1e} /m vs today's K = {K0:.1e} /m"
               f" (ratio {lP/cell**2/K0:.0e})")
Rmin = np.sqrt(lP/K0)
out += [f"   K^(-1/2) is undefined where K fluctuates through zero -> Claim 1 must read K averaged over regions R with l_P/R^2 < K0:",
        f"   R > sqrt(l_P/K0) = {Rmin:.1e} m (~{Rmin*1e3:.2f} mm; the familiar dark-energy length, geometric mean of Planck length and",
        f"   Hubble radius -- a restatement of scales, not a new prediction).",
        "   Ripples shorter than R do not change the averaged K at linear order; their second-order effect averages incoherently.",
        "   -> loop leakage suppressed by powers of (cell/R) ~ 1e-23 or smaller: c ~ 0 for all practical purposes.", ""]
# consequence for N window
out += ["Consequence: the graviton-speed bound no longer forces N to ~1e15. Allowed window returns to N = 2.9e7 .. 5.2e15,",
        f"cell size {0.49*np.sqrt(2.9e7)*lP:.1e} .. 5.7e-28 m. The 'cells right at the gamma-ray-burst limit' near-prediction is LOST."]
txt = "\n".join(out); print(txt); open("iter48_leakage.txt", "w").write(txt + "\n")
