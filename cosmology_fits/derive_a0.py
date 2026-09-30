"""
Deriving a0 from the grid's cosmic tension. Principle fixed BEFORE comparing:
  the grid is slack where the energy stored in the local pull, g^2/(8 pi G), falls below the cosmic background tension,
  i.e. the dark-energy density rho_DE c^2 (dark energy's negative pressure = tension).  ->  g* = c sqrt(8 pi G rho_DE)
Also shown, for context, other 'cosmic acceleration' scales people use (horizon/Unruh matching: c*H0; c*H_Lambda).
Measured a0: 1.15e-10 (our 164-galaxy fit), 0.95e-10 (blind half-split), 1.2e-10 (standard literature value).
Then the distinctive prediction: if a0 is set by the cosmic tension, it changes with cosmic time as the tension does.
"""
import numpy as np
from scipy.integrate import solve_ivp
G = 6.674e-11; c = 2.998e8; Mpc = 3.0857e22
for lab, H0, OL in (("dark-energy law fit (H0 67.2, Om 0.31)", 67.2, 0.69), ("Planck LCDM (H0 67.4)", 67.4, 0.685), ("SH0ES H0 73", 73.0, 0.70)):
    H = H0*1e3/Mpc; rho_c = 3*H**2/(8*np.pi*G); rho_de = OL*rho_c
    g_star = c*np.sqrt(8*np.pi*G*rho_de)
    print(f"{lab}:  derived g* = c*sqrt(8 pi G rho_DE) = {g_star:.2e} m/s^2   (= sqrt(3 Omega_L) c H0 = {np.sqrt(3*OL)*c*H:.2e})")
    print(f"      other scales: c*H0 = {c*H:.2e}, c*H0*sqrt(Omega_L) = {c*H*np.sqrt(OL):.2e}")
    for a0, al in ((1.15e-10, "our fit"), (0.95e-10, "blind"), (1.2e-10, "literature")):
        print(f"      measured a0 ({al}) {a0:.2e}: derived / measured = {g_star/a0:.1f}")
print("\nSo the principle lands at the right order of magnitude (10^-10 out of a possible range spanning ~100 powers of ten)")
print("but about 8x too high. The famous empirical match a0 ~ c H0 / 2 pi needs a factor 1/2pi that the principle does not supply.")
# how the transition point relates to a0 in the interpolation used on the galaxies
nu = lambda y: 1/(1 - np.exp(-np.sqrt(y)))
for boost in (1.1, 1.5, 2.0):
    ys = np.geomspace(1e-3, 1e3, 200000); y = ys[np.argmin(abs(nu(ys) - boost))]
    print(f"   in the fitted law, the extra pull reaches {boost:.1f}x at g_N = {y:.2f} a0")
# distinctive prediction: a0(z) tracks the cosmic tension
print("\nPrediction: a0 changes with cosmic time as the background tension does (baryonic Tully-Fisher: V^4 = G M a0)")
Om = 0.31; orad = 9.1e-5
def rhoDE_law(z):
    f = lambda x, y: [0.5*(Om*np.exp(-3*x)/2 + orad*np.exp(-4*x) - np.exp(y[0]))/(Om*np.exp(-3*x) + orad*np.exp(-4*x) + np.exp(y[0]))]
    x_end = -np.log(1 + z); s = solve_ivp(f, [0, x_end], [np.log(1 - Om)], rtol=1e-9); return np.exp(s.y[0][-1])/(1 - Om)
for z in (0.5, 1.0, 2.0):
    E = np.sqrt(Om*(1+z)**3 + (1-Om))
    r_law = np.sqrt(rhoDE_law(z)); r_H = E
    print(f"   z = {z}: a0 set by dark-energy tension (our law): x{r_law:.3f} -> rotation speed at fixed mass x{r_law**0.25:.3f} ({np.log10(r_law**0.25):+.3f} dex);"
          f"  if set by c*H(z) instead: x{r_H:.2f} -> x{r_H**0.25:.3f} ({np.log10(r_H**0.25):+.3f} dex); Lambda: no change")
