"""ITERATION 71 (pre-registered in PREREG_71.md)."""
import numpy as np
from scipy.integrate import solve_ivp, quad
G, c, Gpc = 6.674e-11, 2.998e8, 3.0857e25
H0 = 67.4e3/3.0857e22; OL = 0.685
k = np.sqrt(4*np.pi/3)*16/ (np.sqrt(4*np.pi/3)*16) * 16.37    # coefficient c^4/(k G R^2): 2 x 4 x sqrt(4 pi/3) = 16.37
k = 2*4*np.sqrt(4*np.pi/3)
ch = np.sqrt(8*np.pi/(3*k))
Om = 0.315
E = lambda a: np.sqrt(Om/a**3 + (1 - Om))
Reh = c/H0*quad(lambda x: 1/(np.exp(x)*E(np.exp(x))), 0, 60, limit=400)[0]
rho_pred = c**2/(k*G*Reh**2)/c**2*c**2                         # kg/m^3 equivalent: c^4/(k G R^2) / c^2
rho_pred = c**2/(k*G*Reh**2)
rho_meas = OL*3*H0**2/(8*np.pi*G)
out = ["ITERATION 71: dark-energy horizon treated like a black-hole horizon (expectations pre-registered)", "",
       f"coefficient k = {k:.2f}; holographic c_h = {ch:.3f} (derived)",
       f"H1: with the LCDM event horizon R = {Reh/Gpc:.2f} Gpc: rho_DE = {rho_pred:.2e} kg/m^3 = {rho_pred/rho_meas:.2f} x measured"]
# self-consistent holographic dark energy: dOmega/dlna = Omega (1 - Omega)(1 + 2 sqrt(Omega)/c_h); w = -1/3 - 2 sqrt(Omega)/(3 c_h)
def shape(Ode0, cc=ch):
    sol = solve_ivp(lambda x, O: O*(1 - O)*(1 + 2*np.sqrt(np.clip(O, 0, 1))/cc), [0, -np.log(1e4)], [Ode0], dense_output=True, rtol=1e-10, atol=1e-12)
    z = np.concatenate([np.linspace(0, 3, 600), np.geomspace(3.01, 1e4, 600)]); x = -np.log(1 + z); O = sol.sol(x)[0]
    H2 = (1 - Ode0)*(1 + z)**3/(1 - O)                             # matter only (radiation negligible for the shape)
    return z, H2*O/Ode0, -1/3 - 2*np.sqrt(O)/(3*cc)
z, sh, w = shape(OL)
np.savetxt("iter71_hde_shape.txt", np.column_stack([z, sh]), header="z rho_DE/rho_DE0 holographic c_h=0.716, Ode0=0.685")
out.append(f"H3 shape: w today = {w[0]:.3f}, at z = 1: {np.interp(1, z, w):.3f} (never above -1: {np.all(w < -1)})")
out.append("H2: self-consistent solutions exist for ANY starting share Omega_DE0 -> the size is not pinned by this route")
txt = "\n".join(out); print(txt); open("iter71_horizon.txt", "w").write(txt + "\n")
