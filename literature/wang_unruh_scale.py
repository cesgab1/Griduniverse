# Wang-Unruh: averaged expansion H ~ Lam * exp(-beta * sqrt(G) * Lam). How large must the cutoff be, and how tuned is it?
import numpy as np
EP = 1.221e28            # Planck energy, eV (sqrt(G) = 1/E_P in natural units)
H0 = 67.5e3/3.0857e22*6.582e-16   # eV
for beta in (0.5, 1, 2):
    # solve L exp(-beta L/EP) = H0  -> x = L/EP: ln x + ln EP - beta x = ln H0
    from scipy.optimize import brentq
    x = brentq(lambda x: np.log(x) + np.log(EP) - beta*x - np.log(H0), 1, 1e4)
    dlnH_dlnL = 1 - beta*x
    print(f"beta={beta}: cutoff = {x:.0f} x Planck energy;  a 1% change in the cutoff changes H by a factor {np.exp(abs(dlnH_dlnL)*0.01):.2f}"
          f" (dark energy by {np.exp(2*abs(dlnH_dlnL)*0.01):.2f})")
