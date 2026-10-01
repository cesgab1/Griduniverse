# Total early-universe amplification of AeST's vector-sector mode: N = Int (growth rate / aH) dln a, from the full 4x4 block
import numpy as np
exec(open("fast_modes.py").read().split('print(f"check')[0])
lna = np.linspace(np.log(1e-8), 0, 1500)
for KB in (0.5, 0.2, 0.1, 1e-3, 2.5e-5):
    g = np.array([max(J(np.exp(x), 0.1, KB)[0], 0) for x in lna])
    N = np.trapezoid(g, lna)
    print(f"K_B = {KB:<8g} amplification e^{N:7.1f}   (analytic early-time estimate rate ~ sqrt(3 Omega_fld / K_B) aH)")
