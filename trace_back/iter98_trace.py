import numpy as np
from scipy.optimize import brentq
G, MSUN, KPC, GYR = 6.674e-11, 1.989e30, 3.086e19, 3.156e16
t0 = 13.8*GYR
def timing(r_kpc, v_kms):
    r, v = r_kpc*KPC, v_kms*1e3          # v < 0 approaching
    # r/t and v: v t / r = sin(eta)(eta - sin eta)/(1 - cos eta)^2
    f = lambda e: np.sin(e)*(e - np.sin(e))/(1 - np.cos(e))**2 - v*t0/r
    eta = brentq(f, np.pi + 1e-6, 2*np.pi - 1e-6)
    B = t0/(eta - np.sin(eta)); A = r/(1 - np.cos(eta)); M = A**3/(G*B**2)
    t_turn = B*np.pi/GYR; r_max = 2*A/KPC
    return eta, M/MSUN, t_turn, r_max
out = []
e, M, tt, rm = timing(770, -110)
out += [f"Milky Way + Andromeda: started together moving apart, max separation {rm:.0f} kpc at {tt:.1f} Gyr, now falling in",
        f"   mass needed {M:.2e} Msun vs visible ~1.5e11 -> x{M/1.5e11:.0f}"]
# Virgo zero-velocity surface: galaxies at R0 have v = 0 now -> eta = pi, turnaround now
R0 = 7.2e3*KPC; M = np.pi**2*R0**3/(8*G*t0**2)/MSUN
out += [f"Virgo: galaxies at 7.2 Mpc are turning around today (eta = pi); everything inside has already started falling in",
        f"   mass needed {M:.2e} Msun (no Lambda; Kashibadze+2020 with Lambda: 7.4e14) vs visible ~1e14 -> x{M/1e14:.1f}"]
print("\n".join(out)); open("iter98_trace.txt", "w").write("\n".join(out) + "\n")
