"""
ITERATION 55 (pre-registered in PREREG_55.md): numbers for the two-level checks.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
H0 = 67.4e3/3.0857e22; Om = 0.315; Or = 9.1e-5; OL = 1 - Om - Or; Gyr = 3.156e16
E = lambda a: np.sqrt(Om/a**3 + Or/a**4 + OL)
out = ["ITERATION 55: two-level checks (expectations pre-registered)", ""]
# Route Q, GR counterpart: (1/4) <rho_m>_4-volume over a book that ends when the a^3-weighted time tally reaches T x today's
I = lambda a1: quad(lambda a: a**2/(H0*E(a)), 1e-10, a1, limit=200)[0]       # integral a^3 dt
t = lambda a1: quad(lambda a: 1/(a*H0*E(a)), 1e-10, a1, limit=200)[0]
past = I(1.0)
out.append(" T (book time tally / so far)   ends in [Gyr]   (1/4)<rho_m> / rho_DE(measured)")
for T in (1.0, 1.5, 2.0, 3.0, 5.0, 10.0):
    af = brentq(lambda a1: I(a1)/past - T, 0.5, 50) if T > 1 else 1.0
    avg_rho_m = Om*t(af)/I(af)                         # rho_m a^3 = Om (critical units, a0 = 1): <rho_m> = Om * t_end / INT a^3 dt
    out.append(f"   {T:5.1f}                         {(t(af) - t(1.0))/Gyr:6.1f}          {0.25*avg_rho_m/OL:7.3f}")
out += ["   (space factor S cancels for a homogeneous universe; sign of the sequestered term not checked here)", ""]
# Unification: dark energy absent before the handover (z_s = 124)
def E1(zz):
    a = 1/(1 + zz); g = lambda le: np.exp(2*le) - Om*a**-3 - Or*a**-4 - OL*(a*np.exp(le))**-0.5
    return np.exp(brentq(g, -10, 40))
for zz in (124, 300, 1090):
    e = E1(zz); share = OL*((1/(1+zz))*e)**-0.5/e**2
    out.append(f"dark-energy share of the total at z = {zz:5d} (Claim 1): {share:.1e}  -> switching it off before z = 124 changes H by < {share/2:.0e}")
txt = "\n".join(out); print(txt); open("iter55_two_levels.txt", "w").write(txt + "\n")
