import numpy as np
from scipy.integrate import cumulative_trapezoid as ctz
Om = 0.31
z = np.concatenate([np.linspace(0, 3, 601)[:-1], np.geomspace(3, 1e4, 1500)])
H = np.sqrt(Om*(1+z)**3 + 9e-5*(1+z)**4 + 1 - Om)
psi = 0.015*(1+z)**2.7/(1 + ((1+z)/2.9)**5.6)
dtdz = 1/((1+z)*H)
def cum_from_inf(f):            # integral from z to infinity
    c = ctz(f[::-1], -z[::-1], initial=0)[::-1]
    return c
s1 = cum_from_inf(psi*dtdz); s2 = cum_from_inf(psi*(1+z)**3*dtdz)
for name, s in (("S1_stars", s1), ("S2_coupled_bh", s2)):
    sh = s/s[0]
    np.savetxt(f"{name}.txt", np.column_stack([z, sh]), header=f"z rho_DE(z)/rho_DE(0) {name}")
    print(name, " ".join(f"z={zz}:{np.interp(zz, z, sh):.3f}" for zz in (0.5, 1, 2, 3, 5)))
