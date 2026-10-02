"""
The grid dark-energy law's own signature, compared with the newest published constraints (Oct 2026).
Law: rho_DE ~ adot^(-1/2)  ->  d ln rho_DE / d ln a = (rho_m/2 + rho_r - rho_DE)/rho_tot,   w = -1 - (1/3) d ln rho_DE/d ln a.
Signature: w < -1 (phantom) while the expansion decelerates, w > -1 once it accelerates; the crossing is EXACTLY at the
deceleration->acceleration moment. No free dark-energy parameter.
Outputs: w(z), crossing redshift, the CPL (w0, wa) that best mimics it over 0 < z < 2.5, and the growth of structure
(sigma8 / S8 change vs LCDM with the same early universe).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
Or = 9.1e-5
def law(Om):
    OL = 1 - Om - Or
    def rhs(lna, y):
        a = np.exp(lna); r = np.exp(y[0]); rm, rr = Om*a**-3, Or*a**-4
        return [(rm/2 + rr - r)/(rm + rr + r)]
    lna = np.linspace(0, -np.log(31), 3001)
    s = solve_ivp(rhs, [0, lna[-1]], [np.log(OL)], t_eval=lna, rtol=1e-10, atol=1e-12)
    a = np.exp(lna); rde = np.exp(s.y[0]); rm, rr = Om*a**-3, Or*a**-4
    q = (rm/2 + rr - rde)/(rm + rr + rde); w = -1 - q/3
    return a[::-1], rde[::-1], w[::-1], (rm + rr + rde)[::-1]
def growth(a, E2, Om):
    # linear growth D(a): D'' + (3/a + dlnE/da) D' = 1.5 Om a^-5 /E^2 D
    lnE = 0.5*np.log(E2); dlnE = np.gradient(lnE, a)
    f = lambda x, y: [y[1], -(3/x + np.interp(x, a, dlnE))*y[1] + 1.5*Om/(x**5*np.interp(x, a, E2))*y[0]]
    s = solve_ivp(f, [a[0], 1], [a[0], 1], rtol=1e-9, atol=1e-12); return s.y[0, -1]
print("Om    crossing z   z(accel starts)  w0       CPL fit w0, wa (0<z<2.5)   sigma8 law/LCDM")
for Om in (0.30, 0.31, 0.32):
    a, rde, w, E2 = law(Om); z = 1/a - 1
    i = np.where(np.diff(np.sign(w + 1)))[0][0]; zc = z[i]
    rm = Om*a**-3; acc = np.where(np.diff(np.sign(rm/2 + Or*a**-4*2 - rde)))[0]; za = z[acc[0]] if len(acc) else np.nan
    m = z < 2.5
    fit = least_squares(lambda p: (p[0] + p[1]*(1 - a[m])) - w[m], [-0.9, -0.5]).x
    E2L = Om*a**-3 + Or*a**-4 + (1 - Om - Or)
    print(f"{Om:.2f}  {zc:9.3f}   {za:12.3f}     {w[-1]:+.3f}   {fit[0]:+.3f}, {fit[1]:+.3f}                {growth(a, E2, Om)/growth(a, E2L, Om):.4f}")
a, rde, w, E2 = law(0.31); z = 1/a - 1
print("\nw(z) for Om = 0.31:", "  ".join(f"z={zz}: {np.interp(zz, z[::-1], w[::-1]):+.3f}" for zz in (0, 0.3, 0.5, 0.65, 1, 1.5, 2, 3)))
