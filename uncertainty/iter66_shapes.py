"""ITERATION 66: build the fixed dark-energy shapes rho(z)/rho(0) for the light tally (sqrt(eta)) and clock tally (sqrt(t)),
self-consistently with the expansion they cause (Om = 0.31, h = 0.68)."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
Om, h = 0.31, 0.68; Or = 4.15e-5/h**2; OL = 1 - Om - Or
a = np.geomspace(1e-9, 1.0, 40000); z = 1/a - 1
for name, kind in (("U3_light_sqrt_eta", "eta"), ("U4_clock_sqrt_t", "t")):
    sh = np.ones_like(a)
    for _ in range(30):
        E = np.sqrt(Om*a**-3 + Or*a**-4 + OL*sh)
        if kind == "eta": X = cumulative_trapezoid(1/(a**2*E), a, initial=0) + a[0]/(a[0]**2*E[0])*a[0]
        else: X = cumulative_trapezoid(1/(a*E), a, initial=0) + 0.5*a[0]/(a[0]*E[0])
        new = np.sqrt(X/X[-1])
        if np.max(np.abs(new - sh)) < 1e-10: break
        sh = new
    m = z <= 1e4
    zz, ss = z[m][::-1], sh[m][::-1]
    np.savetxt(f"{name}.txt", np.column_stack([zz, ss]), header=f"z rho_DE/rho_DE0 ({name}, self-consistent, Om={Om})")
    w = -1 - np.gradient(np.log(sh), np.log(a))/3
    print(f"{name}: w today = {w[-1]:.3f}; w at z=0.5 {np.interp(0.5, zz, w[m][::-1]):.3f}; rho(z=1)/rho0 = {np.interp(1, zz, ss):.3f}")
