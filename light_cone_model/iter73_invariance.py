"""Iteration 73: is matter x horizon^2 at the crossing independent of D0? See PREREG_73.md."""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import cumulative_trapezoid as ctz
att = 3/(8*np.pi)*25/16
def run(ODE, Ok, fading=True, Or=0.0):
    Om = 1 - ODE - Ok - Or
    x = np.linspace(-8, 45, 40001)
    def E1(xi):
        z1 = np.exp(-xi); base = Om*z1**3 + Or*z1**4 + Ok*z1**2
        if not fading: return np.sqrt(base + ODE)
        return np.exp(brentq(lambda l: np.exp(2*l) - base - ODE*(np.exp(l)*np.exp(xi))**-0.5, -80, 80))
    E = np.array([E1(xi) for xi in x])
    rDE = ODE*(E*np.exp(x))**-0.5 if fading else ODE + 0*x
    rM = Om*np.exp(-3*x)
    tail = ctz((np.exp(-x)/E)[::-1], -x[::-1], initial=0)[::-1]
    Reh = np.exp(x)*tail                                  # in c/H0
    ic = np.where(np.diff(np.sign(rDE - rM)))[0][0]
    return 3/(8*np.pi)*rM[ic]*Reh[ic]**2, np.exp(-x[ic]) - 1
out = ["Iteration 73: matter x event-horizon^2 at the crossing, for different dial settings (Omega_DE today)", ""]
for fading in (True, False):
    for Ok in (0.0, 0.0023):
        row = []
        for ODE in (0.3, 0.5, 0.685, 0.8, 0.9):
            kc, zc = run(ODE, Ok, fading)
            row.append(f"{ODE}: {kc/att:.3f}x (z_cross {zc:.2f})")
        out.append(f"{'fading' if fading else 'Lambda'}, Omega_k = {Ok}: " + " | ".join(row))
txt = "\n".join(out); print(txt); open("iter73_invariance.txt", "w").write(txt + "\n")
