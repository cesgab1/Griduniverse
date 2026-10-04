import numpy as np
from scipy.optimize import brentq
from scipy.integrate import cumulative_trapezoid as ctz
D = 0.685*np.sqrt(1/np.sqrt(0.0023))      # D0 in units of today's critical density = 3.13
tH = 977.8/67.4                           # 1/H0 in Gyr
A = np.exp(np.linspace(np.log(1e-4), np.log(50), 200001))   # curvature radius in c/H0
def E(a):  # H in H0 units: E^2 = 1/a^2 + D (a E)^(-1/2)
    return np.exp(brentq(lambda l: np.exp(2*l) - 1/a**2 - D*(a*np.exp(l))**-0.5, -50, 50))
Ev = np.array([E(a) for a in A])
t = ctz(1/(A*Ev), A, initial=0)*tH + A[0]*tH          # early: a = c t
de = D*(A*Ev)**-0.5; share = de/Ev**2
q = -1 - np.gradient(np.log(Ev), np.log(A))
for s in (0.1, 0.5, 0.9, 0.99):
    i = np.argmax(share >= s); print(f"dark energy share {s:.0%} at t = {t[i]:.2f} Gyr")
i = np.argmax(q < 0); print(f"acceleration starts at t = {t[i]:.2f} Gyr; share then {share[i]:.2f}")
print("1/sqrt(8piG D0/3) =", tH/np.sqrt(D), "Gyr")
