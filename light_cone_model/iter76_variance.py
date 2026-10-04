"""Iteration 76: dark-energy-free universe and the variance test. See PREREG_76.md."""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import cumulative_trapezoid as ctz
Ok, Or0, ODE = 0.0023, 9.1e-5, 0.685
Om0 = 1 - ODE - Ok - Or0
tH = 977.8/67.4                                    # Gyr
A = np.exp(np.linspace(np.log(1e-9), np.log(3e3), 120001))
def history(sM=1.0, sD=1.0, dark=True):
    Om = Om0*sM; D = ODE*sD
    def E(a):
        base = Om/a**3 + Or0/a**4 + Ok/a**2
        if not dark: return np.sqrt(base)
        return np.exp(brentq(lambda l: np.exp(2*l) - base - D*(a*np.exp(l))**-0.5, -60, 80))
    Ev = np.array([E(a) for a in A])
    t = (ctz(1/(A*Ev), A, initial=0) + A[0]/(2*A[0]*Ev[0]))*tH
    rde = D*(A*Ev)**-0.5 if dark else 0*A
    q = -1 - np.gradient(np.log(Ev), np.log(A))
    return dict(t=t, E=Ev, rm=Om/A**3, rk=Ok/A**2, rde=rde, q=q)
out = ["Iteration 76: dark-energy-free universe and variance test. PREREG_76.md", ""]
h = history(dark=False)
sk = h["rk"]/h["E"]**2
out.append("A. Dark-energy-free (same matter, radiation, curvature):")
i13 = np.argmin(abs(h["t"] - 13.8))
out.append(f"   at t = 13.8 Gyr: size = {A[i13]:.2f} x today's, curvature share {sk[i13]:.4f}, decelerating (q = {h['q'][i13]:.2f})")
for s in (0.1, 0.5, 0.9):
    i = np.argmax(sk >= s); out.append(f"   empty light cone (curvature) reaches {s:.0%} at t = {h['t'][i]:.0f} Gyr (size {A[i]:.0f} x today's)")
out.append(f"   ever accelerates? {'YES' if (h['q'] < -1e-6).any() else 'NO'}")
G_SI, rhoc0 = 6.674e-11, 3*(67.4e3/3.0857e22)**2/(8*np.pi*6.674e-11)
for tt in (1, 13.8):
    i = np.argmin(abs(h['t']-tt)); rho = h['rm'][i]*rhoc0; pred = 1/(6*np.pi*G_SI*(tt*3.156e16)**2)
    out.append(f"   matter density at {tt} Gyr: {rho:.2e} kg/m^3 vs 1/(6 pi G t^2) = {pred:.2e} (ratio {rho/pred:.3f})")
out.append("")
out.append("B. Variance (full model):")
def marks(hh):
    ia = np.argmax(hh["q"] < 0); ic = np.argmax(hh["rde"] > hh["rm"]); return hh["t"][ia], hh["t"][ic]
base_on, base_cr = marks(history())
for sM in (0.25, 0.5, 1, 2, 4):
    on, cr = marks(history(sM=sM))
    out.append(f"   matter x {sM:<4}: acceleration starts {on:5.2f} Gyr ({on/base_on:.3f}x) | dark energy passes matter {cr:5.2f} Gyr ({cr/base_cr:.3f}x)")
for sD in (0.5, 2):
    on, cr = marks(history(sD=sD))
    out.append(f"   dial   x {sD:<4}: acceleration starts {on:5.2f} Gyr ({on/base_on:.3f}x; D0^-1/2 would be {sD**-0.5:.3f}x) | passes matter {cr:5.2f} Gyr")
D0 = ODE/np.sqrt(np.sqrt(Ok))  # = 3.13 critical
tD = tH/np.sqrt(D0)
out.append("")
out.append(f"C. Dial clock 1/sqrt(8 pi G D0/3) = {tD:.2f} Gyr; acceleration starts at {base_on:.2f} Gyr = {base_on/tD:.3f} x dial clock.")
txt = "\n".join(out); print(txt); open("iter76_variance.txt", "w").write(txt + "\n")
