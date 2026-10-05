"""The cosmic horizon as dark energy's 'container': sizes, reachable fraction, future. Constant vs fading (beta = 1/2)."""
import numpy as np, warnings; warnings.filterwarnings("ignore")
from scipy.optimize import brentq
from scipy.integrate import cumulative_trapezoid as ctz
Om, Or = 0.311, 9.1e-5; ODE = 1 - Om - Or; c = 299792.458; H0 = 67.4; tH = 977.8/H0
x = np.linspace(-12, 40, 80001)
lP = 1.616e-35; Mpc = 3.0857e22
out = ["The cosmic horizon as a container (constant vs fading dark energy)", ""]
for name, fading in (("constant", False), ("fading", True)):
    def E1(xi):
        base = Om*np.exp(-3*xi) + Or*np.exp(-4*xi)
        if not fading: return np.sqrt(base + ODE)
        return np.exp(brentq(lambda l: np.exp(2*l) - base - ODE*(np.exp(l)*np.exp(xi))**-0.5, -80, 80))
    E = np.array([E1(v) for v in x]); t = ctz(1/E, x, initial=0)*tH; t -= np.interp(0, x, t)
    integ = np.exp(-x)/E
    tail = ctz(integ[::-1], -x[::-1], initial=0)[::-1]           # comoving event horizon in c/H0 units
    past = ctz(integ, x, initial=0) + 0                          # comoving particle horizon (approx; early part negligible offset)
    i0 = np.argmin(abs(x))
    chi_eh0 = tail[i0]*c/H0/1000; chi_obs = past[i0]*c/H0/1000
    Reh_phys = np.exp(x)*tail*c/H0                               # Mpc
    S = np.pi*(Reh_phys*Mpc)**2/lP**2/np.log(2)                  # bits on the wall (one per 4 ln2 Planck areas)
    T = 1.0546e-34*2.998e8/(2*np.pi*1.381e-23*Reh_phys*Mpc)
    out.append(f"{name}: wall (event horizon) today {chi_eh0:.2f} Gpc; observable universe {chi_obs:.1f} Gpc; still-reachable share of "
               f"what we see {100*(chi_eh0/chi_obs)**3:.1f}%; wall holds {S[i0]:.1e} bits, glows at {T[i0]:.1e} K")
    for dt in (10, 50, 100, 200):
        j = np.argmin(abs(t - dt))
        out.append(f"   +{dt} Gyr: wall radius {Reh_phys[j]/1000:.2f} Gpc, share of today.s reachable galaxies still reachable {100*(tail[j]/tail[i0])**3:.2f}%, "
                   f"bits {S[j]:.1e}, glow {T[j]:.1e} K")
txt = "\n".join(out); print(txt); open("container.txt", "w").write(txt + "\n")
