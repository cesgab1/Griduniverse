"""
Fate of the universe under the grid dark-energy law rho_DE ~ adot^-beta (fit: beta = 1/2; free fit 0.63 +/- 0.22).
Friedmann: E^2 = Om a^-3 + Or a^-4 + Ok a^-2 + Ode (a E)^-beta,  E = H/H0, normalised so E(1) = 1.
Analytic late time (matter gone): E^(2+beta) ~ a^-beta  ->  a ~ t^((2+beta)/beta)  (t^5 for beta = 1/2): accelerates forever.
Recollapse would need adot -> 0 first, but rho_DE ~ adot^-beta -> infinity there: the slower the grid stretches, the harder the
tension pushes. So no turnaround for any beta > 0, even in a closed universe (curvature ~ a^-2 always loses to a^-beta/(2+beta)*...).
Also computed: CMB temperature, horizon temperature T_h = hbar H / (2 pi k_B), and when structure growth freezes; vs LCDM.
"""
import numpy as np
np.seterr(over="ignore")
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
Om, Or, H0 = 0.31, 9.1e-5, 68/3.086e19                        # H0 in 1/s
yr = 3.156e7; hbar, kB = 1.0546e-34, 1.3807e-23
def Efun(a, beta, Ok=0.0):
    Ode = 1 - Om - Or - Ok
    if beta == 0: return np.sqrt(Om/a**3 + Or/a**4 + Ok/a**2 + Ode)
    f = lambda E: E**2 - Om/a**3 - Or/a**4 - Ok/a**2 - Ode*(a*E)**-beta
    return brentq(f, 1e-12, 1e8)
def history(beta, Ok=0.0, lnamax=40):
    # integrate t(a) and linear growth D(a) (dD/dlna = f D, matter only drives growth)
    def rhs(lna, y):
        a = np.exp(lna); E = Efun(a, beta, Ok); t, D, Dp = y
        dl = 1e-4; dlnE = (np.log(Efun(a*np.exp(dl), beta, Ok)) - np.log(Efun(a*np.exp(-dl), beta, Ok)))/(2*dl)
        Omz = Om/a**3/E**2
        return [1/(H0*E), Dp, -(2 + dlnE)*Dp + 1.5*Omz*D]
    a0 = 1e-3; t0 = 2/(3*H0*np.sqrt(Om))*a0**1.5
    s = solve_ivp(rhs, [np.log(a0), lnamax], [t0, a0, a0], dense_output=True, rtol=1e-8, atol=1e-30, max_step=0.2)
    return s
for name, beta in [("LCDM (Lambda)", 0.0), ("grid law beta = 0.5", 0.5), ("beta = 0.85 (fit +1 sigma)", 0.85), ("beta = 0.41 (fit -1 sigma)", 0.41)]:
    s = history(beta)
    lna = np.linspace(-6.9, 40, 6000); t, D, _ = s.sol(lna); a = np.exp(lna)
    i0 = np.argmin(abs(lna)); age = t[i0]/yr; Dinf = D[-1]
    tf = t[np.argmax(D > 0.99*Dinf)]/yr
    out = f"{name}: age today {age/1e9:.2f} Gyr; growth 99% frozen at {tf/1e9:.0f} Gyr"
    if beta: out += f"; late-time a ~ t^{(2+beta)/beta:.2f}"
    print(out)
    for T_yr in (1e11, 1e12, 1e14, 1e20, 1e40):
        j = np.argmax(t/yr > T_yr)
        if t[-1]/yr < T_yr:                                   # extrapolate with the asymptotic law
            p = (2+beta)/beta if beta else None
            if beta: aa = a[-1]*(T_yr*yr/t[-1])**p; HH = p/(T_yr*yr)
            else: HH = H0*np.sqrt(1-Om-Or); aa = a[-1]*np.exp(HH*(T_yr*yr - t[-1]))
        else: aa = a[j]; HH = H0*Efun(aa, beta)
        Th = hbar*HH/(2*np.pi*kB)
        Tc = 2.7255/aa if np.isfinite(aa) else 0.0
        print(f"   t = {T_yr:.0e} yr: size x {aa:.2e}, CMB {Tc:.1e} K, horizon temperature {Th:.1e} K")
# Closed universe check: does a curved universe ever turn around?
for Ok in (-0.01, -0.05):
    s = history(0.5, Ok, lnamax=20); lna = np.linspace(-6.9, 20, 2000)
    print(f"closed universe Omega_k = {Ok}: keeps expanding to size x {np.exp(lna[-1]):.1e} (t = {s.sol(lna[-1])[0]/yr:.1e} yr) -- no turnaround")
