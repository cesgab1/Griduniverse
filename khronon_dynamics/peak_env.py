"""
(1) Realistic surroundings: mean linear density around a galaxy-forming peak (Gaussian field, conditional mean):
      Delta(R) = delta_c/D(z_c) * sigma^2(R, R_L) / sigma^2(R_L),  sigma^2(R, R_L) = int P(k) W(kR) W(kR_L) k^2 dk / 2 pi^2
    for a 1.5e12 Msun region (R_L) collapsing at z_c = 1 (a ~1.5-sigma peak). This replaces the crude R^-1.2 profile.
(2) MOND halo target with the external field effect (1-D QUMOND estimate):
      g = nu(|gN + gNe|/a0)(gN + gNe) - nu(gNe/a0) gNe,  gNe chosen so that nu(gNe/a0) gNe = g_e (external field g_e)
    phantom(<r) = g r^2/G - Mb, for g_e = 0, 0.01, 0.025, 0.05 a0.
Then runs crossover.galaxy with these initial conditions: CDM control and growing crossovers.
"""
import numpy as np, sys, os
from scipy.integrate import quad
from scipy.optimize import brentq
import crossover as C
h, Om = C.h, C.Om
def T(k): q = k/(Om*h**2); return np.log(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**-0.25
W = lambda x: 3*(np.sin(x)-x*np.cos(x))/x**3
ks = np.logspace(-4, 2.5, 4000)                                          # 1/Mpc comoving
Pk = ks**0.965*T(ks)**2
s8n = np.trapezoid(Pk*W(ks*8/h)**2*ks**2, ks)/(2*np.pi**2); Pk *= 0.81**2/s8n
def sig2(R1, R2): return np.trapezoid(Pk*W(ks*R1)*W(ks*R2)*ks**2, ks)/(2*np.pi**2)
def Dgr(z):
    g = lambda a: quad(lambda x: 1/(x*np.sqrt(Om/x**3 + 1 - Om))**3, 1e-6, a)[0]*2.5*Om*np.sqrt(Om/a**3 + 1 - Om)
    return g(1/(1+z))/g(1.0)
Mh, zc = 1.5e12, 1.0
RL = (3*Mh/(4*np.pi*Om*C.rhoc))**(1/3)/1e3                               # Mpc comoving
D30, Dc = Dgr(30), Dgr(zc)
def Delta_peak(Rkpc):
    R = np.atleast_1d(Rkpc)/1e3
    return np.array([1.686/Dc*sig2(x, RL)/sig2(RL, RL) for x in R])*D30
nu_pk = 1.686/Dc/np.sqrt(sig2(RL, RL))
if __name__ == "__main__" and os.environ.get("TARGET"):
    Mb = 1e11; G = C.G; a0 = C.a0
    print(f"peak height nu = {nu_pk:.2f}; R_L = {RL:.2f} Mpc; Delta(z=30) at R = 0.5, 1, 2, 4, 8 R_L: " +
          ", ".join(f"{Delta_peak(f*RL*1e3)[0]:.4f}" for f in (0.5, 1, 2, 4, 8)))
    print("MOND halo target for Mb = 1e11 (phantom mass within r):")
    for ge in (0.0, 0.01, 0.025, 0.05):
        gNe = brentq(lambda x: C.nu(x/a0)*x - ge*a0, 1e-12, 10*a0) if ge else 0.0
        row = []
        for r in (30, 100, 300, 1000):
            gN = G*Mb/r**2; g = C.nu((gN + gNe)/a0)*(gN + gNe) - (C.nu(gNe/a0)*gNe if ge else 0)
            row.append(g*r**2/G - Mb)
        print(f"   g_e = {ge:5.3f} a0: " + "  ".join(f"{r}: {x:.2e}" for r, x in zip((30, 100, 300, 1000), row)))
    raise SystemExit
if __name__ == "__main__":
    l3, l0 = float(sys.argv[1]), float(sys.argv[2])
    import types
    def galaxy_peak(l3, l0, Mb=1e11, N=300, zi=30.0):
        R = np.geomspace(20, 12000, N); edges = np.r_[0, R]
        m_f = C.Oc*C.rhoc*4/3*np.pi*(edges[1:]**3 - edges[:-1]**3)
        Delta = Delta_peak(R)
        ai = 1/(1+zi); r = ai*R*(1 - Delta/3); v = C.Hof(ai)*r*(1 - Delta/3)
        return C.run_shells(R, m_f, r, v, l3, l0, Mb, zi)
    res = galaxy_peak(l3, l0)
    for lab in ("z=0.9", "z=0"):
        rs, Mf, av = res[lab]
        print(f"lam {l3:g}->{l0:g} {lab}: " + "  ".join(f"{R_}: {np.interp(R_, rs, Mf) - 4/3*np.pi*R_**3*C.Oc*C.rhoc/av**3:.2e}" for R_ in (30, 100, 300, 1000)))
