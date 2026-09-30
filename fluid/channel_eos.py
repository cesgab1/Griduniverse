"""
Step 1: does fluid confined to narrow channels (1D) give P ∝ n^3, the superfluid equation of state that yields MOND?
Exact 1D Bose gas with contact repulsion (Lieb-Liniger, Bethe ansatz). Coupling gamma = m g /(hbar^2 n):
dilute or strongly repelling -> gamma large.  Energy per length E/L = hbar^2 n^3 e(gamma) / 2m.
Solve  G(x) - (1/2pi) Int_-1^1 2 lam G(y)/(lam^2+(x-y)^2) dy = 1/2pi ;  gamma = lam / Int G ;  e = gamma^3/lam^3 Int x^2 G
Then P = n mu - E/L and the effective exponent Gamma = dlnP/dln n at fixed g (MOND superfluid needs Gamma = 3).
Also coarse-grain: a random 3D network of such channels with fixed channel density has P_3D ∝ n_3D^Gamma (same exponent).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
x, w = leggauss(1200)
def solve(lam):
    K = (1/(2*np.pi))*2*lam/(lam**2 + (x[:, None]-x[None, :])**2)*w[None, :]
    G = np.linalg.solve(np.eye(len(x)) - K, np.full(len(x), 1/(2*np.pi)))
    gam = lam/np.sum(w*G); e = gam**3/lam**3*np.sum(w*x**2*G)
    return gam, e
lams = np.logspace(-0.7, 3.5, 250)
ge = np.array([solve(l) for l in lams]); gam, e = ge[:, 0], ge[:, 1]
# at fixed g: n ∝ 1/gamma. E/L ∝ n^3 e(gamma). Use u = ln n = -ln gamma (+const)
lnn = -np.log(gam); lnE = 3*lnn + np.log(e)
E = np.exp(lnE); n = np.exp(lnn)
mu = np.gradient(E, n); P = n*mu - E
Gam = np.gradient(np.log(P), lnn)
print(" gamma      e(gamma)    e/(pi^2/3)   pressure exponent Gamma (MOND superfluid: 3; ordinary BEC: 2)")
for gv in (0.1, 0.3, 1, 3, 10, 30, 100, 1000):
    i = np.argmin(abs(np.log(gam)-np.log(gv)))
    print(f" {gam[i]:8.3f}   {e[i]:9.4f}    {e[i]/(np.pi**2/3):7.4f}     {Gam[i]:.3f}")
# where is Gamma within 5% / 1% of 3?
for tol in (0.05, 0.01):
    ok = gam[np.abs(Gam-3)/3 < tol]
    print(f"Gamma within {tol*100:.0f}% of 3 for gamma > {ok.min():.1f}")
