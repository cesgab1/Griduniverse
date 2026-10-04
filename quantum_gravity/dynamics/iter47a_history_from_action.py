"""
ITERATION 47a (pre-registered in PREREG_47_49.md): the whole cosmic history found by making the grid's own action stationary.
No Friedmann equation is given to the computer. Units: H0 = 1, 8 pi G = 1 (critical density today = 3).
"""
import numpy as np
from scipy.optimize import brentq
from scipy.spatial import cKDTree
Om, Or = 0.315, 9.1e-5; Ode = 1 - Om - Or
rm, rr, C = 3*Om, 3*Or, 3*Ode
rng = np.random.default_rng(47)
# --- grid reader (iteration 45): K_ij fitted from link-length changes on random points, for a uniform stretch rate r = 1
X = rng.uniform(0, 12, (1728, 3)); tree = cKDTree(X, boxsize=12); pr = tree.query_pairs(2.0, output_type="ndarray")
D = X[pr[:, 1]] - X[pr[:, 0]]; D -= 12*np.round(D/12); u = D/np.linalg.norm(D, axis=1)[:, None]
rows = np.stack([u[:, 0]**2, u[:, 1]**2, u[:, 2]**2, 2*u[:, 0]*u[:, 1], 2*u[:, 0]*u[:, 2], 2*u[:, 1]*u[:, 2]], 1)
k6 = np.linalg.lstsq(rows, np.ones(len(rows)), rcond=None)[0]           # rate of every link = 1 x r
K = np.array([[k6[0], k6[3], k6[4]], [k6[3], k6[1], k6[5]], [k6[4], k6[5], k6[2]]])
T = np.trace(K); Q = np.sum(K*K) - T**2                                  # grid's K and K_ijK^ij - K^2 per unit rate^2
def step_action(N, l0, l1):
    ab = np.exp(0.5*(l0 + l1)); r = (l1 - l0)/N                          # stretch rate read on the grid = r
    return N*ab**3*(0.5*Q*r**2 - rm*ab**-3 - rr*ab**-4 - (2/3)*C*(ab*T*r/3)**-0.5)
def dS_dN(N, l0, l1, h=1e-7): return (step_action(N*(1+h), l0, l1) - step_action(N*(1-h), l0, l1))/(2*N*h)
def solve(M, a0=1/1101, a1=1.0):
    la = np.linspace(np.log(a0), np.log(a1), M + 1); N = np.empty(M)
    for n in range(M):
        d = la[n+1] - la[n]
        N[n] = brentq(dS_dN, d/1e6, d*1e3, args=(la[n], la[n+1]), xtol=1e-15, rtol=1e-14)   # computer finds the step
    return la, N
def E_exact(a):
    return brentq(lambda E: E**2 - Om*a**-3 - Or*a**-4 - Ode*(a*E)**-0.5, 1e-6, 1e12)
out = ["ITERATION 47a: history from making the grid action stationary (expectations pre-registered)", "",
       f"grid reader: K = {T:.12f} x rate (exact 3), K_ijK^ij - K^2 = {Q:.12f} x rate^2 (exact -6)", ""]
for M in (250, 500, 1000):
    la, N = solve(M); lm = 0.5*(la[1:] + la[:-1]); H = np.diff(la)/N
    Eex = np.array([E_exact(np.exp(l)) for l in lm]); err = np.max(np.abs(H/Eex - 1))
    # dark energy on the found history, its w, crossing and today's value
    rde = C*(np.exp(lm)*H)**-0.5; w = -1 - np.gradient(np.log(rde), lm)/3; z = np.exp(-lm) - 1
    zc = np.interp(0, (w + 1)[::-1], z[::-1]) if np.any(w > -1) else np.nan
    i = np.where((w[:-1] + 1)*(w[1:] + 1) < 0)[0]; zc = z[i[-1]] + (z[i[-1]+1] - z[i[-1]])*(-(w[i[-1]]+1))/(w[i[-1]+1] - w[i[-1]])
    # the a-equations (NOT imposed): dS/d(ln a_n) at interior points, relative to the size of the separate terms
    res = []
    for n in range(1, M, max(1, M//50)):
        h = 1e-6; f = lambda l: step_action(N[n-1], la[n-1], l) + step_action(N[n], l, la[n+1])
        g = lambda l: abs(step_action(N[n-1], la[n-1], l)) + abs(step_action(N[n], l, la[n+1]))
        res.append(abs(f(la[n] + h) - f(la[n] - h))/(2*h)/g(la[n]))
    t_age = np.sum(N)                                                     # in units of 1/H0
    out.append(f"steps {M:5d}: max |H/H_Friedmann - 1| = {err:.2e}; crossing z = {zc:.3f}; w0 = {w[-1]:.4f}; "
               f"time since z=1100 = {t_age:.4f}/H0; unimposed a-equation residual (median) = {np.median(res):.2e}")
txt = "\n".join(out); print(txt); open("iter47a_history_from_action.txt", "w").write(txt + "\n")
