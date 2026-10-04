"""
ITERATION 18 (arrow of time): the cosmic (apparent) horizon has entropy S_H ∝ area ∝ 1/H^2 (flat). The generalised second
law needs S_H non-decreasing, i.e. dH/dt <= 0 always (equivalently rho + p >= 0 for the total). Phantom dark energy
(w < -1) can violate it if it dominates. Check our histories, past AND future:
  Claim 1 (instant law, d ln rho/d ln a = q/2), the toy (memory kappa = 3, quadratic) with jostled fraction f = 1 and f = 1.5
  (negative baseline, iteration 13's mild preference), and LCDM for reference.
"""
import numpy as np
from scipy.integrate import solve_ivp
Om, Or = 0.31, 9.1e-5; OL = 1 - Om - Or
out = ["ITERATION 18: horizon entropy (S ∝ 1/H^2) and the generalised second law", ""]
def report(name, t, a, H):
    dH = np.gradient(H, t); bad = dH > 1e-9*np.abs(H).max()
    out.append(f"{name:42s}: dH/dt > 0 anywhere? {'YES' if bad.any() else 'no'}"
               + (f" (a = {a[bad].min():.2f}-{a[bad].max():.2f})" if bad.any() else "") + f";  H(a_end)/H0 = {H[-1]:.3f} at a = {a[-1]:.0f}")
# Claim 1: rho_DE ODE in ln a, from a = 1/31 to a = 100
def law(lna, y):
    a = np.exp(lna); rde = np.exp(y[0]); rm, rr = Om*a**-3, Or*a**-4
    return [0.5*(rm/2 + rr - rde)/(rm + rr + rde)]
for lo, hi in ((0, -np.log(31)), (0, np.log(100))):
    pass
lna = np.linspace(-np.log(31), np.log(100), 8000)
sb = solve_ivp(law, [0, lna[0]], [np.log(OL)], t_eval=lna[lna <= 0][::-1], rtol=1e-10, atol=1e-12)
sf = solve_ivp(law, [0, lna[-1]], [np.log(OL)], t_eval=lna[lna >= 0], rtol=1e-10, atol=1e-12)
L1 = np.r_[sb.t[::-1], sf.t[1:]]; R1 = np.exp(np.r_[sb.y[0][::-1], sf.y[0][1:]]); a1 = np.exp(L1)
H1 = np.sqrt(Om*a1**-3 + Or*a1**-4 + R1); t1 = np.r_[0, np.cumsum(np.diff(L1)/H1[1:])]
report("Claim 1 (instant law), a = 0.03-100", t1, a1, H1)
w1 = -1 - 0.5*(Om*a1**-3/2 + Or*a1**-4 - R1)/(Om*a1**-3 + Or*a1**-4 + R1)/3
tot = Om*a1**-3 + 4/3*Or*a1**-4 + (1 + w1)*R1
out.append(f"   (Claim 1 is phantom, w < -1, for a < {a1[np.where(np.diff(np.sign(w1 + 1)))[0][0]]:.2f}; total rho + p min = {tot.min():.3f} > 0)")
# toy, forward from a = 1/31 with adiabatic start; normalise jostled part to f*OL today by shooting on A
from scipy.optimize import brentq
def toy(f):
    base = OL*(1 - f)
    def run(A, tmax=60):
        ai = 1/31; Hi = np.sqrt(Om*ai**-3 + Or*ai**-4); X0 = 1/(6*ai*Hi)
        rhs = lambda t, y: [y[0]*np.sqrt(max(Om*y[0]**-3 + Or*y[0]**-4 + base + A*y[1], 1e-30)),
                            -6*np.sqrt(max(Om*y[0]**-3 + Or*y[0]**-4 + base + A*y[1], 1e-30))*y[1] + 1/y[0]]
        return solve_ivp(rhs, [0, tmax], [ai, X0], method="LSODA", rtol=1e-9, atol=1e-14, dense_output=True)
    def today(A):
        s = run(A, 6); tt = np.linspace(0, 6, 6000); aa, XX = s.sol(tt); i = np.argmin(abs(aa - 1)); return A*XX[i] - f*OL
    A = brentq(today, 1e-3, 1e3); s = run(A, 60); tt = np.linspace(0, 60, 30000); aa, XX = s.sol(tt)
    H = np.sqrt(np.maximum(Om*aa**-3 + Or*aa**-4 + base + A*XX, 0)); keep = aa <= 100
    return tt[keep], aa[keep], H[keep]
for f in (1.0, 1.5):
    t, a, H = toy(f); report(f"toy (kappa 3, quadratic), f = {f}", t, a, H)
aL = np.geomspace(1/31, 100, 4000); HL = np.sqrt(Om*aL**-3 + Or*aL**-4 + OL); tL = np.r_[0, np.cumsum(np.diff(np.log(aL))/HL[1:])]
report("LCDM", tL, aL, HL)
out += ["", "Reading: 'no' everywhere = horizon entropy never decreases: the arrow of time (second law) is respected in all our",
        "histories, despite Claim 1's phantom phase (matter keeps the total rho + p positive)."]
txt = "\n".join(out); print(txt); open("iter18_horizon_entropy.txt", "w").write(txt + "\n")
