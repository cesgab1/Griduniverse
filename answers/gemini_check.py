"""
Check of an outside (Gemini) reading of our equations. Each item recomputed from our own equations.
Planck-like background: Om = 0.315, Or = 9.1e-5, flat.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
Om, Or = 0.315, 9.1e-5; Ode = 1 - Om - Or
out = []
# --- Claim 1: rho_DE = C (aH)^(-1/2), solved self-consistently: E^2 = Om a^-3 + Or a^-4 + Ode (a E)^(-1/2)
def E_of_a(a):
    f = lambda E: E**2 - Om*a**-3 - Or*a**-4 - Ode*(a*E)**-0.5
    return brentq(f, 1e-6, 1e12)
def w_claim1(a, h=1e-5):
    lnrho = lambda la: np.log(Ode*(np.exp(la)*E_of_a(np.exp(la)))**-0.5)
    return -1 - (lnrho(np.log(a)+h) - lnrho(np.log(a)-h))/(2*h)/3
zc = brentq(lambda z: w_claim1(1/(1+z)) + 1, 0.1, 2)
out += ["CLAIM 1  rho_DE = C (aH)^(-1/2)",
        f"  w today {w_claim1(1.0):.3f}; crossing w = -1 at z = {zc:.2f}; matter era (z=50) w = {w_claim1(1/51):.3f} (formula -13/12 = {-13/12:.3f})",
        f"  far future (a = 1e6): w = {w_claim1(1e6):.4f} (attractor a ~ t^5: -1 + 4/30 = {-1+4/30:.4f}); de Sitter (w = -5/6) is NOT reached"]
# recombination ratio
zr = 1090; ar = 1/(1+zr); Er = E_of_a(ar)
ratio = ar*Er/1.0
out += [f"  recombination: (aH)_rec/(aH)_0 = {ratio:.1f} (Gemini 33, matter-only); rho_DE(rec)/rho_DE(0) = {ratio**-0.5:.3f} (Gemini 0.174)",
        f"  dark energy share at recombination: {Ode*ratio**-0.5/Er**2:.2e} of the total (negligible: no early-dark-energy problem)"]
# --- memory version E9: dX/dlna = -6X + 1/(aH), rho_DE = A X ; A fixed so rho_DE(today) = Ode
def run(A, a0=1e-3, a1=1e4):
    def E(a, X): return np.sqrt(Om*a**-3 + Or*a**-4 + A*X)
    def rhs(la, y):
        a = np.exp(la); return [-6*y[0] + 1/(a*E(a, y[0]))]
    a = a0; X0 = 1/(a*np.sqrt(Om*a**-3+Or*a**-4))/(6.5)   # start on the tracking solution
    s = solve_ivp(rhs, [np.log(a0), np.log(a1)], [X0], dense_output=True, rtol=1e-10, atol=1e-14)
    return s, E
def today(A): s, E = run(A, a1=1.0); return A*s.sol(0.0)[0] - Ode
A = brentq(today, 1e-3, 1e4); s, E = run(A)
def w_mem(a, h=1e-4):
    lr = lambda la: np.log(s.sol(la)[0]); return -1 - (lr(np.log(a)+h) - lr(np.log(a)-h))/(2*h)/3
zc2 = brentq(lambda z: w_mem(1/(1+z)) + 1, 0.05, 2)
afar = 1e4; aHX = afar*E(afar, s.sol(np.log(afar))[0])*s.sol(np.log(afar))[0]
out += ["", "MEMORY VERSION (E9)",
        f"  w today {w_mem(1.0):.3f}; crossing z = {zc2:.2f}; matter era (z=200) w = {w_mem(1/201):.3f} (Gemini -7/6 = {-7/6:.3f})",
        f"  far future (a = 1e4): w = {w_mem(afar):.4f}, aH*X = {aHX:.4f}  (analytic: a ~ t^3, w = -7/9 = {-7/9:.4f}, aHX = 3/16 = {3/16:.4f})"]
# --- Claim 2 sound horizon: r_d scales roughly as 1/m_e (recombination happens earlier, at higher T ~ m_e)
out += ["", "CLAIM 2  m_e(rec)/m_e(0) = 1.01",
        f"  simple scaling r_d ~ 1/m_e: 147.09 -> {147.09/1.01:.1f} Mpc (Gemini 143-144 would need m_e ratio ~ {147.09/143.5:.3f})",
        "  our full pipeline (predictions/README.md): H0 67.5 -> ~69.5; SH0ES 73.0: tension eased, NOT solved"]
# --- size: rho_P / sqrt(N4) with N4 = (t/t_P)^4 cells in the visible past
tP = 5.39e-44; t0 = 13.8e9*3.156e7; rhoP = 5.16e96   # kg/m^3
rho_crit = 3*(67.5e3/3.086e22)**2/(8*np.pi*6.674e-11)
out += ["", "SIZE  rho_P / sqrt(N4), N4 = (t/t_P)^4",
        f"  today: {rhoP/(t0/tP)**2:.2e} kg/m^3 vs critical density {rho_crit:.2e} (ratio {rhoP/(t0/tP)**2/rho_crit:.2f})",
        f"  at t = 1 s: {rhoP/(1/tP)**2:.2e} vs total density then ~ 3/(32 pi G t^2) = {3/(32*np.pi*6.674e-11):.2e}",
        "  -> it tracks the TOTAL density at every epoch (rho_P t_P^2 / t^2 = 1/(G t^2)): automatic, not a prediction"]
txt = "\n".join(out); print(txt); open("gemini_check.txt", "w").write(txt + "\n")
