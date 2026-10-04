"""
ITERATION 14: does a NEGATIVE baseline (jostled fraction f > 1, mildly preferred by data in iteration 13) make the
universe turn around and recollapse - i.e. allow a future bounce and a cycle?
Forward integration in TIME (H can reach 0): dX/dt = -6 H X + 1/a (memory kappa = 3, iteration 4), jostled energy
rho_j = A X (quadratic), A fixed so rho_j(today) = f * OL; baseline OL (1 - f) (negative for f > 1).
H^2 = Om a^-3 + Or a^-4 + OL(1 - f) + A X ; stop when H^2 <= 0 (turnaround). Units: H0 = 1, t in Hubble times
(1 Hubble time = 14.5 Gyr for h = 0.674).
Note: the jostled part's memory grows as H -> 0 (memory 1/(3H)); the toy equation handles that non-adiabatically.
"""
import numpy as np
from scipy.integrate import solve_ivp
Om, Or, h = 0.31, 9.1e-5, 0.674
OL = 1 - Om - Or
src = open("iter1_memory.py").read().split("# ---------- Part A")[0]
import os; g = {"__name__": "x", "__file__": os.path.abspath("iter1_memory.py")}; exec(src, g)
a_past = np.array([1.0]); X_today_shape = None
out = ["ITERATION 14: future of the universe with baseline (1 - f) OL + jostled part f OL (today)", ""]
for f in (1.0, 1.25, 1.5, 2.0, 3.0):
    # today's X from the self-consistent past history (shape only matters through dX/dt today); reuse rho_hist to get X0
    # here: start in the adiabatic state today, X0 = 1/(6 a H) (kappa = 3), and set A so rho_j(1) = f OL
    X0 = 1/(6*1.0*1.0); A = f*OL/X0
    def rhs(t, y):
        a, X = y; H2 = Om*a**-3 + Or*a**-4 + OL*(1 - f) + A*X
        H = np.sqrt(max(H2, 0.0)); return [a*H, -6*H*X + 1/a]
    def turn(t, y):
        a, X = y; return Om*a**-3 + Or*a**-4 + OL*(1 - f) + A*X
    turn.terminal = True; turn.direction = -1
    s = solve_ivp(rhs, [0, 200], [1.0, X0], events=turn, rtol=1e-9, atol=1e-12, max_step=0.05)
    if s.t_events[0].size:
        tt = s.t_events[0][0]; at = s.y_events[0][0][0]
        out.append(f"f = {f:4.2f} (baseline {OL*(1-f):+.2f}): TURNAROUND after {tt:5.2f} Hubble times ({tt*14.5:5.0f} Gyr), at a = {at:.2f}; then recollapse -> next bounce")
    else:
        out.append(f"f = {f:4.2f} (baseline {OL*(1-f):+.2f}): no turnaround within 200 Hubble times; a grows to {s.y[0][-1]:.1e}")
out += ["", "Data (iteration 13), Delta chi2 vs Lambda (Pantheon+ / DES-Dovekie / Union3):",
        "  f = 1.00: -7.1 / -9.2 / -8.5 ;  f = 1.25: -7.4 / -9.9 / -9.8 ;  f = 1.5: -7.3 / -10.2 / -10.9 ;  f = 2: -5.7 / -9.5 / -12.2 ;  f = 3: +1.6 / -3.7 / -12.2"]
txt = "\n".join(out); print(txt); open("iter14_future.txt", "w").write(txt + "\n")
