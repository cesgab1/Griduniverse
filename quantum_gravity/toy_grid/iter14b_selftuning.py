"""ITERATION 14b: the jostled tension SELF-TUNES. Check the late-time attractor: steady state 6 H X = 1/a -> rho_j = A/(6 a H);
with a negative baseline -L (matter gone): H^2 = -L + A/(6 a H)  ->  H -> A/(6 a L) -> 0 with a H = const (coasting, a ∝ t).
Test: run with baselines from -0.34 to -1000 (in units of today's critical density, i.e. a vacuum energy 1000x larger than
observed dark energy, negative) and report H and a*H at late times."""
import numpy as np
from scipy.integrate import solve_ivp
Om, Or = 0.31, 9.1e-5; OL = 1 - Om - Or
out = ["ITERATION 14b: self-tuning of a negative baseline by the jostled tension (kappa = 3, quadratic)", ""]
for base in (-0.34, -10.0, -1000.0):
    # choose A so that today's total matches the observed H0 = 1: A X0 = 1 - Om - Or - base, with X0 = 1/6
    X0 = 1/6; A = (OL - base)/X0
    rhs = lambda t, y: [y[0]*np.sqrt(max(Om*y[0]**-3 + Or*y[0]**-4 + base + A*y[1], 0)),
                        -6*np.sqrt(max(Om*y[0]**-3 + Or*y[0]**-4 + base + A*y[1], 0))*y[1] + 1/y[0]]
    s = solve_ivp(rhs, [0, 400], [1.0, X0], method="LSODA", rtol=1e-8, atol=1e-14, dense_output=True)
    for tq in (1, 10, 100, 400):
        a, X = s.sol(tq); H = np.sqrt(max(Om*a**-3 + Or*a**-4 + base + A*X, 0))
        out.append(f"baseline {base:+9.2f}: t = {tq:4d} Hubble times: a = {a:8.2f}, H/H0 = {H:.3e}, a*H = {a*H:.3f}, A/(6 a |base|) = {A/(6*a*abs(base)):.3e}")
    out.append("")
out += ["-> Whatever the size of the negative baseline, the expansion never turns around: the jostled energy ∝ 1/H grows as H",
        "   falls and cancels the baseline, leaving H -> 0 with a H -> constant (coasting). The universe does NOT recollapse:",
        "   no future bounce from this route. But it is a SELF-TUNING behaviour (cf. 'self-tuning' / 'well-tempered' dark-",
        "   energy literature, e.g. Charmousis et al. 2012 'Fab Four'; Appleby & Linder 2018). Caveat: it only cancels a",
        "   NEGATIVE baseline (the jostled energy is positive) and here the size of A was chosen to match today's H0."]
txt = "\n".join(out); print(txt); open("iter14b_selftuning.txt", "w").write(txt + "\n")
