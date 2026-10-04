"""
ITERATION 29: can the 'book' end in a recollapse (crunch) in our model?  Sequestering needs a finite total lifetime.
Expectation written before running (analytic): with Claim 1, rho_DE = C (aH)^(-1/2) (s = 1). Friedmann at fixed a:
   f(H) = H^2 - C a^(-1/2) H^(-1/2) - Om a^-3 - Or a^-4 + B = 0,  B >= 0 a negative 'binding' constant
 f -> -inf as H -> 0+, f -> +inf as H -> inf, and f'(H) = 2H + (C/2) a^(-1/2) H^(-3/2) > 0, so exactly ONE positive H exists
 for every a: H can never reach 0, so there is NO turnaround for any B. Late times: H^(1/2) -> C a^(-1/2)/B, i.e. aH -> constant
 (coasting), the same self-tuning found for the toy in iteration 14.
Numerical check below (forward in a, solving f(H) = 0); C fixed so today's dark energy (C - B) = 0.69.
"""
import numpy as np
from scipy.optimize import brentq
Om, Or = 0.31, 9.1e-5; ODE = 1 - Om - Or
out = ["ITERATION 29: Claim 1 + negative constant (the book's 'binding') -- can expansion turn around?", ""]
for B in (0.0, 0.1, 0.69, 3.0, 100.0):
    C = ODE + B
    rows = []
    for a in (1, 3, 10, 100, 1e4, 1e6):
        f = lambda H: H*H - C*a**-0.5*H**-0.5 - Om*a**-3 - Or*a**-4 + B
        H = brentq(f, 1e-30, 1e10)
        rows.append(f"a={a:.0e}: H={H:.2e}, aH={a*H:.3f}")
    out.append(f"B = {B:6.2f}: " + "; ".join(rows))
out += ["", "Result: H > 0 at every a for every B (no turnaround); with B > 0, aH -> constant (coasting forever).",
        "Together with iteration 14 (toy with memory, baselines up to -1000x): our dark-energy mechanism FORBIDS a crunch.",
        "Consequence for the book: if time has a finite total, the end cannot be a recollapse visible in advance. It would have to",
        "be the grid simply stopping (no precursor, nothing to measure today). The book's size rule then makes no testable",
        "prediction beyond the size itself -> PARKED as untestable within this model (not adopted; P4 stays as registered).",
        "What WOULD reopen it: data showing dark energy heading below zero (a crunch-type w(z), as in Andrei, Ijjas & Steinhardt,",
        "PNAS 119, e2200539119 (2022), arXiv:2201.07704), which would also contradict Claim 1 and the toy."]
txt = "\n".join(out); print(txt); open("iter29_can_the_book_close.txt", "w").write(txt + "\n")
