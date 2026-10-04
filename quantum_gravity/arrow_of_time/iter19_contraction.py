"""
ITERATION 19 (time asymmetry of the mechanism): run the jostled-tension equation in a CONTRACTING universe (H < 0).
d<tau^2>/dt = -6 H <tau^2> + Theta_c/(2 pi a): in expansion the -6H term is friction (memory, iteration 4); in contraction
it is ANTI-friction. Compare growth with matter (a^-3), radiation (a^-4) and anisotropic shear = Weyl-type energy (a^-6,
the BKL chaos Penrose associates with a high-entropy crunch). Contracting matter-dominated background from a = 1 to 1e-3.
"""
import numpy as np
from scipy.integrate import solve_ivp
Om = 1.0
def rhs(lna, y):                       # ln a decreasing; dX/dlna = (dX/dt)/H with H = -sqrt(Om a^-3)
    a = np.exp(lna); H = -np.sqrt(Om*a**-3); return [(-6*H*y[0] + 1/a)/H]
s = solve_ivp(rhs, [0, np.log(1e-3)], [1/6], t_eval=np.log([1, 1e-1, 1e-2, 1e-3]), rtol=1e-10, atol=1e-30)
out = ["ITERATION 19: the jostled tension in a CONTRACTING universe", "", "   a        tension energy (rel.)   matter a^-3   radiation a^-4   shear/Weyl a^-6"]
for lna, X in zip(s.t, s.y[0]):
    a = np.exp(lna); out.append(f"   {a:7.0e}   {X/(1/6):14.3e}        {a**-3:9.1e}     {a**-4:9.1e}       {a**-6:9.1e}")
slope = np.polyfit(s.t[1:], np.log(s.y[0][1:]), 1)[0]
out += ["", f"-> in contraction the tension energy grows as a^{slope:.1f}: anti-friction turns the memory into AMPLIFICATION, like",
        "   stiff matter (w = +1). It grows exactly as fast as anisotropic shear (a^-6), so it cannot smooth a collapse",
        "   (ekpyrotic smoothing needs faster growth, w >> 1). The mechanism is intrinsically time-asymmetric: it calms",
        "   the tension only while space expands. A collapse would end with high Weyl curvature, unlike our smooth start:",
        "   Penrose's asymmetry between Big Bang and Big Crunch, here built into the dark-energy mechanism."]
txt = "\n".join(out); print(txt); open("iter19_contraction.txt", "w").write(txt + "\n")
