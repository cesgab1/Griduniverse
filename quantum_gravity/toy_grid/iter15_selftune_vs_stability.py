"""
ITERATION 15 (iterations 10 + 14 combined): can self-tuning explain the SIZE once the coupling obeys iteration 10's
stability/featherweight bound g <= (3H)^2 ?
Self-tuned steady state: jostled energy cancels the negative baseline L:  N g^2 T_s / (24 pi H) ≈ L (+ observed DE)
-> H = N g^2 T_s / (24 pi L). With g <= 9 H^2:  H <= N 81 H^4 T_s / (24 pi L)  ->  H^3 >= 24 pi L / (81 N T_s).
N T_s is capped by Delta N_eff (iteration 10: max 2.2e7 eV at N = 1.2e15).
"""
import numpy as np
H0 = 1.44e-33; NT = 2.2e7; rho_obs = 2.54e-11
out = ["ITERATION 15: self-tuning vs stability (eV units)", ""]
for lab, L in (("baseline as small as observed dark energy", rho_obs), ("baseline ~ (1 meV)^4", 1e-12),
               ("baseline ~ (1 TeV)^4 (particle-physics vacuum)", 1e48)):
    Hmin = (24*np.pi*L/(81*NT))**(1/3)
    out.append(f"{lab:48s}: stable self-tuning needs H >= {Hmin:.1e} eV  = {Hmin/H0:.0e} x today's H0")
out += ["", "-> With a stable coupling, self-tuning would hold the expansion rate at >= 1e26-1e46 times today's value: it cannot produce",
        "   our small H. The self-tuning attractor is real in the toy, but it does not solve the size problem; it moves it.",
        "   (Without the stability bound it would 'work' for any baseline, which is why the bound matters.)"]
txt = "\n".join(out); print(txt); open("iter15_selftune_vs_stability.txt", "w").write(txt + "\n")
