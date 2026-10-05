import numpy as np
hbar, c, G = 1.0546e-34, 2.998e8, 6.674e-11
lP = np.sqrt(hbar*G/c**3); T0, kB = 2.7255, 8.617e-14
def gs(T): return 106.75 if T > 160 else 61.75 if T > 0.15 else 10.75 if T > 5e-4 else 3.91
a_of_T = lambda T: (T0*kB/T)*(3.91/gs(T))**(1/3)
rDE = 5.26e-10; RH = c/(67.4e3/3.0857e22)
LDE = (hbar*c/rDE)**0.25
kap = lambda l: rDE/(hbar*c/l**4*np.sqrt(RH/l))
out = ["Iteration 86: Planck-sized cells born at a moment, stretched to today. PREREG_86.md", "",
       f"dark energy's own length (hbar c/rho)^(1/4) = {LDE*1e3:.3f} mm; averaging scale (iteration 48) 0.03 mm; Eot-Wash tested to 0.052 mm"]
for name, T in (("Planck time", 1.2209e19), ("reheating 1e15 GeV", 1e15), ("electroweak 160 GeV", 160)):
    a = a_of_T(T); l0 = lP/a
    out.append(f"born at {name}: stretch x{1/a:.2e} -> cell today {l0:.2e} m | kappa needed {kap(l0):.1e} | vs DE length x{l0/LDE:.1f} | "
               f"vs photon limit x{l0/5.7e-28:.0e}")
out.append(f"with 60 e-folds of inflation before a Planck birth: x{np.exp(60):.1e} more -> {lP/a_of_T(1.2209e19)*np.exp(60):.1e} m")
txt = "\n".join(out); print(txt); open("iter86_grid_history.txt", "w").write(txt + "\n")
