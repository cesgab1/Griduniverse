import numpy as np
hbar, c, G = 1.0546e-34, 2.998e8, 6.674e-11
H0 = 67.4e3/3.0857e22; RH = c/H0
rDE = 5.85e-27*c**2                       # J/m^3
lP = np.sqrt(hbar*G/c**3)
f = lambda l, k: k*hbar*c/l**4*np.sqrt(RH/l)
out = ["Iteration 85: block size from the tally rule rho = kappa (hbar c/l^4) sqrt(N), N = (c/H0)/l. PREREG_85.md", ""]
for k in (1.0, 0.5):
    l = (k*hbar*c*np.sqrt(RH)/rDE)**(1/4.5)
    out.append(f"kappa = {k}: block size l = {l:.2e} m ({l/lP:.1e} Planck lengths); blocks across horizon N = {RH/l:.1e}")
for name, l in (("window top 5.7e-28 m", 5.7e-28), ("window bottom 4.3e-32 m", 4.3e-32), ("Planck length", lP)):
    out.append(f"kappa needed for l = {name}: {rDE/(hbar*c/l**4*np.sqrt(RH/l)):.1e}")
out.append("")
out.append("Verdict: natural kappa gives centimetre-to-metre blocks (excluded); blocks inside the allowed window need kappa ~ 1e-120..1e-140.")
txt = "\n".join(out); print(txt); open("iter85_block_size.txt", "w").write(txt + "\n")
