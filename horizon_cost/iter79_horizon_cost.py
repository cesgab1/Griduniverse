"""Iteration 79: horizon cost rules, circularity, and a star clock. See PREREG_79.md."""
import numpy as np
G, c, hbar, kB = 6.674e-11, 2.998e8, 1.0546e-34, 1.381e-23
H0 = 67.4e3/3.0857e22; Gpc = 3.0857e25; Gyr = 3.156e16
rhoc = 3*H0**2/(8*np.pi*G); rDE = 0.685*rhoc               # kg/m^3
lP2 = hbar*G/c**3
out = ["Iteration 79: what does a horizon cost? PREREG_79.md", "", f"measured dark energy {rDE:.2e} kg/m^3", "Step 1:"]
for hname, R in (("event horizon 6.24 Gpc", 6.24*Gpc), ("Hubble radius 4.45 Gpc", c/H0)):
    h1 = 3*c**2/(8*np.pi*G*R**2)
    h2 = hbar*c/(2*2.68*lP2*R**2)/c**2
    T = hbar*c/(2*np.pi*kB*R); h3 = (np.pi**2/15)*(kB*T)**4/(hbar*c)**3/c**2
    out.append(f"  {hname}: H1 {h1/rDE:.2f}x | H2 {h2/rDE:.2f}x | H3 {h3/rDE:.1e}x  (glow T = {T:.1e} K)")
out.append(f"  Circularity: Einstein's equation gives total density x (Hubble radius)^2 = 3c^2/(8 pi G) exactly, for any contents;"
           f" H1 with the Hubble radius = total density = 1/0.685 = {1/0.685:.2f} x dark energy. The horizon is set BY the contents.")
out.append("Step 2 (outside clock: stars):")
alpha, mp, me = 1/137.036, 1.6726e-27, 9.109e-31
mPl = np.sqrt(hbar*c/G)
t_cr = alpha**2*(mp/me)**2*(mPl/mp)**2*hbar/(mp*c**2)
for name, t in (("particle constants (Carr-Rees)", t_cr), ("Sun's measured main-sequence life 10 Gyr", 10*Gyr),
                ("dial clock 8.2 Gyr (reference)", 8.2*Gyr)):
    rho = 3/(8*np.pi*G*t**2)
    out.append(f"  {name}: t = {t/Gyr:.2f} Gyr -> rho = {rho:.2e} kg/m^3 = {rho/rDE:.2f} x measured")
naive = c**5/(hbar*G**2)
out.append(f"  For scale: the naive quantum estimate (Planck density) is {naive/rDE:.1e} x measured.")
txt = "\n".join(out); print(txt); open("iter79_horizon_cost.txt", "w").write(txt + "\n")
