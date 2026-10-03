"""How big must each 'exchange kick' be for the observed dark-energy density? (shape of Claim 1 fixed by counting; size not)."""
import numpy as np
c, G, hb = 2.998e8, 6.674e-11, 1.0546e-34
H0 = 67.5e3/3.0857e22; lP = np.sqrt(hb*G/c**3); rhoP = c**5/(hb*G**2)
rho_DE = 0.69*3*H0**2/(8*np.pi*G)
out = [f"rho_DE today = {rho_DE:.2e} kg/m^3 = {rho_DE/rhoP:.1e} Planck densities"]
for lab, l in (("Planck length", lP), ("max allowed cell 5.7e-28 m", 5.7e-28)):
    N = c/(2*H0*l); out.append(f"link = {lab:28s}: exchanges per link per stretch time N = {N:.1e}; sqrt(N) = {np.sqrt(N):.1e};"
                               f" kick needed per exchange = {rho_DE/np.sqrt(N)/rhoP:.1e} Planck densities")
out.append("A natural kick would be ~1 Planck density (or ~1 cell energy); the counting therefore does NOT explain the size of dark energy.")
out.append("Breakdown: rho ∝ adot^-1/2 diverges where adot = 0 (the Bounce). There the stretch time 1/H is infinite and the counting")
out.append("(N crossings per stretch time) is meaningless: the law is a late-universe law; it must be cut off near the Bounce.")
out.append("Links stretching since z_f: allowed stretch per link <= 3.5e7, so z_f <= 3.5e7 if links began at the Planck length;")
out.append("dark energy only matters at z < ~5, so the s = 1 epoch needs z_f > ~5 only: no conflict with the cell-creation requirement.")
txt = "\n".join(out); print(txt); open("magnitude.txt", "w").write(txt + "\n")
