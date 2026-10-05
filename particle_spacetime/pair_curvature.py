"""How much does space-time bend when light turns into an electron-positron pair? Standard physics numbers."""
import numpy as np
G, c, hbar, me = 6.674e-11, 2.998e8, 1.0546e-34, 9.109e-31
E = 2*me*c**2                          # pair threshold energy
lc = hbar/(me*c)                       # reduced Compton wavelength: size of the region
M = E/c**2
phi = G*M/(lc*c**2)                    # dimensionless bending (potential / c^2)
rs = 2*G*me/c**2
out = [f"Pair threshold energy {E/1.602e-13:.3f} MeV, region ~ reduced Compton wavelength {lc:.2e} m",
       f"Space-time bending at the region's edge (G M / r c^2): {phi:.1e}",
       f"Electron: gravity radius 2Gm/c^2 = {rs:.1e} m vs its quantum size {lc:.1e} m -> ratio {rs/lc:.1e}",
       f"Gravity vs electric force between two electrons: {G*me**2/(8.988e9*(1.602e-19)**2):.1e}",
       f"Best detectors: LIGO strain ~1e-21 -> short by x{1e-21/phi:.0e}; smallest mass whose gravity was measured ~9e-5 kg (90 mg) -> x{9e-5/M:.0e} heavier",
       f"Grid: region = {lc/2.5e-3:.1e} of a 2.5 mm cell; = {lc/1.616e-35:.1e} Planck lengths across"]
txt = "\n".join(out); print(txt); open("pair_curvature.txt", "w").write(txt + "\n")
