"""
Which physical quantity could switch the dark fluid between 'clumps' and 'stays out'?
It must clump in: CMB-era perturbations, Lyman-alpha forest structures (z~3, ~1e9-1e11 Msun), today's large-scale structure, clusters.
It must stay out of: galaxy halos (1e11-1e12.5 Msun, forming z~1).
A simple threshold works only if galaxies sit at an EXTREME of the quantity (above or below every 'clump' case).
Order-of-magnitude values (physical units, typical systems):
  rate     = 1/(collapse or growth time)   [1/Gyr]
  speed    = typical velocity (dispersion or bulk flow)   [km/s]
  stress   = fluid density x speed^2   [rho_mean_today x (km/s)^2]
  mass     = mass scale involved   [Msun]
  depth    = sqrt(|potential|)   [km/s]
"""
import numpy as np
H0inv = 14.4
H = lambda z: np.sqrt(0.31*(1+z)**3 + 0.69)
cases = {
  # name: (z, rate, speed, density/mean_today, mass, depth)
  "CMB-era perturbations (must clump)":   (1000, H(1000)/H0inv, 10.0, 0.31**-1*0+ (1001)**3, 1e15, 950.0),
  "Lyman-alpha structures z~3 (clump)":    (3, H(3)/H0inv, 100.0, 64*3, 1e10, 100.0),
  "large-scale structure today (clump)":   (0, 1/H0inv, 300.0, 1.0, 1e15, 950.0),
  "clusters at formation (clump)":         (0.4, 1/1.3, 1000.0, 990, 1e15, 2000.0),
  "GALAXY halos at formation (stay out)":  (1.0, 1/0.9, 200.0, 2045, 1e12, 500.0),
}
names = ["rate 1/Gyr", "speed km/s", "stress", "mass Msun", "depth km/s"]
vals = {k: [v[1], v[2], v[3]*v[2]**2, v[4], v[5]] for k, v in cases.items()}
print(f"{'':40s}" + "".join(f"{n:>13s}" for n in names))
for k, v in vals.items(): print(f"{k:40s}" + "".join(f"{x:13.2e}" for x in v))
g = vals["GALAXY halos at formation (stay out)"]
print("\nIs the galaxy case at an extreme (needed for a single threshold)?")
for i, n in enumerate(names):
    others = [v[i] for k, v in vals.items() if "GALAXY" not in k]
    ext = "yes (lowest)" if g[i] < min(others) else ("yes (highest)" if g[i] > max(others) else "NO")
    lo = [k.split(" (")[0] for k, v in vals.items() if "GALAXY" not in k and v[i] < g[i]]
    hi = [k.split(" (")[0] for k, v in vals.items() if "GALAXY" not in k and v[i] > g[i]]
    print(f"   {n:12s}: {ext:14s}  below galaxies: {', '.join(lo) or '-'} | above: {', '.join(hi) or '-'}")
