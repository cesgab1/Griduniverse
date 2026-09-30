"""
Step 3: collapse triggered by grid STRAIN. A superposition of an object in two places puts two different dimples on the grid.
Tension picture: the extra strain energy of 'two dimples at once' is  dE = (strain energy of the difference of the two dimple patterns),
with the grid's stiffness c^4/G. This is exactly the Diosi-Penrose gravitational self-energy; collapse time t = hbar/dE.
The one free choice: how finely the grid resolves each atom's dimple (smearing size R0).
   - grid-sized R0 (0.6 Planck lengths)  = the natural tension-grid value
   - Gran Sasso underground germanium detector (Donadi et al. 2021): spontaneous radiation implies R0 >= 0.54e-10 m
"""
import numpy as np
G = 6.674e-11; hbar = 1.0546e-34; lP = 1.616e-35
mA = 72.6*1.6605e-27                    # germanium atom
def dE(N, R0): return N*G*mA**2/R0*1.2   # displaced by more than R0: each atom's self-energy counted twice minus overlap (~6/5 factor, uniform ball)
print("collapse time of a germanium grain displaced into two places (t = hbar/dE):")
for R0, lab in ((0.603*lP, "grid-sized dimples (0.6 Planck lengths)"), (1e-15, "nucleus-sized"), (0.54e-10, "Gran Sasso minimum allowed")):
    row = []
    for mass in (1e-18, 1e-12, 1e-6, 1e-3):
        N = mass/mA; t = hbar/dE(N, R0); row.append(f"{mass:.0e} kg: {t:.1e} s")
    print(f"   {lab:40s} " + " | ".join(row))
# what the Gran Sasso bound says about grid-sized strain
print(f"\n Gran Sasso excludes R0 below 0.54e-10 m; the grid-sized value is {0.54e-10/(0.603*lP):.0e} times smaller -> excluded")
print(" (even the heaviest spread-out objects measured so far -- molecules of 25,000 atomic mass units -- would collapse in"
      f" {hbar/(G*(25000*1.66e-27)**2/(0.603*lP)):.0e} s with grid-sized strain, yet interfere for milliseconds)")
