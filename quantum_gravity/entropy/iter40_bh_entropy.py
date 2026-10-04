"""
ITERATION 40 (QG): black-hole entropy coefficient 1/4 from the grid.
Ingredients: (1) our measured entanglement entropy per unit wall area of ONE grid field: s0 = 0.0603 per cell area
(area_law.txt; S ~ R^2.00); (2) N layers ('species') of that field; (3) INDUCED gravity (Sakharov): Newton's constant is produced
by the same N fields. BORROWED theorem (Jacobson 1994; Susskind & Uglum 1994): for scalar fields with a common cutoff, the
entanglement entropy equals A/(4 G_induced) automatically -- the 1/4 is not tuned.
What the grid adds: matching S = N s0 A / a^2 to A / (4 l_P^2) fixes the cell size a in terms of N:
     a = 2 sqrt(s0 N) l_P      (the species bound guessed a = sqrt(N) l_P with an unknown O(1) factor)
Criterion written before computing: the factor 2 sqrt(s0) must be O(1) (0.1-10) for the picture to be consistent.
Caveats: gauge fields / gravitons add 'contact terms' (Kabat 1995); our tension field is a scalar, so the clean case applies.
"""
import numpy as np
s0 = 0.0603; lP = 1.616e-35
f = 2*np.sqrt(s0)
a_grb = 5.7e-28                                      # largest cell allowed by gamma-ray-burst timing (earlier work)
Nmax = (a_grb/(f*lP))**2; Nmin = 2.9e7               # smoothness bound from real Planck low-l data
out = ["ITERATION 40: black-hole entropy coefficient from the grid (induced gravity)", "",
       f"entanglement per cell area, one field: s0 = {s0}",
       f"cell size from S = A/4: a = {f:.3f} x sqrt(N) x l_P   (factor {f:.2f}: O(1) -> consistent; species-bound guess was 1)",
       f"window: N = {Nmin:.1e} (smoothness) to {Nmax:.1e} (gamma-ray bursts) -> cell size {f*np.sqrt(Nmin)*lP:.1e} to {a_grb:.1e} m",
       f"gravity's true cutoff (energy where the grid shows): {1.22e19/(f*np.sqrt(Nmax)):.1e} to {1.22e19/(f*np.sqrt(Nmin)):.1e} GeV",
       "", "Reading: with gravity induced by the layers, Bekenstein-Hawking's 1/4 comes out automatically (borrowed theorem), and",
       "the grid's own entanglement fixes the previously unknown factor linking cell size to the number of layers (0.49).",
       "Not new physics by itself: it ties three of our numbers (layers, cell size, entropy) into one consistent relation."]
txt = "\n".join(out); print(txt); open("iter40_bh_entropy.txt", "w").write(txt + "\n")
