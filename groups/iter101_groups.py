import numpy as np, glob
G = 4.30e-6
v2 = []
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(fn, comments="#", ndmin=2); V = d[:, 1][d[:, 1] > 0]
    if len(V) >= 6: v2.append(V.max()**2)
v2 = np.array(v2); lg = np.log10(np.median(v2)); top = np.log10(np.percentile(v2, 90))
lgrp = np.log10(2*128**2); lcl = np.log10(2*1000**2)
n = np.log10(5.3)/(lcl - lg)
k_pred = 10**(n*(lgrp - lg))
out = [f"SPARC median depth v^2 = {10**lg:.2e} (km/s)^2; 90th percentile {10**top:.2e}",
       f"groups depth 2 sigma^2 (sigma 128) = {10**lgrp:.2e}; clusters (sigma 1000) = {10**lcl:.2e}",
       f"depth-law exponent needed for clusters x5.3: n = {n:.3f}",
       f"predicted group boost k = {k_pred:.2f};  observed k = 0.72/(0.6-1.0) = {0.72/1.0:.2f}-{0.72/0.6:.2f}",
       f"predicted boost for the deepest 10% of SPARC galaxies: {10**(n*(top-lg)):.2f} (iteration 96 robust: ~1.0-1.1 above the lowest quarter)"]
print("\n".join(out)); open("iter101_groups.txt", "w").write("\n".join(out) + "\n")
