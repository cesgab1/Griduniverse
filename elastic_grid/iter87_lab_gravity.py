import numpy as np
rDE, rW = 5.85e-27, 1.93e4           # kg/m^3: dark energy (mass-equivalent); tungsten test masses
lP = 1.616e-35; T0kB = 2.7255*8.617e-14
out = ["Iteration 87: millimetre tension cells vs short-range gravity (Eot-Wash 2020). PREREG_87.md", ""]
out.append(f"V1 (a): cell-structure gravity / test-mass gravity <= rho_DE/rho_W = {rDE/rW:.1e}; Eot-Wash sensitivity ~1e-2 of gravity near 52 um")
out.append("   -> undetectable by ~29 orders. (b) no extra force (tension field does not propagate). V1 prediction: NULL. Data: null -> consistent, not a test.")
lmax = 38.6e-6
a_min = lP/lmax
T_max = T0kB/a_min*(3.91/106.75)**(1/3)
out.append(f"V2: gravity on 2.5 mm cells would break Newton's law below ~2.5 mm; verified to 52 um -> EXCLUDED (x{2.53e-3/lmax:.0f} too big).")
out.append(f"   Gravity's cells must be <= 38.6 um -> if born Planck-sized, birth at a >= {a_min:.1e}, i.e. T <= {T_max:.1e} GeV (Planck: 1.2e19).")
txt = "\n".join(out); print(txt); open("iter87_lab_gravity.txt", "w").write(txt + "\n")
