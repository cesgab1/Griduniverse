import numpy as np
G, c, hbar = 6.674e-11, 2.998e8, 1.0546e-34; Ms = 1.989e30
lP = np.sqrt(hbar*G/c**3)
rs = 2*G*Ms/c**2; A = 4*np.pi*rs**2
S_BH = A/(4*lP**2)
out = [f"Solar-mass black hole: horizon radius {rs/1e3:.2f} km, area {A:.2e} m^2, entropy {S_BH:.2e} (k_B)", ""]
for name, l in (("Planck-born stretched tension block 2.5 mm", 2.5e-3), ("tension block at window top 5.7e-28 m", 5.7e-28)):
    S = np.log(2)*A/l**2
    out.append(f"B1 {name}: one bit per block -> {S:.1e}; short by x{S_BH/S:.0e}")
for s0, nm in ((np.log(2), "one bit"), (1.0, "one nat")):
    out.append(f"B2 required block size for {nm} per block: l = {2*np.sqrt(s0):.3f} l_P = {2*np.sqrt(s0)*lP:.2e} m")
out.append(f"Curiosity (not evidence): far-future value 1.637 l_P vs one-bit 1.665 l_P (l^2: 2.68 vs 2.77, 3% apart)")
txt = "\n".join(out); print(txt); open("iter89_bh_blocks.txt", "w").write(txt + "\n")
