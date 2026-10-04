"""GR track G0b: the two places GR stops (time zero, overall shape) and how the grid's answers must connect. Numbers only."""
import numpy as np
lP, TP_K, T0 = 1.616e-35, 1.417e32, 2.7255
Gpc = 3.0857e25
stretch = TP_K/T0                                   # how much space has stretched since Planck temperature (radiation era T ~ 1/a)
out = ["GR track G0b: connecting the two gaps", "",
       f"stretch since the Planck temperature: {stretch:.1e}",
       f"a Planck-length cell from then is today {lP*stretch:.1e} m across"]
# piece 1 -> piece 2: a universe BORN small (quantum creation, compact shape) with no later inflation
for L0_lP in (1, 1e10, 1e20):
    out.append(f"  born as a wrap-around space {L0_lP:.0e} Planck lengths across -> wraps today at {L0_lP*lP*stretch:.1e} m")
for Lt in (27.7, 30, 60):
    L0 = Lt*Gpc/stretch
    out.append(f"  wrap today {Lt:5.1f} Gpc needs a starting size {L0:.1e} m = {L0/lP:.1e} Planck lengths -> {(L0/lP)**3:.1e} starting cells")
n_gamma = 4.11e8
out.append(f"cells (stretched Planck cells) per CMB photon today: {1/(n_gamma*(lP*stretch)**3):.1f}  (automatic: both scale with temperature^3)")
# grid density cap: one cell-quantum of energy per cell, cells a = 0.49 sqrt(N) l_P
for N in (4, 41, 163):
    out.append(f"grid's maximum density with N = {N:3d} layers: rho_P / {0.058*N**2:.0f}  (cell {0.49*np.sqrt(N):.1f} l_P)")
txt = "\n".join(out); print(txt); open("G0b_two_gaps.txt", "w").write(txt + "\n")
