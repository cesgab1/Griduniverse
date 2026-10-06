"""Solver sanity: kappa from the starting guess (point-mass potential) vs the converged A2 -> proves the solver moved."""
import os, numpy as np
here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, "iter114_bullet_eg.py")).read()
src = src.split("# ---------- fine grid")[0] + src.split("# ---------- A2: full non-linear solution on a coarse grid ----------")[1].split("t = time.time()")[0]
exec(src)
lap = sum(np.gradient(np.gradient(phi0, hc, axis=i), hc, axis=i) for i in range(3))
Sig = (lap/(4*np.pi*Gk))[:, :, 1:-1].sum(axis=2)*hc
l, k = score("STARTING GUESS (single point mass at the centre, not a solution)", cx, cy, cz, hc, Sig_ph=Sig)
open(os.path.join(here, "iter114_solvercheck.txt"), "w").write("\n".join(l) + "\n")
