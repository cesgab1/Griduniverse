"""Iteration 62 validation: transfer-function normalisation; large-L torus -> isotropic."""
import numpy as np
from scipy.spatial.transform import Rotation
from iter62_lib import *
out = ["ITERATION 62 validation", ""]
for kc in (35/CHI, 60/CHI):
    cl = iso_cl_from_transfer(kc)
    out.append(f"k cut = {kc*CHI:.0f}/chi: C_l(transfer)/C_l(CAMB) for l=2..12: " + " ".join(f"{cl[l]/cl_camb[l]:.3f}" for l in range(2, LMAX_T + 1)))
kc = 35/CHI; cl = iso_cl_from_transfer(kc); norm = np.mean([cl_camb[l]/cl[l] for l in range(2, 13)])
out.append(f"normalisation factor applied to torus sums (CAMB/transfer, mean l=2..12): {norm:.4f}")
C = torus_cov(100000.0, Rotation.identity(), kc, 1.0)
diag = np.array([np.mean(np.diag(C)[LL == l]) for l in range(2, 13)])
out.append("L = 100 Gpc torus, mean diagonal / isotropic (same k cut), l=2..12: " + " ".join(f"{diag[l-2]/cl[l]:.3f}" for l in range(2, 13)))
off = np.abs(C - np.diag(np.diag(C))).max()/np.diag(C).max()
out.append(f"   largest off-diagonal / largest diagonal: {off:.3f} (should be small)")
txt = "\n".join(out); print(txt); open("iter62_validate.txt", "w").write(txt + "\n")
