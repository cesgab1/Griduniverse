"""
ITERATION 17: Weyl curvature on the grid. Ricci-like = cells squeezed evenly; Weyl-like = cells DISTORTED coherently
(stretched one way, squeezed another). Per-cell shape anisotropy = traceless part of each cell's second-moment tensor;
signal = mean of (xx - yy) component; averaged over patches of n cells.
Two ways a tidal strain e can act:
  (A) GEOMETRIC: the existing cells are stretched (link lengths change): shapes carry the distortion.
  (B) RE-TESSELLATED: cell centres move with the strain but cells are rebuilt by the mosaic rule (as when cells are
      added/renewed): how much distortion survives in the shapes?
"""
import numpy as np
from scipy.spatial import Voronoi
rng = np.random.default_rng(2); L, n = 16.0, 4096
sh = np.array([(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)])*L
def cells(P):
    Q = np.vstack([P + s for s in sh]); home = 13*len(P); vor = Voronoi(Q); out = []
    for i in range(len(P)):
        reg = vor.regions[vor.point_region[home + i]]
        out.append(None if (-1 in reg or not reg) else vor.vertices[reg] - P[i])
    return out
def signal(Vs, S=np.eye(3)):
    s = []
    for V in Vs:
        if V is None: s.append(np.nan); continue
        W = V @ S; M = W.T @ W/len(W); M = M/np.trace(M); s.append(M[0, 0] - M[1, 1])
    return np.array(s)
P0 = rng.uniform(0, L, (n, 3)); V0 = cells(P0); base = signal(V0)
def patches(sig, P, m):
    ok = ~np.isnan(sig); k = int(round((ok.sum()/m)**(1/3))); idx = np.floor((P[ok] % L)/(L/k)).astype(int).clip(0, k - 1)
    lab = idx[:, 0]*k*k + idx[:, 1]*k + idx[:, 2]; return np.array([sig[ok][lab == j].mean() for j in np.unique(lab)])
out = ["ITERATION 17: coherent distortion (Weyl-like) vs random irregularity on a random 3-D mosaic (4096 cells)", "",
       f"no strain: single-cell (xx - yy) spread {np.nanstd(base):.3f}; patch spread: " +
       ", ".join(f"~{m} cells {patches(base, P0, m).std():.4f}" for m in (8, 64, 512))]
for e in (0.02, 0.05, 0.10):
    S = np.diag([1 + e, 1 - e, 1.0])
    gA = signal(V0, S); PB = ((P0 - L/2) @ S + L/2); VB = cells(PB % L); gB = signal(VB)
    out.append(f"strain e = {e:4.2f}: (A) cells stretched: mean shift {np.nanmean(gA) - np.nanmean(base):+.4f}   "
               f"(B) cells rebuilt: mean shift {np.nanmean(gB) - np.nanmean(base):+.4f}   -> retained fraction {(np.nanmean(gB) - np.nanmean(base))/(np.nanmean(gA) - np.nanmean(base)):.2f}")
out += ["", "-> Random irregularity averages away over patches (about 1/sqrt(n)): it is not Weyl curvature, so a random mosaic is",
        "   compatible with Penrose's 'no Weyl at the start'. A coherent strain survives averaging (A).",
        "-> When cells are REBUILT by the mosaic rule (B), most of the distortion is erased from the shapes: tidal distortion",
        "   cannot be stored in cell shapes of a grid that renews its cells; it must live in the LINK LENGTHS (curvature as",
        "   deficit angles, quantum_gravity/regge_gauge.py). Gravitational entropy = coherent link-length distortion."]
txt = "\n".join(out); print(txt); open("iter17_grid_weyl.txt", "w").write(txt + "\n")
