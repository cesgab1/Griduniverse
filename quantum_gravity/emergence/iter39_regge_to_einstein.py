"""
ITERATION 39 (QG): does the grid's curvature action become EINSTEIN'S action at large scales?
Put a smooth curved geometry on our random Delaunay grid (the Mosaic) purely through the LINK LENGTHS:
   metric g = e^{2 phi} delta,  phi = A exp(-r^2 / 2 s^2)  (a gentle bump in the middle of the box; flat near the edges)
   l_e = integral of e^{phi} along the link (Simpson, 7 points)
Grid action (Regge 1961): S_grid = sum over interior links of  l_e x deficit_e     (curvature lives on links)
Einstein action for the same space (3-D): S_E = (1/2) INT R sqrt(g) d^3x,  R sqrt(g) = -e^{phi}(4 lap(phi) + 2 |grad phi|^2)
Expectations written BEFORE running (BORROWED theorem: Cheeger, Mueller & Schrader 1984 -- Regge converges to Einstein
in the mean as the grid is refined):
 R1 ratio S_grid / S_E -> 1 as the number of cells grows; within ~10% by a few thousand cells (averaged over random grids).
 R2 the scatter between random grids shrinks with refinement.
 R3 a jittered (shaken regular) grid converges at least as well as the fully random one.
 Kill: no trend toward 1 (e.g. ratio stuck away from 1 by > 20% at the finest grid) would mean our random grid does NOT give
 Einstein gravity at large scales.
"""
import numpy as np
from scipy.spatial import Delaunay, ConvexHull
A, s = 0.06, 0.12; c0 = np.array([0.5, 0.5, 0.5])
def phi(x): r2 = ((x - c0)**2).sum(-1); return A*np.exp(-r2/(2*s*s))
def S_einstein(n=160):
    g = (np.arange(n) + 0.5)/n; X, Y, Z = np.meshgrid(g, g, g, indexing="ij"); P = np.stack([X, Y, Z], -1)
    ph = phi(P); r2 = ((P - c0)**2).sum(-1)
    grad2 = ph**2*r2/s**4; lap = ph*(r2/s**4 - 3/s**2)
    return 0.5*np.sum(-np.exp(ph)*(4*lap + 2*grad2))/n**3
pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
def S_grid(X):
    tri = Delaunay(X); T = tri.simplices
    E = np.unique(np.sort(np.vstack([T[:, list(p)] for p in pairs]), axis=1), axis=0)
    key = E[:, 0]*len(X) + E[:, 1]; order = np.argsort(key); ks = key[order]
    def eidx(a, b):
        lo, hi = np.minimum(a, b), np.maximum(a, b); return order[np.searchsorted(ks, lo*len(X) + hi)]
    TE = np.stack([eidx(T[:, a], T[:, b]) for a, b in pairs], 1)
    hull = ConvexHull(X); hf = hull.simplices
    hull_e = np.unique(np.concatenate([eidx(hf[:, a], hf[:, b]) for a, b in ((0, 1), (0, 2), (1, 2))]))
    interior = np.setdiff1d(np.arange(len(E)), hull_e)
    t = np.linspace(0, 1, 7); w = np.array([1, 4, 2, 4, 2, 4, 1])/18.0   # Simpson
    P0, P1 = X[E[:, 0]], X[E[:, 1]]
    l = sum(wi*np.exp(phi(P0 + ti*(P1 - P0))) for ti, wi in zip(t, w))*np.linalg.norm(P1 - P0, axis=1)
    L2 = l[TE]**2; pid = {p: k for k, p in enumerate(pairs)}; pid.update({(b, a): k for (a, b), k in pid.items()})
    th = np.zeros(L2.shape)
    for k, (a, b) in enumerate(pairs):
        cc, d = [v for v in range(4) if v not in (a, b)]
        uu = L2[:, pid[(a, b)]]; vv = L2[:, pid[(a, cc)]]; ww = L2[:, pid[(a, d)]]
        uv = (uu + vv - L2[:, pid[(b, cc)]])/2; uw = (uu + ww - L2[:, pid[(b, d)]])/2; vw = (vv + ww - L2[:, pid[(cc, d)]])/2
        th[:, k] = np.arccos(np.clip((uu*vw - uw*uv)/np.sqrt((uu*vv - uv**2)*(uu*ww - uw**2)), -1, 1))
    tot = np.zeros(len(E)); np.add.at(tot, TE.ravel(), th.ravel()); defi = 2*np.pi - tot
    return np.sum(l[interior]*defi[interior])
SE = S_einstein()
out = ["ITERATION 39: grid (Regge) action vs Einstein action for the same gentle curvature (expectations committed first)", "",
       f"Einstein action S_E = {SE:.6f}", "", " grid        cells     ratio S_grid/S_E (mean over random grids)   scatter"]
rng = np.random.default_rng(39)
for kind in ("random", "jittered"):
    for n in (1000, 4000, 16000, 64000):
        rs = []
        for rep in range(4 if n <= 16000 else 2):
            if kind == "random": X = rng.uniform(0, 1, (n, 3))
            else:
                m = int(round(n**(1/3))); g = (np.arange(m) + 0.5)/m
                X = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T + rng.uniform(-0.3, 0.3, (m**3, 3))/m
            rs.append(S_grid(X)/SE)
        out.append(f" {kind:9s} {n:7d}         {np.mean(rs):.4f}                                  {np.std(rs):.4f}")
txt = "\n".join(out); print(txt); open("iter39_regge_to_einstein.txt", "w").write(txt + "\n")
