"""
ITERATION 39b (declared revision of 39, same criteria R1-R3): the first run failed by DESIGN, not necessarily physics:
 (i) the bump (width 0.12) spanned only ~1-5 cells; (ii) random Delaunay grids contain 'sliver' tetrahedra, which the
 convergence theorem (Cheeger-Mueller-Schrader) excludes (it needs well-shaped cells).
Revision: wider bump (s = 0.2), finer grids, and a well-shaped grid (each cube split into 6 tetrahedra: Freudenthal/Kuhn) next
to the random one. Same comparison: S_grid / S_E -> 1?
"""
import numpy as np, itertools
src = open("iter39_regge_to_einstein.py").read().split("SE = S_einstein()")[0]
g = {}; exec(src, g)
g["s"] = 0.2
def S_grid_tets(X, T):
    pairs = g["pairs"]
    E = np.unique(np.sort(np.vstack([T[:, list(p)] for p in pairs]), axis=1), axis=0)
    key = E[:, 0]*len(X) + E[:, 1]; order = np.argsort(key); ks = key[order]
    def eidx(a, b):
        lo, hi = np.minimum(a, b), np.maximum(a, b); return order[np.searchsorted(ks, lo*len(X) + hi)]
    TE = np.stack([eidx(T[:, a], T[:, b]) for a, b in pairs], 1)
    # interior edges: both endpoints strictly inside the box (boundary of the cube grid = outer layer of nodes)
    inside = np.all((X > 1e-9) & (X < 1 - 1e-9), axis=1); interior = np.where(inside[E[:, 0]] & inside[E[:, 1]])[0]
    t = np.linspace(0, 1, 7); w = np.array([1, 4, 2, 4, 2, 4, 1])/18.0
    P0, P1 = X[E[:, 0]], X[E[:, 1]]
    l = sum(wi*np.exp(g["phi"](P0 + ti*(P1 - P0))) for ti, wi in zip(t, w))*np.linalg.norm(P1 - P0, axis=1)
    L2 = l[TE]**2; pid = {p: k for k, p in enumerate(pairs)}; pid.update({(b, a): k for (a, b), k in pid.items()})
    th = np.zeros(L2.shape)
    for k, (a, b) in enumerate(pairs):
        cc, d = [v for v in range(4) if v not in (a, b)]
        uu = L2[:, pid[(a, b)]]; vv = L2[:, pid[(a, cc)]]; ww = L2[:, pid[(a, d)]]
        uv = (uu + vv - L2[:, pid[(b, cc)]])/2; uw = (uu + ww - L2[:, pid[(b, d)]])/2; vw = (vv + ww - L2[:, pid[(cc, d)]])/2
        th[:, k] = np.arccos(np.clip((uu*vw - uw*uv)/np.sqrt((uu*vv - uv**2)*(uu*ww - uw**2)), -1, 1))
    tot = np.zeros(len(E)); np.add.at(tot, TE.ravel(), th.ravel()); defi = 2*np.pi - tot
    return np.sum(l[interior]*defi[interior]), np.abs(defi[interior]).max()
def kuhn(m):
    g1 = np.arange(m + 1)/m; X = np.array(np.meshgrid(g1, g1, g1, indexing="ij")).reshape(3, -1).T
    idx = lambda i, j, k: (i*(m + 1) + j)*(m + 1) + k
    T = []
    for i, j, k in itertools.product(range(m), repeat=3):
        base = np.array([i, j, k])
        for perm in itertools.permutations(range(3)):
            v = base.copy(); tet = [idx(*v)]
            for ax in perm: v = v.copy(); v[ax] += 1; tet.append(idx(*v))
            T.append(tet)
    return X, np.array(T)
SE = g["S_einstein"](200)
out = ["ITERATION 39b: revised Regge-to-Einstein test (wider bump, well-shaped grid added)", "", f"S_E = {SE:.6f}", "",
       " grid                      cells      ratio S_grid/S_E"]
for m in (8, 16, 24, 32, 40):
    X, T = kuhn(m); r, _ = S_grid_tets(X, T)
    out.append(f" well-shaped (cube->6 tets) {6*m**3:7d}     {r/SE:.4f}")
rng = np.random.default_rng(7)
from scipy.spatial import Delaunay
for n in (4000, 16000, 64000):
    rs = []
    for rep in range(3):
        X = np.vstack([rng.uniform(0, 1, (n, 3)), np.array(list(itertools.product([0, 1], repeat=3)), float)])
        rs.append(g["S_grid"](X)/SE)
    out.append(f" random Delaunay           {n:7d}     {np.mean(rs):.4f} +/- {np.std(rs):.4f}")
txt = "\n".join(out); print(txt); open("iter39b_regge_revised.txt", "w").write(txt + "\n")
