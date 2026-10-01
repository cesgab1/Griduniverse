"""
Four grid geometries in the same ball (radius R, ~same number of nodes), node 0 = the centre where the mass sits.
  cubic   : symmetric cubic lattice (spacing 1)
  mosaic  : random points, Delaunay tetrahedra edges (isotropic random tiling)
  web     : 3-D spider web: radial spokes along evenly spread directions + rings joining neighbouring spokes at each radius
  fungal  : grown network: hyphal tips walk with persistence, branch, and fuse (anastomose) when they meet another hypha;
            started from many spores, so it is one connected, irregular, low-connectivity network
Edge conductance (static) = 1/length (unit cross-section).
"""
import numpy as np
from scipy.spatial import Delaunay, cKDTree
R = 16.0
def _ball_mask(P, R): return np.linalg.norm(P, axis=1) <= R
def cubic(seed=0):
    g = np.arange(-int(R), int(R) + 1); X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    P = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1).astype(float); P = P[_ball_mask(P, R)]
    P = P[np.argsort(np.linalg.norm(P, axis=1))]                       # centre first
    t = cKDTree(P); E = np.array(sorted(t.query_pairs(1.01)))
    return P, E
def mosaic(seed=0):
    rng = np.random.default_rng(seed); n = int(4/3*np.pi*R**3)          # same density as the cubic lattice
    P = rng.uniform(-R, R, (int(n*6/np.pi*1.05), 3)); P = P[_ball_mask(P, R)][:n]
    P = np.vstack([[0, 0, 0], P])
    tri = Delaunay(P); S = tri.simplices
    E = np.unique(np.sort(np.vstack([S[:, [i, j]] for i in range(4) for j in range(i+1, 4)]), axis=1), axis=0)
    return P, E
def web(seed=0, nspoke=None):
    n_target = int(4/3*np.pi*R**3)
    nr = int(R); nspoke = nspoke or n_target//nr
    i = np.arange(nspoke) + 0.5; th = np.arccos(1 - 2*i/nspoke); ph = np.pi*(1 + 5**0.5)*i
    U = np.stack([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)], 1)
    P = [np.zeros(3)]; E = []
    for k in range(1, nr + 1):
        P.extend(U*k)
    P = np.array(P)
    sid = lambda s, k: 1 + (k - 1)*nspoke + s
    for s in range(nspoke):                                             # spokes
        E.append((0, sid(s, 1)))
        for k in range(1, nr): E.append((sid(s, k), sid(s, k + 1)))
    nb = cKDTree(U).query(U, 7)[1][:, 1:]                               # rings: join each spoke to its ~6 angular neighbours
    for k in range(1, nr + 1):
        for s in range(nspoke):
            for t2 in nb[s]:
                if s < t2: E.append((sid(s, k), sid(t2, k)))
    return P, np.unique(np.sort(np.array(E), axis=1), axis=0)
def fungal(seed=0, branch=0.08, persist=0.85, fuse=0.6):
    rng = np.random.default_rng(seed); n_target = int(4/3*np.pi*R**3)
    P = [np.zeros(3)]; E = []; tips = []
    spores = rng.uniform(-R, R, (400, 3)); spores = spores[_ball_mask(spores, R)][:150]
    for s in spores:
        P.append(s); d = rng.normal(size=3); tips.append([len(P)-1, d/np.linalg.norm(d)])
    d = rng.normal(size=(6, 3))
    for v in d: tips.append([0, v/np.linalg.norm(v)])                  # hyphae from the centre too
    tree = None; rebuild = 0
    while len(P) < n_target and tips:
        if tree is None or len(P) - rebuild > 300: tree = cKDTree(np.array(P)); rebuild = len(P)
        new_tips = []
        for node, dirn in tips:
            dirn = persist*dirn + (1 - persist)*rng.normal(size=3); dirn /= np.linalg.norm(dirn)
            x = P[node] + dirn
            if np.linalg.norm(x) > R: continue
            near = tree.query_ball_point(x, fuse)
            near = [j for j in near if j != node]
            if near:                                                    # anastomosis: fuse into the existing hypha
                E.append((node, near[0])); continue
            P.append(x); j = len(P) - 1; E.append((node, j)); new_tips.append([j, dirn])
            if rng.random() < branch:
                b = np.cross(dirn, rng.normal(size=3)); b /= np.linalg.norm(b)
                nd = dirn + b; new_tips.append([j, nd/np.linalg.norm(nd)])
        tips = new_tips
        if not tips:                                                    # re-sprout from random existing nodes
            for j in rng.integers(0, len(P), 30):
                v = rng.normal(size=3); tips.append([j, v/np.linalg.norm(v)])
    P = np.array(P); E = np.unique(np.sort(np.array(E), axis=1), axis=0); E = E[E[:, 0] != E[:, 1]]
    return P, E
BUILDERS = {"cubic": cubic, "mosaic": mosaic, "web": web, "fungal": fungal}
