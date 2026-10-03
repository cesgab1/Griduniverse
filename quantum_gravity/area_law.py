"""
Area-law test on the random grid (quantum version of the tension links).
Treat each node as a quantum oscillator coupled by the tension links (K = weighted graph Laplacian, w_ij = A_ij/d_ij,
plus a tiny mass to remove the zero mode). Compute the ground-state entanglement entropy S of a ball of cells with the rest.
Question: does S follow the VOLUME of the ball (cells inside), the WALL AREA crossing its boundary, or the NUMBER of cut links?
Black-hole entropy scales with area; this is Srednicki (1993) / Bombelli et al. (1986) done on our random Mosaic.
"""
import numpy as np, os
from scipy.spatial import Voronoi, ConvexHull
rng = np.random.default_rng(11); HERE = os.path.dirname(os.path.abspath(__file__))
L, n = 14.0, 2744; P = rng.uniform(0, L, (n, 3))
sh = np.array([(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)])*L
Q = np.vstack([P + s for s in sh]); home = 13*n; vor = Voronoi(Q)
K = np.zeros((n, n)); links = []
for (a, b), rv in zip(vor.ridge_points, vor.ridge_vertices):
    if -1 in rv: continue
    ia, ib = a - home, b - home
    if not (0 <= ia < n or 0 <= ib < n): continue
    i, j = a % n, b % n                       # periodic identification (copy index = k*n + i)
    if i == j: continue
    if not (0 <= ia < n): continue            # count each wall once from the home copy of a
    V = vor.vertices[rv]; c = V.mean(0); area = 0.5*np.linalg.norm(sum(np.cross(V[k]-c, V[(k+1) % len(V)]-c) for k in range(len(V))))
    d = np.linalg.norm(Q[a] - Q[b]); w = area/d
    if b - home >= 0 and b - home < n and b < a: continue   # both in home: keep one
    K[i, j] -= w; K[j, i] -= w; K[i, i] += w; K[j, j] += w; links.append((i, j, area, Q[a], Q[b]))
K += 1e-4*np.eye(n)
ev, U = np.linalg.eigh(K); om = np.sqrt(ev)
X = (U/om) @ U.T/2; Pm = (U*om) @ U.T/2
def S_of(idx):
    M = X[np.ix_(idx, idx)] @ Pm[np.ix_(idx, idx)]; nu = np.sqrt(np.clip(np.linalg.eigvals(M).real, 0.25, None))
    a, b = nu + .5, nu - .5; return float(np.sum(a*np.log(a) - np.where(b > 1e-12, b*np.log(np.where(b > 1e-12, b, 1)), 0)))
ctr = np.full(3, L/2); rows = []
for R in np.linspace(1.6, 5.6, 11):
    dvec = (P - ctr + L/2) % L - L/2; inside = np.linalg.norm(dvec, axis=1) < R; idx = np.where(inside)[0]
    cutA = sum(A for i, j, A, qa, qb in links if inside[i] != inside[j]); cutN = sum(1 for i, j, *_ in links if inside[i] != inside[j])
    rows.append((R, len(idx), cutA, cutN, S_of(idx)))
rows = np.array(rows); out = []
out.append(f"random 3-D Mosaic, {n} cells, {len(links)} tension links; ground state of coupled oscillators\n")
out.append(" R     cells_in   wall_area_cut   links_cut    S")
for r in rows: out.append(f"{r[0]:4.1f}  {int(r[1]):8d}   {r[2]:12.1f}   {int(r[3]):8d}   {r[4]:7.2f}")
def fit(x, y):
    A = np.c_[x, np.ones_like(x)]; c, res, *_ = np.linalg.lstsq(A, y, rcond=None); yhat = A @ c
    return c, 1 - np.sum((y - yhat)**2)/np.sum((y - y.mean())**2)
S = rows[:, 4]
for name, x in (("volume (cells inside)", rows[:, 1]), ("wall area cut", rows[:, 2]), ("links cut", rows[:, 3])):
    c, r2 = fit(x, S); out.append(f"S vs {name:22s}: slope {c[0]:.4f}, R^2 = {r2:.4f}")
sl = np.polyfit(np.log(rows[:, 0]), np.log(S), 1)[0]; out.append(f"S ~ R^{sl:.2f}   (area law: 2, volume law: 3)")
c, _ = fit(rows[:, 2], S); out.append(f"entropy per unit wall area = {c[0]:.4f} (in cell units); Bekenstein-Hawking needs 1/(4 l_P^2) = 0.25 per Planck area")
txt = "\n".join(out); print(txt); open(os.path.join(HERE, "area_law.txt"), "w").write(txt + "\n")
