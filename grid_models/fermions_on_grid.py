"""
Weak point 2: spin-1/2 matter (electrons, quarks) on the random grid.
The known obstacle: on a regular lattice the simplest Dirac operator gives 2^d copies of every fermion ("doublers",
Nielsen-Ninomiya). Test: does the random grid (piece 0) give ONE fermion, or doublers too?

Operator (Christ-Friedberg-Lee style, 2D Euclidean, periodic box L x L, density 1):
  (D psi)_i = (1/V_i) sum_j (A_ij / 2) gamma.n_ij psi_j      A_ij = shared Voronoi edge length, n_ij = unit vector i->j, V_i = cell area
  This is exact for constant and linear fields (closed Voronoi cells), so it is a consistent discretisation of gamma.grad.
Continuum answer: eigenvalues +-|k|, k = 2 pi (n1, n2)/L, one species -> count N(|lambda| < Lam) = 2 * (# k-points with |k| < Lam).
Regular square lattice, same operator: 4 species in 2D (doublers at k = (pi,0), (0,pi), (pi,pi)).
Also measured: how spread out each eigenmode is (participation ratio); physical modes fill the box, localised ones do not.
"""
import numpy as np, json
from scipy.spatial import Voronoi
from scipy.linalg import eigh
rng = np.random.default_rng(5)
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]])

def poly_area(V):
    c = V.mean(0); a = np.arctan2(V[:, 1] - c[1], V[:, 0] - c[0]); V = V[np.argsort(a)]
    return 0.5 * abs(np.dot(V[:, 0], np.roll(V[:, 1], -1)) - np.dot(V[:, 1], np.roll(V[:, 0], -1)))

def random_grid(L):
    n = rng.poisson(L * L); P = rng.uniform(0, L, (n, 2))
    shifts = np.array([(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)]) * L
    Q = np.concatenate([P + s for s in shifts]); home = np.tile(np.arange(n), 9); centre = 4 * n
    vor = Voronoi(Q)
    links = {}
    for (a, b), rv in zip(vor.ridge_points, vor.ridge_vertices):
        if -1 in rv: continue
        if not (centre <= a < centre + n): a, b = b, a
        if not (centre <= a < centre + n): continue
        A = np.linalg.norm(vor.vertices[rv[0]] - vor.vertices[rv[1]])
        d = Q[b] - Q[a]; links[(home[a], home[b], round(d[0], 9), round(d[1], 9))] = (A, d / np.linalg.norm(d))
        if centre <= b < centre + n:                     # both ends in the home box: this ridge is listed once, add the reverse link
            links[(home[b], home[a], round(-d[0], 9), round(-d[1], 9))] = (A, -d / np.linalg.norm(d))
    Vc = np.array([poly_area(vor.vertices[vor.regions[vor.point_region[centre + i]]]) for i in range(n)])
    return n, links, Vc

def square_grid(L):
    n = L * L; links = {}
    for i in range(L):
        for j in range(L):
            a = i * L + j
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                b = ((i + di) % L) * L + (j + dj) % L
                links[(a, b, di, dj)] = (1.0, np.array([di, dj], float))
    return n, links, np.ones(n)

def spectrum(n, links, Vc):
    H = np.zeros((2 * n, 2 * n), complex)
    for (a, b, *_), (A, nv) in links.items():
        H[2*a:2*a+2, 2*b:2*b+2] += 1j * (A / 2) * (nv[0] * sx + nv[1] * sy)     # i * (anti-Hermitian) = Hermitian
    H = (H + H.conj().T) / 2
    M = np.repeat(Vc, 2)
    lam, U = eigh(H, np.diag(M))
    w = (np.abs(U)**2 * M[:, None]); w = w.reshape(n, 2, -1).sum(1); w /= w.sum(0)
    pr = 1 / (w**2).sum(0) / n                                   # fraction of the box a mode occupies (1 = fully spread)
    return lam, pr

def continuum_count(L, Lam):
    m = int(Lam * L / (2 * np.pi)) + 2; k = 2 * np.pi / L * np.array([(a, b) for a in range(-m, m + 1) for b in range(-m, m + 1)])
    return 2 * np.sum(np.linalg.norm(k, axis=1) < Lam)

L = 36
out = {}
for name, grid in [("square lattice", square_grid), ("random grid", random_grid)]:
    reps = 1 if name == "square lattice" else 3
    counts = {Lam: [] for Lam in (0.3, 0.5, 0.7)}; prs_low = []; first = None
    for r in range(reps):
        n, links, Vc = grid(L)
        lam, pr = spectrum(n, links, Vc)
        for Lam in counts: counts[Lam].append(np.sum(np.abs(lam) < Lam))
        low = np.abs(lam) < 0.5; prs_low.append(np.median(pr[low]))
        if first is None: first = np.sort(np.abs(lam))[:24]
    out[name] = dict(counts={str(k): float(np.mean(v)) for k, v in counts.items()}, pr=float(np.mean(prs_low)), low=first.tolist())
    print(f"\n{name}  ({n} points)")
    for Lam, v in counts.items():
        cc = continuum_count(L, Lam); print(f"  modes with |energy| < {Lam}: {np.mean(v):6.1f}   continuum (one fermion) {cc:4d}   ->  {np.mean(v)/cc:.2f} species")
    print(f"  median spread of low modes: {np.mean(prs_low):.2f} of the box")
    print("  lowest |energies|:", np.round(first[:16], 3))
kc = np.sort(np.linalg.norm(2*np.pi/L*np.array([(a, b) for a in range(-3, 4) for b in range(-3, 4)]), axis=1))
print("\ncontinuum lowest |k| (each x2 for spin):", np.round(np.repeat(kc[:8], 2), 3))
json.dump(out, open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/fermion_results.json", "w"), indent=1)
