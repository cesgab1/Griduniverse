"""
Weak point 1: electromagnetism directly on the random grid (no extra dimension).
Grid: Poisson-random points, density 1 per l_d^3 (spatial slice of the random grid; Christ-Friedberg-Lee random lattice).
Links: Voronoi neighbours. Field: a U(1) phase theta_ij on every link (compact: theta ~ theta + 2 pi).
Hamiltonian lattice gauge theory: H = (g^2/2) sum_links E_ij^2 / w_ij + (1/g^2) sum_plaquettes (1 - cos theta_P)
  E_ij = integer (conjugate to a compact phase)  ->  Gauss law  sum_j E_ij = integer charge at each point.
Static charges: the electric field minimising the energy with Gauss's law obeys the discrete Poisson equation
  sum_j w_ij (phi_i - phi_j) = q_i ,  w_ij = (Voronoi face area)/(link length)   (the grid's own geometry, no tuning)
Tests
  1. Coulomb: charge at the centre of a grounded sphere -> continuum answer phi = (1/4pi)(1/r - 1/R). Fit coefficient + power law.
  2. Isotropy: scatter of phi at fixed r (does the random grid pick directions?).
  3. Light speed vs wavelength: plane-wave dispersion omega^2(k) on the same links -> frame-dependent (Lorentz-violating) corrections.
"""
import numpy as np, json
from scipy.spatial import Voronoi
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
rng = np.random.default_rng(7)

def build(Rb=12.0, pad=3.0):
    Lh = Rb + pad
    n = rng.poisson((2 * Lh)**3)
    P = rng.uniform(-Lh, Lh, (n, 3))
    vor = Voronoi(P)
    I, J, W = [], [], []
    for (a, b), rv in zip(vor.ridge_points, vor.ridge_vertices):
        if -1 in rv or len(rv) < 3: continue
        V = vor.vertices[rv]; nrm = P[b] - P[a]; L = np.linalg.norm(nrm); nrm /= L
        cen = V.mean(0); u = V[0] - cen; u -= u.dot(nrm) * nrm; u /= np.linalg.norm(u); v = np.cross(nrm, u)
        ang = np.arctan2((V - cen) @ v, (V - cen) @ u); V = V[np.argsort(ang)]
        area = 0.5 * np.linalg.norm(np.cross(V - cen, np.roll(V, -1, 0) - cen).sum(0))
        I.append(a); J.append(b); W.append(area / L)
    return P, np.array(I), np.array(J), np.array(W)

def coulomb(P, I, J, W, Rb):
    c = np.argmin(np.linalg.norm(P, axis=1)); X = P - P[c]          # charge on the point nearest the middle
    r = np.linalg.norm(X, axis=1); inside = r < Rb
    idx = -np.ones(len(P), int); idx[inside] = np.arange(inside.sum())
    n = inside.sum()
    diag = np.zeros(n); rows, cols, vals = [], [], []
    for a, b, w in zip(I, J, W):
        for s, t in ((a, b), (b, a)):
            if inside[s]:
                diag[idx[s]] += w
                if inside[t]: rows.append(idx[s]); cols.append(idx[t]); vals.append(-w)   # outside = grounded (phi = 0)
    A = coo_matrix((np.r_[vals, diag], (np.r_[rows, np.arange(n)], np.r_[cols, np.arange(n)])), (n, n)).tocsr()
    q = np.zeros(n); q[idx[c]] = 1.0
    phi = spsolve(A, q)
    rr = r[inside]
    sel = (rr > 2.5) & (rr < Rb - 1.5)
    x = 1 / rr[sel] - 1 / Rb
    k = (phi[sel] @ x) / (x @ x)                                   # phi = k (1/r - 1/R); continuum k = 1/(4 pi)
    # power law of the field: fit phi + k/R = k' r^-p
    p = -np.polyfit(np.log(rr[sel]), np.log(phi[sel] + k / Rb), 1)[0]
    bins = np.arange(2.5, Rb - 1.5, 1.0); scat = []
    for lo in bins:
        m = sel & (rr >= lo) & (rr < lo + 1)
        if m.sum() > 20: scat.append((lo + 0.5, float(np.std(phi[m] - k * (1 / rr[m] - 1 / Rb)) / (k * (1 / (lo + .5) - 1 / Rb)))))
    return k, p, scat

def dispersion(P, I, J, W, Rb):
    r = np.linalg.norm(P, axis=1); inb = r < Rb
    wt = (inb[I].astype(float) + inb[J]) / 2                       # half weight for links crossing the edge
    D = P[J] - P[I]; out = []
    for kmag in [0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.3]:
        vals = []
        for _ in range(12):
            nhat = rng.normal(size=3); nhat /= np.linalg.norm(nhat); kv = kmag * nhat
            num = (wt * W * 2 * (1 - np.cos(D @ kv))).sum(); den = inb.sum()     # mean Voronoi volume = 1
            vals.append(num / den / kmag**2)
        out.append((kmag, float(np.mean(vals)), float(np.std(vals))))
    return out

res = []
for trial in range(4):
    P, I, J, W = build()
    k, p, scat = coulomb(P, I, J, W, 12.0)
    res.append(dict(k=k, p=p, scat=scat))
    print(f"grid {trial}: {len(P)} points, {len(W)} links | Coulomb coefficient k = {k:.5f} (continuum 1/4pi = {1/(4*np.pi):.5f}, ratio {k*4*np.pi:.3f}) | power {p:.3f}")
ks = np.array([r_["k"] for r_ in res]) * 4 * np.pi
print(f"\nCoulomb strength / continuum = {ks.mean():.3f} ± {ks.std()/np.sqrt(len(ks)):.3f}")
print("direction scatter of the potential at distance r (fraction):", [(round(a, 1), round(b, 3)) for a, b in res[0]["scat"]])
disp = dispersion(P, I, J, W, 12.0)
print("\nlight speed^2 vs wavenumber (units of the grid spacing), averaged over directions:")
for kk, m, s in disp: print(f"  k l_d = {kk:4.1f}:  (omega/k)^2 = {m:.4f} ± {s:.4f} (direction spread)")
ks_ = np.array([d[0] for d in disp]); ms = np.array([d[1] for d in disp])
a2 = np.polyfit(ks_[:4]**2, ms[:4], 1)
print(f"fit: (omega/k)^2 = {a2[1]:.3f} + ({a2[0]:.3f}) (k l_d)^2")
json.dump(dict(coulomb=res, disp=disp, fit=a2.tolist()), open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/u1_results.json", "w"), indent=1)
