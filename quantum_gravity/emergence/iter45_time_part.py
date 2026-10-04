"""
ITERATION 45 (QG): the TIME part of Einstein's action on the random grid.
With the preferred time slicing, Einstein's action needs the slice's stretching rate K_ij (extrinsic curvature):
    kinetic density = K_ij K^ij - lambda K^2   (GR: lambda = 1).
On a random grid only LINK LENGTHS between ticks are available: d ln(l)/dt = n.K.n for a link along unit direction n. So in each
small neighbourhood fit the six numbers of K_ij to the length changes of all links there (least squares = averaging).
Known test geometries on a random grid (periodic box, ~1 point per unit volume):
  (H)  uniform expansion g_ij = a(t)^2 delta_ij, H = 0.01       -> K^i_j = H delta; K_ij K^ij - K^2 = -6 H^2 (Friedmann term)
  (GW) gravitational wave h_xx = -h_yy = A cos(k z - w t), w = k   -> traceless K = h_dot/2; mean K_ij K^ij = A^2 w^2 / 4 (graviton
       kinetic term), K = 0
Expectations written BEFORE running:
 W1 (H): recovered K = 3H within 1%; kinetic term -6H^2 within 2%.
 W2 (GW): mean K_ij K^ij within 10% of A^2 w^2/4 (small bias from the wave varying across a neighbourhood, ~ (k r)^2); mean K ~ 0.
 W3 both pieces (K_ij K^ij and K^2) come out separately, so any lambda can be read off -- the time part is complete on the grid.
"""
import numpy as np
from scipy.spatial import cKDTree
rng = np.random.default_rng(45)
L = 40.0; n = int(L**3); X = rng.uniform(0, L, (n, 3)); tree = cKDTree(X, boxsize=L)
rfit = 2.0; dt = 1e-3
pairs = tree.query_pairs(rfit, output_type="ndarray")
D = X[pairs[:, 1]] - X[pairs[:, 0]]; D -= L*np.round(D/L); M = X[pairs[:, 0]] + D/2
def lengths(t, case, A=0.01, k=2*np.pi/20, H=0.01):
    if case == "H":
        a = np.exp(H*t); return a*np.linalg.norm(D, axis=1)
    hp = A*np.cos(k*M[:, 2] - k*t)
    return np.sqrt(D[:, 0]**2*(1 + hp) + D[:, 1]**2*(1 - hp) + D[:, 2]**2)
def fit_K(case, **kw):
    rate = (np.log(lengths(dt, case, **kw)) - np.log(lengths(-dt, case, **kw)))/(2*dt)
    u = D/np.linalg.norm(D, axis=1)[:, None]
    rows = np.stack([u[:, 0]**2, u[:, 1]**2, u[:, 2]**2, 2*u[:, 0]*u[:, 1], 2*u[:, 0]*u[:, 2], 2*u[:, 1]*u[:, 2]], 1)
    # normal equations per point (each link belongs to both its end points)
    AtA = np.zeros((n, 6, 6)); Atb = np.zeros((n, 6))
    for end in (0, 1):
        np.add.at(AtA, pairs[:, end], rows[:, :, None]*rows[:, None, :]); np.add.at(Atb, pairs[:, end], rows*rate[:, None])
    ok = np.linalg.cond(AtA) < 1e6
    k6 = np.linalg.solve(AtA[ok], Atb[ok][:, :, None])[:, :, 0]
    K = np.zeros((ok.sum(), 3, 3)); K[:, 0, 0], K[:, 1, 1], K[:, 2, 2] = k6[:, 0], k6[:, 1], k6[:, 2]
    K[:, 0, 1] = K[:, 1, 0] = k6[:, 3]; K[:, 0, 2] = K[:, 2, 0] = k6[:, 4]; K[:, 1, 2] = K[:, 2, 1] = k6[:, 5]
    tr = np.trace(K, axis1=1, axis2=2); KK = np.einsum("nij,nij->n", K, K)
    return tr, KK, ok.sum()
H = 0.01; A = 0.01; k = 2*np.pi/20
tr, KK, m = fit_K("H", H=H)
out = ["ITERATION 45: time part (slice stretching K_ij) read from link-length changes on the random grid (expectations first)", "",
       f"grid: {n} random points, links within {rfit} spacings (~{2*len(pairs)/n:.0f} links per point), points fitted: {m}", "",
       f"(H) uniform expansion H = {H}: mean K = {tr.mean():.6f} (exact {3*H:.6f}, ratio {tr.mean()/(3*H):.4f});",
       f"    K_ij K^ij - K^2 = {np.mean(KK - tr**2):.3e} (exact {-6*H**2:.3e}, ratio {np.mean(KK - tr**2)/(-6*H**2):.4f})"]
tr, KK, m = fit_K("GW", A=A, k=k)
ex = A**2*k**2/4
out += [f"(GW) wave A = {A}, wavelength 20: mean K_ij K^ij = {KK.mean():.3e} (exact {ex:.3e}, ratio {KK.mean()/ex:.4f});",
        f"    mean K = {tr.mean():+.2e} (exact 0), mean K^2 = {np.mean(tr**2):.2e} (vs K_ij K^ij {KK.mean():.2e})"]
txt = "\n".join(out); print(txt); open("iter45_time_part.txt", "w").write(txt + "\n")
