"""
Can the random grid carry a graviton (a spin-2 wave), and does it keep the symmetry that leaves only 2 polarisations?

Part A (counting). What each kind of grid variable can carry, by helicity (spin along the direction of travel):
  one number per link/node (tension, potential)      -> helicity 0 only            (no gravitational waves)
  a displacement per node (springs, 'draw-in')        -> 0, ±1                      (sound + shear; still no graviton)
  a symmetric tensor per cell (cell SHAPE, from the
  lengths of its links: |d|^2 = d.g.d fitted per cell)-> 0,0,±1,±2                  (contains the graviton, ±2)
  GR keeps ONLY ±2 because of gauge symmetry: h -> h + D_i xi_j + D_j xi_i (moving the nodes around) changes nothing.
  Note: moving the nodes = pure gauge. Gravity proper is when the link lengths cannot be fitted by any node placement
  (curvature). This is Regge calculus (1961); linearised Regge gravitons: Rocek & Williams 1981 (BORROWED).

Part B (computed). Build derivatives on the random Mosaic (least-squares gradient, exact for linear fields) and on a
cubic lattice (central differences). Apply the linearised Einstein operator E(h) (static spatial part) to
  (1) a pure-gauge wave h = D_i xi_j + D_j xi_i  -> continuum: E = 0 exactly
  (2) a transverse-traceless (TT, helicity ±2) wave -> continuum: E = -k^2 h
Gauge residual r = |E(h_gauge)| / (k^2 |h_gauge|). If r -> 0 as k -> 0, the symmetry (and so 2 polarisations) emerges at
long wavelengths. Also check that a constant h costs nothing (graviton massless).

RESULT (graviton_test.txt): cubic lattice keeps the symmetry exactly (residual 1e-16). Random Mosaic with derivatives
does NOT: residual 0.9 at the longest wavelength and growing (k^-0.8); even the TT wave is off by 15-25%. Cause: the
method (derivatives of derivatives amplify cell-to-cell noise on an irregular grid), not physics. Conclusion: on a random
grid gravity cannot be written with derivatives; it must be written geometrically, in link lengths (regge_gauge.py).
"""
import numpy as np, os, scipy.sparse as sp
from scipy.spatial import Delaunay
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(7)
L, N = 24.0, 13824
def random_D():
    P = rng.uniform(0, L, (N, 3)); m = 2.6; pts, idx = [P], [np.arange(N)]
    for s in [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1) if (a, b, c) != (0, 0, 0)]:
        Q = P + np.array(s)*L; keep = np.all((Q > -m) & (Q < L + m), axis=1); pts.append(Q[keep]); idx.append(np.where(keep)[0])
    X = np.vstack(pts); I = np.concatenate(idx); tri = Delaunay(X)
    e = np.vstack([tri.simplices[:, [a, b]] for a in range(4) for b in range(4) if a < b]); e = np.vstack([e, e[:, ::-1]])
    e = e[e[:, 0] < N]; e = np.unique(e, axis=0)                       # edges from home points
    i, jx = e[:, 0], e[:, 1]; d = X[jx] - X[i]; j = I[jx]
    pair = np.unique(np.c_[i, j, np.round(d, 6)], axis=0); i, j, d = pair[:, 0].astype(int), pair[:, 1].astype(int), pair[:, 2:]
    M = np.zeros((N, 3, 3)); np.add.at(M, i, d[:, :, None]*d[:, None, :]); Mi = np.linalg.inv(M)
    c = np.einsum("nab,nb->na", Mi[i], d); D = []
    for a in range(3):
        A = sp.csr_matrix((c[:, a], (i, j)), shape=(N, N)); D.append(A - sp.diags(np.asarray(A.sum(1)).ravel()))
    deg = np.bincount(i, minlength=N).mean(); return P, D, deg
def cubic_D():
    n = int(round(L)); g = np.arange(n); x, y, z = np.meshgrid(g, g, g, indexing="ij"); P = np.c_[x.ravel(), y.ravel(), z.ravel()].astype(float)
    def id_(x, y, z): return ((x % n)*n + (y % n))*n + (z % n)
    D = []
    for a in range(3):
        s = np.zeros(3, int); s[a] = 1; r = id_(x, y, z).ravel(); p = id_(x + s[0], y + s[1], z + s[2]).ravel(); q = id_(x - s[0], y - s[1], z - s[2]).ravel()
        D.append(sp.csr_matrix((np.r_[np.full(n**3, .5), np.full(n**3, -.5)], (np.r_[r, r], np.r_[p, q])), shape=(n**3, n**3)))
    return P, D, 6
def E_op(h, D):
    dd = lambda a, b, f: D[a] @ (D[b] @ f); tr = h[0][0] + h[1][1] + h[2][2]
    lap = lambda f: sum(dd(k, k, f) for k in range(3)); s = sum(dd(k, l, h[k][l]) for k in range(3) for l in range(3)); ltr = lap(tr)
    return [[lap(h[i][j]) - sum(dd(i, k, h[k][j]) for k in range(3)) - sum(dd(j, k, h[k][i]) for k in range(3)) + dd(i, j, tr)
             + (s - ltr if i == j else 0) for j in range(3)] for i in range(3)]
nrm = lambda h: np.sqrt(sum(np.vdot(h[i][j], h[i][j]).real for i in range(3) for j in range(3)))
def run(name, P, D, out):
    out.append(f"\n{name}"); out.append("  k*spacing   gauge residual r   TT: E/(-k^2 h)   TT leak out of TT")
    const = [[np.ones(len(P)) if (i, j) in ((0, 1), (1, 0)) else np.zeros(len(P)) for j in range(3)] for i in range(3)]
    out.append(f"  constant h (graviton mass test): |E| = {nrm(E_op(const, D)):.1e}  (0 = massless)")
    rows = []
    for nvec in ([1, 0, 0], [1, 1, 0], [2, 1, 0], [2, 2, 1], [3, 2, 0], [4, 2, 1], [5, 3, 2], [7, 4, 2], [9, 6, 3]):
        k = 2*np.pi/L*np.array(nvec, float); kk = np.linalg.norm(k); kh = k/kk
        e1 = np.cross(kh, [0.3, 0.5, 0.81]); e1 /= np.linalg.norm(e1); e2 = np.cross(kh, e1); ph = np.exp(1j*P @ k)
        xi = [rng.normal()*ph for _ in range(3)]
        hg = [[D[i] @ xi[j] + D[j] @ xi[i] for j in range(3)] for i in range(3)]
        r = nrm(E_op(hg, D))/(kk**2*nrm(hg))
        T = np.outer(e1, e1) - np.outer(e2, e2); htt = [[T[i, j]*ph for j in range(3)] for i in range(3)]; Et = E_op(htt, D)
        proj = sum(np.vdot(htt[i][j], Et[i][j]) for i in range(3) for j in range(3)).real/(-kk**2*nrm(htt)**2)
        res = [[Et[i][j] + proj*kk**2*htt[i][j] for j in range(3)] for i in range(3)]; leak = nrm(res)/(kk**2*nrm(htt))
        rows.append((kk, r, proj, leak)); out.append(f"  {kk:8.3f}     {r:12.2e}      {proj:8.4f}        {leak:8.3f}")
    R = np.array(rows); lo = R[:, 0] < 1.2
    if np.all(R[:, 1] > 1e-10):
        p = np.polyfit(np.log(R[lo, 0]), np.log(R[lo, 1]), 1)[0]; out.append(f"  gauge residual scales as (k*spacing)^{p:.2f} at long wavelength")
    return R
out = ["Graviton test: linearised Einstein operator on lattice derivatives (spacing ~1 = one cell)"]
P, D, deg = random_D(); out.append(f"random Mosaic: {N} cells, mean neighbours {deg:.2f}")
Rr = run("RANDOM MOSAIC (least-squares derivatives)", P, D, out)
P, D, _ = cubic_D(); Rc = run("CUBIC LATTICE (central differences; for comparison)", P, D, out)
txt = "\n".join(out); print(txt); open(os.path.join(HERE, "graviton_test.txt"), "w").write(txt + "\n")
np.save(os.path.join(HERE, "graviton_rows.npy"), {"random": Rr, "cubic": Rc}, allow_pickle=True)
