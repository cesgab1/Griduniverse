"""
Fungus / slime-mould type grid (Coalesce, Oct 2026): links (tubes of fluid, dead space between) thicken where more flux flows
through them and thin where less does -- the Tero et al. (Physarum) adaptive-network rule, steady state D = F(|Q|).
Continuum steady state: flux Q = D grad p, so div(D(|grad p|) grad p) = source  -- the AQUAL (MOND) field equation.
Deep MOND (D ~ g/a0) needs reinforcement D ~ |Q|^(1/2); saturation of link thickness (D -> 1) gives Newton. Then a0 = the flux
per unit area at which links are fully built.
Test on a 3-D cubic lattice (point source of flux 4 pi G M at the centre, p = 0 on the outer sphere):
  adapt D = min( sqrt(|Q|/a0), 1 ) (hard cap) or smooth D = u/(1+u), u(1+u) = |Q|/a0  ... iterate to steady state;
  measure g(r) = flux per area along axes and diagonals (isotropy), compare with MOND for the same interpolation, and check
  v^4 ~ M (baryonic Tully-Fisher) across M.
Lattice units: spacing 1, a0 = 1, G = 1.
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
L = 41; c = L//2; R_out = c - 1
idx = np.arange(L**3).reshape(L, L, L)
X, Y, Z = np.meshgrid(*(np.arange(L) - c,)*3, indexing="ij"); Rn = np.sqrt(X**2 + Y**2 + Z**2)
inside = Rn <= R_out
edges = []
for ax in range(3):
    sl_a = [slice(None)]*3; sl_b = [slice(None)]*3; sl_a[ax] = slice(0, -1); sl_b[ax] = slice(1, None)
    a = idx[tuple(sl_a)].ravel(); b = idx[tuple(sl_b)].ravel()
    keep = inside.ravel()[a] & inside.ravel()[b]; edges.append(np.stack([a[keep], b[keep]], 1))
E = np.concatenate(edges); nE = len(E); ins = np.flatnonzero(inside.ravel())
pos = -np.ones(L**3, int); pos[ins] = np.arange(len(ins))
bnd = (Rn.ravel()[ins] > R_out - 1)                                       # Dirichlet p = 0 on the outer shell
B = sp.csr_matrix((np.r_[np.ones(nE), -np.ones(nE)], (np.r_[np.arange(nE), np.arange(nE)], np.r_[pos[E[:, 0]], pos[E[:, 1]]])),
                  shape=(nE, len(ins)))
free = ~bnd
_p0 = {}
def solve(D, M):
    Lap = (B.T @ sp.diags(D) @ B).tocsr()
    s = np.zeros(len(ins)); s[pos[idx[c, c, c]]] = 4*np.pi*M
    p = np.zeros(len(ins)); A = Lap[free][:, free]
    dg = A.diagonal(); Pre = sp.diags(1/dg)
    x, info = spl.cg(A, s[free], x0=_p0.get(M), rtol=1e-10, maxiter=20000, M=Pre)
    _p0[M] = x; p[free] = x; return p
def Dof(q, mode):
    if mode == "hard": return np.minimum(np.sqrt(q), 1.0)
    sq = np.sqrt(q); return np.maximum(sq/(1 + sq), 1e-6)            # D ~ sqrt(Q) when thin, -> 1 when fully built
def run(M, mode, it=40):
    D = np.ones(nE)
    for _ in range(it):
        p = solve(D, M); Q = np.abs(D*(B @ p)); Dn = Dof(Q, mode)
        if np.max(np.abs(Dn - D)) < 1e-6: break
        D = 0.5*D + 0.5*Dn
    return p, D
def g_profile(p):
    pf = np.zeros(L**3); pf[ins] = p; pf = pf.reshape(L, L, L)
    out = {}
    for name, dirn in (("axis", (1, 0, 0)), ("diag", (1, 1, 1))):
        dirn = np.array(dirn)/np.linalg.norm(dirn); rs = np.arange(3, R_out - 3)
        vals = []
        for r in rs:
            p1 = trilin(pf, c + (r - 0.5)*dirn); p2 = trilin(pf, c + (r + 0.5)*dirn); vals.append(p1 - p2)
        out[name] = (rs, np.array(vals))
    return out
def trilin(f, x):
    i = np.floor(x).astype(int); d = x - i; v = 0
    for a in (0, 1):
        for b in (0, 1):
            for cc in (0, 1):
                w = (d[0] if a else 1 - d[0])*(d[1] if b else 1 - d[1])*(d[2] if cc else 1 - d[2])
                v += w*f[i[0] + a, i[1] + b, i[2] + cc]
    return v
def mond_g(gN, mode):
    if mode == "hard": return np.where(gN >= 1, gN, np.sqrt(gN))
    # smooth: D = u/(1+u), u(1+u) = Q = D g -> solve D g = gN for g
    from scipy.optimize import brentq
    return np.array([brentq(lambda g: Dsteady(g)*g - x, 1e-12, 1e8) for x in gN])
def Dsteady(g):                                     # smooth rule at steady state: D = u/(1+u) with u^2 = D g -> u = g/(1+u)...
    u = (-1 + np.sqrt(1 + 4*g))/2                  # u(1+u) = g  =>  D = u/(1+u), and Q = D g = u^2  (consistent)
    return u/(1 + u)
if __name__ == "__main__":
    for mode in ("hard", "smooth"):
        print(f"\n== link rule: {mode} saturation ==")
        flats = []
        for M in (5, 10, 20, 40):
            p, D = run(M, mode); prof = g_profile(p)
            rs, ga = prof["axis"]; _, gd = prof["diag"]
            gN = M/rs**2                                                      # Newtonian (G = 1): flux 4 pi M over 4 pi r^2
            gm = mond_g(gN, mode)
            rM = np.sqrt(M)
            print(f"  M = {M:4d} (MOND radius {rM:.1f}):  r | g axis | g diag | MOND prediction | Newton")
            for j in range(0, len(rs), 2):
                print(f"      {rs[j]:3d} | {ga[j]:7.3f} | {gd[j]:7.3f} | {gm[j]:7.3f} | {gN[j]:7.3f}")
            print(f'      isotropy (axis/diag) at r = 8-12: {np.median(ga[(rs>=8)&(rs<=12)]/gd[(rs>=8)&(rs<=12)]):.2f}')
            outer = (rs > 1.5*rM) & (rs <= 12)
            if outer.sum() > 2: flats.append((M, np.median((ga*rs)[outer])**2))   # v^2 = g r -> v^4 = (g r)^2
        if len(flats) > 1:
            Ms, v4 = np.array(flats).T; sl = np.polyfit(np.log(Ms), np.log(v4), 1)[0]
            print(f"  outer v^4 vs M slope (MOND: 1, i.e. v^4 = G M a0): {sl:.2f}; v^4/M = " + ", ".join(f"{a/b:.2f}" for a, b in zip(v4, Ms)))
