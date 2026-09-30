"""
Tension in space AND time: does a principle fix how much space stretches relative to time (gamma), and so the bending of light?
Grid deformation: a symmetric 4x4 field h_mn (h_00 = stretch in time, h_ij = stretch in space).
Most general tension (quadratic, rotation-invariant) energy:
   L = c1 d_l h_mn d_l h_mn + c2 d_m h_mn d_l h_ln + c3 d_m h_mn d_n h + c4 d_l h d_l h        (h = trace)
PRINCIPLE (your picture made precise): sliding/relabelling the grid points without stretching them, h_mn -> h_mn + d_m xi_n + d_n xi_m,
must cost no energy. Step A: find which (c1..c4) obey it; solve for a static mass; read off gamma = (space stretch)/(time stretch).
Step B: on the random 3D grid, solve for the dimple of a point mass and bend a light ray past it: time-only tension vs the principle.
Light bending in the weak field: deflection = (1 + gamma)/2 x 4GM/(c^2 b).  Measured at the Sun (VLBI): 1.7512 arcsec, gamma - 1 = (2 +- 2) x 1e-5 (Cassini).
"""
import numpy as np, sympy as sp
print("=== A. Which tension laws respect 'relabelling the grid costs nothing'? ===")
c1, c2, c3, c4, k = sp.symbols("c1 c2 c3 c4 k", real=True)
eta = sp.diag(-1, 1, 1, 1)
kv = sp.Matrix(sp.symbols("k0:4", real=True))
# quadratic form in momentum space: S = sum h_mn(-k) M^{mn,ab}(k) h_ab(k); build M by differentiating the Lagrangian with d -> i k
H = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f"h{min(i,j)}{max(i,j)}"))
up = lambda i: eta[i, i]
def L_of(Hm):
    tr = sum(up(m) * Hm[m, m] for m in range(4))
    t1 = sum(up(l)*kv[l]**2 * sum(up(m)*up(n)*Hm[m, n]**2 for m in range(4) for n in range(4)) for l in range(4))
    div = [sum(up(m)*kv[m]*Hm[m, n] for m in range(4)) for n in range(4)]
    t2 = sum(up(n)*div[n]**2 for n in range(4))
    t3 = sum(up(n)*div[n]*kv[n] for n in range(4)) * tr
    t4 = sum(up(l)*kv[l]**2 for l in range(4)) * tr**2
    return c1*t1 + c2*t2 + c3*t3 + c4*t4
L = sp.expand(L_of(H))
syms = sorted({s for s in H}, key=lambda s: s.name)
EOM = [sp.diff(L, s) for s in syms]                       # field equations (up to factor 2)
xi = sp.Matrix(sp.symbols("x0:4", real=True))
gauge = {sp.Symbol(f"h{min(i,j)}{max(i,j)}"): kv[i]*xi[j]*up(j) + kv[j]*xi[i]*up(i) for i in range(4) for j in range(4)}
conds = set()
for e in EOM:
    ge = sp.expand(e.subs(gauge, simultaneous=True))
    for mono, coef in sp.Poly(ge, *kv, *xi).terms():
        if coef != 0: conds.add(sp.factor(coef))
sol = sp.solve(list(conds), [c2, c3, c4], dict=True)
print("   relabelling-invariance forces:", sol)
if not sol: raise SystemExit("no invariant solution")
s = sol[0]
# static point source: k = (0, k, 0, 0); ansatz h00 = A, hij = B delta_ij, h0i = 0 ; source T00 = rho (loads time only)
A, B, rho = sp.symbols("A B rho", real=True)
sub = {sp.Symbol(f"h{i}{j}"): 0 for i in range(4) for j in range(i, 4)}
sub.update({sp.Symbol("h00"): A, sp.Symbol("h11"): B, sp.Symbol("h22"): B, sp.Symbol("h33"): B})
kstatic = {kv[0]: 0, kv[1]: k, kv[2]: 0, kv[3]: 0}
E = {sy.name: sp.simplify(sp.expand(e.subs(s)).subs(kstatic).subs(sub)) for sy, e in zip(syms, EOM)}
eqs = [sp.Eq(E["h00"], rho), sp.Eq(E["h22"], 0)]          # h22 equation (transverse space) sourced by nothing
solAB = sp.solve(eqs, [A, B], dict=True)
print("   static mass, principle imposed:", solAB)
gam = sp.simplify(solAB[0][B] / solAB[0][A])
print(f"   => space stretch / time stretch = gamma = {gam}   (general relativity: gamma = 1)")
# without the principle: arbitrary tension laws
rng = np.random.default_rng(0); gs = []
for _ in range(2000):
    vals = {c1: 1.0, c2: rng.normal(0, 2), c3: rng.normal(0, 2), c4: rng.normal(0, 2)}
    Eg = {sy.name: sp.expand(e).subs(vals).subs(kstatic).subs(sub) for sy, e in zip(syms, EOM) if sy.name in ("h00", "h22")}
    try:
        sAB = sp.solve([sp.Eq(Eg["h00"], 1), sp.Eq(Eg["h22"], 0)], [A, B], dict=True)
        if sAB: gs.append(float(sAB[0][B] / sAB[0][A]))
    except Exception: pass
    if len(gs) >= 300: break
gs = np.array(gs)
print(f"   WITHOUT the principle (random tension laws): gamma ranges {np.percentile(gs,5):.2f} to {np.percentile(gs,95):.2f} (5-95%), within 1e-3 of 1 in {np.mean(abs(gs-1)<1e-3)*100:.1f}% of cases")

print("\n=== B. On the random grid: dimple of a point mass, and light bending ===")
from scipy.spatial import Voronoi
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.interpolate import griddata
g = np.random.default_rng(3); Rb, pad = 14.0, 3.0; Lh = Rb + pad
n = g.poisson((2*Lh)**3); P = g.uniform(-Lh, Lh, (n, 3)); vor = Voronoi(P)
I, J = vor.ridge_points[:, 0], vor.ridge_points[:, 1]; ell = np.linalg.norm(P[I] - P[J], axis=1)
r = np.linalg.norm(P, axis=1); inside = r < Rb; idx = -np.ones(n, int); idx[inside] = np.arange(inside.sum()); N = inside.sum()
rows, cols, vals = [], [], []; diag = np.zeros(N)
for a, b, w in zip(I, J, 1/ell):
    for s_, t in ((a, b), (b, a)):
        if inside[s_]:
            diag[idx[s_]] += w
            if inside[t]: rows.append(idx[s_]); cols.append(idx[t]); vals.append(-w)
Lap = coo_matrix((np.r_[vals, diag], (np.r_[rows, np.arange(N)], np.r_[cols, np.arange(N)])), (N, N)).tocsr()
src = np.zeros(N); near = np.argsort(r[inside])[:30]; src[near] = 1/30          # a small ball of mass at the centre
u = spsolve(Lap, src)                                          # the grid's dimple profile (shape ~ 1/r)
pts = P[inside]; ug = u
def deflection(bimp, gamma):
    zs = np.linspace(-Rb*0.8, Rb*0.8, 400)
    line = np.c_[np.full_like(zs, bimp), np.zeros_like(zs), zs]
    eps = 0.3
    up_ = griddata(pts, ug, line + [eps, 0, 0], method="linear"); dn_ = griddata(pts, ug, line - [eps, 0, 0], method="linear")
    grad = (up_ - dn_) / (2*eps)
    return -(1 + gamma) / 2 * 2 * np.trapezoid(grad, zs)       # time part contributes 1, space part gamma
for b_ in (4.0, 6.0):
    d_time = deflection(b_, 0.0); d_full = deflection(b_, float(gam))
    print(f"   ray passing at {b_} grid spacings: bending with time-only tension {d_time:.4f}, with the principle {d_full:.4f}  -> ratio {d_full/d_time:.3f}")
full = 4*6.6743e-11*1.98847e30/(2.99792458e8**2*6.957e8)*206265
print(f"   scaled to the Sun's edge: time-only {full/2:.3f} arcsec, principle {full*(1+float(gam))/2:.3f} arcsec, measured 1.7512 arcsec")
