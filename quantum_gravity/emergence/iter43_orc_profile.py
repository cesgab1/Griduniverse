"""
ITERATION 43 (QG): build Einstein's action from neighbourhood-averaged curvature on the RANDOM grid.
Measure R_eps(r) = (30/eps^2) <kappa_curved - kappa_flat> in shells at r = 0, 0.2, 0.35, 0.5 around the bump, with common random
points (as 41b) and random pair directions (the shell average of Ric(v,v) over directions is R/3). True values (conformally flat
bump A = 0.3, s = 0.25):  R(0) = 31.6, R(0.2) = 20.7, R(0.35) = 5.3, R(0.5) = -2.6  (note the SIGN CHANGE at the edge).
Then the action over the ball r < 0.5: S = (1/2) INT R sqrt(g) dV from the measured profile vs the true profile.
Expectations written BEFORE running:
 P1 each shell within ~30% (or 2 sigma) of the true R, including the NEGATIVE sign at r = 0.5.
 P2 the action from the measured profile within ~25% of the true action.
 Kill: wrong sign at r = 0.5, or the profile shape clearly off (e.g. no fall from centre to edge).
"""
import numpy as np, ot, sys
from scipy.spatial import cKDTree
from scipy.integrate import cumulative_trapezoid
A0, s = 0.3, 0.25
def phi(x, A): return A*np.exp(-(x**2).sum(-1)/(2*s*s))
def plen(P, Q, A):
    t = np.linspace(0, 1, 5); w = np.array([1, 4, 2, 4, 1])/12.0
    return sum(wi*np.exp(phi(P + ti*(Q - P), A)) for ti, wi in zip(t, w))*np.linalg.norm(Q - P, axis=-1)
def Rtrue(r):
    ph = A0*np.exp(-r**2/(2*s*s)); lap = ph*(r**2/s**4 - 3/s**2); g2 = ph**2*r**2/s**4
    return -np.exp(-2*ph)*(4*lap + 2*g2)
rg = np.linspace(0, 1.0, 100001)
Vc = cumulative_trapezoid(4*np.pi*rg**2*np.exp(3*A0*np.exp(-rg**2/(2*s*s))), rg, initial=0)
def to_curved(X):
    r = np.linalg.norm(X, axis=1); rn = np.interp(4/3*np.pi*r**3, Vc, rg); return X*(rn/np.maximum(r, 1e-12))[:, None]
eps = 0.07; rho = float(sys.argv[1]) if len(sys.argv) > 1 else 1.2e6; npairs = int(sys.argv[2]) if len(sys.argv) > 2 else 30
rng = np.random.default_rng(43)
Rb = 0.65; n = rng.poisson(rho*4/3*np.pi*Rb**3); u = rng.normal(size=(n, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
Xf = u*(Rb*rng.uniform(0, 1, n)**(1/3))[:, None]; Xc = to_curved(Xf); tf, tc = cKDTree(Xf), cKDTree(Xc)
def kap(X, A, i, j, tree):
    def nb(k):
        c = np.array(tree.query_ball_point(X[k], eps*1.05)); d = plen(X[k][None], X[c], A); return c[(d <= eps) & (c != k)]
    a, b = nb(i), nb(j); C = plen(X[a][:, None, :], X[b][None, :, :], A)
    return 1 - ot.emd2(np.full(len(a), 1/len(a)), np.full(len(b), 1/len(b)), C, numItermax=5000000)/plen(X[i], X[j], A), len(a)
rf = np.linalg.norm(Xf, axis=1); out = []
for r0 in (0.0, 0.2, 0.35, 0.5):
    # shells are defined on the CURVED radius; flat radius of a curved shell from the volume map
    rc = np.linalg.norm(Xc, axis=1); sel = np.where(np.abs(rc - r0) < 0.012)[0] if r0 > 0 else np.where(rc < 0.02)[0]
    ds = []; m = 0
    for i in rng.permutation(sel):
        c = np.array(tf.query_ball_point(Xf[i], 0.55*eps)); c = c[np.linalg.norm(Xf[c] - Xf[i], axis=1) > 0.45*eps]
        if len(c) == 0: continue
        j = rng.choice(c)
        kf, m = kap(Xf, 0.0, i, j, tf); kc, _ = kap(Xc, A0, i, j, tc); ds.append(kc - kf)
        if len(ds) >= npairs: break
    ds = np.array(ds); Rm, Re = 30/eps**2*ds.mean(), 30/eps**2*ds.std()/np.sqrt(len(ds))
    line = f"r = {r0:4.2f}: measured R = {Rm:+7.2f} +/- {Re:5.2f}   true R = {Rtrue(r0):+7.2f}   (pairs {len(ds)}, ~{m} points/neighbourhood)"
    out.append(line); print(line, flush=True); open("iter43_orc_profile.txt", "a").write(line + "\n")
